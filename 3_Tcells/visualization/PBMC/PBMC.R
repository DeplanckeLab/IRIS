#Author: Joern Pezoldt
#Date init: 12.01.2024
#Scopes:
## 1) Embed all CD8 ER data processed with IRIS
## 2) Superimposte image data classification for DEG analysis and visualization

#libraries
library(dplyr)
library(plyr)
library(tidyr)
library(patchwork)
library(Seurat)
library(SeuratData)
library(UpSetR)
library(grid)
library(ggplot2)
library(topGO)
 
###
#Todos
# 1) check stripy in lymphocytes
# 2) nuclear invaginations for myeloids and lymphoids
# 3) Export respective annotation tables with n_cell per ROI
# 4) Clean script with new annotation exports
# 5) Color UMAP, DEG heatmap to 5 DEGs


#load custom functions around scRNA-seq with IRIS, imaging and tools used for analysis
source("/home/pezoldt/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/seuratV5.R")
source("/home/pezoldt/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/support_transcriptome_integration.R")

#####
#Set PATHs
#####
#UserID
userID <- "pezoldt"
analysisID <- "2024_PBMCs"
#Samples to integrate
sample_ids_IRIS <- c("JP278", "NG039", "NG042", "NG050")
sample_ids <- c(sample_ids_IRIS)
#Define paths
PATH_input_IRIS_sequencing <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis")
PATH_input_IRIS_imaging <- paste0("/home/",userID,"/updepla/projects/iris/6_imaging_analysis")

PATH_output <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis/0_Combined/",analysisID,"/results")
PATH_output_figures <- paste0(PATH_output,"/plots")
PATH_output_objects <- paste0(PATH_output,"/objects")
PATH_output_tables <- paste0(PATH_output,"/tables")

#####
#Seurat object load
#####
seurat.object <- readRDS(paste0(PATH_output_objects,"/",analysisID,"_manuscript_nuclearshapeannotation.Rds"))

########
#Signatures
########
#Path signatures
PATH_signature <- "/home/pezoldt/NAS2/iris/1_scripts/iris_scRNAseq/3_signatures_genes"
#Signatures PBMCs
signatures.PBMCs <- read.delim(paste0(PATH_signature, "/Tim_extended_PBMC_sig.txt"), sep = "\t")
# Signatures cell cycle
signature_cell_cycle_human_mouse <- readRDS(paste0(PATH_signature,
                                                   "/signature_cell_cycle_human_mouse.Rds"))

#####
#Parameters to be modified
#####
#Variables
species = "human"
# Batch over "CC" or "expID"
integrate_over = "expID"

#Threshing
min_cells_percent = 0.005
min_gene_number = 1000
mito_cutoff = 15
min_nUMI = 3000

#Set starting dimensions
dims_use = 30
#common variables to regress are: "n_UMI", "percent.mt", "Phase"
vars_to_regress <- c("n_UMI", "percent.mt")

set.seed(300)

#UMAP resolution
resolution <- 0.4

#######################################
#Image data
#######################################
#####
#Doublet data
#####
#Annotate doublets and store them in a list
l.image.doublet.ROI.annotation <- annotate_doublets(PATH_input_IRIS_imaging, 
                                                    sample_ids_IRIS)
l_doublets <- lapply(l.image.doublet.ROI.annotation, function(x){
  FCI = x$Full_Cell_ID
  cell_count = x$cell_count
  data_frame_i = data.frame(Full_Cell_ID = FCI,
                            cell_count = cell_count)
  data_frame_i
})
t_doublets <- do.call(rbind, l_doublets)

#####
#Stripy annotation data
#####
#Manual curation + cells per ROI
l_image_label_stripie <- list()
sample_ids_IRIS_ER_Nucleus <- c("JP278", "NG039", "NG042", "NG050")

for(i in 1:length(sample_ids_IRIS_ER_Nucleus)){
  #i = 1
  #Define path imaging data
  sample_ids_IRIS_i <- sample_ids_IRIS_ER_Nucleus[i]
  PATH_input_IRIS_imaging_i <- paste0(PATH_input_IRIS_imaging,"/",sample_ids_IRIS_i,"/image_annotation/manual/",sample_ids_IRIS_i,"_quantification.csv")
    #Read tables
  t_curation_stripie_state <- read.delim(PATH_input_IRIS_imaging_i)
  nrow(t_curation_stripie_state)
  head(t_curation_stripie_state)
  #Get ID
  WC <- paste0("WC",(as.numeric(unlist(lapply(strsplit(t_curation_stripie_state$Cell, "_"), function(x){x[[1]]}))))+1)
  Partial_ID <- unlist(lapply(strsplit(t_curation_stripie_state$Path, "/"), function(x){x[[9]]}))
  Exp_ID <-  unlist(lapply(strsplit(Partial_ID, "_"), function(x){x[[1]]}))
  PC <-  unlist(lapply(strsplit(Partial_ID, "_"), function(x){x[[2]]}))
  CC <- unlist(lapply(strsplit(Partial_ID, "_"), function(x){x[[3]]}))
  Full_Cell_ID <- paste(Exp_ID, PC, WC, CC ,sep="_")
  t_curation_stripie_state$Full_Cell_ID <- Full_Cell_ID
  
  n_count_ROI <- unlist(lapply(split(t_curation_stripie_state, t_curation_stripie_state$Full_Cell_ID), nrow))
  t_count <- data.frame(Full_Cell_ID = names(n_count_ROI), n_cell_ROI = n_count_ROI)
  nrow(t_count)
  
  #Eliminate duplicates in imaging data, aka multiplet encapsulations
  t_curation_stripie_state <- t_curation_stripie_state[!(duplicated(t_curation_stripie_state$Full_Cell_ID)),]
  nrow(t_curation_stripie_state)
  t_curation_stripie_state <- merge(t_curation_stripie_state, t_count, by = "Full_Cell_ID")
  t_curation_stripie_state <- t_curation_stripie_state[,c("Full_Cell_ID","Stripie","n_cell_ROI")]
  colnames(t_curation_stripie_state) <- c("Full_Cell_ID","Stripie_manual","n_cell_ROI")
  l_image_label_stripie[[i]] <- t_curation_stripie_state
}
names(l_image_label_stripie) <- sample_ids_IRIS_ER_Nucleus
lapply(l_image_label_stripie, head)
unlist(lapply(l_image_label_stripie, function(x){
  lapply((split(x, x$Stripie_manual)),nrow)}))
t_stripy_annotation <- do.call(rbind, l_image_label_stripie)

#Cleanup stripy annotation
t_stripy_annotation$ExpID <- unlist(lapply(strsplit(t_stripy_annotation$Full_Cell_ID, "_"), function(x){x[[1]]}))
t_stripy_annotation$Stripie_manual <- str_replace(t_stripy_annotation$Stripie_manual, "stripy_nucleus", "nucleus_indent")
t_stripy_annotation$Stripie_manual <- str_replace(t_stripy_annotation$Stripie_manual, "normi_circle", "nucleus_ball")
t_stripy_annotation$Stripie_manual <- str_replace(t_stripy_annotation$Stripie_manual, "normi_demicircle", "Undefined")
t_stripy_annotation$Stripie_manual <- str_replace(t_stripy_annotation$Stripie_manual, "stripy_center", "Undefined")
colnames(t_stripy_annotation) <- c("Full_Cell_ID","Annotation_Nuclear_Shape","n_cell_ROI","ExpID")
rownames(t_stripy_annotation) <- NULL

#####
#Load Count matrix IRIS platform
#####
#LiquideDrop data
l_pivot_table <- list()
l_UMI_all_cells_SYMB <- list()
l_SM_annotation <- list()

for(i in 1:length(sample_ids_IRIS)){
  #i = 1
  print(i)
  print(sample_ids_IRIS[i])
  
  #Load files
  load(paste0(PATH_input_IRIS_sequencing,"/",sample_ids_IRIS[i],"/pipeline_output/",species,"/analysis/objects/objects_plotting.RData"), envir = parent.frame(), verbose = FALSE)
  #Select cells of species choice
  SM_annotation_i <- read.delim(paste0(PATH_input_IRIS_sequencing,"/",sample_ids_IRIS[i],"/results/tables/",sample_ids_IRIS[i],"_transcriptome_species_annotation.txt"))
  Cells_wanted_i <- subset(SM_annotation_i, species_annotation == species)$Full_Cell_ID
  print(paste0(nrow(pivot_table), "| cells detected"))
  #Pivot table
  pivot_table_i <- subset(pivot_table, Full_Cell_ID %in% Cells_wanted_i)
  print(paste0(nrow(pivot_table_i), "| human cells"))
  #Select cells with min_UMI
  pivot_table_i <- subset(pivot_table_i, n_UMI > min_nUMI)
  print(paste0(nrow(pivot_table_i), "| cells >3000UMI"))
  
  #Select encapsulation with only one cell in ROI
  t_cells_per_ROI_i <- l_doublets[[i]]
  t_cells_per_ROI_i <- subset(t_cells_per_ROI_i, Full_Cell_ID %in% pivot_table_i$Full_Cell_ID)
  Cells_wanted_i <- subset(t_cells_per_ROI_i, cell_count == "1")$Full_Cell_ID
  #Pivot table
  pivot_table_i <- subset(pivot_table_i, Full_Cell_ID %in% Cells_wanted_i)
  print(paste0(nrow(pivot_table_i), "| cells single encapsulations"))
  
  #Count matrix
  #ncol(UMI_all_cells_SYMB)
  UMI_all_cells_SYMB_i <- UMI_all_cells_SYMB[,colnames(UMI_all_cells_SYMB) %in% pivot_table_i$Full_Cell_ID]
  #ncol(UMI_all_cells_SYMB_i)                     
  #Store data
  l_SM_annotation[[i]] <- SM_annotation_i
  l_pivot_table[[i]] <- pivot_table_i
  l_UMI_all_cells_SYMB[[i]] <- UMI_all_cells_SYMB_i
}
#names to lists
names(l_SM_annotation) <- sample_ids_IRIS
names(l_pivot_table) <- sample_ids_IRIS
names(l_UMI_all_cells_SYMB) <- sample_ids_IRIS
lapply(l_pivot_table, nrow)
lapply(l_pivot_table, ncol)
lapply(l_UMI_all_cells_SYMB, nrow)
lapply(l_UMI_all_cells_SYMB, ncol)
pivot_all <- do.call(rbind, l_pivot_table)

#QC on a glance
print("meadian(mappedRead)")
lapply(l_pivot_table, function(x){median(x$n_mapped_reads)})
print("meadian(nUMI)")
lapply(l_pivot_table, function(x){median(x$n_UMI)})
print("mean(nUMI)")
lapply(l_pivot_table, function(x){mean(x$n_UMI)})
print("sd(nUMI)")
lapply(l_pivot_table, function(x){sd(x$n_UMI)})
print("median(nGenes)")
lapply(l_pivot_table, function(x){median(x$n_genes)})
print("nCells")
unlist(lapply(l_pivot_table, nrow))


#####
#Prep meta data
#####
#make meta data
lapply(l_pivot_table, nrow)
meta_rna <- do.call(rbind, l_pivot_table)
meta_image_stripy <- t_stripy_annotation

#meta_image_ROI <- t_count_cells_ROI
meta_data <- merge(meta_rna, meta_image_stripy, by = "Full_Cell_ID",
                   all.x = TRUE)
rownames(meta_data) <- meta_data$Full_Cell_ID
#Replace non annotated cells
meta_data <- as.data.frame(meta_data)
meta_data$Stripie_manual <- replace_na(meta_data$Stripie_manual, "Undefined")
meta_data$n_cell_ROI <- replace_na(meta_data$n_cell_ROI, 1)

names(meta_data)[names(meta_data) == "Sample_Plate"] <- "expID"
levels(as.factor(meta_data$expID))

#####
#Prep matrix
#####
matrix <- l_UMI_all_cells_SYMB[[1]]
for(i in 1:(length(l_UMI_all_cells_SYMB)-1)){
  print(names(l_UMI_all_cells_SYMB)[i])
  print(names(l_UMI_all_cells_SYMB)[i+1])
  matrix <- merge(matrix, l_UMI_all_cells_SYMB[[i+1]], by=0, all=TRUE)
  rownames(matrix) <- matrix$Row.names
  matrix <- matrix[,2:ncol(matrix)]
  rownames(matrix)[grep("^MT.", rownames(matrix))]
}
sum(unlist(lapply(l_UMI_all_cells_SYMB, ncol)))
ncol(matrix)
nrow(matrix)
#Fill NAs
matrix[is.na(matrix)] <- 0
matrix <- as.matrix(matrix)

#Fix "MT-"
rownames(matrix)[grep("^MT\\.", rownames(matrix))]
rownames(matrix) <- gsub("^MT\\.", "MT-", rownames(matrix))
rownames(matrix)[grep("^MT-", rownames(matrix))]

#Select genes with defined threshold
nrow(matrix)
#matrix <- matrix[rowSums(matrix)>200,]
ncol(matrix)

#####
#Batch integration with Seurat
#####
#Seurat
seurat.object <- seurat_integrate(
  signature_cell_cycle_human_mouse,
  species,
  matrix,
  meta_data,
  integrate_over,
  min_cells_percent,
  min_gene_number,
  mito_cutoff,
  dims_use,
  resolution,
  vars_to_regress)

 #Replace NAs in imaging data
seurat.object@meta.data$Annotation_Nuclear_Shape <- sub("^$", "Undefined", seurat.object@meta.data$Annotation_Nuclear_Shape)
seurat.object@meta.data$Annotation_Nuclear_Shape <- seurat.object@meta.data$Annotation_Nuclear_Shape %>% replace_na('Undefined')


######
#Cluster annotation
######
seurat.object <- RenameIdents(object = seurat.object,
                              `0` = "CD8_CMlike",
                              `1` = "CD4_naive",
                              `2` = "CD4_mem",
                              `3` = "CD14_mono",
                              `4` = "CD8_EM",
                              `5` = "CD8_naive",
                              `6` = "NK",
                              `7` = "BC",
                              `8` = "FCGR3A_mono",
                              `9` = "DC")
seurat.object@meta.data$Cluster_names <- as.character(seurat.object@active.ident)

#####
#UMAP
#####
DimPlot(seurat.object,
        reduction = "umap.cca",
        combine = FALSE, label.size = 5,
        label = TRUE)
DimPlot(seurat.object,
        reduction = "umap.cca",
        combine = FALSE, label.size = 5,
        label = TRUE,
        group.by = "Phase",
        cols = c("red","green","orange"))
DimPlot(seurat.object,
        reduction = "umap.cca",
        combine = FALSE, label.size = 5,
        label = TRUE,
        cols = c("#91f089","#b3cde0","#6497b1",
                 "purple","#48bf53","#11823b",
                 "grey50","black","#D6BEF1",
                 "#B0ABE7"))

DimPlot(seurat.object,
        reduction = "umap.cca",
        label = FALSE,
        pt.size = 0.3,
        split.by = "expID",
        cols = c("#91f089","#b3cde0","#6497b1",
                 "purple","#48bf53","#11823b",
                 "grey50","black","#D6BEF1",
                 "#B0ABE7"),
        ncol = 2)

VlnPlot(seurat.object, features = c("nFeature_RNA"),
        pt.size = 0,
        #flip = FALSE,
        split.by = "expID")
FeaturePlot(seurat.object, features = c("nCount_RNA", "percent.mt","percent.Hsp","nFeature_RNA"))
FeaturePlot(seurat.object, features = c("CD8A", "SELL","CCR7","IL7R","CD4","FCGR3A"))
FeaturePlot(seurat.object, features = c("TOP2A"))
FeaturePlot(seurat.object, features = c("IL7R", "LTB","CD247","IL32",
                                        "LYZ","STAT4","BACH2","RORA",
                                        "IFITM3"))

DimPlot(seurat.object,
        reduction = "umap.cca",
        #label = TRUE,
        pt.size = 0.5,
        group.by = "Annotation_Nuclear_Shape",
        order = c("nucleus_indent","nucleus_ball","Undefined"),
        cols = c("grey95","blue","magenta"))

#####
#Summarize distribution states per cluster
#####
meta <- seurat.object@meta.data
library(plyr)
counts <- ddply(meta, .(meta$Cluster_names, meta$Annotation_Nuclear_Shape), nrow)
names(counts) <- c("Cell_Type", "Annotation_Nuclear_Shape", "Freq")

#Plot
ggplot(counts, aes(x = Cell_Type, y = Freq, fill = Annotation_Nuclear_Shape)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  #geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("nucleus_ball" = "deepskyblue","nucleus_indent" = "orange", "Undefined" = "grey70"))
#Plot
counts_tidy <- subset(counts, Annotation_Nuclear_Shape %in% c("nucleus_ball","nucleus_indent"))
#counts_tidy$Cell_Type <- as.factor(counts_tidy$Cell_Type)
counts_tidy$Cell_Type <- factor(counts_tidy$Cell_Type, levels=c("CD4_mem","CD4_naive","CD8_CMlike","CD8_EM","CD8_naive","NK","BC","FCGR3A_mono","CD14_mono","DC"))

ggplot(counts_tidy, aes(x = Cell_Type, y = Freq, fill = Annotation_Nuclear_Shape)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  #geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_fill_manual(values = c("nucleus_ball" = "blue","nucleus_indent" = "magenta")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(c("CD4_mem","CD4_naive","CD8_CMlike","CD8_EM","CD8_naive","NK","BC","FCGR3A_mono","CD14_mono","DC"))
#Lymphocytes vs myeloids
counts_tidy_lympho_ball <- subset(counts_tidy,
                             Cell_Type %in% c("BC","CD4_mem","CD4_naive","CD8_CMlike","CD8_EM","CD8_naive","NK") &
                             Annotation_Nuclear_Shape == "nucleus_ball")
counts_tidy_lympho_indent <- subset(counts_tidy,
                                  Cell_Type %in% c("BC","CD4_mem","CD4_naive","CD8_CMlike","CD8_EM","CD8_naive","NK") &
                                    Annotation_Nuclear_Shape == "nucleus_indent")
counts_tidy_mye_ball <- subset(counts_tidy,
                                  Cell_Type %in% c("FCGR3A_mono","CD14_mono","DC") &
                                    Annotation_Nuclear_Shape == "nucleus_ball")
counts_tidy_mye_indent <- subset(counts_tidy,
                                    Cell_Type %in% c("FCGR3A_mono","CD14_mono","DC") &
                                      Annotation_Nuclear_Shape == "nucleus_indent")
counts_mainsubsets <- data.frame(Cell_Type = c("Lymphocytes cells","Lymphocytes cells",
                                               "Myeloid cells","Myeloid cells"),
                                 Annotation_Nuclear_Shape = c("nucleus_ball","nucleus_indent",
                                                              "nucleus_ball","nucleus_indent"),
                                 Freq = c(sum(counts_tidy_lympho_ball$Freq), sum(counts_tidy_lympho_indent$Freq),
                                          sum(counts_tidy_mye_ball$Freq), sum(counts_tidy_mye_indent$Freq)))
#Plot
ggplot(counts_tidy, aes(x = Cell_Type, y = Freq, fill = Annotation_Nuclear_Shape)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  #geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_fill_manual(values = c("nucleus_ball" = "blue","nucleus_indent" = "magenta")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
