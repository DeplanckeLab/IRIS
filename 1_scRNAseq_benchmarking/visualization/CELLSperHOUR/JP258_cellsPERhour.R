#Author: Joern Pezoldt
#Purpose: Annotate mouse and human cells in mixed datasets

##############
#Figure 2
##############

#load custom functions around scRNA-seq with IRIS, imaging and tools used for analysis
source("/home/pezoldt/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/seuratV5.R")
source("/home/pezoldt/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/support_transcriptome_integration.R")

#libraries
library(dplyr)
library(patchwork)
library(Seurat)
library(UpSetR)
library(grid)

#####
#Set PATHs
#####
#UserID
userID <- "pezoldt"
#Samples to integrate
sample_ids_IRIS <- c("JP258")
sample_ids <- c(sample_ids_IRIS)
#Define paths sequencing data
PATH_input <- paste0("/home/",userID,"/updepla/projects/iris/")
PATH_input_IRIS_sequencing <- paste0(PATH_input,"/4_sequencing_analysis/",sample_ids_IRIS)
PATH_output <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis/",sample_ids_IRIS,"/results")
PATH_output_figures <- paste0(PATH_output,"/plots")
PATH_output_objects <- paste0(PATH_output,"/objects")
PATH_output_tables <- paste0(PATH_output,"/tables")

#Path signatures
PATH_signature <- "/home/pezoldt/NAS2/iris/1_scripts/iris_scRNAseq/3_signatures_genes"

#####
#Parameters Seurat
#####
#Variables
# Batch over "CC" or "expID"
species = "human"
integrate_over = "CC"

#Threshing
min_cells_percent = 0.005
min_gene_number = 2000
mito_cutoff = 15
min_nUMI = 3000
#Seurat
dims_use = 25
resolution = 0.4

#common variables to regress are: "n_UMI", "percent.mt", "Phase"
vars_to_regress <- c()

#####
#Precomputed
#####
#Seurat object
#seurat.object <- readRDS(paste0(PATH_output_objects,"/NG027_all_CD8.Rds"))
#Signatures
#  Cell cycle
signature_cell_cycle_human_mouse <- readRDS(paste0(PATH_signature,"/signature_cell_cycle_human_mouse.Rds"))
if(species == "human"){
  signature_cell_cycle <- readRDS(paste0(PATH_signature,"/signature_cell_cycle_",species,".Rds"))
  mito_genes = "^MT-"
  ribo_genes = "^RPL"
  heatshock_genes = "^HSP"
}
if(species == "mouse"){
  signature_cell_cycle <- readRDS(paste0(PATH_signature,"/signature_cell_cycle_",species,".Rds"))
  mito_genes = "^mt-"
  ribo_genes = "^Rpl"
  heatshock_genes = "^Hsp"
}

#####
#Load Count matrix IRIS platform
#####
#LiquideDrop data
l_pivot_table <- list()
l_UMI_all_cells_SYMB <- list()
for(i in 1:length(sample_ids_IRIS)){
  #i = 1
  print(i)
  print(sample_ids_IRIS[i])
  
  #Load files
  load(paste0(PATH_input_IRIS_sequencing,"/pipeline_output/",species,"/analysis/objects/objects_plotting.RData"), envir = parent.frame(), verbose = FALSE)

  ## Image annotation
  #Select cells with min_UMI
  pivot_table_i <- subset(pivot_table, n_UMI > min_nUMI)
  nrow(pivot_table_i)
  
  #Count matrix
  #ncol(UMI_all_cells_SYMB)
  UMI_all_cells_SYMB_i <- UMI_all_cells_SYMB[,colnames(UMI_all_cells_SYMB) %in% pivot_table_i$Full_Cell_ID]
  #ncol(UMI_all_cells_SYMB_i)                     
  #Store data
  print(nrow(pivot_table_i))
  l_pivot_table[[i]] <- pivot_table_i
  l_UMI_all_cells_SYMB[[i]] <- UMI_all_cells_SYMB_i
}
#names to lists
names(l_pivot_table) <- sample_ids_IRIS
names(l_UMI_all_cells_SYMB) <- sample_ids_IRIS
lapply(l_pivot_table, nrow)
lapply(l_pivot_table, ncol)
lapply(l_UMI_all_cells_SYMB, nrow)
lapply(l_UMI_all_cells_SYMB, ncol)
pivot_table_all <- do.call(rbind, l_pivot_table)

#QC on a glance
print("meadian(mappedRead)")
lapply(l_pivot_table, function(x){median(x$n_mapped_reads)})
print("meadian(nUMI)")
lapply(l_pivot_table, function(x){mean(x$n_mapped_reads)})
print("mean(nUMI)")
lapply(l_pivot_table, function(x){median(x$n_UMI)})
print("s(nUMI)")
lapply(l_pivot_table, function(x){sd(x$n_UMI)})
print("median(nGenes)")
lapply(l_pivot_table, function(x){median(x$n_genes)})
print("nCells")
unlist(lapply(l_pivot_table, nrow))


#####
#Stats
#####
meta <- l_pivot_table[[1]]
#UMI per cell
mean_UMI <- unlist(lapply(split(meta, meta$CC), function(x){mean(x$n_UMI)}))
mean(mean_UMI)
sd(mean_UMI)
mean_genes <- unlist(lapply(split(meta, meta$CC), function(x){mean(x$n_genes)}))
mean(mean_genes)
sd(mean_genes)

#####
#Prep meta data
#####
#make meta data
meta_data <- l_pivot_table[[1]]
rownames(meta_data) <- meta_data$Full_Cell_ID
nrow(meta_data)

#####
#Prep matrix
#####
matrix <- l_UMI_all_cells_SYMB[[1]]
ncol(matrix)

#####
#Integration
#####
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


#####
#Plot along time
#####
#Add time
meta <- seurat.object@meta.data
#Order meta
head(meta[with(meta, order(CC)), ])
tail(meta[with(meta, order(CC)), ])
meta <- meta[with(meta, order(CC)), ]

time_per_cell <- round(4000/nrow(meta),2)
#calculate time
l_time_of_cell <- list()
time_i = 0
for(i in 1:(nrow(meta))){
  time_i = i*time_per_cell
  l_time_of_cell[[i]] <- time_i
}
tail(unlist(l_time_of_cell))
meta$time_of_cell <- round(unlist(l_time_of_cell),1)
nrow(subset(meta, time_of_cell <= 3600))

#Mitochondria
ggplot(meta, aes(x = time_of_cell, y = percent.mt)) +
  geom_point(aes(color = CC), size = 0.1, ) +
  geom_smooth(method = "loess", n = 10, colour = "deepskyblue4")  +
  scale_color_manual(values = c("CC13" = "grey60", "CC14" = "grey60")) +
  scale_x_continuous(breaks = c(600,1200,1800,
                                2400,3000,3600)) +
  theme_bw() + guides(color = 'none') +
  ylim(0,15)

#Heatshock proteins
ggplot(meta, aes(x = time_of_cell, y = percent.Hsp)) +
  geom_point(aes(color = CC), size = 0.2, ) +
  geom_smooth(method = "loess", n = 10, colour = "deepskyblue4")  +
  scale_color_manual(values = c("CC13" = "grey60", "CC14" = "grey60")) +
  scale_x_continuous(breaks = c(600,1200,1800,
                                2400,3000,3600)) +
  theme_bw() + guides(color = 'none') +
  ylim(0,4)

#Ribosomal protein
ggplot(meta, aes(x = time_of_cell, y = percent.Rpl)) +
  geom_point(aes(color = CC), size = 0.2, ) +
  geom_smooth(method = "loess", n = 10, colour = "black")  +
  scale_color_manual(values = c("CC13" = "grey60", "CC14" = "grey60")) +
  scale_x_continuous(breaks = c(600,1200,1800,
                                2400,3000,3600)) +
  theme_bw() + guides(color = 'none') +
  ylim(0,15)

#Number Genes
ggplot(meta, aes(x = time_of_cell, y = nFeature_RNA)) +
  geom_point(aes(color = CC), size = 0.2, ) +
  geom_smooth(method = "loess", n = 10, colour = "black")  +
  scale_color_manual(values = c("CC13" = "deepskyblue", "CC14" = "deepskyblue4")) +
  scale_x_continuous(breaks = c(600,1200,1800,
                                2400,3000,3600)) +
  scale_y_log10(limits = c(2000,10000),
                breaks = c(2000,7000,10000)) +
  theme_bw() + guides(color = 'none') 

#Number UMI
ggplot(meta, aes(x = time_of_cell, y = n_UMI)) +
  geom_point(aes(color = CC), size = 0.2, ) +
  geom_smooth(method = "loess", n = 10, colour = "black")  +
  scale_color_manual(values = c("CC13" = "deepskyblue", "CC14" = "deepskyblue4")) +
  scale_x_continuous(breaks = c(600,1200,1800,
                                2400,3000,3600)) +
  scale_y_log10(limits = c(1000,1e5),
                breaks = c(1000,10000,50000)) +
  theme_bw() + guides(color = 'none')

#####
#Scatter UMI|nGenes per mapped read
#####
ggplot(meta, aes(x=n_mapped_reads, y=n_all_reads, color=CC)) +
  scale_color_manual(values=c("deepskyblue","deepskyblue4")) +
  geom_point(size = 1.0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  #labs(title=title_i) + 
  geom_smooth(method = "lm",
              formula = y ~ x,
              se = FALSE)
ggplot(meta, aes(x=n_mapped_reads, y=n_UMI, color=CC)) +
  scale_color_manual(values=c("deepskyblue","deepskyblue4")) +
  geom_point(size = 1.0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  #labs(title=title_i) + 
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = FALSE)
ggplot(meta, aes(x=n_mapped_reads, y=n_genes, color=CC)) +
  scale_color_manual(values=c("deepskyblue","deepskyblue4")) +
  geom_point(size = 1.0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  #labs(title=title_i) + 
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = FALSE)
ggplot(meta, aes(x=n_mapped_reads, y=ratio_UMI_to_mapped_reads, color=CC)) +
  scale_color_manual(values=c("deepskyblue","deepskyblue4")) +
  geom_point(size = 1.0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  #labs(title=title_i) + 
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = FALSE)
