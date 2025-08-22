
#Author: Joern Pezoldt
#Date: 14.03.2023
#Function: Annotate mouse and human cells in mixed datasets



#libararies
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)

#####
#Global Variables
#####
## Cutoffs
mapped_genome <- "human_mouse"
#options: "human", "mouse", "droso"
#Mote: Only a vector with two elements permitted, second item is spike-in
species_in_Exp <- c("human","mouse")
min_UMI <- 5000
#Percent of species mixing allowed to be considered a "clean" cell
thresh_sm <- 7.5

## PATHs
#General experimental ID used throughout analysis
expID <- "JP310"
userID <- "pezoldt"
#FolderName in experimental folder
simple_ID <- "JP310"

#Output name of scRNA-seq pipeline
#Automatic paths
PATH_input <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis/",expID,"/pipeline_output/",mapped_genome)
PATH_output <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis/",expID,"/results")
PATH_output_figures <- paste0(PATH_output,"/plots")
PATH_output_objects <- paste0(PATH_output,"/objects_transcriptome")
PATH_output_tables <- paste0(PATH_output,"/tables")

#Imaging data
PATH_imaging <- paste0("/home/",userID,"/updepla/projects/iris/6_imaging_analysis/",expID)
PATH_imaging_yolo <- paste0(PATH_imaging,"/image_annotation/",expID,"_quantification.csv")

#####
#Imaging data
#####
# YOLO quantification doublets
t_imaging_yolo <- read.table(PATH_imaging_yolo, sep = ",", header = TRUE)
#Get ID
WC <- unlist(lapply(strsplit(t_imaging_yolo$cell_path, "/"), function(x){x[[9]]}))
WC <- as.numeric(lapply(strsplit(WC, "_"), function(x){x[[1]]}))
WC <- paste0("WC", WC+1)
Partial_ID <- unlist(lapply(strsplit(t_imaging_yolo$cell_path, "/"), function(x){x[[8]]}))
Exp_ID <-  unlist(lapply(strsplit(Partial_ID, "_"), function(x){x[[1]]}))
PC <-  unlist(lapply(strsplit(Partial_ID, "_"), function(x){x[[2]]}))
CC <- unlist(lapply(strsplit(Partial_ID, "_"), function(x){x[[3]]}))
Full_Cell_ID <- paste(Exp_ID, PC, WC, CC ,sep="_")
t_imaging_yolo$Full_Cell_ID <- Full_Cell_ID

#Select singlets and multiplets
t_imaging_yolo_doublets <- t_imaging_yolo[(duplicated(t_imaging_yolo$Full_Cell_ID)),]
t_imaging_yolo_doublets <- data.frame(Full_Cell_ID = t_imaging_yolo_doublets$Full_Cell_ID,
                                      Fluo488 = t_imaging_yolo_doublets$measurement_mean_intensity_seg_488_mean,
                                      Fluo630 = t_imaging_yolo_doublets$measurement_mean_intensity_seg_638_mean,
                                      n_cells = rep(2, nrow(t_imaging_yolo_doublets)))
t_imaging_yolo_doublets <- t_imaging_yolo_doublets[!(duplicated(t_imaging_yolo_doublets$Full_Cell_ID)),]
subset(t_imaging_yolo_doublets, Full_Cell_ID == "JP308_P3_WC38_CC11")

t_imaging_yolo_singlets <- t_imaging_yolo[!(duplicated(t_imaging_yolo$Full_Cell_ID)),]
subset(t_imaging_yolo_singlets, Full_Cell_ID == "JP308_P3_WC38_CC11")
t_imaging_yolo_singlets <- subset(t_imaging_yolo_singlets, !(Full_Cell_ID %in% t_imaging_yolo_doublets$Full_Cell_ID))

t_imaging_yolo_singlets <- data.frame(Full_Cell_ID = t_imaging_yolo_singlets$Full_Cell_ID,
                                      Fluo488 = t_imaging_yolo_singlets$measurement_mean_intensity_seg_488_mean,
                                      Fluo630 = t_imaging_yolo_singlets$measurement_mean_intensity_seg_638_mean,
                                      n_cells = 1)
t_imaging_counts <- rbind(t_imaging_yolo_doublets, t_imaging_yolo_singlets)
#Calculate plotting score
score_Fluo <- t_imaging_counts$Fluo630 - t_imaging_counts$Fluo488
hist(score_Fluo)

#####
#Load sequencing data
#####
load(paste0(PATH_input,"/analysis/objects/objects_plotting.RData"), envir = parent.frame(), verbose = FALSE)
ucm_lid <- UMI_all_cells_ENSG
colnames(ucm_lid) <- colnames(UMI_all_cells_SYMB)
sum(as.numeric(colSums(ucm_lid) > min_UMI))
ucm_lid <- ucm_lid[,colSums(ucm_lid) > min_UMI]
ncol(ucm_lid)
mean(colSums(ucm_lid))
median(colSums(ucm_lid))
print("CCs detected:")
levels(as.factor(pivot_table$CC))

#####
#Compile Metadata for species mixing
#####
#Split Table according to species
genes_human <- ucm_lid[grepl("^ENSG",rownames(ucm_lid)),]
nrow(genes_human)
genes_mouse <- ucm_lid[grepl("^ENSM",rownames(ucm_lid)),]
nrow(genes_mouse)
#Count UMIs per Species
nUMI_mouse <- colSums(genes_mouse)
nUMI_human <- colSums(genes_human)

#Count nGene per Species
nGene_mouse <- colSums(genes_mouse > 0)
nGene_human <- colSums(genes_human > 0)

#build table with UMI and Full_Cell_ID data
to_pivot_ID_UMI_SM <- data.frame(Full_Cell_ID = names(nUMI_mouse),
                                 nUMI_mouse = nUMI_mouse, 
                                 nUMI_human = nUMI_human, 
                                 nGene_mouse = nGene_mouse, 
                                 nGene_human = nGene_human)
pivot_UMI_tresh_double <- merge(pivot_table, to_pivot_ID_UMI_SM, by = "Full_Cell_ID")

#Merge with cell number data
pivot_UMI_tresh_double <- merge(pivot_UMI_tresh_double, t_imaging_counts,
              by = "Full_Cell_ID")
pivot_UMI_tresh_double$n_cells <- as.factor(pivot_UMI_tresh_double$n_cells)
nrow(subset(pivot_UMI_tresh_double, n_cells == 2))

#####
#Calculate scores and categorization for SM
#####
l_sm_percent <- list()
l_sm_score <- list()
k = 2
for(i in 1:length(species_in_Exp)){
  #i = 1
  #Calculate extend of species mixing
  species_i <- species_in_Exp[i]
  sm_percent_i <- round(pivot_UMI_tresh_double[,paste0("nUMI_",species_i)] / pivot_UMI_tresh_double$n_UMI * 100, 1)
  l_sm_percent[[i]] <- sm_percent_i
  #Score which species or if mixed
  species_score_i <- as.numeric(sm_percent_i > (100 - thresh_sm))
  species_score_i <- replace(species_score_i, species_score_i==1, k)
  l_sm_score[[i]] <- species_score_i
  k = k-1
}
names(l_sm_percent) <- species_in_Exp
names(l_sm_score) <- species_in_Exp
sm_score <- l_sm_score[[1]] + l_sm_score[[2]]

#label with character vector
sm_score <- replace(sm_score, sm_score==0, "sm")
sm_score <- replace(sm_score, sm_score==2, species_in_Exp[1])
sm_score <- replace(sm_score, sm_score==1, species_in_Exp[2])

pivot_UMI_tresh_double$sm_annotation <- sm_score

#####
#Calculate Fluo Score
#####
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
pivot_UMI_tresh_double$Fluo488_center <- rescale(pivot_UMI_tresh_double$Fluo488)
pivot_UMI_tresh_double$Fluo630_center <- rescale(pivot_UMI_tresh_double$Fluo630)
pivot_UMI_tresh_double$Fluo488_630_centered <- pivot_UMI_tresh_double$Fluo488_center - pivot_UMI_tresh_double$Fluo630_center


#####
#Scatter Transcriptome
#####
theme_set(theme_bw())
pivot_UMI_tresh_double_doublet <- subset(pivot_UMI_tresh_double, n_cells == 2)

ggplot(pivot_UMI_tresh_double,
       aes_string(x = paste0("nUMI_",species_in_Exp[1]),
                  y = paste0("nUMI_",species_in_Exp[2]))) + 
  geom_point(aes(col = Fluo488_630_centered),size=2) +  # draw points
  scale_colour_gradient2(low = "deepskyblue3",
                        midpoint = 0, mid = "white",
                        high="green") + 
  #xlim(c(0, 20000)) + 
  #ylim(c(0, 180000)) +  
  labs(x = paste0("nUMI_",species_in_Exp[1]),
       y = paste0("nUMI_",species_in_Exp[2]),
       title = paste0(simple_ID))

#pivot_UMI_tresh_double_doubletsPlot <- pivot_UMI_tresh_double[order(pivot_UMI_tresh_double$n_cells),]
ggplot(pivot_UMI_tresh_double,
       aes_string(x = paste0("nUMI_",species_in_Exp[1]),
                  y = paste0("nUMI_",species_in_Exp[2]))) + 
  geom_point(aes(col = Fluo488_630_centered),size=2) +  # draw points
  scale_color_gradientn(colors = c("deepskyblue3","deepskyblue2","lightblue", "white","green","#8DC71E", "#69B41E")) +
  geom_point(data = pivot_UMI_tresh_double_doublet, 
             aes_string(x = paste0("nUMI_",species_in_Exp[1]),
                        y = paste0("nUMI_",species_in_Exp[2])),
             color = "#8300c4", size = 2) + 
  #xlim(c(0, 20000)) + 
  #ylim(c(0, 180000)) +  
  labs(x = paste0("nUMI_",species_in_Exp[1]),
       y = paste0("nUMI_",species_in_Exp[2]),
       title = paste0(simple_ID))

#####
#Percentage mixing
#####
percent_human <- round(pivot_UMI_tresh_double$nUMI_human / pivot_UMI_tresh_double$n_UMI * 100, 1)
percent_mouse <- round(pivot_UMI_tresh_double$nUMI_mouse / pivot_UMI_tresh_double$n_UMI * 100, 1)
t_percent <- data.frame(percent_human = percent_human, percent_mouse = percent_mouse)

pivot_UMI_tresh_double$percent_UMI_human <- percent_human
pivot_UMI_tresh_double$percent_UMI_mouse <- percent_mouse

sm_percent_human <- mean(subset(pivot_UMI_tresh_double, sm_annotation == "human")$percent_UMI_human)
sm_percent_mouse <- mean(subset(pivot_UMI_tresh_double, sm_annotation == "mouse")$percent_UMI_mouse)
sm_percent_sm <- mean(subset(pivot_UMI_tresh_double, sm_annotation == "sm")$percent_UMI_mouse)
sm_percent_human
sm_percent_mouse
sm_percent_sm

#####
#Save data
#####
#Save
write.table(pivot_UMI_tresh_double, paste0(PATH_output_tables,"/",simple_ID,"_manuscript_transcriptome_image_annotation_all_toCamille.txt"), 
            sep = "\t",
            row.names = FALSE)




