#Author: Joern Pezoldt
#Purpose: Annotate mouse and human cells in mixed datasets

##############
#Figure 2
##############


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
#Mote: Only a vector with two elements permitted, second item is spike-in
species_in_Exp <- c("human","mouse")
min_UMI <- 5000
#Percent of species mixing allowed to be considered a "clean" cell
thresh_sm <- 10

## PATHs
#General experimental ID used throughout analysis
expID <- "JP298"
userID <- "pezoldt"
#FolderName in experimental folder
simple_ID <- "JP298"

#Output name of scRNA-seq pipeline
#Automatic paths
PATH_input <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis/",expID,"/pipeline_output/",mapped_genome)
PATH_output <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis/",expID,"/results")
PATH_output_figures <- paste0(PATH_output,"/plots")
PATH_output_objects <- paste0(PATH_output,"/objects_transcriptome")
PATH_output_tables <- paste0(PATH_output,"/tables")

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
#Percentage mixing
#####
pivot_UMI_tresh_double <- subset(pivot_UMI_tresh_double, Description == "JP298_P2")

percent_human <- round(pivot_UMI_tresh_double$nUMI_human / pivot_UMI_tresh_double$n_UMI * 100, 1)
percent_mouse <- round(pivot_UMI_tresh_double$nUMI_mouse / pivot_UMI_tresh_double$n_UMI * 100, 1)
t_percent <- data.frame(percent_human = percent_human, percent_mouse = percent_mouse)

pivot_UMI_tresh_double$percent_UMI_human <- percent_human
pivot_UMI_tresh_double$percent_UMI_mouse <- percent_mouse

sm_percent_human <- round(mean(subset(pivot_UMI_tresh_double, sm_annotation == "human")$percent_UMI_human), 2)
sm_percent_mouse <- round(mean(subset(pivot_UMI_tresh_double, sm_annotation == "mouse")$percent_UMI_mouse), 2)
sm_percent_sm <- mean(subset(pivot_UMI_tresh_double, sm_annotation == "sm")$percent_UMI_mouse)
sm_percent_human
sm_percent_mouse
sm_percent_sm
nrow(subset(pivot_UMI_tresh_double, sm_annotation == "human"))
nrow(subset(pivot_UMI_tresh_double, sm_annotation == "mouse"))


#####
#Scatter Transcriptome
#####
#For experiment
theme_set(theme_bw())

ggplot(pivot_UMI_tresh_double,
       aes_string(x = paste0("nUMI_",species_in_Exp[1]),
                  y = paste0("nUMI_",species_in_Exp[2]))) + 
  geom_point(aes(col = sm_annotation),size=2) +   # draw points
  #xlim(c(0, 20000)) + 
  #ylim(c(0, 180000)) +  
  scale_color_manual(values = c("brown2","darkolivegreen3","grey40")) +
  labs(x = paste0("nUMI_",species_in_Exp[1]),
       y = paste0("nUMI_",species_in_Exp[2]),
       title = paste0(simple_ID)) +
  xlim(c(0, max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[1])]))) + 
  ylim(c(0, max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[2])]))) +
  annotate(geom="text",
           size = 6,
           x = 10000,
           y = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[2])]) - 1500, 
           color="green",
           label=(round(sum(pivot_UMI_tresh_double$sm_annotation == species_in_Exp[2]) / nrow(pivot_UMI_tresh_double) *100,1))) +
  annotate(geom="text",
           size = 6,
           x = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[1])]) - 5000,
           y = 8000, 
           color="deepskyblue3",
           label=(round(sum(pivot_UMI_tresh_double$sm_annotation == species_in_Exp[1]) / nrow(pivot_UMI_tresh_double) *100,1))) +
  annotate(geom="text",
           size = 6,
           x = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[1])]) - 1000,
           y = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[2])]) - 1000,
           color="grey40",
           label=(round(sum(pivot_UMI_tresh_double$sm_annotation == "sm") / nrow(pivot_UMI_tresh_double) *100,1))) +
  theme(legend.position = "none")

