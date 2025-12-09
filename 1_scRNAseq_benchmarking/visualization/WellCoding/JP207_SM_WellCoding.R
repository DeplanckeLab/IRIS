
#Author: Joern Pezoldt
#Date: 14.03.2023
#Function: Annotate mouse and human cells in mixed datasets



#libararies
library(dplyr)
library(patchwork)

#####
#Global Variables
#####
## Cutoffs
#Mote: Only a vector with two elements permitted, second item is spike-in
species_in_Exp <- c("human","mouse")
min_UMI <- 3000
#Percent of species mixing allowed to be considered a "clean" cell
thresh_sm <- 10

## PATHs
#General experimental ID used throughout analysis
expID <- "JP207"
#FolderName in experimental folder
simple_ID <- "JP207"

#Output name of scRNA-seq pipeline
#Automatic paths
#Note: Get count matrix and QC files from GEO and unzip
#Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9072151
PATH_countmatrix <- "/YOURPATH/GSM9072151_JP207_human_mouse_UMI_all_cells_ENSG.txt"
PATH_pivot <- "/YOURPATH/GSM9072151_JP207_human_mouse_QC.txt"
#Output
PATH_output <- paste0("/YOURPATH")
PATH_output_figures <- paste0(PATH_output,"/plots")
PATH_output_objects <- paste0(PATH_output,"/objects_transcriptome")
PATH_output_tables <- paste0(PATH_output,"/tables")

#####
#Load sequencing data
#####
UMI_all_cells_SYMB <- read.delim(PATH_countmatrix)
pivot_table <- read.delim(PATH_pivot)
ucm_lid <- UMI_all_cells_SYMB
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
genes_droso <- ucm_lid[grepl("^FBgn",rownames(ucm_lid)),]
nrow(genes_droso)
#Count UMIs per Species
nUMI_mouse <- colSums(genes_mouse)
nUMI_human <- colSums(genes_human)
nUMI_droso <- colSums(genes_droso)
#plot(nUMI_mouse, nUMI_human)
#plot(nUMI_droso, nUMI_human)
#Quick QC
#plot(log10(nUMI_droso), log10(nUMI_human))

#Count nGene per Species
nGene_mouse <- colSums(genes_mouse > 0)
nGene_human <- colSums(genes_human > 0)
nGene_droso <- colSums(genes_droso > 0)
#Quick QC
#plot(nGene_droso, nGene_human)

#build table with UMI and Full_Cell_ID data
to_pivot_ID_UMI_SM <- data.frame(Full_Cell_ID = names(nUMI_mouse),
                                 nUMI_mouse = nUMI_mouse, 
                                 nUMI_human = nUMI_human, 
                                 nUMI_droso = nUMI_droso,
                                 nGene_mouse = nGene_mouse, 
                                 nGene_human = nGene_human,
                                 nGene_droso = nGene_droso)
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
       title = paste0(simple_ID," | cutoff=",thresh_sm, "%")) +
  xlim(c(0, max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[1])]))) + 
  ylim(c(0, max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[2])]))) +
  annotate(geom="text",
           size = 6,
           x = 2500,
           y = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[2])]) - 1500, 
           color="brown2",
           label=(round(sum(pivot_UMI_tresh_double$sm_annotation == species_in_Exp[2]) / nrow(pivot_UMI_tresh_double) *100,1))) +
  annotate(geom="text",
           size = 6,
           x = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[1])]) - 1500,
           y = 2500, 
           color="darkolivegreen3",
           label=(round(sum(pivot_UMI_tresh_double$sm_annotation == species_in_Exp[1]) / nrow(pivot_UMI_tresh_double) *100,1))) +
  annotate(geom="text",
           size = 6,
           x = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[1])]) - 1500,
           y = max(pivot_UMI_tresh_double[,paste0("nUMI_",species_in_Exp[2])]) - 1500,
           color="grey40",
           label=(round(sum(pivot_UMI_tresh_double$sm_annotation == "sm") / nrow(pivot_UMI_tresh_double) *100,1))) +
  theme(legend.position = "none")


#Per plate
l_plot_barnyard <- list()
conditions <- levels(as.factor(pivot_UMI_tresh_double$Description))
for(i in 1:length(conditions)){
  #i = 1
  condition_i <-conditions[i]
  print(i)
  print(condition_i)
  One_Plate <- subset(pivot_UMI_tresh_double, Description %in% condition_i)
  plot_i <- ggplot(One_Plate,
                   aes_string(x = paste0("nUMI_",species_in_Exp[1]),
                              y = paste0("nUMI_",species_in_Exp[2]))) + 
    geom_point(aes(col = sm_annotation),size=2) +   # draw points
    #xlim(c(0, 20000)) + 
    #ylim(c(0, 180000)) +  
    scale_color_manual(values = c("deepskyblue3","green","grey40")) +
    labs(x = paste0("nUMI_",species_in_Exp[1]),
         y = paste0("nUMI_",species_in_Exp[2]),
         title = paste0(condition_i," | cutoff=",thresh_sm, "%")) +
    xlim(c(0, max(One_Plate[,paste0("nUMI_",species_in_Exp[1])]))) + 
    ylim(c(0, max(One_Plate[,paste0("nUMI_",species_in_Exp[2])]))) +
    annotate(geom="text",
             size = 6,
             x = 2500,
             y = max(One_Plate[,paste0("nUMI_",species_in_Exp[2])]) - 1500, 
             color="green",
             label=(round(sum(One_Plate$sm_annotation == species_in_Exp[2]) / nrow(One_Plate) *100,1))) +
    annotate(geom="text",
             size = 6,
             x = max(One_Plate[,paste0("nUMI_",species_in_Exp[1])]) - 1500,
             y = 2500, 
             color="deepskyblue3",
             label=(round(sum(One_Plate$sm_annotation == species_in_Exp[1]) / nrow(One_Plate) *100,1))) +
    annotate(geom="text",
             size = 6,
             x = max(One_Plate[,paste0("nUMI_",species_in_Exp[1])]) - 1500,
             y = max(One_Plate[,paste0("nUMI_",species_in_Exp[2])]) - 1500,
             color="grey40",
             label=(round(sum(One_Plate$sm_annotation == "sm") / nrow(One_Plate) *100,1))) +
    theme(legend.position = "none")
  l_plot_barnyard[[i]] <- plot_i
}
grid.arrange(grobs = l_plot_barnyard,
             ncol = 2,
             nrow = 2)

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

#per plate
plates <- levels(as.factor(pivot_UMI_tresh_double$Description))
for(i in 1:length(plates)){
  #i = 1
  plate_i <- plates[i]
  print(paste0("PlateID = ", plate_i))
  pivot_i <- subset(pivot_UMI_tresh_double, Description == plate_i)
  #Compute
  percent_human <- round(pivot_i$nUMI_human / pivot_i$n_UMI * 100, 1)
  percent_mouse <- round(pivot_i$nUMI_mouse / pivot_i$n_UMI * 100, 1)
  t_percent <- data.frame(percent_human = percent_human, percent_mouse = percent_mouse)
  
  pivot_i$percent_UMI_human <- percent_human
  pivot_i$percent_UMI_mouse <- percent_mouse
  
  sm_percent_human <- mean(subset(pivot_i, sm_annotation == "human")$percent_UMI_human)
  sm_percent_mouse <- mean(subset(pivot_i, sm_annotation == "mouse")$percent_UMI_mouse)
  sm_percent_sm <- mean(subset(pivot_i, sm_annotation == "sm")$percent_UMI_mouse)
  print(paste0("% human in human = ",sm_percent_human))
  print(paste0("% mouse in mouse = ",sm_percent_mouse))
  print(paste0("% sm in sm = ",sm_percent_sm))
}

