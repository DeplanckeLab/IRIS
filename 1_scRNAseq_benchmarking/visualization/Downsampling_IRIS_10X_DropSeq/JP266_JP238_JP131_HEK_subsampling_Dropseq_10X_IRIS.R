# Author: Joern Pezoldt
# Date: 08.08.2023
# Description: Comparative QC of 10X data to IRIS biochem
##1) Read distribution
##2) Mapped reads vs. UMI

#libraries
library(Seurat)
library(tidyr)
library(rtracklayer)


#####
#Directories and global variables
#####
#PATH
expID_IRIS <- c("JP266")
expID_10X <- c("JP131")
expID_Dropseq <- c("JP238")
genome_10X <- c("hg")
species_in_Exp <- c("human")
genome_IRIS <- c("human")

PATH_input_10X_LT <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP131/results/tables"
PATH_input_IRIS <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis"
PATH_input_Dropseq <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/results/tables"

PATH_output <- "/home/pezoldt/updepla/projects/0_Experiments/226_10X_comparison/results"
PATH_human_genes <- "/home/pezoldt/updepla/projects/iris/2_annotation_genomes/genomes/human/GRCh38.90/gtf/genesymbol_ensg_GRCh38.90.txt"
  
#Globals
min_UMI <- 3000
min_mapped_reads <- 100000
options(scipen = 999)

###########################################################
#################### IRIS #################################
###########################################################
#####
#Load Count matrix IRIS platform
#####
#LiquideDrop data
l_pivot_table <- list()
l_pivot_UMI_tresh <- list()
l_rcm_all_cells <- list()
l_rcm_mapped_threshed <- list()
l_rcm_UMI_threshed <- list()
l_UMI_all_cells_ENSG <- list()
l_UMI_all_cells_SYMB <- list()
l_unknown_barcodes <- list()
l_SIS_wellStats <- list()

for(i in 1:length(expID_IRIS)){
  #i = 1
  print(i)
  print(expID_IRIS[i])
  load(paste0(PATH_input_IRIS,"/",expID_IRIS[i],"/pipeline_output/human/analysis/objects/objects_plotting.RData"), envir = parent.frame(), verbose = FALSE)
  l_pivot_table[[i]] <- pivot_table
  l_UMI_all_cells_ENSG[[i]] <- UMI_all_cells_ENSG
  l_UMI_all_cells_SYMB[[i]] <- UMI_all_cells_SYMB
  l_SIS_wellStats[[i]] <- readRDS(paste0(PATH_input_IRIS,"/",expID_IRIS[i],"/pipeline_output/human/analysis/objects/SIS_wellStats.Rds"))
}
#names to lists
names(l_pivot_table) <- expID_IRIS
pivot_all <- do.call(rbind, l_pivot_table)

names(l_UMI_all_cells_ENSG) <- expID_IRIS
names(l_UMI_all_cells_SYMB) <- expID_IRIS
lapply(l_pivot_table, head)
lapply(l_UMI_all_cells_ENSG, ncol)
#QC on a glance
print("meadian(mappedRead)")
lapply(l_pivot_table, function(x){median(x$n_mapped_reads)})
print("meadian(nUMI)")
lapply(l_pivot_table, function(x){median(x$n_UMI)})
print("s(nUMI)")
lapply(l_pivot_table, function(x){sd(x$n_UMI)})
print("median(nGenes)")
lapply(l_pivot_table, function(x){median(x$n_genes)})
print("nCells")
unlist(lapply(l_UMI_all_cells_SYMB, ncol))

#################
## Read stats ###
#################
i = 1
#summarize per plate
t_SIS_wellStats <- l_SIS_wellStats[[i]]
t_SIS_wellStats <- t_SIS_wellStats[complete.cases(t_SIS_wellStats[,c("n_discarded_at_CC","n_Unmapped","n_noFeature","n_ambigFeature","n_UMI","n_saturated_UMI")]),]
nrow(t_SIS_wellStats)
t_SIS_wellStats_sum <- aggregate(t_SIS_wellStats[, c("n_discarded_at_CC","n_Unmapped","n_noFeature","n_ambigFeature","n_UMI","n_saturated_UMI")], list(t_SIS_wellStats$Description), sum)
colnames(t_SIS_wellStats_sum) <- c("Sample_ID","n_discarded_at_CC","n_Unmapped","n_noFeature","n_ambigFeature","n_UMI","n_saturated_UMI")

#################
## UMI vs Reads #
#################
expID_i <- expID_IRIS[[i]]
t_IRIS_pivot <- l_pivot_table[[i]]
t_IRIS_UMI_reads <- t_IRIS_pivot[,c("Full_Cell_ID","Sample_ID","n_UMI","n_genes","n_all_reads","n_mapped_reads")]
colnames(t_IRIS_UMI_reads) <- c("Cell_ID","orig.ident","nCount_RNA","nFeature_RNA","n_all_reads","n_mapped_reads")
t_IRIS_UMI_reads$orig.ident <- rep(expID_i, nrow(t_IRIS_UMI_reads))

#################
## Subsampling ##
#################
# 10X LT HEK
pivot_subsamp_10X_LT <- read.delim(paste0(PATH_input_10X_LT,"/JP131_subsampling_100000.txt"))
# Dropseq 20min HEK
pivot_subsamp_Dropseq <- read.delim(paste0(PATH_input_Dropseq,"/JP238_subsampling.txt"))

# IRIS----------------
load(paste0(PATH_input_IRIS,"/",expID_IRIS,"/pipeline_output/human/analysis/objects/subsampling_objects_plotting.RData"), envir = parent.frame(), verbose = FALSE)
t_IRIS_subsamp_gene <- gene_sum_table
t_IRIS_subsamp_UMI <- UMI_sum_table
# pivot data
## gene
pivot_t_IRIS_subsamp_gene <- t_IRIS_subsamp_gene %>%
  pivot_longer(
    cols = starts_with("CC"),
    names_to = "CC",
    names_prefix = "",
    values_to = "nGene",
    values_drop_na = TRUE)
pivot_t_IRIS_subsamp_gene$Full_Cell_ID <- paste0(pivot_t_IRIS_subsamp_gene$Sample_well_id,"_",pivot_t_IRIS_subsamp_gene$CC)
## nUMI
pivot_t_IRIS_subsamp_UMI <- t_IRIS_subsamp_UMI %>%
  pivot_longer(
    cols = starts_with("CC"),
    names_to = "CC",
    names_prefix = "",
    values_to = "nUMI",
    values_drop_na = TRUE)
pivot_t_IRIS_subsamp_UMI$Full_Cell_ID <- paste0(pivot_t_IRIS_subsamp_UMI$Sample_well_id,"_",pivot_t_IRIS_subsamp_UMI$CC)
## Merge by Full_Cell_ID
pivot_IRIS_subsamp <- merge(pivot_t_IRIS_subsamp_gene, pivot_t_IRIS_subsamp_UMI,
                            by = c("Full_Cell_ID","subsample_size","CC","Sample_well_id","Description"), 
                            no.dups = FALSE)

#Combine subsampling and all information
pivot_all_condense <- pivot_all[,c("Full_Cell_ID",
                                   "n_all_reads","n_mapped_reads",
                                   "n_UMI","n_genes")]
colnames(pivot_all_condense) <- c("Full_Cell_ID",
                                  "n_all_reads","n_mapped_reads",
                                  "n_all_UMI","n_all_genes")
pivot_IRIS_subsamp_condense <- pivot_IRIS_subsamp[,c("Full_Cell_ID",
                                                     "subsample_size",
                                                     "nUMI","nGene")]
colnames(pivot_IRIS_subsamp_condense) <- c("Full_Cell_ID",
                                  "subsample_size",
                                  "n_subsamp_UMI","n_subsamp_genes")
#Merge
pivot_subsamp_IRIS <- merge(pivot_all_condense, pivot_IRIS_subsamp_condense, by = "Full_Cell_ID")
pivot_subsamp_IRIS$expID <- rep("JP266", nrow(pivot_subsamp_IRIS))

# combine tables
pivot_subsamp_all <- rbind(pivot_subsamp_10X_LT,
                           pivot_subsamp_Dropseq,
                           pivot_subsamp_IRIS)

#Force integer
levels(as.factor(pivot_subsamp_all$subsample_size))
as.numeric(pivot_subsamp_all$subsample_size)
levels(as.factor(as.numeric(pivot_subsamp_all$subsample_size)))
pivot_subsamp_all$subsample_size <- as.character(as.numeric(pivot_subsamp_all$subsample_size))
#eliminate duplicates
pivot_subsamp_all <- unique(pivot_subsamp_all)
#Adapt to minimum mapped reads chosen for cells sub-sampled for public 10X data
pivot_subsamp_all <- subset(pivot_subsamp_all, n_mapped_reads > 40000)
#Add Column for biochem annotation
expID <- pivot_subsamp_all$expID
Biochem <- gsub("JP131", "10X_V3", expID)
Biochem <- gsub("JP238", "Drop-seq", Biochem)
Biochem <- gsub("JP266", "IRIS", Biochem)
pivot_subsamp_all$Biochem <- Biochem

#####
#Distribution of cell sequencing depth
#####
#JP131 10X V3
FCI_JP131 <- subset(pivot_subsamp_all,  n_all_UMI > 1000 & expID == "JP131")$Full_Cell_ID
t_mapped_JP131 <- t_mappedReads <- read.table(paste0("/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP131/results/tables/JP131_mappedreads_per_barcode.txt"),
                                              header = FALSE, skip = 1)
colnames(t_mapped_JP131) <- c("n_mapped_reads", "Cell_ID")
t_mapped_JP131 <- subset(t_mapped_JP131, Cell_ID %in%  FCI_JP131)
hist(log10(t_mapped_JP131$n_mapped_reads))
t_mapped_JP131$n_mapped_reads_scaled <- scale(t_mapped_JP131$n_mapped_reads)
t_mapped_JP131$Biochem <- rep("10X_V3", nrow(t_mapped_JP131))
hist(t_mapped_JP131$n_mapped_reads_scaled)

#JP238 Drop-seq
FCI_JP238 <- subset(pivot_subsamp_all,  n_all_UMI > 1000 & expID == "JP238")$Full_Cell_ID
t_mapped_JP238 <- t_mappedReads <- read.table(paste0("/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/tables/JP238_20min_mappedreads_per_barcode.txt"),
                                              header = FALSE, skip = 1)
colnames(t_mapped_JP238) <- c("n_mapped_reads", "Cell_ID")
t_mapped_JP238 <- subset(t_mapped_JP238, Cell_ID %in%  FCI_JP238)
hist(log10(t_mapped_JP238$n_mapped_reads))
t_mapped_JP238$n_mapped_reads_scaled <- scale(t_mapped_JP238$n_mapped_reads)
t_mapped_JP238$Biochem <- rep("Drop-seq", nrow(t_mapped_JP238))
hist(t_mapped_JP238$n_mapped_reads_scaled)

#JP266 IRIS
t_mapped_JP266 <- subset(pivot_subsamp_all,  n_all_UMI > 1000 & expID == "JP266")
t_mapped_JP266 <- t_mapped_JP266[,c("n_mapped_reads", "Full_Cell_ID")]
colnames(t_mapped_JP266) <- c("n_mapped_reads", "Cell_ID")
hist(log10(t_mapped_JP266$n_mapped_reads))
t_mapped_JP266$n_mapped_reads_scaled <- scale(t_mapped_JP266$n_mapped_reads)
t_mapped_JP266$Biochem <- rep("IRIS", nrow(t_mapped_JP266))
hist(t_mapped_JP266$n_mapped_reads_scaled)

#Combine data
t_mapped_scaled <- rbind(t_mapped_JP131,t_mapped_JP238,t_mapped_JP266)

#Histogram overlay
# Change histogram plot fill colors by groups
ggplot(t_mapped_scaled, aes(x=n_mapped_reads_scaled, fill=Biochem, color=Biochem)) +
  geom_density(alpha=0.5) +
  xlim(-3.5,6) +
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  scale_fill_manual(values=c("deepskyblue", "lightblue","blue")) +
  geom_histogram(position="identity",
              alpha=0.2,
             aes(y=..density..),
            binwidth = 0.15) +
  labs(x = "mean-centered(mapped reads)",
       y = "proportion cells")

ggplot(t_mapped_scaled, aes(x=n_mapped_reads, fill=Biochem, color=Biochem)) +
  geom_density(alpha=0.5) +
  #xlim(-3.5,6) +
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  scale_fill_manual(values=c("deepskyblue", "lightblue","blue")) +
  geom_histogram(position="identity",
                 alpha=0.2,
                 aes(y=..density..)) +
  scale_x_continuous(trans = "log10") +
  labs(x = "mean-centered(mapped reads)",
       y = "proportion cells")

#10X_V3
t_mapped_scaled_10X <- subset(t_mapped_scaled, Biochem == "10X_V3")
ggplot(t_mapped_scaled_10X, aes(x=n_mapped_reads_scaled, fill=Biochem, color=Biochem)) +
  geom_density(alpha=0.5) +
  xlim(-3.5,6) +
  scale_color_manual(values=c("deepskyblue")) +
  scale_fill_manual(values=c("deepskyblue")) +
  geom_histogram(position="identity",
                 alpha=0.2,
                 aes(y=..density..),
                 binwidth = 0.3) +
  labs(x = "mean-centered(mapped reads)",
       y = "proportion cells") +
  theme_classic()

#Drop-seq
t_mapped_scaled_Dropseq <- subset(t_mapped_scaled, Biochem == "Drop-seq")
ggplot(t_mapped_scaled_Dropseq, aes(x=n_mapped_reads_scaled, fill=Biochem, color=Biochem)) +
  geom_density(alpha=0.5) +
  xlim(-3.5,6) +
  scale_color_manual(values=c("lightblue")) +
  scale_fill_manual(values=c("lightblue")) +
  geom_histogram(position="identity",
                 alpha=0.2,
                 aes(y=..density..),
                 binwidth = 0.3) +
  labs(x = "mean-centered(mapped reads)",
       y = "proportion cells") +
  theme_classic()

#IRIS
t_mapped_scaled_IRIS <- subset(t_mapped_scaled, Biochem == "IRIS")
ggplot(t_mapped_scaled_IRIS, aes(x=n_mapped_reads_scaled, fill=Biochem, color=Biochem)) +
  geom_density(alpha=0.5) +
  xlim(-3.5,6) +
  scale_color_manual(values=c("blue")) +
  scale_fill_manual(values=c("blue")) +
  geom_histogram(position="identity",
                 alpha=0.2,
                 aes(y=..density..),
                 binwidth = 0.3) +
  labs(x = "mean-centered(mapped reads)",
       y = "proportion cells") +
  theme_classic()

##################################
# Plotting
##################################
#####
#full range 200000
#####
level_order <- c("10000","20000","50000","100000")
pivot_subsamp_all_select <- subset(pivot_subsamp_all, subsample_size %in% level_order)
pivot_subsamp_all_select_JP266 <- subset(pivot_subsamp_all_select, expID == "JP266")
hist(log10(pivot_subsamp_all_select_JP266$n_all_genes))
pivot_subsamp_all_select_JP131 <- subset(pivot_subsamp_all_select, expID == "JP131")
hist(log10(pivot_subsamp_all_select_JP131$n_all_genes))

#Boxplot
## Genes
ggplot(pivot_subsamp_all_select, aes(x=factor(subsample_size, level = level_order), 
                                     y=n_subsamp_genes ,
                                     fill=Biochem)) +
  geom_boxplot(outliers = FALSE,
               notch = TRUE,
               notchwidth = 0.2,
               staplewidth = 0.3,
               varwidth = FALSE) +
  labs(x = "Subsample MappedReads",
       y = "n_genes",
       colour = "Biochemistry") +
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  scale_fill_manual(values=c("deepskyblue", "lightblue","blue")) +
  theme_classic()
#Genes comparison mean
pivot_100000 <- subset(pivot_subsamp_all_select, subsample_size == 100000)
nrow(pivot_100000)
genes_Dropseq <- mean(subset(pivot_100000, Biochem == "Drop-seq")$n_subsamp_genes)
nrow(subset(pivot_100000, Biochem == "Drop-seq"))
genes_10X <- mean(subset(pivot_100000, Biochem == "10X_V3")$n_subsamp_genes)
nrow(subset(pivot_100000, Biochem == "10X_V3"))
genes_IRIS <- mean(subset(pivot_100000, Biochem == "IRIS")$n_subsamp_genes)
nrow(subset(pivot_100000, Biochem == "IRIS"))
genes_IRIS/genes_Dropseq
genes_IRIS/genes_10X *100

## UMI
ggplot(pivot_subsamp_all_select, aes(x=factor(subsample_size, level = level_order), 
                                     y=n_subsamp_UMI,
                                     fill=Biochem)) +
  geom_boxplot(outliers = FALSE,
               notch = TRUE,
               notchwidth = 0.2,
               staplewidth = 0.3,
               varwidth = FALSE) +
  labs(x = "Subsample MappedReads",
       y = "n_UMI",
       colour = "Biochemistry") +
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  scale_fill_manual(values=c("deepskyblue", "lightblue","blue")) +
  theme_classic()

#####
#Plot UMI vs mapped reads
#####
pivot_subsamp_all_select_nondup <- subset(pivot_subsamp_all_select, subsample_size == 100000)
#n Genes
ggplot(pivot_subsamp_all_select_nondup, aes(x=n_mapped_reads, y=n_all_genes, color=Biochem)) + #n_mapped_reads
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  geom_point(size = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  labs(title="mapped_reads vs. UMI") + 
  xlim(c(0,500000)) +
  #ylim(c(0,200000)) +
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = FALSE)
#n UMI
ggplot(pivot_subsamp_all_select, aes(x=n_mapped_reads, y=n_all_UMI, color=Biochem)) + #n_mapped_reads
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  geom_point(size = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  labs(title="mapped_reads vs. UMI") + 
  xlim(c(0,500000)) +
  #ylim(c(0,200000)) +
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = FALSE)

#n Genes
ggplot(pivot_subsamp_all_select, aes(x=n_all_UMI, y=n_all_genes, color=Biochem)) + #n_mapped_reads
  scale_color_manual(values=c("deepskyblue", "lightblue","blue")) +
  geom_point(size = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted')) +
  labs(title="mapped_reads vs. UMI") + 
  xlim(c(0,100000)) +
  #ylim(c(0,200000)) +
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = FALSE)

