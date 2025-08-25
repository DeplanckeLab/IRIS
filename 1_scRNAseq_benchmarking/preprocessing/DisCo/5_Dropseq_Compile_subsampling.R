#Author: Joern Pezoldt

#Purpose: Compile Subsampling

#####
#Libraries
#####
library(plyr)

#####
#Paths & Globals
#####
expID <- "JP238"
PATH_all <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/results"
PATH_10X <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/objects/Solo.out/Gene/filtered"
PATH_subsampling <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/objects/subsampling/merged"

PATH_output_tables <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/results/tables"

#Depth
subsampling_depth <- c(5000,10000,15000,20000,50000,100000,200000)
subsampling_depth_character <- c("5000","10000","15000","20000","50000","100000","200000")
subsampling_folders <- paste0("merged_",subsampling_depth_character)

options(scipen = n)

#####
#Count matrix all data
#####
#From all_reads data
## All reads
t_allReads <- read.table(paste0("/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/tables/JP238_20min_reads_per_barcode.txt"),
                         header = FALSE, skip = 1)
colnames(t_allReads) <- c("n_all_reads", "Cell_ID")
## Mapped Reads
t_mappedReads <- read.table(paste0("/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP238/tables/JP238_20min_mappedreads_per_barcode.txt"),
                           header = FALSE, skip = 1)
colnames(t_mappedReads) <- c("n_mapped_reads", "Cell_ID")
## Count matrix
t_counts <- as.data.frame(Read10X(PATH_10X))
t_UMI <- data.frame(Cell_ID = colnames(t_counts),
                    n_UMI = colSums(t_counts))
nonzero <- function(x) sum(x != 0)
v_Genes <- numcolwise(nonzero)(t_counts)
rownames(v_Genes) <- NULL
t_Genes <- data.frame(Cell_ID = colnames(v_Genes),
                      n_Genes = as.numeric(v_Genes[1,]))
# Merge Intel
t_all <- merge(t_allReads, t_mappedReads, by = "Cell_ID")
t_all <- merge(t_all, t_Genes, by = "Cell_ID",
               all.y = TRUE)
t_all <- merge(t_all, t_UMI, by = "Cell_ID",
               all.y = TRUE)
colnames(t_all) <- c("Full_Cell_ID","n_all_reads","n_mapped_reads","n_all_genes","n_all_UMI")
t_all <- t_all[,c("Full_Cell_ID","n_all_reads","n_mapped_reads","n_all_UMI","n_all_genes")]

#####
#Count matrices subsampling
#####
#Loop over subsampling
#Select barcodes with sufficient sequencing depth
#Compute genes and UMI per barcode
l_subsamp <- list()
for(i in 1:length(subsampling_folders)){
  #i = 1
  print(i)
  subsampling_folders_i <- subsampling_folders[i]
  PATH_10_subsamp_i <- paste0(PATH_subsampling,"/",
                              subsampling_folders_i,
                              "/Solo.out/Gene/filtered")
  ## Count matrix
  t_counts_i <- as.data.frame(Read10X(PATH_10_subsamp_i))
  t_UMI_i <- data.frame(Cell_ID = colnames(t_counts_i),
                      n_UMI = colSums(t_counts_i))
  t_UMI_i$subsample_size <- rep(subsampling_depth[i], nrow(t_UMI_i))
  nonzero <- function(x) sum(x != 0)
  v_Genes_i <- numcolwise(nonzero)(t_counts_i)
  rownames(v_Genes_i) <- NULL
  t_Genes_i <- data.frame(Cell_ID = colnames(v_Genes_i),
                        n_Genes = as.numeric(v_Genes_i[1,]))
  t_subsamp_i <- merge(t_UMI_i, t_Genes_i)
  colnames(t_subsamp_i) <- c("Full_Cell_ID","n_subsamp_UMI","subsample_size", "n_subsamp_genes")
  t_subsamp_i <- t_subsamp_i[,c("Full_Cell_ID","subsample_size","n_subsamp_UMI","n_subsamp_genes")]
  l_subsamp[[i]] <- t_subsamp_i
}
lapply(l_subsamp, head)
t_subsamp_all <- do.call(rbind, l_subsamp)

#Merge
pivot_subsamp_IRIS <- merge(t_all, t_subsamp_all, by = "Full_Cell_ID")
pivot_subsamp_IRIS$expID <- rep(expID, nrow(pivot_subsamp_IRIS))

#save
write.table(pivot_subsamp_IRIS, paste0(PATH_output_tables,"/",expID,"_subsampling.txt"), sep="\t")

