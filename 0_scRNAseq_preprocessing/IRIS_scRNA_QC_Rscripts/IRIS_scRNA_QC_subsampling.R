# Author: Caro
# Date: 17.08.2021
# Description: This script extracts quality control statitistics from the subsampled reads produced by the IRIS_scRNAseq_QC_analysis pipeline


# Print status of script 
cat("Succesfully entered subsampling script\n")

#############################################
########## Load required packages ###########
#############################################
  
## Load packages and supress loading messages to be printed to console
suppressPackageStartupMessages({
library(dplyr)
library(tidyverse)
library(Matrix)
library(data.table)
library(stringr)
library(readr)
library(tibble)
})

# Print status of script 
cat("Loaded packages\n")


#############################################
## Assign variables coming from config.txt ##
#############################################

# Print status of script
cat("Entered subsampling Rscript\n")

# Read variables from command line (T = drop internal arguments)
args <- commandArgs(T)

work_dir <- args[1]
#work_dir <- "/home/pezoldt/NAS2/iris/experiments/pezoldt_JP188_iris12_manual"

expID <- args[2]
#expID <- "JP188"

barcodes_brb_file <- args[3]
#barcodes_brb_file <- "/home/pezoldt/NAS2/iris/genomes/barcodes/iris_genomics/20210618_scRNAseq_cellcode_BRBtoolbox_CC_7bp.txt"

sample_info_sheet <- args[4]
#sample_info_sheet <- "/home/pezoldt/NAS2/iris/sample_information_sheets/sis_iris12__mapping_QC_pipeline.csv"

# Print status of script 
cat("Read arguments from command line:\n", args)

#############################################
######### Create file directories  ##########
#############################################

# Set fixed parts of directories 
data <- "data"
analysis <- "analysis"
rcm <- "count_matrices"
plots <- "plots"
objects <- "objects"

# Set file paths 
rcm_path <- file.path(work_dir, data, rcm)
plots_path <- file.path(work_dir, analysis, plots)
objects_path <- file.path(work_dir, analysis, objects)

# Print status of script
cat("Set file names\n")

#############################################
############## Load files  ##################
#############################################

# Load brb barcodes file
barcodes_brb <- read.delim(barcodes_brb_file)

# Load sample information sheet 
SIS <- read.csv(sample_info_sheet, skip = 1, header = T)

# Load pivot table of experiment
pivot_table <- readRDS(file.path(objects_path, "pivot_table.Rds"))

# Print status of script
cat("Loaded files\n")

#############################################
############## Variables ####################
#############################################

## Define variables needed across the script

# Set subsample sizes 
subsample_sizes <- c("200000", "100000", "50000", "20000", "15000", "10000", "5000")

subsample_folder <- paste0("subsample_", subsample_sizes, "_reads")

# Extract information of specific experiment from SIS
SIS_exp <- SIS[which(str_detect(SIS$Sample_ID, expID)),]

# Get barcodes used in experiment 
barcodes_used <- str_split(SIS_exp$CC_used[1], "_") %>% unlist()

# Split SIS by Sample ID (i.e. by well)
SIS_list <- split(SIS_exp, SIS_exp$Sample_ID)

# Create SIS with one row per cell
SIS_list <- lapply(SIS_list, function(x){
  
  # Determine used barcodes
  barcodes_used <- str_split(x$"CC_used", "_") %>% unlist()
  
  # Repeat each well information by the amount of cells per well 
  n_cells <- length(barcodes_used)
  well <- x[rep(seq_len(nrow(x)), each = n_cells), ]
  
  # Add Lid_ID and Full_Cell_ID (= well ID + Lid ID) as columns to SIS 
  well$CC <- barcodes_used
  well <- well %>% unite("Full_Cell_ID", c(Sample_ID, CC), remove = FALSE)
  
  return(well)
})

# Combine information per cell to one data frame 
SIS_exp_perCell <- do.call(rbind, SIS_list)
rownames(SIS_exp_perCell) <- NULL


# Load brb barcodes 
barcodes_brb <- read.delim(barcodes_brb_file)


#############################################
############## PreVariables #################
#############################################

## Define variables needed across the script

# Get names of well folders from rcm directory
well_folders <- list.dirs(rcm_path, recursive = F)

# Get sample/well names 
well_names <- gsub(".*/", "", well_folders)

# Print status of script
cat("listed well names\n")

#############################################
############## STAR solo data present########
#############################################
## Check if file "Log.final.out" exists 
well_mapped <- lapply(well_folders, function(x){
  log_file <- file.path(x, "Log.final.out")
  return(file.exists(log_file))
})

#Eliminate well-folders without "Log.final.out" file
well_folders <- well_folders[unlist(well_mapped)]
well_names <- well_names[unlist(well_mapped)]
length(well_names)
#Adapt SIS
SIS_exp <- subset(SIS_exp, Sample_ID %in% well_names)

## Check if there are more than 100 UMIs detected in total
UMI_detected <- lapply(well_folders, function(x){
  #x <- well_folders[[1]]
  UMI_count <- read.table(paste0(x,"/Solo.out/Gene/UMIperCellSorted.txt"))
  colSums(UMI_count) > 100
})
UMI_detected <- unlist(UMI_detected)
names(UMI_detected) <- NULL
#Eliminate well-folders without sufficient UMI counts per well
well_folders <- well_folders[UMI_detected]
well_names <- well_names[UMI_detected]
#Adapt SIS
SIS_exp <- subset(SIS_exp, Sample_ID %in% well_names)

#############################################
####### Re-define Variables #################
#############################################

## Re-Define variables needed across the script

# Get names of well folders from rcm directory
well_folders <- paste0(rcm_path,"/",SIS_exp$Sample_ID)

# Get sample/well names 
well_names <- gsub(".*/", "", well_folders)

#############################################
####### Read in STARsolo UMI matrices #######
#############################################

# Get file names of UMI matrices per subsample size 
subsample_files <- lapply(subsample_folder, function(x) {
  
  lapply(well_folders, function(y){
    file.path(y, x, "Solo.out/Gene/filtered")
    })
})

names(subsample_files) <- subsample_folder

# Extract UMI matrix per well and per subsample size
UMI_per_subsample_size <- lapply(seq_along(subsample_files), function(x){
  
  # Determine which barcodes were used 
  barcodes_used <- str_split(SIS_exp[x, "CC_used"], "_") %>% unlist()
  
  lapply(subsample_files[[x]], function(y){
    
    # Load UMI data per well i.e. barcodes, features and expression matrix
    barcodes <- read.delim(file.path(y, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
    features <- read.delim(file.path(y, "features.tsv"), header = FALSE, stringsAsFactors = FALSE)
    matrix <- readMM(file.path(y, "matrix.mtx"))
    
    # Set gene symbols as rownames of expression matrix
    rownames(matrix) = features$V2
    
    # Match and replace the barcode sequences with the Lid ID 
    barcodes_detected <- merge(barcodes_brb, barcodes, by.x = "B1", by.y = "V1")
    colnames(matrix) <- barcodes_detected$Name 
    
    # Only keep the used barcodes
    matrix <- as.matrix(matrix)
    matrix <- matrix[,(colnames(matrix) %in% barcodes_used)]
    
    return(matrix)
    
  })
})

names(UMI_per_subsample_size) <- subsample_folder

# Set names for each well in the subsample lists
UMI_per_subsample_size <- lapply(UMI_per_subsample_size, function(x){
  
  names(x) <- well_names
  return(x)
  
})

# Print status of script
cat("extracted UMI per subsample size\n")


#############################################
###### Calculate UMI per subsample size #####
#############################################

#Sum UMIs per CC for each subsample for each condition
l_well <- list()
l_subsample <- list()
for(i in 1:length(UMI_per_subsample_size)){
  #i = 1
  subsample_i <- UMI_per_subsample_size[[i]]
  for(j in 1:length(subsample_i)){
    #j = 1
    well_j <- subsample_i[[j]]
    well_name_j <- subsample_i[j]
    well_UMI_j <- colSums(well_j)
    well_UMI_j <- well_UMI_j[order(names(well_UMI_j))]
    l_well[[j]] <- well_UMI_j
  }
  names(l_well) <- names(subsample_i)
  l_subsample[[i]] <- l_well
  l_well <- list()
}
names(l_subsample) <- names(UMI_per_subsample_size)
UMI_per_CC <- l_subsample

# Bind sums of UMI per subsample size in one dataframe 
UMI <- lapply(UMI_per_CC, function(x){
  
  per_well_sub <- do.call(rbind, x)
  
})
lapply(UMI, nrow)

## Split data by conditions defined in pivot table 

# Split pivot table by conditions 
conditions_list <- split(pivot_table, pivot_table$Description) 

# Extract sample names related to conditions
conditions_list <- lapply(conditions_list, function(x){x$Sample_ID})

# Set condition description as names
condition_names <- names(conditions_list)

# Get UMI sum tables per condition 
UMI_per_condition <- lapply(seq_along(condition_names), function(x){
  
  filtered_UMI_subsamp <- lapply(UMI, function(y){
    
    table <- as.data.frame(y)
    table <- table[rownames(table) %in% conditions_list[[x]], ]
    return(table)
  })
  return(filtered_UMI_subsamp)
})
names(UMI_per_condition) <- condition_names

# Print status of script
cat("Calculated UMI per subsample size\n")

# Make table of UMIs with all subsamples and respective CCs
l_condition_tables <- list()
for(i in 1:length(UMI_per_condition)){
  #i = 1
  l_condition <- UMI_per_condition[[i]]
  condition_i <- names(UMI_per_condition)[i]
  l_subsample <- list()
  for(j in 1:length(l_condition)){
    #j = 1
    sub_size_j <- unlist(strsplit(names(l_condition)[j], "_"))[2]
    subsample_j <- l_condition[[j]]
    subsample_j$subsample_size <- rep(sub_size_j, nrow(subsample_j))
    subsample_j$Sample_well_id <- rownames(subsample_j)
    l_subsample[[j]] <- subsample_j
  }
  t_condition <- do.call(rbind, l_subsample)
  t_condition$Description <- rep(condition_i, nrow(t_condition))
  #t_condition$Sample_well_id <- rownames(t_condition)
  l_condition_tables[[i]] <- t_condition
}
UMI_sum_table <- do.call(rbind, l_condition_tables)

# Print status of script
cat("Compiled UMI summary\n")

#############################################
##### Calculate genes per subsample size ####
#############################################

#Sum Genes per CC for each subsample for each condition
l_subsample <- list()
for(i in 1:length(UMI_per_subsample_size)){
  #i = 1
  subsample_i <- UMI_per_subsample_size[[i]]
  l_well <- list()
  for(j in 1:length(subsample_i)){
    #j = 1
    well_j <- subsample_i[[j]]
    well_name_j <- subsample_i[j]
    well_gene_j <- colSums(well_j[,colnames(well_j)] > 0)
    well_gene_j <- well_gene_j[order(names(well_gene_j))]
    l_well[[j]] <- well_gene_j
  }
  names(l_well) <- names(subsample_i)
  l_subsample[[i]] <- l_well
}
names(l_subsample) <- names(UMI_per_subsample_size)
gene_per_CC <- l_subsample

# Bind sums of UMI per subsample size in one dataframe 
GENE <- lapply(gene_per_CC, function(x){
  
  per_well_sub <- do.call(rbind, x)
  
})
lapply(GENE, nrow)

# Get gene sum tables per condition 
genes_per_condition <- lapply(seq_along(condition_names), function(x){
  
  filtered_gene_subsamp <- lapply(GENE, function(y){
    
    table <- as.data.frame(y)
    table <- table[rownames(table) %in% conditions_list[[x]], ]
    return(table)
  })
  return(filtered_gene_subsamp)
})
names(genes_per_condition) <- condition_names

# Print status of script
cat("Calculates genes per condition\n")

# Make table of UMIs with all subsamples and respective CCs
l_condition_tables <- list()
for(i in 1:length(genes_per_condition)){
  #i = 1
  l_condition <- genes_per_condition[[i]]
  condition_i <- names(genes_per_condition)[i]
  l_subsample <- list()
  for(j in 1:length(l_condition)){
    #j = 1
    sub_size_j <- unlist(strsplit(names(l_condition)[j], "_"))[2]
    subsample_j <- l_condition[[j]]
    subsample_j$subsample_size <- rep(sub_size_j, nrow(subsample_j))
    subsample_j$Sample_well_id <- rownames(subsample_j)
    l_subsample[[j]] <- subsample_j
  }
  t_condition <- do.call(rbind, l_subsample)
  t_condition$Description <- rep(condition_i, nrow(t_condition))
  l_condition_tables[[i]] <- t_condition
}
gene_sum_table <- do.call(rbind, l_condition_tables)

# Print status of script
cat("Calculated gene sums\n")

#######################################
#### Save data needed for plotting ####
#######################################

# Set file name 
saved_objects_path <- file.path(objects_path, "subsampling_objects_plotting.RData")

# Save UMI and gene table as one object
save(UMI_sum_table, gene_sum_table, file = saved_objects_path)

# Print status of script
cat("Saved objects for plotting\n")
