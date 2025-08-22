# Author: Caro
# Date: 17.08.2021
# Description: This script extracts quality control statistics from STARsolo output produced by the IRIS_scRNAseq_QC_analysis pipeline


# Print status of script 
cat("Succesfully entered RCMtoDEG script\n")

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

# Read variables from command line (T = drop internal arguments)
args <- commandArgs(T)

# Assign args 
work_dir <- args[1]
#work_dir <- "/home/pezoldt/NAS2/iris/experiments/ngrennin_JP280_aris03_human_mouse_droso"

expID <- args[2]
#expID <- "JP280"

ERCC <- args[3]
#ERCC <- "no"

barcodes_wl <- args[4]
#barcodes_wl <- "/home/pezoldt/NAS2/iris/4_genomes/barcodes/iris_genomics/20210618_scRNAseq_cellcode_whitelist_CC_7bp.txt"

barcodes_file <- args[5]
#barcodes_file <- "/home/pezoldt/NAS2/iris/4_genomes/barcodes/iris_genomics/20210618_scRNAseq_cellcode_BRBtoolbox_CC_7bp.txt"

sample_info_sheet <- args[6]
#sample_info_sheet <- "/home/pezoldt/NAS2/iris/0_seq_run_processing/c_SIS/sis_aris03_mapping_QC_pipeline_v2.csv"

ERCC_concentrations <- args[7]
#ERCC_concentrations <- "/home/pezoldt/NAS2/iris/4_genomes/genomes/ERCC/20201212_ERCC_concentrations.txt"

ERCC_gtf <- args[8]
#ERCC_gtf <- "/home/pezoldt/NAS2/iris/4_genomes/genomes/ERCC/ERCC92.GeneAdded.mCherry.EGFP.gtf"

# Print status of script 
cat("Read arguments from command line:\n")

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
cat("Created file directories\n")

#############################################
######## Global variables & cut-offs ########
#############################################
min_UMI <- 1000

#############################################
############## Load files  ##################
#############################################

# Load barcodes to CC file
t_barcodes_file <- read.delim(barcodes_file)

## Load ERCC information files 

# Load files with ERCC information
ERCC_concentrations <- read.delim(ERCC_concentrations)
ERCC_gtf <- read.delim(ERCC_gtf, header = FALSE)

# Set colanmes 
colnames(ERCC_gtf) <- c("ID", "content", "region", "chr", "length",
                        "k1", "k2", "k3", "k4")

# Merge file information 
ERCC_intel <- merge(ERCC_concentrations, ERCC_gtf, by.x = "ERCC.ID", by.y = "ID")


## Add number of ERCC molecules 

# Volume of ... per droplet in nl
volume <- 2.5

# Dilution (in 1000ul, used 1ul from 1:10 diluted stock)
dilution <- 10

# Molecules per attomol
attomol <- 6.022*10^23 * 10^(-18)

# Calculate ERCC molecules per droplet for two mixes 
ERCC_intel$Molecules_Droplet_Mix1 <- round((ERCC_intel$concentration.in.Mix.1..attomoles.ul. * attomol / dilution) / (1000 * 1000 * volume),3)
ERCC_intel$Molecules_Droplet_Mix2 <- round((ERCC_intel$concentration.in.Mix.2..attomoles.ul. * attomol / dilution) / (1000 * 1000 * volume),3)


## Create sample information with one row per cell 

# Load sample information sheet 
SIS <- read.csv(sample_info_sheet, skip = 1, header = T, stringsAsFactors = FALSE)
nrow(SIS)
#head(SIS)
#tail(SIS)
# Extract information of specific experiment from SIS
SIS_exp <- SIS[which(str_detect(SIS$Sample_ID, expID)),]
SIS_exp$Sample_ID_2 <- paste0(SIS_exp$Description,"_",SIS_exp$I5_Index_ID)

# Split SIS by Sample ID (i.e. by well)
SIS_list <- split(SIS_exp, SIS_exp$Sample_ID_2)

#####################################
# Create SIS with one row per cell ##
#####################################
#####
#One CC per cell
#####
for(i in 1:length(SIS_list)){
  #i = 1
  x = SIS_list[[i]]
  # Determine used barcodes
  barcodes_used <- str_split(x$"CC_used", "_") %>% unlist()
  
  # Repeat each well information by the amount of cells per well 
  n_cells <- length(barcodes_used)
  well <- x[rep(seq_len(nrow(x)), each = n_cells), ]
  
  # Add Lid_ID and Full_Cell_ID (= well ID + Lid ID) as columns to SIS 
  well$CC <- barcodes_used
  well <- unite(well, "Full_Cell_ID", c(Description, I5_Index_ID, CC), sep="_", remove = FALSE)
  SIS_list[[i]] <- well
}


# Combine information per cell to one data frame 
SIS_exp_perCell <- do.call(rbind, SIS_list)
rownames(SIS_exp_perCell) <- NULL

#####
#Allocate Dropouts per Plate
#####
#Check if DropOut annotation needs to be performed
#If more than one condition is present Dropouts have been annotated in SIS
#Select only plates with a dropout type
SIS_exp_perCell_Dropout <- subset(SIS_exp_perCell, Dropout_Cell_Type != "no")
SIS_exp_perCell_noDropout <- subset(SIS_exp_perCell, Dropout_Cell_Type == "no")

if(nrow(SIS_exp_perCell_Dropout) > 0){
  #List of plates dropping out
  l_plate <- split(SIS_exp_perCell_Dropout, SIS_exp_perCell_Dropout$Description)
  #Per plate
  l_Dropout_plate <- list()
  for(i in 1:length(l_plate)){
      #i = 1
      print("i")
      print(i)
      t_plate <- l_plate[[i]]
      #nrow(t_plate)
      Dropout_Type <- unlist(strsplit(t_plate$Dropout_Cell_Type[[1]], "__"))
      Dropout_CC <- unlist(strsplit(t_plate$Dropout_Cell_CC[[1]], "__"))
      Dropout_WC <- unlist(strsplit(t_plate$Dropout_Cell_WC[[1]], "__"))
      #Note: For each Dropout type there is a matching CC and a set of wells
      l_Dropout <- list()
      for(j in 1:length(Dropout_Type)){
        print("j")
        print(j)
        #j = 1
        Dropout_WC_j <- unlist(strsplit(Dropout_WC[j],"_"))
        Dropout_Type_j <- Dropout_Type[j]
        Dropout_CC_j <- Dropout_CC[j]
        Dropout_WC_j <- paste0("WC",
                             Dropout_WC_j[1]:Dropout_WC_j[2])
        Full_Cell_ID <- paste(t_plate$Description[1], Dropout_WC_j, Dropout_CC_j ,sep="_")
        t_Dropout_j <- data.frame(Full_Cell_ID = Full_Cell_ID,
                                Dropout_Type = rep(Dropout_Type_j, length(Full_Cell_ID)))
        l_Dropout[[j]] <- t_Dropout_j
      }
    t_Dropout_toMerge <- do.call(rbind, l_Dropout)
    #Merge with plate
    l_Dropout_plate[[i]] <- merge(t_plate, t_Dropout_toMerge, by = "Full_Cell_ID", all=TRUE)
  }
  SIS_exp_perCell_Dropout <- do.call(rbind, l_Dropout_plate)
  SIS_exp_perCell_Dropout <- replace(SIS_exp_perCell_Dropout, is.na(SIS_exp_perCell_Dropout), "no")
  SIS_exp_perCell_noDropout$Dropout_Type <- rep("no", nrow(SIS_exp_perCell_noDropout))
  #Combine
  SIS_exp_perCell <- rbind(SIS_exp_perCell_Dropout, SIS_exp_perCell_noDropout)
  SIS_exp_perCell <- SIS_exp_perCell[,!(colnames(SIS_exp_perCell) %in% c("Dropout_Cell_Type", "Dropout_Cell_CC", "Dropout_Cell_WC"))]
}else{
#If no Dropouts
SIS_exp_perCell$Dropout_Type <- rep("no", nrow(SIS_exp_perCell))
SIS_exp_perCell <- SIS_exp_perCell[,!(colnames(SIS_exp_perCell) %in% c("Dropout_Cell_Type", "Dropout_Cell_CC", "Dropout_Cell_WC"))]
}
#############################################
############## PreVariables #################
#############################################

## Define variables needed across the script

# Get names of well folders from rcm directory
well_folders <- paste0(rcm_path,"/",SIS_exp$Sample_ID_2)

# Get sample/well names 
well_names <- gsub(".*/", "", well_folders)


#############################################
############## QC for excluding wells #######
#############################################
cat("Number of processed wells\n")
print(nrow(SIS_exp))

#####
## Check if file "Log.final.out" exists 
#####
well_mapped_STARlogfile <- lapply(well_folders, function(x){
  log_file <- file.path(x, "Log.final.out")
  return(file.exists(log_file))
})
sum(as.numeric(well_mapped_STARlogfile))
#####
## Check if there is a "UMIperCellSorted.txt" present
#####
suppressWarnings(
  well_mapped_UMPperCell <- lapply(well_folders, function(x){
    #x = well_folders[[1]]
    #print(x)
    status_table_x <- try(read.table(paste0(x,"/Solo.out/GeneFull/UMIperCellSorted.txt")))
    if(inherits(status_table_x, "try-error"))
      return(FALSE)
    else
      return(TRUE)
  })
)
sum(as.numeric(well_mapped_UMPperCell))

#Combine dropout wells
n_assessed_qc = 2
dropouts <- as.numeric(well_mapped_STARlogfile) + as.numeric(well_mapped_UMPperCell)
#Assemble vector to define wells that passed the mapping...
v_pass_mapping <- replace(dropouts, dropouts==n_assessed_qc, "TRUE")
v_pass_mapping <- replace(v_pass_mapping, dropouts==1, "FALSE")
v_pass_mapping <- replace(v_pass_mapping, dropouts==0, "FALSE")
v_pass_mapping <- as.logical(v_pass_mapping)

#Build SIS store
SIS_store <- SIS_exp
#Mapping
SIS_store$Well_Star_mapped <- replace(SIS_store$Well_Star_mapped, as.numeric(well_mapped_STARlogfile)==1, "mapped")
SIS_store$Well_Star_mapped <- replace(SIS_store$Well_Star_mapped, as.numeric(well_mapped_STARlogfile)==0, "not_mapped")
#UMI
SIS_store$UMI_detected <- replace(SIS_store$UMI_detected, as.numeric(well_mapped_UMPperCell)==1, "yes")
SIS_store$UMI_detected <- replace(SIS_store$UMI_detected, as.numeric(well_mapped_UMPperCell)==0, "no")

#Adapt SIS
SIS_exp <- subset(SIS_exp, Sample_ID_2 %in% well_names[v_pass_mapping])
cat("Number of processed wells minus failed mapping and/or no UMIs\n")
print(nrow(SIS_exp))

#############################################
####### Re-define Variables #################
#############################################

## Re-Define variables needed across the script

# Get names of well folders from rcm directory
well_folders <- paste0(rcm_path,"/",SIS_exp$Sample_ID_2)

# Get sample/well names 
well_names <- gsub(".*/", "", well_folders)


#############################################
######### STARsolo stats per well ###########
#############################################

## Process stats of STARsolo Log.final.out (stats per well)

# Extract statistic information per well from STARsolos Log.final.out 
STARsolo_well_stats <- lapply(well_folders, function(x){
  log_file <- file.path(x, "Log.final.out")
  log_table <- read.delim(log_file, header = F)  %>%
    filter(grepl("Number of input reads", V1) |     
             grepl("Uniquely mapped reads", V1) |
             grepl("Uniquely mapped reads %", V1)) %>%
    set_names("variables", "values")
  return(log_table)
  
})

# Rename stats variables and set values as numeric  
STARsolo_well_stats <- lapply(STARsolo_well_stats, function(x){
  
  variables <- c("n_Reads_input", "n_Reads_uniqueMapped", "%_Reads_mapped")    
  values <- map(x$values, function(x){gsub("%", "", x) %>% as.numeric(x)}) %>% unlist()
  table <- data.frame(variables = variables, values = values)
  return(table)
  
})

# Set well names as names of stats list 
names(STARsolo_well_stats) <- well_names

# Transpose STARsolo alignment stats per well 
STARsolo_well_stats <- lapply(seq_along(STARsolo_well_stats), function(x){
  
  align_per_row <- data.frame(Fastq_ID = well_names[x],
                              STARsolo_i5_n_reads_sequenced = STARsolo_well_stats[[x]]$value[1],
                              STARsolo_i5_n_reads_unique_mapped = STARsolo_well_stats[[x]]$value[2],
                              STARsolo_i5_percent_reads_unique_mapped = STARsolo_well_stats[[x]]$value[3])
  return(align_per_row)
  
})

# Create one data frame with alignment STARsolo stats per well 
data_STARsolo_well_stats <- do.call(rbind, STARsolo_well_stats)

cat("Compiled STARsolo stats\n")

#############################################
#### Mapping stats averaged per cell/well ###
#############################################
l_mapping_statsPERwell <- list()
mappingWELLcell_folder <- lapply(well_folders, function(x) {file.path(x, "Solo.out")})
for(i in 1:length(mappingWELLcell_folder)){
  #i = 19
  #print(i)
  mappingWELLcell_folder[[i]]
  #UMI stats
  CELLcount_i <- read.delim(paste0(mappingWELLcell_folder[[i]],"/GeneFull/UMIperCellSorted.txt"), header = FALSE)
  CELLcount_i <- sum(as.numeric(CELLcount_i$V1 > min_UMI))
  #cells above threshold min_UMI
  t_n_cells_i <- data.frame(Feature = "n_cells", n_reads = CELLcount_i)
  #barcode stats
  n_discarded_barcode_i <- read.table(paste0(mappingWELLcell_folder[[i]],"/Barcodes.stats"), header = FALSE)
  n_all_reads <- sum(n_discarded_barcode_i$V2)
  n_discarded_barcode_i <- n_discarded_barcode_i[n_discarded_barcode_i$V1 %in% c("nNinUMI","nUMIhomopolymer","nNoMatch"),]
  n_discarded_barcode_i <- sum(n_discarded_barcode_i$V2)
  #Reads that are duplicated/triplicated relating for a UMI at transcript
  t_discarded_barcode_i <- data.frame(Feature = "n_discarded_at_CC", n_reads = n_discarded_barcode_i)
  #Total Reads
  t_all_reads_i <- data.frame(Feature = "n_all_reads", n_reads = n_all_reads)
  
  #feature-stats
  features_i <- read.table(paste0(mappingWELLcell_folder[[i]],"/GeneFull/Features.stats"), header = FALSE)
  features_summary_i <- features_i[features_i$V1 %in% c("nUnmapped","nNoFeature","nAmbigFeature","nUMIs","nMatchUnique"),]
  colnames(features_summary_i) <- c("Feature", "n_reads")
  duplicate_UMI_i <- subset(features_summary_i, Feature == "nMatchUnique")$n_reads - subset(features_summary_i, Feature == "nUMIs")$n_reads
  #Bug in STAR solo for low CC matches, give then nUMI of 1B or higher
  if(duplicate_UMI_i > 100000000){
    duplicate_UMI_i = 0
  }
  #Reads that are duplicated/triplicated relating for a UMI at transcript
  t_duplicate_UMI_i <- data.frame(Feature = "n_saturated_UMI", n_reads = duplicate_UMI_i)
  #Put table together
  features_summary_i <- rbind(t_n_cells_i,t_all_reads_i,t_discarded_barcode_i,features_summary_i, t_duplicate_UMI_i)
  features_summary_i$Feature <- c("n_cells","n_all_reads","n_discarded_at_CC","n_Unmapped", "n_noFeature","n_ambigFeature","n_read_unique","n_UMI","n_saturated_UMI")
  features_summary_i <- features_summary_i[features_summary_i$Feature %in% c("n_cells","n_all_reads","n_discarded_at_CC","n_Unmapped", "n_noFeature","n_ambigFeature","n_UMI","n_saturated_UMI"),]
  #Store data
  l_mapping_statsPERwell[[i]] <- features_summary_i
}
names(l_mapping_statsPERwell) <- well_names

#Make table to merge with SIS_store
l_pivot_mapping_statsPERwell <- list()
for(i in 1:length(l_mapping_statsPERwell)){
  #i = 1
  mapping_statsPERwell_i <- l_mapping_statsPERwell[[i]]
  mapping_statsPERwell_i <- pivot_wider(mapping_statsPERwell_i, names_from = c("Feature"), values_from = "n_reads")
  l_pivot_mapping_statsPERwell[[i]] <- mapping_statsPERwell_i
  
}
t_pivot_mapping_statsPERwell <- as.data.frame(do.call(rbind, l_pivot_mapping_statsPERwell))
t_pivot_mapping_statsPERwell$Sample_ID <- well_names

cat("Compiled all and mapped reads per well\n")

#############################################
####### All & mapped Reads per cell #########
#############################################
# Load folder for read count stats 
READ_count_files <- lapply(well_folders, function(x) {file.path(x, "Read_counts")})

#Store read counts per well
l_reads_all_mapped <- list()### Per well
for(i in 1:length(READ_count_files)){
  #i = 1
  print(i)
  #Get WellName for Full_Cell_ID
  ID_i <- unlist(strsplit(READ_count_files[[i]],"/"))[[length(unlist(strsplit(well_folders[[i]],"/")))]]
  # Determine which barcodes were used 
  barcodes_used <- str_split(SIS_exp[i, "CC_used"], "_") %>% unlist()
  #Table all reads
  reads_noCC_present <- read.table(paste0(READ_count_files[[i]],"/CC_all_reads.txt"), header = FALSE, fill=TRUE)
  #Check if CC is present in first line
  CC_present <- nchar(reads_noCC_present[1,2]) > 1
  if(CC_present == TRUE){
    t_all_i <- read.table(paste0(READ_count_files[[i]],"/CC_all_reads.txt"), header = FALSE)
  }else{
    t_all_i <- read.table(paste0(READ_count_files[[i]],"/CC_all_reads.txt"), header = FALSE, skip = 1)
  }
  colnames(t_all_i) <- c("n_all_reads","B1")
  #Table mapped reads
  reads_noCC_present <- read.table(paste0(READ_count_files[[i]],"/CC_mapped_reads.txt"), header = FALSE, fill=TRUE)
  #Check if CC is present in first line
  CC_present <- nchar(reads_noCC_present[1,2]) > 1
  if(CC_present == TRUE){
    t_mapped_i <- read.table(paste0(READ_count_files[[i]],"/CC_mapped_reads.txt"), header = FALSE)
  }else{
    t_mapped_i <- read.table(paste0(READ_count_files[[i]],"/CC_mapped_reads.txt"), header = FALSE, skip = 1)
  }
  colnames(t_mapped_i) <- c("n_mapped_reads","B1")

  #CCs present
  t_barcodes_file_i <- t_barcodes_file[t_barcodes_file$Name %in% barcodes_used,]
  #Merge read coutns with barcode files
  t_CC_all_i <- merge(t_barcodes_file_i, t_all_i, by = "B1")
  t_CC_all_read_i <- merge(t_CC_all_i, t_mapped_i, by = "B1")
  t_CC_all_read_i <- t_CC_all_read_i[,c("Name","n_all_reads","n_mapped_reads")]
  t_CC_all_read_i$Name <- paste0(ID_i,"_",t_CC_all_read_i$Name)
  colnames(t_CC_all_read_i) <- c("Full_Cell_ID", "n_all_reads", "n_mapped_reads")
  #Store in list
  l_reads_all_mapped[[i]] <- t_CC_all_read_i
}
t_reads_all_mapped <- as.data.frame(do.call(rbind,l_reads_all_mapped))
rownames(t_reads_all_mapped) <- NULL

cat("Compiled reads per CellCode\n")

#############################################
######### UMI count matrices  ###############
#############################################

## UMI read count matrices from STARsolo "Solo.out/Gene/filtered"
# Load UMI stats files 
STARsolo_UMI_files <- lapply(well_folders, function(x) {file.path(x, "Solo.out/GeneFull/filtered")})

## Extract UMI matrices
# Create expression matrices with rownames = ENSG id OR gene symbols and colnames = Full_Cell_ID
UMI_SYMB_list <- list()
UMI_ENSG_list <- list()
UMI_allDetected_CCs_SYMB_list <- list()
#l_store_UMI_allDetected_CCs_ENSG -> list()
v_store_CC_state <- c()
#UMI_mats <- lapply(seq_along(STARsolo_UMI_files), function(x){
for(i in 1:length(STARsolo_UMI_files)){
  #i = 723
  print(i)
  x = STARsolo_UMI_files[[i]]
  #print(x)
  # Determine which barcodes were used 
  barcodes_used <- str_split(SIS_exp[i, "CC_used"], "_") %>% unlist()
  
  # Load UMI data per well i.e. barcodes, features and expression matrix
  barcodes <- read.delim(file.path(STARsolo_UMI_files[[i]], "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
  features <- read.delim(file.path(STARsolo_UMI_files[[i]], "features.tsv"), header = FALSE, stringsAsFactors = FALSE)
  matrix <- readMM(file.path(STARsolo_UMI_files[[i]], "matrix.mtx"))
  matrix <- as.data.frame(as.matrix(matrix))
  
  # Extract ENSG gene ids and gene symbols 
  ENSGids <- features$V1
  GeneSymbols <- features$V2
  
  # Match and replace the barcode sequences with the Lid ID 
  barcodes_detected <- merge(t_barcodes_file, barcodes, by.x = "B1", by.y = "V1")
  colnames_matrix <- barcodes_detected$Name 
  colnames(matrix) <- colnames_matrix
  ## Store all detected CCs
  matrix_allCCs <- matrix
  colnames(matrix_allCCs) <- paste0(well_names[i], "_", colnames(matrix_allCCs))
  UMI_allDetected_CCs_SYMB_list[[i]] <- matrix_allCCs
  
  # Only keep and store cells of the used barcodes
  ## If there is no match between CellCodes identified and CellCodes noted in SIS
  ## Keep CC with highest UMI count
  # If there is no match all columns will be set to FALSE and thus the sum is Zero
  if(!(sum(as.numeric(colnames(matrix) %in% barcodes_used)) >= 1)){
    max(colSums(matrix))
    #matrix <- matrix[,sort(colSums(matrix))]
    colCount = colSums(matrix)
    topColSums = order(colCount,decreasing=TRUE)[1]
    matrix = matrix[1:nrow(matrix),topColSums, drop = FALSE]
    CC_state <- "Non_of_expected_CCs"
    colnames(matrix) <- paste0(well_names[i], "_", colnames(matrix))
  }else{
    matrix <- matrix[which(colnames(matrix) %in% barcodes_used)]
    head(matrix)
    colnames(matrix) <- paste0(well_names[i], "_", colnames(matrix))
    CC_state <- "Expected_CCs"
  }
  # Set ENSG gene ids and GeneSymbols as rownames of expression matrix
  matrix_ensg <- matrix 
  rownames(matrix_ensg) <- ENSGids
  matrix_symb <- as.matrix(matrix) # matrix allows for duplicated rownames 
  rownames(matrix_symb) <- GeneSymbols
  
  # Return two matrices, one with ENSGids and one with GeneSymbols
  #return(list(matrix_ensg, matrix_symb))
  UMI_SYMB_list[[i]] <- matrix_symb
  UMI_ENSG_list[[i]] <- matrix_ensg
  v_store_CC_state[i] <- CC_state
}


names(UMI_SYMB_list) <- well_names
names(UMI_ENSG_list) <- well_names
names(UMI_allDetected_CCs_SYMB_list) <- well_names
t_store_CC_state <- data.frame(Sample_ID = well_names, CC_state = v_store_CC_state)

## Bind expression data of all cells with rownames = ENSGid 
# Bind expression data of all cells to one expression matrix
UMI_all_cells_ENSG <- as.matrix(do.call(cbind, UMI_ENSG_list))
# Set GeneSymbols as rownames
ENSGids <- rownames(UMI_ENSG_list[[1]])
rownames(UMI_all_cells_ENSG) <- ENSGids
# Filter out genes that are not expressed in any cell  
UMI_all_cells_ENSG <- UMI_all_cells_ENSG[which(rowSums(UMI_all_cells_ENSG) > 0), ]

## Bind expression data of all cells with rownames = GeneSymbol
# Bind expression data of all cells to one expression matrix
UMI_all_cells_SYMB <- as.matrix(do.call(cbind, UMI_SYMB_list))
# Set GeneSymbols as rownames
GeneSymbols <- rownames(UMI_SYMB_list[[1]])
rownames(UMI_all_cells_SYMB) <- GeneSymbols
# Filter out genes that are not expressed in any cell  
UMI_all_cells_SYMB <- UMI_all_cells_SYMB[which(rowSums(UMI_all_cells_SYMB) > 0), ]

## Bind expression data of all detected cell codes with rownames = GeneSymbol
# Bind expression data of all cells to one expression matrix
UMI_all_detected_CellCodes_SYMB <- as.matrix(do.call(cbind, UMI_allDetected_CCs_SYMB_list))
# Set GeneSymbols as rownames
GeneSymbols <- rownames(UMI_SYMB_list[[1]])
rownames(UMI_all_detected_CellCodes_SYMB) <- GeneSymbols
# Filter out genes that are not expressed in any cell  
UMI_all_detected_CellCodes_SYMB <- UMI_all_detected_CellCodes_SYMB[which(rowSums(UMI_all_detected_CellCodes_SYMB) > 0), ]

cat("Compiled UMI count matrix\n")

#####
## Extract overall UMI stats per cells
#####
# Extract UMI stats under regard of CC number per well
UMI_stats <- lapply(seq_along(UMI_SYMB_list), function(x){
  # Loop through list 
  UMI_table <-  UMI_SYMB_list[[x]]
  
  # Check if table only has only one column = one CC per weill
  if(ncol(UMI_table) == 1){
    
    # Get number of ERCC per cell 
    n_UMI_ERCC <- as.vector(sum(UMI_table[grepl("ERCC-00", rownames(UMI_table)),]))
    
    # Get number of mitochondrial UMI per cell 
    #n_UMI_mito <- as.vector(sum(UMI_table[grepl("^MT-", rownames(UMI_table)),]))
    
    # Get number of UMI per cell (= sum all counts per cell)
    n_UMI <- as.vector(sum(UMI_table))
    
    # Get number of detected genes per cell (i.e. amount of genes with a count > 0)
    n_genes <- apply(UMI_table, 2, function(x){sum(x > 0)})
    
    # Get number of detected mitochondrial genes per cell 
    #n_genes_mito <- as.data.frame(UMI_table[which(grepl("MT-", rownames(UMI_table))), ])
    #n_genes_mito <- apply(n_genes_mito, 2, function(x){sum(x > 0)})
    
    # Percentage of mito UMI 
    #UMI_mito_percent <- n_UMI_mito/n_UMI
    
    #Percentage of mito genes 
    #mito_percentage <- n_genes_mito/n_genes
    
    # more than one CC per well  
  }else{
    
    # Get number of ERCC per cell 
    n_UMI_ERCC <- colSums(UMI_table[grepl("ERCC-00", rownames(UMI_table)),])
    
    # Get number of mitochondrial UMI per cell 
    #n_UMI_mito <- colSums(UMI_table[grepl("^MT-", rownames(UMI_table)),])
    
    # Get number of UMI per cell (= sum all counts per cell)
    n_UMI <- colSums(UMI_table)
    
    # Get number of detected genes per cell (i.e. amount of genes with a count > 0)
    n_genes <- apply(UMI_table, 2, function(x){sum(x > 0)})
    
    # Get number of detected mitochondrial genes per cell 
    #n_genes_mito <- UMI_table[which(grepl("^MT-", rownames(UMI_table))), ]
    #n_genes_mito <- apply(n_genes_mito, 2, function(x){sum(x > 0)})
    
    # Percentage of mito UMI 
    #UMI_mito_percent <- n_UMI_mito/n_UMI
    
    #Percentage of mito genes 
    #mito_percentage <- n_genes_mito/n_genes
    
  }
  
  # Create dataframe of stats 
  if(ERCC  == "yes"){
    
    temp_umi <- data.frame(n_UMI = n_UMI,
                           n_UMI_ERCC = n_UMI_ERCC,
                           #n_UMI_mito = n_UMI_mito,
                           n_genes = n_genes,
                           #n_genes_mito = n_genes_mito,
                           #percent_UMI_mito = UMI_mito_percent,
                           #percent_genes_mito = mito_percentage,
                           Full_Cell_ID = paste0(names(n_genes)))
    
    return(temp_umi)
    
  }else{
    
    temp_umi <- data.frame(n_UMI = n_UMI,
                           #n_UMI_mito = n_UMI_mito,
                           n_genes = n_genes,
                           #n_genes_mito = n_genes_mito,
                           #percent_UMI_mito = UMI_mito_percent,
                           #percent_genes_mito = mito_percentage,
                           Full_Cell_ID = paste0(names(n_genes)))
    
    return(temp_umi)
    
  }
  
})

# Create one data frame of all UMI stats per cell
UMI_stats <- do.call(rbind, UMI_stats)
rownames(UMI_stats) <- NULL

#####
## Extract overall UMI stats of all detected CellCodes
#####
# Extract UMI stats under regard of CC number per well
UMI_stats_detected_CellCodes <- lapply(seq_along(UMI_allDetected_CCs_SYMB_list), function(x){
  # Loop through list 
  UMI_table <-  UMI_allDetected_CCs_SYMB_list[[x]]
  
  # Check if table only has only one column = one CC per weill
  if(ncol(UMI_table) == 1){
    
    # Get number of ERCC per cell 
    n_UMI_ERCC <- as.vector(sum(UMI_table[grepl("ERCC-00", rownames(UMI_table)),]))
    
    # Get number of mitochondrial UMI per cell 
    #n_UMI_mito <- as.vector(sum(UMI_table[grepl("^MT-", rownames(UMI_table)),]))
    
    # Get number of UMI per cell (= sum all counts per cell)
    n_UMI <- as.vector(sum(UMI_table))
    
    # Get number of detected genes per cell (i.e. amount of genes with a count > 0)
    n_genes <- apply(UMI_table, 2, function(x){sum(x > 0)})
    
    # Get number of detected mitochondrial genes per cell 
    #n_genes_mito <- as.data.frame(UMI_table[which(grepl("MT-", rownames(UMI_table))), ])
    #n_genes_mito <- apply(n_genes_mito, 2, function(x){sum(x > 0)})
    
    # Percentage of mito UMI 
    #UMI_mito_percent <- n_UMI_mito/n_UMI
    
    #Percentage of mito genes 
    #mito_percentage <- n_genes_mito/n_genes
    
    # more than one CC per well  
  }else{
    
    # Get number of ERCC per cell 
    n_UMI_ERCC <- colSums(UMI_table[grepl("ERCC-00", rownames(UMI_table)),])
    
    # Get number of mitochondrial UMI per cell 
    #n_UMI_mito <- colSums(UMI_table[grepl("^MT-", rownames(UMI_table)),])
    
    # Get number of UMI per cell (= sum all counts per cell)
    n_UMI <- colSums(UMI_table)
    
    # Get number of detected genes per cell (i.e. amount of genes with a count > 0)
    n_genes <- apply(UMI_table, 2, function(x){sum(x > 0)})
    
    # Get number of detected mitochondrial genes per cell 
    #n_genes_mito <- UMI_table[which(grepl("^MT-", rownames(UMI_table))), ]
    #n_genes_mito <- apply(n_genes_mito, 2, function(x){sum(x > 0)})
    
    # Percentage of mito UMI 
    #UMI_mito_percent <- n_UMI_mito/n_UMI
    
    #Percentage of mito genes 
    #mito_percentage <- n_genes_mito/n_genes
    
  }
  
  # Create dataframe of stats 
  if(ERCC  == "yes"){
    
    temp_umi <- data.frame(n_UMI = n_UMI,
                           n_UMI_ERCC = n_UMI_ERCC,
                           #n_UMI_mito = n_UMI_mito,
                           n_genes = n_genes,
                           #n_genes_mito = n_genes_mito,
                           #percent_UMI_mito = UMI_mito_percent,
                           #percent_genes_mito = mito_percentage,
                           Full_Cell_ID = paste0(names(n_genes)))
    
    return(temp_umi)
    
  }else{
    
    temp_umi <- data.frame(n_UMI = n_UMI,
                           #n_UMI_mito = n_UMI_mito,
                           n_genes = n_genes,
                           #n_genes_mito = n_genes_mito,
                           #percent_UMI_mito = UMI_mito_percent,
                           #percent_genes_mito = mito_percentage,
                           Full_Cell_ID = paste0(names(n_genes)))
    
    return(temp_umi)
    
  }
  
})

# Create one data frame of all UMI stats per cell
UMI_stats_detected_CellCodes <- do.call(rbind, UMI_stats_detected_CellCodes)
rownames(UMI_stats_detected_CellCodes) <- NULL

cat("Compiled complete UMI stats\n")


#############################################
########## Create pivot table  ##############
#############################################

## Merge all stats per cell to one data frame 

# Merge stats from STARsolo 
merged_stats <- inner_join(UMI_stats, t_reads_all_mapped, by = "Full_Cell_ID")

# Merge stats and SIS 
pivot_table <- inner_join(SIS_exp_perCell, merged_stats, by = "Full_Cell_ID")


# Add ratio UMI to mapped reads
pivot_table$ratio_UMI_to_mapped_reads <- pivot_table$n_UMI / pivot_table$n_mapped_reads

# Add percentage of ERCC
if(ERCC == "yes"){
  pivot_table$percent_UMI_ERCC <- pivot_table$n_UMI_ERCC/pivot_table$n_UMI * 100
}

cat("Pivot table created\n")

#######################################
####Add Condition per CC ##############
#######################################
meta_lid <- pivot_table
#Split by plate
l_meta_lid <- split(meta_lid, meta_lid$Description)
l_meta_lid_addon <- list()
for(i in 1:length(l_meta_lid)){
  #i = 1
  meta_lid_i <- l_meta_lid[[i]]
  CCs_i <- unlist(strsplit(meta_lid_i$CC_used[1], "_"))
  Condition_per_CC_i <- unlist(strsplit(meta_lid_i$Condition_per_CC[1], "_"))
  names(CCs_i) <- Condition_per_CC_i
  CC_condition <- c()
  for(j in 1:nrow(meta_lid_i)){
    #print(nrow(meta_lid_i))
    #j = 1
    CC_j <- meta_lid_i[j,]$CC
    CC_condition_j <- names(CCs_i[CCs_i %in% CC_j])
    CC_condition[j] <- CC_condition_j
  }
  meta_lid_i$CellType_Condition <- CC_condition
  l_meta_lid_addon[[i]] <- meta_lid_i
}
lapply(l_meta_lid_addon, head)
lapply(l_meta_lid_addon, nrow)
#Pivot Table
meta_lid <- do.call(rbind, l_meta_lid_addon)
pivot_table <- meta_lid

#######################################
#### Data preparation for plotting ####
#######################################

## Filtered read count/UMI matrices 

# Apply threshhold for UMI and create separated ERCC objects
#Note: If ERCC section of the pipeline is to be used it is required to pull data from the STARsolo output
#Note: The BRBseqtools output is obsolete
#if(ERCC == "yes"){
  
  # Filter out cells with less than 1000 UMI 
  #pivot_UMI_tresh <- subset(pivot_table, n_UMI > min_UMI)
  
  # Create filtered RCM of mapped reads (filtered by UMI threshhold)
  #rcm_mapped_threshed <- rcm_all_cells[,colnames(rcm_all_cells) %in% pivot_UMI_tresh$Full_Cell_ID]
  
  # Exclude ERCCs from mapped RCM and save them separately 
  #rcm_mapped_ERCC <- rcm_mapped_threshed[grepl("ERCC-", rownames(rcm_mapped_threshed)),]
  #rcm_mapped_endogenous <- rcm_mapped_threshed[!(grepl("ERCC-",rownames(rcm_mapped_threshed))),]
  
  # Create filtered RCM of UMI 
  #rcm_UMI_threshed <- UMI_all_cells_SYMB[,colnames(UMI_all_cells_SYMB) %in% pivot_UMI_tresh$Full_Cell_ID]
  
  # Exclude ERCCs filtered UMI matrix and save them separately 
  #rcm_UMI_ERCC <- rcm_UMI_threshed[grepl("ERCC-", rownames(rcm_UMI_threshed)),]
  #rcm_UMI_endogenous <- rcm_UMI_threshed[!(grepl("ERCC-", rownames(rcm_UMI_threshed))),]
  
#}else{
  
  # Filter out cells with less than 1000 UMI 
  #pivot_UMI_tresh <- subset(pivot_table, n_UMI > min_UMI)
  
  # Create filtered RCM of mapped reads (filtered by UMI threshhold)
  #rcm_mapped_threshed <- rcm_all_cells[,colnames(rcm_all_cells) %in% pivot_UMI_tresh$Full_Cell_ID]
  
  # Create filtered RCM of UMI 
  #rcm_UMI_threshed <- UMI_all_cells_SYMB[,colnames(UMI_all_cells_SYMB) %in% pivot_UMI_tresh$Full_Cell_ID]
#}


#######################################
##### Additional ERCC stats ###########
#######################################

## ERCC per cell according to abundance

#if(ERCC == "yes"){
  
  # Get order of ERCCs by abundance of each ERCC across all cells 
  #ERCC_order <- names(sort(rowSums(rcm_UMI_ERCC)))
  
  # Order ERCC by abundance 
  #rcm_UMI_ERCC_order <- rcm_UMI_ERCC[ERCC_order,]
  #pivot_ERCC <- as.data.frame(rcm_UMI_ERCC_order)
  
  # Add ERCC ID as column
  #pivot_ERCC$ERCC_ID <- rownames(pivot_ERCC)
  
  # Create overview of ERCC count per cell 
  #pivot_ERCC <- pivot_longer(pivot_ERCC, cols = 1:ncol(rcm_UMI_ERCC_order), names_to = "Cell", values_to = "UMI_ERCC")
  
  # Add log10 of ERCC_UMI
  #pivot_ERCC$log10_UMI_ERCC <- log10(pivot_ERCC$UMI_ERCC + 1)
  
  # Change ERCC_ID column to factor 
  #pivot_ERCC$ERCC_ID <- as.factor(pivot_ERCC$ERCC_ID)
  
  # Save ERCC pivot table as data frame
  #pivot_ERCC <- as.data.frame(pivot_ERCC)
  
  # Get mean per ERCC across all cells 
  #ERCC_UMI_mean <- data.frame(ERCC.ID = rownames(rcm_UMI_ERCC), mean_UMI_ERCC_detected = rowMeans(rcm_UMI_ERCC))
  
  # Merge ERCC information with mean value 
  #ERCC_intel_UMI_mean <- merge(ERCC_intel, ERCC_UMI_mean, by = "ERCC.ID")
#}

#######################################
##### Additional SIS intel ############
#######################################
SIS_store <- merge(SIS_store, t_store_CC_state, by.x = "Sample_ID_2", by.y = "Sample_ID", all=TRUE)
nrow(SIS_store)
SIS_store$CC_state[is.na(SIS_store$CC_state)] <- "notMapped_noUMI"
#Merge with read data
SIS_store <- merge(SIS_store, t_pivot_mapping_statsPERwell, by.x = "Sample_ID_2", by="Sample_ID", all = TRUE)
SIS_store <- SIS_store[,!colnames(SIS_store) %in% c("Sample_ID")]
colnames(SIS_store)[1] <- "Sample_ID"

#######################################
#### Save data needed for plotting ####
#######################################

## Save object needed for plotting in one file

#if(ERCC == "yes"){
  
 # save(pivot_table, rcm_all_cells, UMI_all_cells_SYMB, UMI_all_cells_ENSG, UMI_stats_detected_CellCodes, unknown_barcodes, pivot_UMI_tresh, rcm_mapped_threshed, rcm_UMI_threshed, 
  #     rcm_mapped_endogenous, rcm_UMI_endogenous, rcm_mapped_ERCC, rcm_UMI_ERCC, pivot_ERCC, ERCC_intel_UMI_mean, file = file.path(objects_path, "objects_plotting.RData"))
  
#}else{
  
  save(pivot_table, UMI_all_cells_SYMB, UMI_all_cells_ENSG, UMI_stats_detected_CellCodes, file = file.path(objects_path, "objects_plotting.RData"))
#}


## Save the pivot table seperately 
saveRDS(pivot_table, file = file.path(objects_path, "pivot_table.Rds"))
## Save SIS with mapping information per well
saveRDS(SIS_store, file = file.path(objects_path, "SIS_wellStats.Rds"))

## Save conventional tables
write.table(pivot_table, paste0(objects_path,"/pivot_table.txt"), sep="\t")
write.table(UMI_all_cells_SYMB, paste0(objects_path,"/UMI_all_cells_SYMB.txt"), sep="\t")
write.table(UMI_all_cells_ENSG, paste0(objects_path,"/UMI_all_cells_ENSG.txt"), sep="\t")

## Save MTX data
dir.create(paste0(objects_path,"/mtx"))
gbm <- UMI_all_cells_SYMB
# save sparse matrix
sparse.gbm <- Matrix(gbm , sparse = T )
## Market Exchange Format (MEX) format
writeMM(obj = sparse.gbm, paste0(objects_path,"/mtx/matrix.mtx"))
# save genes and cells names
write(x = rownames(gbm), paste0(objects_path,"/mtx/genes.tsv"))
write(x = colnames(gbm), paste0(objects_path,"/mtx/barcodes.tsv"))


# Print status of script
cat("Saved objects for plotting\n")

