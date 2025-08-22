suppressPackageStartupMessages({
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(readr)
library(tibble)
})
  
# Read variables from command line (T = drop internal arguments)
args <- commandArgs(T)

expID <- args[1]
barcodes_brb_file <- args[2]
sample_info_sheet <- args[3]
temp_dir <- args[4]

# Load brb barcodes 
barcodes_brb <- read.delim(barcodes_brb_file)

# Load sample information sheet 
SIS <- read.csv(sample_info_sheet, skip = 1, header = T, stringsAsFactors = FALSE)

# Extract information of specific experiment from SIS
SIS_exp <- SIS[which(str_detect(SIS$Sample_ID, expID)),]

# Determine which barcodes were used 
barcodes_used <- str_split(SIS_exp$CC_used, "_") %>% unlist() %>% unique()

# Create file with subset of barcodes based on used barcodes
barcode_subset <- barcodes_brb[which(barcodes_brb$Name %in% barcodes_used),]

# Write subset to file
write.table(barcode_subset, file = file.path(temp_dir,"barcodes_subset.txt"), sep = "\t", row.names = F, quote = F)

