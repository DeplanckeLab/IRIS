#!/bin/bash

# Author: Joern Pezoldt
# Date: Summer 2023
# Description: Process 10X fastq to QC and count matrices

# Path to conda environments
source /opt/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate sc_rna_seq_v2

# Define directory paths 
#Adapt
PATH_output="/YOURPATH"
#GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148091
PATH_fastq="/YOURPATH/fastq"
PATH_starindexdir="/data/genome/homo_sapiens/GRCh38.104/STAR_Index"
read2_1=$PATH_fastq/"HEK_10X_S6_R2_001.fastq.gz"
read1_1=$PATH_fastq/"HEK_10X_S6_R1_001.fastq.gz"
#read2_2=$PATH_fastq/"hgmm_1k_v3_S1_L002_R2_001.fastq.gz"
#read1_2=$PATH_fastq/"hgmm_1k_v3_S1_L002_R1_001.fastq.gz"


#Permanent
cellcodes_wl="YOURPATH/10X_genomics/low_throughput_9K-green19-LT-dev.txt"

	## Run STARsolo to map and demultiplex 

	# See https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf for explanation of parameters 
	# --outFileNamePrefix : Output directory 
	# --readFilesCommand  : zcat - to uncompress .gz files (execute for every file)
	# --readFilesIn : For paired-end reads use comma separated list for read1 followed by space followed by comma separated list for read2 (!!! 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read !!!)
	# --soloType : type of single-cell RNA-seq (CB_UMI_Simple: one UMI and one Cell Barcode of fixed length in read1; CB_samTagOut: output cell barcode as CR and/or CB SAM tag, requires -outSAMtype BAM Unsorted and/or SortedByCoordiante)  
	# --soloCBstart : cell barcode start base
	# --soloCBlen : cell barcode length
	# --soloUMIstart : UMI start base 
	# --soloUMIlen : UMI length
	# --soloCBmatchWLtype : matching the Cell Barcodes to the WhiteList (1MM: only one match in whitelist with 1 mismatched base allowed)
	# --soloCBwhitelist : file(s) with whitelist(s) of cell barcodes
	# --soloUMIdedup : type of UMI deduplication (collapsing) algorithm (1MM_All: all UMIs with 1 mismatch distance to each other are collapsed (i.e.counted once))
	# --outSAMtype : output of SAM/BAM file (first word: BAM, SAM or NONE; second word: Unsorted (needed for soloType), SortedByCoordinate)
	# --outSAMattributes :  a string of desired SAM/BAM attributes
	# --outSAMunmapped : output of unmapped reads in the SAM format 
	# --limitBAMsortRAM : >=0 maximum available RAM (bytes) for sorting BAM; must be > 0 when --genomeLoad is LoadAndKeep

STAR --genomeDir $PATH_starindexdir \
	 --runThreadN 10 \
	 --outFileNamePrefix $PATH_output \
	 --readFilesCommand zcat \
	 --readFilesIn $read2_1 $read1_1 \
	 --soloType CB_UMI_Simple \
	 --outFilterScoreMin 30 \
	 --soloCBstart 1 \
	 --soloCBlen 16 \
	 --soloUMIstart 17 \
	 --soloUMIlen 12 \
	 --soloCBwhitelist $cellcodes_wl \
	 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
	 --soloUMIfiltering MultiGeneUMI_CR \
	 --soloUMIdedup 1MM_CR \
	 --clipAdapterType CellRanger4 \
	 --soloCellFilter EmptyDrops_CR \
	 --outSAMtype BAM SortedByCoordinate \
	 --soloFeatures Gene GeneFull SJ Velocyto \
	 --soloMultiMappers Uniform PropUnique EM Rescue \
	 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
	 --outSAMunmapped Within \
	 --limitBAMsortRAM 10000000000






