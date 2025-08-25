#!/bin/bash

# Author: Joern Pezoldt
# Date: Summer 2023
# Description: Count all reads per CB

source /opt/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate sc_rna_seq_v2

# Define directory paths 
# Adapt
PATH_bam="/data/4_sequencing_analysis/JP131/results/human/STARsoloAligned.sortedByCoord.out.bam"
PATH_output_mapped="/data/4_sequencing_analysis/JP131/results/mapped_supi.bam"
PATH_output_CB="/data/4_sequencing_analysis/JP131/results/human/reads_per_barcode.txt"
PATH_output_CBmapped="/data/4_sequencing_analysis/JP131/results/human/mappedreads_per_barcode.txt"

#Get mapped reads
samtools view -b -F 4 $PATH_bam > $PATH_output_mapped

#Get total reads per CB
samtools view $PATH_bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $PATH_output_CB

#Get totamapped reads per CB
samtools view $PATH_output_mapped | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $PATH_output_CBmapped

