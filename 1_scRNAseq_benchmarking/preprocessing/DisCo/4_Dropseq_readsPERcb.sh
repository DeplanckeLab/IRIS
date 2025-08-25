#!/bin/bash

# Author: Joern Pezoldt
# Date: Summer 2023
# Description: Count all reads per CB

# Path to conda environments
source /opt/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate sc_rna_seq_v2


# Define directory paths 
# Adapt
PATH_output="/data/4_sequencing_analysis/JP238/results/tables/JP238_20min_reads_per_barcode.txt"
PATH_bam="/data/4_sequencing_analysis/JP238/results/objects/Aligned.sortedByCoord.out.bam"

#Get total reads per CB
samtools view $PATH_bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $PATH_output
samtools view possorted_genome_bam.bam | head -n 1000 | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > reads_per_barcode
