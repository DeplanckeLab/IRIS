#!/bin/bash

# Author: Joern Pezoldt
# Date: Summer 2024
# Description: Subsample per cell

# Path to conda environments
source /opt/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate sc_rna_seq_v2

#Dir & Paths
dir_STARIndex="/data/genomes/human/GRCh38.90_GFP_mCherry_ERCC/STAR_index_2_7_9a"
dir_bam="/data/4_sequencing_analysis/JP238/results/objects"
dir_subsampling=$dir_bam"/subsampling"
dir_subsampling_bam=$dir_subsampling"/bam"
dir_subsampling_merged=$dir_subsampling"/merged"
#path_barcodes="/data/4_sequencing_analysis/JP142/results/hg/Solo.out/GeneFull/filtered/barcodes.tsv.gz"
name_bam="Aligned.sortedByCoord.out.bam"

#set subsample sizes
subsample_sizes=(200000 100000 50000 20000 15000 10000 5000)
#(200000 100000 50000 20000 15000 10000 5000)
#Make dirs
mkdir $dir_subsampling
mkdir $dir_subsampling_bam
mkdir $dir_subsampling_merged

#Make test bam

#Load files
## Barcode files
## Note: need a spacer line as header and at the end

# Kick out unmapped reads 
samtools view -@ 3 -b -F 4 $dir_bam/$name_bam > $dir_bam"/mapped.bam"
# Save the header lines
samtools view -H $dir_bam"/mapped.bam" > $dir_subsampling"/BAM_header"

### Solution Caro
# Read barcodes from file 
{
read
while read -r Barcode
do
echo $Barcode
#Filter bam of mapped reads for single barcode
#Barcode="AATCACGTCCGGACGT"
samtools view -@ 3 $dir_bam"/mapped.bam" | grep -F $Barcode > $dir_subsampling_bam/"$Barcode"_subsample.bam
#Combine header with body
cat $dir_subsampling"/BAM_header" $dir_subsampling_bam/"$Barcode"_subsample.bam > $dir_subsampling_bam/"$Barcode"_subsample.sam
#Convert SAM to BAM file 
samtools view -b $dir_subsampling_bam/"$Barcode"_subsample.sam > $dir_subsampling_bam/"$Barcode"_subsample.bam
# Create subsample files for each barcode 
#chmod 777 -R $dir_subsampling_bam
for i in "${subsample_sizes[@]}"; do
	echo $i
	# Calculate fractions of reads to match the desired subsample size
	total=$(samtools view -c $dir_subsampling_bam"/"$Barcode"_subsample.bam")
	echo $total
	frac=$(awk -v var1=$i -v var2=$total 'BEGIN { print  ( var1 / var2 ) }')
	echo $frac
	if [[ $frac > 1 ]]; then
	cp $dir_subsampling_bam"/"$Barcode"_subsample.bam" $dir_subsampling_bam"/"$Barcode"_subsample_"$i".bam"
	else
	samtools view -b -s $frac $dir_subsampling_bam"/"$Barcode"_subsample.bam" -o $dir_subsampling_bam"/"$Barcode"_subsample_"$i".bam"
	samtools view -F 0x04 -c $dir_subsampling_bam"/"$Barcode"_subsample_"$i".bam"
	fi
done
done
} < $dir_subsampling"/barcodes1000UMI.tsv"

# Merge BAM of same sub sample size and copy to well dir foler 
for i in "${subsample_sizes[@]}"; do
	mkdir $dir_subsampling_merged/merged_"$i"
	samtools merge $dir_subsampling_merged/merged_"$i"/merged_"$i".bam $dir_subsampling_bam/*_"$i".bam
done

# Delete temporary subfolder
rm -r $dir_subsampling_bam 

for i in $dir_subsampling_merged/merged_*; do
echo $i
echo $i/*.bam
# Create UMI matrix per subsample BAM
STAR --genomeDir $dir_STARIndex \
	 --genomeLoad LoadAndKeep \
	 --runThreadN 12 \
	 --outFileNamePrefix $i/ \
	 --readFilesIn $i/*.bam \
	 --soloType CB_UMI_Simple \
	 --soloCBstart 1 \
	 --soloCBlen 12 \
	 --soloUMIstart 13 \
	 --soloUMIlen 8 \
	 --soloCBmatchWLtype 1MM \
	 --soloCBwhitelist $dir_subsampling"/barcodes1000UMI.tsv" \
	 --soloUMIdedup 1MM_All \
	 --readFilesType SAM SE \
	 --readFilesCommand samtools view -F 0x100 \
	 --soloInputSAMattrBarcodeSeq CR UR

rm -r $i/*.bam
rm -r $i/*.sam
rm -r $i/Solo.out/Gene/raw
done



