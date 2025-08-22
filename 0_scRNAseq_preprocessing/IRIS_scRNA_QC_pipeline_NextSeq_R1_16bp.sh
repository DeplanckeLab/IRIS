#!/bin/bash

# Author: Caro
# Date: Summer 2021
# Description: This script extracts quality control data from fastq files of IRIS experiments   

# Path to conda environments
source /opt/miniconda3/miniconda3/etc/profile.d/conda.sh

# Access conda environment for mapping
conda activate sc_rna_seq_v2

# Load external packages
BRBseq_toolbox="/opt/miniconda3/miniconda3/envs/sc_rna_seq_v2/bin/BRBseqTools-1.6.jar"
TFseq_tools="/bin/TFseqTools-1.0.jar"

# Define colours used for messages
RED='\033[0;31m' # Red
NC='\033[0m' # No Color

# Help statement for the script 
showHelp() {
cat << EOF # Cat should stop reading when EOF is detected; equivalent to echo but neater when printing out.
Usage: ./load_config -c <config-file> [-h]
Load config file.

-h, -help,         --help                Display help

-c, -config-file,  --config-file         Set location of configuration file

EOF
}


# $@ : all command line parameters are passed to the script
# -o : for short options (like -c)
# -l : for long options with double dash (like --config-file); comma separates different long options
# -a : for long options with single dash (like -config-file)
# :  : mandatory variable (:: optional)

options=$(getopt -l "help,config-file:" -o "hc:" -a -- "$@")

# set --: interprets everything following as argument also when it starts with '-' (good practice to use set -- when 
# you do not know what follows)

eval set -- "$options"

while true # Run until any stop function (exit/break) is called 
do
case $1 in 
-h|--help)
    showHelp
    exit 0
    ;;
-c|--config-file)
    shift # Jump to the variable behind "-c"
    export config_file=$1 # Set config_file as variable in the shell
	break
    ;;
*) # Stop if something unexpected (* matches all patterns) is entered
	echo -e "${RED}Error: Define the location of the config file.${NC}"
	exit 
	;;
esac # end case loop
shift 
done 

## Config file control

# Check if config file exists 
if ! [[ -e $config_file ]]; then 
	echo -e "${RED}Error: "$config_file" does not exist${NC}"
	exit
fi

# Load (variables from) config file
source $config_file

# Define variables from config file 
user=$userID
seqRun=$seqRunID
expID=$experimentID
fastq_dir=$fastq_directory
scripts_path=$scripts_path
cellcodes_wl=$cellcodes_whitelist
cellcodes_brb=$cellcodes_BRBtoolbox
SIS=$sample_info_sheet
genome=$genome
index_version=$STAR_index_version
ERCC=$ERCC
ERCC_conc=$ERCC_concentration
ERCC_gtf=$ERCC_gtf
subsamp=$subsampling
tfseq=$TFseq
fastq_dir_tfseq=$TFseq_fastq_directory
vector_genome=$TFseq_vector_genome
tf_barcodes=$TFseq_barcodes


# Check if any variable was not defined in the config file 
# Array of variables coming from config file that must be assigned 
config_variables=(userID seqRunID experimentID fastq_directory scripts_path cellcodes_whitelist cellcodes_BRBtoolbox 
sample_info_sheet genome STAR_index_version ERCC ERCC_concentration ERCC_gtf subsampling TFseq)

if [[ $tfseq == "yes" ]]; then 
	config_variables+=(TFseq_fastq_directory TFseq_vector_genome TFseq_barcodes)
fi 

# Print name of variables if not defined in config file and exit 
for var_name in "${config_variables[@]}"; 
do
        [ -z "${!var_name}" ] && echo -e "${RED}Error: $var_name is not defined in the config file.${NC}" && var_unset=true
done
        [ -n "$var_unset" ] && exit 1

# Check if directories and files defined in the config file exist
directories=($fastq_directory $scripts_path)
files=($cellcodes_whitelist $cellcodes_BRBtoolbox $sample_info_sheet $ERCC_concentration $ERCC_gtf)

# Extend variable list if TFseq is used
if [[ $tfseq == "yes" ]]; then 
	files+=($TFseq_vector_genome $TFseq_barcodes) 
	directories+=($TFseq_fastq_directory)

fi 

# Check if files exist and exit of not 
for file in "${files[@]}"; do

	if ! [[ -e $file ]]; then 
		echo -e "${RED}Error: "$file" does not exist${NC}"
		exit
	fi
done 

# Check if directories exist and exit if not 
for dir in "${directories[@]}"; do

	if ! [[ -d $dir ]]; then 
		echo -e "${RED}Error: "$dir" does not exist${NC}"
		exit
	fi
done 


# Check if defined genome in config file is valid (i.e. part of available genomes )
if ! [[ "$genome" =~ ^(mouse|mouse_addins|human_mouse|human_ERCC_GFP_mCherry|human_mouse_droso)$ ]]; then
    echo -e "${RED}Error: "$genome" is not part of available genomes${NC}"
    echo "Available genomes: mouse, human_mouse, human_ERCC_GFP_mCherry, human_mouse_droso"
    exit
fi


## Create directory structure 

# Define experiment folder 
exp_folder=""$user"_"$expID"_"$seqRun""

# Define directory paths 
home_dir=/"home"/$userID
global_seq_dir=$home_dir/"sequencing"
final_output=$global_seq_dir/"experiments"/$exp_folder
exp_dir=$home_dir/"experiments"/$exp_folder
analysis_dir=$exp_dir/"analysis"
data_dir=$exp_dir/"data"
temp_dir=$exp_dir/"temp"
rcm_dir=$data_dir/"count_matrices"
tfseq_dir=$data_dir/"TFseq"
plots_dir=$analysis_dir/"plots"
objects_dir=$analysis_dir/"objects"
genomes_dir="/data/genomes"

# Create array of subdirectories 
subdirectories=($analysis_dir $data_dir $temp_dir $rcm_dir $tfseq_dir $plots_dir $objects_dir)

# Create "/experiments" if not existing
if [[ ! -e $home_dir ]]; then 
	mkdir -v $home_dir
fi 

# Promp user input if final output folder already exists
if [[ -e $final_output ]]; then 
	echo -e "${RED}$final_output already exists. By continuing all existing files in $final_output will be overwritten.${NC}"
	# Promp user input 
	while true; do
	    read -p "Do you wish to continue? (y/n)" yn
	    case $yn in
		# Whole directory structure (+ all files) will be removed
		[yYyes]*) rm -rfv $final_output
			  break;;
		# Exit script
		[nNno]*) exit;;
		* ) echo "Please answer y for yes or n for no.";;
	    esac
	done
fi

# Create directory structure if experiment folder does not exist yet
if [[ ! -e $exp_dir ]]; then

	mkdir -v $exp_dir

	for dir in "${subdirectories[@]}"; do
		mkdir $dir 
	done

# Promp user input if experiment folder already exists
else
    echo -e "${RED}$exp_dir exists. By continuing all existing files in $exp_dir will be overwritten.${NC}"
	# Promp user input 
while true; do
    read -p "Do you wish to continue? (y/n)" yn
    case $yn in
		# Whole directory structure (+ all files) will be overwritten with empty folders
        [yYyes]*)	rm -rfv $exp_dir; mkdir -v $exp_dir

				for dir in "${subdirectories[@]}"; do
		
					mkdir -v $dir 

				done
				break;;
		# Exit script
        [nNno]*) exit;;
	* ) echo "Please answer y for yes or n for no.";;
    esac
done
fi

## set directory of genome

# Once there was a human genome ... but than J rn deleted it ... playing god can be so easy :D 
# Mouse genome
if [[ "$genome" == "mouse" ]]
then
	genomedir="$genomes_dir/mouse/m38r9";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/mm_m38r91.gtf;

# Mouse + artifical addin genes
elif [[ "$genome" == "mouse_addins" ]]
then
	genomedir="$genomes_dir/mouse_addins/m38r9";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/mm_m38r91_addins.gtf;

# Human + Mouse genome
elif [[ "$genome" == "human_mouse" ]]
then
	genomedir="$genomes_dir/human_mouse/GRCh38.90_m38r9";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/Mus_musculus.GRCm38.Homo_sapiens.GRCh38.93.mixed.gtf;

# Human + ERCC + GFP + mCherry genome
elif [[ "$genome" == "human_ERCC_GFP_mCherry" ]]
then
	genomedir="$genomes_dir/human/GRCh38.90_GFP_mCherry_ERCC";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/Homo_sapiens.GRCh38.90.dna.primary_assembly_ERCC_GFP_mCherry.gtf;

# Human + Mouse + Drosophila genome
elif [[ "$genome" == "human_mouse_droso" ]]
then
	genomedir="$genomes_dir/human_mouse_droso/GRCh38.90_Drosor6.37FB202006_m38r9";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/GRCh38.90_Drosor6.37FB202006_m38r9.gtf;
fi

# Copy genome files to workstation and set new location
#echo "Copying "$STARIndexdir" and "$inputGTF" to /data/genomes"
#cp -r $STARIndexdir /data/genomes/
#cp -r $inputGTF /data/genomes/
#STARIndexdir=/data/genomes/STAR_index_"$index_version"
#inputGTF=/data/genomes/*.gtf


#Create file that contains name and sequence of used barcodes
if [[ "$subsamp" == "yes" ]]; then 

	# Activate conda right environment to run R script
	conda deactivate
	conda activate r_v4

	#Create file that contains name and sequence of used barcodes
	Rscript $scripts_path/IRIS_scRNA_QC_CCsubset.R $expID $cellcodes_brb $SIS $temp_dir

	# Switch back to mapping conda env 
	conda deactivate 
	conda activate sc_rna_seq_v2	
fi



# Loop through fastq folders (= wells)

for fastq in "$fastq_dir"/*.fastq.gz; 
do 
	# Create folder per well in "/count_matrices" folder and copy perspective fastq files 
	file=${fastq##*/}
	dir_well=$rcm_dir/${file%_*_*_*}
	[[ ! -d "$dir_well" ]] && mkdir -- "$dir_well"
	cp "$fastq" "$dir_well"
done

# Loop through well folders and map reads with STARsolo 

for dir_well in $rcm_dir/*;
do

	# Get read1 and read2 from each well 
	read1=$(ls -d $dir_well/*.fastq.gz | sort -V | head -n 1)
	read2=$(ls -d $dir_well/*.fastq.gz | sort -V | tail -n 1)


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

	STAR --genomeDir $STARIndexdir \
		--genomeLoad LoadAndKeep \
		--runThreadN 18 \
		--outFileNamePrefix $dir_well/ \
		--readFilesCommand zcat \
		--readFilesIn $read2 $read1 \
		--soloType CB_UMI_Simple CB_samTagOut \
		--outFilterScoreMin 30 \
		--soloCBstart 1 \
		--soloCBlen 7 \
		--soloUMIstart 8 \
		--soloUMIlen 9 \
		--soloCBwhitelist $cellcodes_wl \
		--soloCBmatchWLtype 1MM \
		--soloUMIfiltering MultiGeneUMI_CR \
		--soloUMIdedup 1MM_CR \
		--clipAdapterType CellRanger4 \
		--outSAMtype BAM SortedByCoordinate  \
		--soloFeatures Gene GeneFull SJ Velocyto \
		--soloMultiMappers Uniform PropUnique EM Rescue \
		--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
		--outSAMunmapped Within \
		--limitBAMsortRAM 5000000000

	## Generate ReadCountMatrix from STARsolo BAM files with BRBseqTool ExtractReadCountMatrix
	
	# See https://github.com/DeplanckeLab/BRB-seqTools for more information 
	# -b : [Required] Path of STAR-aligned BAM file to analyze.
	# -c : [Required] Path of Barcode/Samplename mapping file.
	# -gtf : [Required] Path of GTF file.
	# -o : Output folder
	# -chunkSize : Maximum number of reads to be stored in RAM (default = 1000000)

	#java	-jar $BRBseq_toolbox ExtractReadCountMatrix \
	#		-b $dir_well/*.bam \
	#		-c $cellcodes_brb \
	#		-gtf $inputGTF \
	#		-o $dir_well


	## Generate Counts of all Reads and mapped Reads per BAM
	#Total Reads and Mapped reads
	dir_read_count=$dir_well"/Read_counts"
	#make folder
	mkdir -v $dir_read_count
	# make mapped reads
	samtools view -b -F 4 $dir_well/*.bam > $dir_read_count"/temp_mapped.bam"
	# Get total reads per CC
	samtools view $dir_well/*.bam | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $dir_read_count"/CC_all_reads.txt"
	# Get mapped reads per CC
	samtools view $dir_read_count"/temp_mapped.bam" | grep CB:Z: | sed 's/.*CB:Z:\([ACGT]*\).*/\1/' | sort | uniq -c > $dir_read_count"/CC_mapped_reads.txt"
	# Delete BAMs
	rm $dir_read_count"/temp_mapped.bam"

	if  [[ "$subsamp" == "yes" ]]
	then 

	# Create temporal subfolder for processing
	subsamp_temp=$dir_well/"temp"
	mkdir $subsamp_temp
	
	# Kick out unmapped reads 
	
	samtools view -b -F 4 $dir_well/Aligned.sortedByCoord.out.bam > $subsamp_temp/mapped.bam

	echo $subsamp_temp
	# Save the header lines
	samtools view -H $subsamp_temp/mapped.bam > $subsamp_temp/SAM_header

	# Set subsample sizes for mapped reads
	subsample_sizes=(100000 50000 20000 15000 10000 5000)

	#Subsample from BAM by barcodes 
	
	# Read barcodes from file 
	{  
	read
	while read -r Name B1
	do  

	
	# Filter bam of mapped reads for single barcode
	samtools view $subsamp_temp/mapped.bam | grep -F "$B1" > $subsamp_temp/"$Name"_subsample.bam
	
	# Combine SAM header with filtered body 
	cat $subsamp_temp/SAM_header $subsamp_temp/"$Name"_subsample.bam > $subsamp_temp/"$Name"_subsample.sam

	# Convert SAM to BAM file 
	samtools view -b $subsamp_temp/"$Name"_subsample.sam > $subsamp_temp/"$Name"_subsample.bam

		# Create subsample files for each barcode 
		for i in "${subsample_sizes[@]}"; do
			echo $i
			
			# Calculate fractions of reads to match the desired subsample size
			total=$(samtools view -c $subsamp_temp/"$Name"_subsample.bam)
			echo $total	
			
			frac=$(awk -v var1=$i -v var2=$total 'BEGIN { print  ( var1 / var2 ) }')
			echo $frac
			
			# Copy all reads to new files if fraction is equal or greater than 1
			if [[ $frac > 1 ]]; then 
				
				cp $subsamp_temp/"$Name"_subsample.bam $subsamp_temp/"$Name"_subsample_"$i".bam
			
			# Copy fraction of reads
			else
				samtools view -bs $frac $subsamp_temp/"$Name"_subsample.bam > $subsamp_temp/"$Name"_subsample_"$i".bam
	
			fi

		done

	done
	} < "$temp_dir/barcodes_subset.txt" 


	# Merge BAM of same sub sample size and copy to well dir foler 
	for i in "${subsample_sizes[@]}"; do 
		
		# Make folders per subsample
		subsamp_folder=$dir_well/"subsample_"$i"_reads"
		mkdir $subsamp_folder
		

		samtools merge $subsamp_folder/merged_"$i".bam $subsamp_temp/*_"$i".bam
    done


	# Delete temporary subfolder
	rm -r $subsamp_temp 

	# Go through subsample folders 
	for sub in $dir_well/subsample*/; do 

	# Create UMI matrix per subsample BAM
	STAR --genomeDir $STARIndexdir \
		 --runThreadN 18 \
		 --outFileNamePrefix $sub \
		 --readFilesIn $sub/*.bam \
		 --soloType CB_UMI_Simple \
		 --soloCBstart 1 \
		 --soloCBlen 7 \
		 --soloUMIstart 8 \
		 --soloUMIlen 9 \
		 --soloCBmatchWLtype 1MM \
		 --soloCBwhitelist $cellcodes_wl \
		 --soloUMIdedup 1MM_All \
		 --readFilesType SAM SE \
		 --readFilesCommand samtools view -F 0x100 \
		 --soloInputSAMattrBarcodeSeq CR UR
	rm -r $sub/*.bam
	rm -r $sub/*.sam
	rm -r $sub/Solo.out/Gene/raw
	done

fi

done

# Remove loaded genome 
STAR --genomeDir $STARIndexdir --genomeLoad Remove

# TFseq processing 
if [[ $tfseq == "yes" ]]; then
	
	# Create folder per well in "/TFseq" folder and copy perspective fastq files 
	for fastq in "$fastq_dir_tfseq"/*.fastq.gz; do 
		file=${fastq##*/} 
		dir_well=$tfseq_dir/${file%_*_*_*}
		[[ ! -d "$dir_well" ]] && mkdir -- "$dir_well"
		cp "$fastq" "$dir_well"
	done

	# Loop through well folders and map read2 with STARsolo 
	for dir_well in $tfseq_dir/*;
	do

		# Get read1 and read2 from each well 
		read1=$(ls -d $dir_well/*.fastq.gz | sort -V | head -n 1)
		read2=$(ls -d $dir_well/*.fastq.gz | sort -V | tail -n 1)

 		# Create mapped read2 BAM file with STARsolo 
		STAR --runMode alignReads \
			--genomeLoad LoadAndKeep \
			--outSAMmapqUnique 60 \
			--runThreadN 8 \
			--genomeDir $vector_genome \
			--outFilterMultimapNmax 1 \
			--readFilesCommand zcat \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix $dir_well/ \
			--readFilesIn $read2

		## Use TFseq_tools to extract RCM and UMI count matrix for TFs
		
		# See https://github.com/DeplanckeLab/TFseqTools for more information 
		# --r1  [Required] Path of R1 FastQ file.
		# --r2  [Required] Path of R2 aligned BAM file [do not need to be sorted or indexed].
		# --tf 	[Required] File containing known TF barcodes
		# -o    Output folder [default = folder of BAM file]
		# --nu  Number of allowed difference (hamming distance) for two UMIs to be counted only once [default = 0].
		# -p 	Cell barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file [default = 'BU', i.e. barcode followed by the UMI]
		# --UMI length of the UMI 
		# --BC 	length of the barcode
		# --log name of detailed log file
		
		java -jar $TFseq_tools Counter \
			--r1 $read1 \
			--r2 $dir_well/Aligned.out.bam \
			--tf $tf_barcodes \
			--nu 1 \
			--UMI 9 \
			--BC 7 \
			--log $dir_well/log_TFseq.txt
	done

STAR --genomeDir $vector_genome --genomeLoad Remove

fi 

## Delete all files and directories not needed 
for well in $rcm_dir/*; do
	rm $well/*.fastq.gz
	#rm $well/*.bam
done

if [[ $tfseq == "yes" ]]; then
	for well in $tfseq_dir/*; do
		rm $well/*.fastq.gz
		rm $well/*.bam
	done
fi

# Delete temporary folder 
rm -r $temp_dir

# Delete genome files 
#rm -r $STARIndexdir
#rm $inputGTF

# Deactivate conda environment for mapping 
conda deactivate 

# Activate conda environment for R statistics
conda activate r_v4

# Run R script for plotting statistics
Rscript $scripts_path/IRIS_scRNA_QC_RCMtoDEG.R $exp_dir $expID $ERCC $cellcodes_wl $cellcodes_brb $SIS $ERCC_conc $ERCC_gtf
#arg1 exp_dir
#arg2 exp_ID
#arg3 ERCC
#arg4 cellcodes_wl
#arg5 ellcodes_brb
#arg6 SIS
#arg7 ERCC_conc
#arg8 ERCC_gtf

# Run R script for subsampling plots 
if  [[ "$subsamp" == "yes" ]]
then 

Rscript $scripts_path/IRIS_scRNA_QC_subsampling.R $exp_dir $expID $cellcodes_brb $SIS

fi 

conda deactivate


# Activate conda environment for R plotting
# Run R script for plotting QC
conda activate r_v4_plot

Rscript $scripts_path/IRIS_scRNA_QC_plotting.R $exp_dir $expID $ERCC $cellcodes_wl $cellcodes_brb $SIS $ERCC_conc $ERCC_gtf

conda deactivate

# Copy folder to sequencing folder
mv $exp_dir $global_seq_dir/"experiments"/

# Delete folder
rm -r $exp_dir
