#!/bin/bash

# Author: Caro
# Date: Summer 2021
# Description: This script extracts quality control data from fastq files of IRIS experiments   

# Path to conda environments
source /opt/miniconda3/miniconda3/etc/profile.d/conda.sh

# Access conda environment for mapping
conda activate sc_rna_seq_v2

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
SIS=$sample_info_sheet
genome=$genome
index_version=$STAR_index_version

# Check if any variable was not defined in the config file 
# Array of variables coming from config file that must be assigned 
config_variables=(userID seqRunID experimentID fastq_directory scripts_path cellcodes_whitelist  
sample_info_sheet genome STAR_index_version)

# Print name of variables if not defined in config file and exit 
for var_name in "${config_variables[@]}"; 
do
        [ -z "${!var_name}" ] && echo -e "${RED}Error: $var_name is not defined in the config file.${NC}" && var_unset=true
done
        [ -n "$var_unset" ] && exit 1

# Check if directories and files defined in the config file exist
directories=($fastq_directory $scripts_path)
files=($cellcodes_whitelist $sample_info_sheet )

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
if ! [[ "$genome" =~ ^(mouse|mouse_addins|human_mouse|human_ERCC_GFP_mCherry|human_mouse_droso|hamster|hamster_droso)$ ]]; then
    echo -e "${RED}Error: "$genome" is not part of available genomes${NC}"
    echo "Available genomes: mouse, human_mouse, human_ERCC_GFP_mCherry, human_mouse_droso, hamster, hamster_droso"
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
plots_dir=$analysis_dir/"plots"
objects_dir=$analysis_dir/"objects"
objects_mtx_dir=$objects_dir/"mtx"
genomes_dir="/data/genomes"

# Create array of subdirectories 
subdirectories=($analysis_dir $data_dir $temp_dir $rcm_dir $plots_dir $objects_dir $objects_mtx_dir)

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
	
# Hamster genome
elif [[ "$genome" == "hamster" ]]
then
	genomedir="$genomes_dir/hamster/CriGri-PICRH-1.0_r113_adds";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/Cricetulus_griseus_picr.CriGri-PICRH-1.0.113_biologics.gtf;
	
# Hamster + Drosophila genome
elif [[ "$genome" == "hamster_droso" ]]
then
	genomedir="$genomes_dir/hamster_droso/CriGri_Drosor6.37FB202006";
	STARIndexdir=$genomedir/STAR_index_"$index_version";
	inputGTF=$genomedir/gtf/CriGri_Drosor6.37FB202006.gtf;
fi


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
done

# Remove loaded genome 
STAR --genomeDir $STARIndexdir --genomeLoad Remove


## Delete all files and directories not needed 
for well in $rcm_dir/*; do
	rm $well/*.fastq.gz
	#rm $well/*.bam
done

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
Rscript $scripts_path/IRIS_scRNA_QC_RCMtoDEG.R $exp_dir $expID $cellcodes_wl $SIS
#arg1 exp_dir
#arg2 exp_ID
#arg3 cellcodes_wl
#arg4 SIS

# Copy folder to sequencing folder
mv $exp_dir $global_seq_dir/"experiments"/

# Delete folder
rm -r $exp_dir
