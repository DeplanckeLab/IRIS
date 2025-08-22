#Author: Joern Pezoldt
#Purpose: Annotate mouse and human cells in mixed datasets

##############
#Figure 2
##############




#libraries
  library(tidyverse)
  library(dplyr)
  library(plyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
  library(data.table)
  library(stringr)
  library(readr)
  library(gridExtra)
  library(pheatmap)
  library(scales)


#############################################
## Assign variables coming from config.txt ##
#############################################
work_dir <- "/home/pezoldt/updepla/projects/iris/4_sequencing_analysis/JP266/pipeline_output/human"
expID <- "JP266"
ERCC <- "no"
barcodes_wl <- "/home/pezoldt/NAS2/iris/4_genomes/barcodes/iris_genomics/20210618_scRNAseq_cellcode_whitelist_CC_7bp.txt"
barcodes_file <- "/home/pezoldt/NAS2/iris/4_genomes/barcodes/iris_genomics/20210618_scRNAseq_cellcode_BRBtoolbox_CC_7bp.txt"
sample_info_sheet <- "/home/pezoldt/NAS2/iris/0_seq_run_processing/c_SIS/sis_aris08_mapping_QC_pipeline_v2.csv"
ERCC_concentrations <- "/home/pezoldt/NAS2/iris/4_genomes/genomes/ERCC/20201212_ERCC_concentrations.txt"

#ERCC_gtf <- args[8]
ERCC_gtf <- "/home/pezoldt/NAS2/iris/4_genomes/genomes/ERCC/ERCC92.GeneAdded.mCherry.EGFP.gtf"

# Print status of script 
cat("Read arguments from command line:\n")

#############################################
############## Global Variables #############
#############################################
#Version Pipeline ID
#Design.Mapping.Compilation.Plotting
Version_Pipeline <- "1.3.5.7"

#Sequencing Run ID
seqRun_ID <- unlist(strsplit(work_dir, "_"))[length(unlist(strsplit(work_dir, "_")))]

#CC IRIS universe
CCs_available <- paste0("CC", 1:14)

#Setup rows of 96 well Plate
rows <- list(c(paste0("WC",1:12)),
             c(paste0("WC",13:24)),
             c(paste0("WC",25:36)),
             c(paste0("WC",37:48)),
             c(paste0("WC",49:60)),
             c(paste0("WC",61:72)),
             c(paste0("WC",73:84)),
             c(paste0("WC",85:96)))

#Color palettes
## CELLCODES
## green scale for supposed to be there cellcodes
cellcode_IN_colors <- c("cadetblue1","cadetblue3","cyan","cyan3","aquamarine","aquamarine3","deepskyblue",
               "deepskyblue3","darkslategray3","darkturquoise","dodgerblue","dodgerblue3","blue3")
## grey scale for not supposed to be there cellcodes
cellcode_notIN_colors <- c("gray100","gray95","gray90","gray85","gray80","gray75","gray70",
                  "gray65","gray60","gray55","gray50","gray45","gray40","gray35")
## Experimental CONDITIONS
#Color palettes
color_condition <-  c("#fff100", "#ff8c00", "#e81123", "#ec008c", "#68217a",
                      "#00188f", "#00bcf2", "#00b294", "#009e49", "#bad80a")
## PLATES
color_plate <- c("#c5f76d", "#96cf06", "#498c02", "#314d1d", "#2b800e",
                 "#63e735", "#95ef75", "#244601", "#6fa909", "#bfdea9")


#############################################
###### Load data & shape#####################
#############################################
## Load ERCC information files  -----------------
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

## Load/Create sample information with one row per cell ----------------------
# Load sample information sheet 
SIS <- read.csv(sample_info_sheet, skip = 1, header = T, stringsAsFactors = FALSE)
nrow(SIS)
#head(SIS)
#tail(SIS)
# Extract information of specific experiment from SIS
SIS_exp <- SIS[which(str_detect(SIS$Sample_ID, expID)),]
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


## Load output files from pipeline----------------
#Pivot table, UMI matrix, CCs detected, read count matrix
load(paste0(work_dir,"/analysis/objects/objects_plotting.RData"), envir = parent.frame(), verbose = FALSE)
#Dropoutwells
SIS_wellStats <- readRDS(paste0(work_dir,"/analysis/objects/SIS_wellStats.Rds"))
#Table matching experimental plate ID (column "Description") with tagmentation plate ID (PCxxx)
t_plateID_tagmentID <- unique(data.frame(experiment_plate_ID = pivot_table$Description,
                                         tagmentation_ID = pivot_table$I7_Index_ID))
#list of pivot tables split per plate
#Add intel to pivot table
##Add Codes
CellCode <- unlist(lapply(strsplit(pivot_table$Full_Cell_ID, "_"), function(x){x[[4]]}))
WellCode <- unlist(lapply(strsplit(pivot_table$Full_Cell_ID, "_"), function(x){x[[3]]}))
PlateCode <- unlist(lapply(strsplit(pivot_table$Full_Cell_ID, "_"), function(x){x[[2]]}))
ExpID <- unlist(lapply(strsplit(pivot_table$Full_Cell_ID, "_"), function(x){x[[1]]}))
pivot_table$CellCode <- CellCode
pivot_table$WellCode <- WellCode
pivot_table$PlateCode <- PlateCode

##Add experimentress/imaging ID
l_pivot_table <- split(pivot_table, pivot_table$PlateCode)

#CCs per plate
l_CC_IN_per_plate <- lapply(l_pivot_table,
                         function(x){
                           levels(as.factor(x$CC))})

##############################################
########## General QC ########################
##############################################
## Proportion of CCs -------------------
#Spot if unused CCs are overproportionally represented
#Output: PieChart with percentages per plate
## CC universe of experiment
#CCs used for experiment per plate
l_CCs_per_plate <-
  lapply(l_pivot_table,
         function(x){levels(as.factor(x$CC))})

## Setup table of detected UMIs
#Setup table
CellCode <- unlist(lapply(strsplit(UMI_stats_detected_CellCodes$Full_Cell_ID, "_"), function(x){x[[4]]}))
WellCode <- unlist(lapply(strsplit(UMI_stats_detected_CellCodes$Full_Cell_ID, "_"), function(x){x[[3]]}))
PlateCode <- unlist(lapply(strsplit(UMI_stats_detected_CellCodes$Full_Cell_ID, "_"), function(x){x[[2]]}))
ExpID <- unlist(lapply(strsplit(UMI_stats_detected_CellCodes$Full_Cell_ID, "_"), function(x){x[[1]]}))
UMI_stats_detected_CellCodes$CellCode <- CellCode
UMI_stats_detected_CellCodes$WellCode <- WellCode
UMI_stats_detected_CellCodes$PlateCode <- PlateCode
UMI_stats_detected_CellCodes$ExpID_PlateCode_WellCode <- paste0(ExpID,"_",PlateCode,"_",WellCode)
#Add experiment ID 
l_UMI_stats_detected_CellCodes <- split(UMI_stats_detected_CellCodes, UMI_stats_detected_CellCodes$PlateCode)

#Condense UMI counts per CC for each experiment_plate_ID
l_UMI_stats_wide <- list()
l_UMI_counts_per_CC <- list()
for(i in 1:length(l_UMI_stats_detected_CellCodes)){
  #i = 1
  UMI_stats_detected_CellCodes_i <- l_UMI_stats_detected_CellCodes[[i]]
  #Collect UMI per CC across all detected CC per per Well
  UMI_stats_wide_i <- UMI_stats_detected_CellCodes_i %>%
    pivot_wider(names_from = CellCode, values_from = n_UMI)
  CCs_i <- colnames(UMI_stats_wide_i)[grepl("CC", colnames(UMI_stats_wide_i))]
  UMI_stats_wide_i <- UMI_stats_wide_i %>%
    mutate_at(CCs_i, ~replace_na(.,0))
  colSums(UMI_stats_wide_i[,CCs_i])
  #Store table
  l_UMI_stats_wide[[i]] <- UMI_stats_wide_i
  #Store UMI counts ordered by number
  UMI_counts_per_CC_i <- data.frame(CCs = names(sort(colSums(UMI_stats_wide_i[,CCs_i]), decreasing = TRUE)),
                                    UMI_total = sort(colSums(UMI_stats_wide_i[,CCs_i]), decreasing = TRUE))
  #Add CCs that were not detected and give "0" UMI
  CCs_to_fill_i <- CCs_available[!(CCs_available %in% names(colSums(UMI_stats_wide_i[,CCs_i])))]
  if(length(CCs_to_fill_i) > 0){
    UMI_counts_per_CC_fill_i <- data.frame(CCs = CCs_to_fill_i, UMI_total = rep(0, length(CCs_to_fill_i)))
    UMI_counts_per_CC_i <- rbind(UMI_counts_per_CC_i, UMI_counts_per_CC_fill_i)
    UMI_counts_per_CC_i$CCs <- factor(UMI_counts_per_CC_i$CCs, levels = as.character(UMI_counts_per_CC_i$CCs))
  }else{
    UMI_counts_per_CC_i <- UMI_counts_per_CC_i
    UMI_counts_per_CC_i$CCs <- factor(UMI_counts_per_CC_i$CCs, levels = as.character(UMI_counts_per_CC_i$CCs))
  }
  rownames(UMI_counts_per_CC_i) <- NULL
  l_UMI_counts_per_CC[[i]] <- UMI_counts_per_CC_i
}
names(l_UMI_stats_wide) <- names(l_UMI_stats_detected_CellCodes)
names(l_UMI_counts_per_CC) <- names(l_UMI_stats_detected_CellCodes)

#assign colors according to CCs in system or not
l_IN_system_colors <- list()
l_notIN_system_colors <- list()
set.seed(123)
for(i in 1:length(l_UMI_counts_per_CC)){
  #i = 1
  set.seed(2)
  #list CCs in or notIn
  CC_detected_per_plate_i <- l_UMI_counts_per_CC[[i]]$CCs
  CC_IN_per_plate_i <- l_CC_IN_per_plate[[i]]
  CC_notIN_per_plate_i <- CC_detected_per_plate_i[!(CC_detected_per_plate_i %in% CC_IN_per_plate_i)]
  #grab colores
  IN_system_colors_i <- cellcode_IN_colors[sample(1:length(cellcode_IN_colors), length(CC_IN_per_plate_i))]
  names(IN_system_colors_i) <- CC_IN_per_plate_i
  l_IN_system_colors[[i]] <- IN_system_colors_i
  notIN_system_colors_i <- cellcode_notIN_colors[sample(1:length(cellcode_notIN_colors), length(CC_notIN_per_plate_i))]
  names(notIN_system_colors_i) <- CC_notIN_per_plate_i
  l_notIN_system_colors[[i]] <- notIN_system_colors_i
}
names(l_IN_system_colors) <- names(l_UMI_stats_detected_CellCodes)
names(l_notIN_system_colors) <- names(l_UMI_stats_detected_CellCodes)

#Plot and store 
#Header for Pie-Chart whether CCs with most UMIs are also the ones used
l_plot_pie_reads_CC_per_plate <- list()
for(i in 1:length(l_UMI_counts_per_CC)){
  #i = 3
  top_CCs_i <- l_UMI_counts_per_CC[[i]][1:length(l_IN_system_colors[[i]]),]$CCs
  Plate_ID_i <- names(l_UMI_counts_per_CC)[i]
  UMI_counts_per_CC_i <- l_UMI_counts_per_CC[[i]]
  if(sum(as.numeric(!(top_CCs_i %in% names(l_IN_system_colors[[i]])))) == 0){
    title <- "CCs match"
    l_plot_pie_reads_CC_per_plate[[i]] <- ggplot(UMI_counts_per_CC_i, aes(x="", y=UMI_total, fill=CCs)) +
      geom_bar(stat="identity", width=1) +
      labs(x = NULL, y = NULL, fill = NULL, 
           title = paste0(Plate_ID_i," | ",title)) +
      guides(fill = guide_legend(reverse = FALSE)) +
      scale_fill_manual(values=c(l_IN_system_colors[[i]],l_notIN_system_colors[[i]])) +
      coord_polar("y", direction = -1) +
      theme(legend.text  = element_text(size = 6),
            legend.key.size = unit(0.8, "lines"),
            plot.title = element_text(size=10),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank())
  }else{
    title <- "ATTENTION: CCs NO-match"
    l_plot_pie_reads_CC_per_plate[[i]] <- ggplot(UMI_counts_per_CC_i, aes(x="", y=UMI_total, fill=CCs)) +
      geom_bar(stat="identity", width=1) +
      labs(x = NULL, y = NULL, fill = NULL, 
           title = paste0(Plate_ID_i," | ",title)) +
      guides(fill = guide_legend(reverse = FALSE)) +
      scale_fill_manual(values=c(l_IN_system_colors[[i]],l_notIN_system_colors[[i]])) +
      coord_polar("y", direction = -1) +
      theme(legend.text  = element_text(size = 6),
            legend.key.size = unit(0.8, "lines"),
            plot.title = element_text(size=10),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank())
  }
}
l_plot_pie_reads_CC_per_plate[[1]]

## Read count table ---------------
#Make table
t_read_summary <- data.frame(Plate_ID = names(l_pivot_table),
                             n_all_reads = unlist(lapply(l_pivot_table,function(x){sum(x$n_all_reads)})),
                             n_mapped_reads = unlist(lapply(l_pivot_table,function(x){sum(x$n_mapped_read)})),
                             n_UMI = unlist(lapply(l_pivot_table,function(x){sum(x$n_UMI)})))
rownames(t_read_summary) <- NULL

## Starsolo output - Reads specs per plate --------------
nrow(SIS_wellStats)

#summarize per plate
t_SIS_wellStats <- SIS_wellStats[complete.cases(SIS_wellStats[,c("n_discarded_at_CC","n_Unmapped","n_noFeature","n_ambigFeature","n_UMI","n_saturated_UMI")]),]
nrow(t_SIS_wellStats)
test <- subset(SIS_wellStats, Description == "JP234_P3")
t_SIS_wellStats_sum <- aggregate(t_SIS_wellStats[, c("n_discarded_at_CC","n_Unmapped","n_noFeature","n_ambigFeature","n_UMI","n_saturated_UMI")], list(t_SIS_wellStats$Description), sum)

#Plot Percent Bargraph
##Pivot longer for reads classifications
t_SIS_wellStats_sum_pivot <- 
  t_SIS_wellStats_sum %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "reads_subclass",
    #names_prefix = "wk",
    values_to = "reads_subclass_n",
    values_drop_na = TRUE)
##Plot
plot_percentage_subclass_reads <- 
  ggplot(t_SIS_wellStats_sum_pivot, aes(fill=factor(reads_subclass), y=reads_subclass_n, x=Group.1)) + 
  #geom_bar(position="fill", stat="identity",
   #        width=0.4) +
  geom_bar(position="fill", stat="identity",
           width = 0.5) +
  scale_x_discrete(expand = expansion(mult = 0, 0.25*nrow(t_read_summary))) +
  scale_fill_manual(values=c('blue', 'red', 'grey', 'deepskyblue3','deepskyblue','black')) + 
  labs(title="Distribution Read Subclasses") + 
  theme(axis.text.x = element_text(angle=90),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"),
        plot.title = element_text(size=10),
        legend.title = element_blank())
plot_percentage_subclass_reads
#Plot Count Bargraphs
##Storage plots
plot_reads_per_plate <-
  ggplot(t_read_summary, aes(x=Plate_ID, y=n_all_reads)) + 
  geom_bar(stat="identity", width=.5, fill="black") + 
  scale_x_discrete(expand = expansion(mult = 0, 0.25*nrow(t_read_summary))) +
  labs(title="All sequenced reads") + 
  theme(axis.text.x = element_text(angle=90),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"),
        plot.title = element_text(size=10))
plot_mapped_per_plate <-
  ggplot(t_read_summary, aes(x=Plate_ID, y=n_mapped_reads)) + 
  geom_bar(stat="identity", width=.5, fill="black") + 
  scale_x_discrete(expand = expansion(mult = 0, 0.25*nrow(t_read_summary))) +
  labs(title="All UMIs") + 
  theme(axis.text.x = element_text(angle=90),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"),
        plot.title = element_text(size=10))
plot_UMI_per_plate <-
  ggplot(t_read_summary, aes(x=Plate_ID, y=n_UMI)) + 
  geom_bar(stat="identity", width=.5, fill="black") + 
  scale_x_discrete(expand = expansion(mult = 0, 0.25*nrow(t_read_summary))) +
  labs(title="All UMIs") + 
  theme(axis.text.x = element_text(angle=90),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"),
        plot.title = element_text(size=10))

#Store plots section
l_plot_bar_reads_per_plate_1 <- list()
l_plot_bar_reads_per_plate_2 <- list()
l_plot_bar_reads_per_plate_1[[1]] <- plot_reads_per_plate
l_plot_bar_reads_per_plate_1[[2]] <- plot_mapped_per_plate
l_plot_bar_reads_per_plate_1[[3]] <- plot_UMI_per_plate
l_plot_bar_reads_per_plate_2[[1]] <- plot_percentage_subclass_reads
names(l_plot_bar_reads_per_plate_1) <- c("reads_all_per_plate","reads_mapped_per_plate","UMIs_per_plate")
names(l_plot_bar_reads_per_plate_2) <- c("percentage_subclass_reads")

## Distribution of UMI per cell per on histogram -----------
#Identify cut-off for high quality cells
#Output:
#1) Histogram (log10(UMI)) for experiment, per plate, and per CC
#Find Plate with most CellCodes
plate_with_most_CCs <- which.max(unlist(lapply(l_IN_system_colors, length)))
names(plate_with_most_CCs) <- NULL
#Generate input data
colors_cellcodes_present <- unname(l_IN_system_colors[[plate_with_most_CCs]])
colors_plates_present <- color_plate[1:length(l_pivot_table)]
colors_conditions_present <- color_condition[1:length(split(pivot_table, pivot_table$CellType_Condition))]

#Plot
## Manual limits
l_histo_UMI_manualLIM <- list()
l_histo_UMI_manualLIM[[1]] <- ggplot(pivot_table, aes(x=n_UMI, color=CellCode, fill=CellCode)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2) +
  geom_density(alpha=0.1, linewidth = 1) +
  scale_fill_discrete(name="") +
  scale_color_manual(values=colors_cellcodes_present, name="") +
  scale_fill_manual(values=colors_cellcodes_present, name="") +
  labs(title="nUMI per CellCode",x="nUMI", y = "Density") +
  theme_classic() +
  scale_x_continuous(trans='log10',
                     limits = c(200, max(pivot_table$n_UMI)+1.0*max(pivot_table$n_UMI)))
#l_histo_UMI_manualLIM[[1]]

l_histo_UMI_manualLIM[[2]] <- ggplot(pivot_table, aes(x=n_UMI, color=Description, fill=Description)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2) +
  geom_density(alpha=0.1, linewidth = 1) +
  scale_color_manual(values=colors_plates_present, name="")+
  scale_fill_manual(values=colors_plates_present, name="")+
  labs(title="nUMI per 96well Plate", x="nUMI", y = "Density")+
  theme_classic() +
  scale_x_continuous(trans='log10',
                     limits = c(200, max(pivot_table$n_UMI)+1.0*max(pivot_table$n_UMI)))

l_histo_UMI_manualLIM[[3]] <- ggplot(pivot_table, aes(x=n_UMI, color=CellType_Condition, fill=CellType_Condition)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2) +
  geom_density(alpha=0.1, linewidth = 1) +
  scale_color_manual(values=colors_conditions_present, name="")+
  scale_fill_manual(values=colors_conditions_present, name="")+
  labs(title="nUMI per Exp. Condition",x="nUMI", y = "Density")+
  theme_classic() +
  scale_x_continuous(trans='log10',
                     limits = c(200, max(pivot_table$n_UMI)+1.0*max(pivot_table$n_UMI)))
names(l_histo_UMI_manualLIM) <- c("Histo_nUMI_CellCode","Histo_nUMI_Plate","Histo_nUMI_Condition_Celltype")               
#l_histo_UMI_manualLIM[[1]]

## Full scale
l_histo_UMI_fullscale <- list()
l_histo_UMI_fullscale[[1]] <- ggplot(pivot_table, aes(x=n_UMI, color=CellCode, fill=CellCode)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2) +
  geom_density(alpha=0.1, linewidth = 1) +
  scale_color_manual(values=colors_cellcodes_present, name="")+
  scale_fill_manual(values=colors_cellcodes_present, name="")+
  labs(title="nUMI per CellCode",x="nUMI", y = "Density")+
  theme_classic() +
  scale_x_continuous(trans='log10',
                     limits = c(1, max(pivot_table$n_UMI)+1.0*max(pivot_table$n_UMI)))

l_histo_UMI_fullscale[[2]] <- ggplot(pivot_table, aes(x=n_UMI, color=Description, fill=Description)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2) +
  geom_density(alpha=0.1, linewidth = 1) +
  scale_color_manual(values=colors_plates_present, name="")+
  scale_fill_manual(values=colors_plates_present, name="")+
  labs(title="Histo nUMI per 96well Plate", x="nUMI", y = "Density")+
  theme_classic() +
  scale_x_continuous(trans='log10',
                     limits = c(1, max(pivot_table$n_UMI)+1.0*max(pivot_table$n_UMI)))

l_histo_UMI_fullscale[[3]] <- ggplot(pivot_table, aes(x=n_UMI, color=CellType_Condition, fill=CellType_Condition)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2) +
  geom_density(alpha=0.1, linewidth = 1) +
  scale_color_manual(values=colors_conditions_present, name="")+
  scale_fill_manual(values=colors_conditions_present, name="")+
  labs(title="nUMI per Exp. Condition",x="nUMI", y = "Density")+
  theme_classic() +
  scale_x_continuous(trans='log10',
                     limits = c(1, max(pivot_table$n_UMI)+1.0*max(pivot_table$n_UMI)))
names(l_histo_UMI_fullscale) <- c("Histo_nUMI_CellCode","Histo_nUMI_Plate","Histo_nUMI_Condition_Celltype")               
#l_histo_UMI_fullscale[[1]]

##############################################
########## QC per plate ######################
##############################################
#####
## Dropout wells in plates ----------
#####
#Output: Heatmap color coded for at mapping, n(UMI), match of expected CCs

#split into plates
l_SIS_wellStats <- split(SIS_wellStats, SIS_wellStats$Description)
#Storage
l_SIS_wellStats_mapped <- list()
l_SIS_wellStats_UMI_detected <- list()
l_SIS_wellStats_CC_state <- list()

for(i in 1:length(l_SIS_wellStats)){
  #i = 1
  SIS_wellStats_i <- l_SIS_wellStats[[i]]
  l_rows_mapped <- list()
  l_rows_UMI_detected <- list()
  l_rows_CC_state <- list()
  for(j in 1:length(rows)){
    #j = 1
    rows_j <- rows[[j]]
    SIS_wellStats_i_j <- subset(SIS_wellStats_i, I5_Index_ID %in% rows_j)
    l_rows_mapped[[j]] <-  unlist(SIS_wellStats_i_j$Well_Star_mapped)
    l_rows_UMI_detected[[j]] <-  unlist(SIS_wellStats_i_j$UMI_detected)
    l_rows_CC_state[[j]] <- unlist(SIS_wellStats_i_j$CC_state)
  }
  #Make tables
  #Mapped
  rows_mapped <- do.call(rbind, l_rows_mapped)
  colnames(rows_mapped) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
  rownames(rows_mapped) <- c("A","B","C","D","E","F","G","H")
  l_SIS_wellStats_mapped[[i]] <- rows_mapped
  #UMI detected
  rows_UMI_detected <- do.call(rbind, l_rows_UMI_detected)
  colnames(rows_UMI_detected) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
  rownames(rows_UMI_detected) <- c("A","B","C","D","E","F","G","H")
  l_SIS_wellStats_UMI_detected[[i]] <- rows_UMI_detected
  #Expected CellCodes
  rows_CC_expected <- do.call(rbind, l_rows_CC_state)
  colnames(rows_CC_expected) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
  rownames(rows_CC_expected) <- c("A","B","C","D","E","F","G","H")
  l_SIS_wellStats_CC_state[[i]] <- rows_CC_expected
  #
}
#Name Plates
names(l_SIS_wellStats_mapped) <- names(l_SIS_wellStats)
names(l_SIS_wellStats_UMI_detected) <- names(l_SIS_wellStats)
names(l_SIS_wellStats_CC_state) <- names(l_SIS_wellStats)

# Plot Binary heatmap
## STAR mapping-------------------------
l_plot_96_STARmapping <- list()
for(i in 1:length(l_SIS_wellStats_mapped)){
  #print(i)
  #i = 3
  plate_id_i <- names(l_SIS_wellStats_mapped)[i]
  SIS_wellStats_i <- as.data.frame(l_SIS_wellStats_mapped[[i]])
  SIS_wellStats_i <- replace(SIS_wellStats_i, SIS_wellStats_i == "mapped", 1)
  SIS_wellStats_i <- replace(SIS_wellStats_i, SIS_wellStats_i == "not_mapped", 0)
  SIS_wellStats_i <- mutate_all(SIS_wellStats_i, function(x) as.numeric(as.character(x)))
  good_or_bad <- sum(SIS_wellStats_i)
  #print(good_or_bad)
  if(good_or_bad < 96){
    x <- 
      pheatmap(SIS_wellStats_i,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              color = c("blue", "yellow"),
              legend = FALSE,
              main = paste0(plate_id_i, " | STAR failed on BLUE wells"))
    pheatmap_i <- x[[4]]
    }
  if(good_or_bad == 96){
    SIS_wellStats_i[4,4] <- 0
    SIS_wellStats_i[5,5] <- 0
    SIS_wellStats_i[6,6] <- 0
    SIS_wellStats_i[5,7] <- 0
    SIS_wellStats_i[4,8] <- 0
    SIS_wellStats_i[3,9] <- 0
    SIS_wellStats_i[2,10] <- 0
    x <- pheatmap(SIS_wellStats_i,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = c("blue", "yellow"),
             legend = FALSE,
             main = paste0(plate_id_i, " | STAR ran on all Wells"))
    pheatmap_i <- x[[4]]
  }
  l_plot_96_STARmapping[[i]] <- pheatmap_i
}
names(l_plot_96_STARmapping) <- names(l_SIS_wellStats_mapped)
#l_plot_96_STARmapping[[1]]
## UMI Detected------------------
l_plot_96_UMI_detected <- list()
for(i in 1:length(l_SIS_wellStats_UMI_detected)){
  #print(i)
  #i = 3
  plate_id_i <- names(l_SIS_wellStats_UMI_detected)[i]
  SIS_wellStats_i <- as.data.frame(l_SIS_wellStats_UMI_detected[[i]])
  SIS_wellStats_i <- replace(SIS_wellStats_i, SIS_wellStats_i == "yes", 1)
  SIS_wellStats_i <- replace(SIS_wellStats_i, SIS_wellStats_i == "no", 0)
  SIS_wellStats_i <- mutate_all(SIS_wellStats_i, function(x) as.numeric(as.character(x)))
  good_or_bad <- sum(SIS_wellStats_i)
  #print(good_or_bad)
  if(good_or_bad < 96){
    x <- pheatmap(SIS_wellStats_i,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = c("deepskyblue2", "orange"),
                           legend = FALSE,
                           main = paste0(plate_id_i, " | No UMI on BLUE wells"))
    pheatmap_i <- x[[4]]
  }
  if(good_or_bad == 96){
    SIS_wellStats_i[4,4] <- 0
    SIS_wellStats_i[5,5] <- 0
    SIS_wellStats_i[6,6] <- 0
    SIS_wellStats_i[5,7] <- 0
    SIS_wellStats_i[4,8] <- 0
    SIS_wellStats_i[3,9] <- 0
    SIS_wellStats_i[2,10] <- 0
    x <- pheatmap(SIS_wellStats_i,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = c("deepskyblue2", "orange"),
                           legend = FALSE,
                           main = paste0(plate_id_i, " | UMI for all Wells"))
    pheatmap_i <- x[[4]]
  }
  l_plot_96_UMI_detected[[i]] <- pheatmap_i
}
names(l_plot_96_UMI_detected) <- names(l_SIS_wellStats_UMI_detected)
## CC matched------------------
l_plot_96_CC_matched <- list()
for(i in 1:length(l_SIS_wellStats_CC_state)){
  #print(i)
  #i = 2
  plate_id_i <- names(l_SIS_wellStats_CC_state)[i]
  SIS_wellStats_i <- as.data.frame(l_SIS_wellStats_CC_state[[i]])
  SIS_wellStats_i <- replace(SIS_wellStats_i, SIS_wellStats_i == "Expected_CCs", 1)
  SIS_wellStats_i <- replace(SIS_wellStats_i, SIS_wellStats_i == "notMapped_noUMI", 0)
  SIS_wellStats_i <- mutate_all(SIS_wellStats_i, function(x) as.numeric(as.character(x)))
  good_or_bad <- sum(SIS_wellStats_i)
  #print(good_or_bad)
  if(good_or_bad < 96){
    x <- pheatmap(SIS_wellStats_i,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = c("deepskyblue2", "orange"),
                           legend = FALSE,
                           main = paste0(plate_id_i, " | No UMI on BLUE wells"))
    pheatmap_i <- x[[4]]
  }
  if(good_or_bad == 96){
    SIS_wellStats_i[4,4] <- 0
    SIS_wellStats_i[5,5] <- 0
    SIS_wellStats_i[6,6] <- 0
    SIS_wellStats_i[5,7] <- 0
    SIS_wellStats_i[4,8] <- 0
    SIS_wellStats_i[3,9] <- 0
    SIS_wellStats_i[2,10] <- 0
    x <- pheatmap(SIS_wellStats_i,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = c("deepskyblue2", "orange"),
                           legend = FALSE,
                           main = paste0(plate_id_i, " | UMI for all Wells"))
    pheatmap_i <- x[[4]]
  }
  l_plot_96_CC_matched[[i]] <- pheatmap_i
}
names(l_plot_96_CC_matched) <- names(l_SIS_wellStats_CC_state)

#####
#Numeric evaluation
#####
l_pivot_table <- split(pivot_table, pivot_table$PlateCode)
#compile intel per plate
#Storage
l_per_well_allReads <- list()
l_per_well_mappedReads <- list()
l_per_well_UMI <- list()

for(i in 1:length(l_pivot_table)){
  #i = 2
  pivot_wellStats_i <- l_pivot_table[[i]]
  l_rows_allReads <- list()
  l_rows_mappedReads <- list()
  l_rows_UMI <- list()
  for(j in 1:length(rows)){
    #j = 2
    rows_j <- rows[[j]]
    pivot_wellStats_i_j <- subset(pivot_wellStats_i, WellCode %in% rows_j)
    #extract all read data
    all_reads_j <- log10(unlist(lapply(split(pivot_wellStats_i_j,
                                       pivot_wellStats_i_j$I5_Index_ID),
                                 function(x){sum(x$n_all_reads)})))
    #fill empty slots
    if(length(all_reads_j) < length(rows_j)){
      rows_missing_wc <- rows_j[!(rows_j %in% names(all_reads_j))]
      rows_missing <- rep(0, length(rows_missing_wc))
      names(rows_missing) <- rows_missing_wc
      all_reads_append_j <- c(all_reads_j,rows_missing)
      all_reads_j <- all_reads_append_j[sort(names(all_reads_append_j))]
    }
    l_rows_allReads[[j]] <- all_reads_j
    
    #extract all mapped data
    mapped_reads_j <- log10(unlist(lapply(split(pivot_wellStats_i_j,
                                       pivot_wellStats_i_j$I5_Index_ID),
                                 function(x){sum(x$n_mapped_reads)})))
    #fill empty slots
    if(length(mapped_reads_j) < length(rows_j)){
      rows_missing_wc <- rows_j[!(rows_j %in% names(mapped_reads_j))]
      rows_missing <- rep(0, length(rows_missing_wc))
      names(rows_missing) <- rows_missing_wc
      mapped_reads_append_j <- c(mapped_reads_j,rows_missing)
      mapped_reads_j <- mapped_reads_append_j[sort(names(mapped_reads_append_j))]
    }
    l_rows_mappedReads[[j]] <- mapped_reads_j
    
    #extract all UMI data
    UMI_reads_j <- log10(unlist(lapply(split(pivot_wellStats_i_j,
                                          pivot_wellStats_i_j$I5_Index_ID),
                                    function(x){sum(x$n_UMI)})))
    #fill empty slots
    if(length(UMI_reads_j) < length(rows_j)){
      rows_missing_wc <- rows_j[!(rows_j %in% names(UMI_reads_j))]
      rows_missing <- rep(0, length(rows_missing_wc))
      names(rows_missing) <- rows_missing_wc
      UMI_reads_append_j <- c(UMI_reads_j,rows_missing)
      UMI_reads_j <- UMI_reads_append_j[sort(names(UMI_reads_append_j))]
    }
    l_rows_UMI[[j]] <- UMI_reads_j
  }
  #Make tables
  #Mapped
  rows_allReads <- do.call(rbind, l_rows_allReads)
  colnames(rows_allReads) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
  rownames(rows_allReads) <- c("A","B","C","D","E","F","G","H")
  l_per_well_allReads[[i]] <- rows_allReads
  #UMI detected
  rows_mappedReads <- do.call(rbind, l_rows_mappedReads)
  colnames(rows_mappedReads) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
  rownames(rows_mappedReads) <- c("A","B","C","D","E","F","G","H")
  l_per_well_mappedReads[[i]] <- rows_mappedReads
  #Expected CellCodes
  rows_UMI <- do.call(rbind, l_rows_UMI)
  colnames(rows_UMI) <- c("1","2","3","4","5","6","7","8","9","10","11","12")
  rownames(rows_UMI) <- c("A","B","C","D","E","F","G","H")
  l_per_well_UMI[[i]] <- rows_UMI
  #
}
#Name Plates
names(l_per_well_allReads) <- names(l_SIS_wellStats)
names(l_per_well_mappedReads) <- names(l_SIS_wellStats)
names(l_per_well_UMI) <- names(l_SIS_wellStats)

## Sequencing depths ----------------
#Output: Heatmap of mapped reads per well
l_plot_96_Reads_all_perWell <- list()
l_plot_96_Reads_mapped_perWell <- list()
l_plot_96_Reads_UMI_perWell <- list()

for(i in 1:length(l_per_well_allReads)){
  #print(i)
  #i = 3
  plate_id_i <- names(l_per_well_allReads)[i]
  #Load Tables per plate
  allReads_i <- as.data.frame(l_per_well_allReads[[i]])
  allReads_i[allReads_i == -Inf] <- 1
  mappedReads_i <- as.data.frame(l_per_well_mappedReads[[i]])
  mappedReads_i[mappedReads_i == -Inf] <- 1
  UMI_i <- as.data.frame(l_per_well_UMI[[i]])
  UMI_i[UMI_i == -Inf] <- 1
  #plot & store
  x_pheatmap_allReads_i <- pheatmap(allReads_i,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         colorRampPalette(c("#040481","deepskyblue","white","#CB1B1B"), space="rgb")(40),
                         font_size = 1,
                         main = paste0(plate_id_i, " | all Reads"))
  x_pheatmap_mappedReads_i <- pheatmap(mappedReads_i,
                                  cluster_rows = FALSE,
                                  cluster_cols = FALSE,
                                  colorRampPalette(c("#040481","deepskyblue","white","#CB1B1B"), space="rgb")(40),
                                  font_size = 1,
                                  main = paste0(plate_id_i, " | mapped Reads"))
  x_pheatmap_UMI_i <- pheatmap(UMI_i,
                                     cluster_rows = FALSE,
                                     cluster_cols = FALSE,
                                     colorRampPalette(c("#040481","deepskyblue","white","#CB1B1B"), space="rgb")(40),
                                     font_size = 1,
                                     main = paste0(plate_id_i, " | UMI"))
  #store
  l_plot_96_Reads_all_perWell[[i]] <- x_pheatmap_allReads_i[[4]]
  l_plot_96_Reads_mapped_perWell[[i]] <- x_pheatmap_mappedReads_i[[4]]
  l_plot_96_Reads_UMI_perWell[[i]] <- x_pheatmap_UMI_i[[4]]
}

l_plot_96_Reads_UMI_perWell[[1]]

names(l_plot_96_Reads_all_perWell) <- names(l_per_well_allReads)
names(l_plot_96_Reads_mapped_perWell) <- names(l_per_well_allReads)
names(l_plot_96_Reads_UMI_perWell) <- names(l_per_well_allReads)

