# Pre-process data
![IRIS workflow scRNA-seq indexing workflow providing CC, WC and PC per cell, also referred to as "Full_Cell_ID"](X_docs/figures/Workflow_Full.png)
## Obtain Fastq files
- From **GEO repository**

## Index genome
- Index genome of interest and directory to genome in script  
  `IRIS_scRNA_QC_pipeline….sh` (for **STAR 2.7.9a**)

---

# Running scRNA-seq pipeline

## Prepare configuration
- Prepare `config_file_ZZZ.txt`
- Example config file: `example_config_file_example.txt`
- Follow instructions in config file

## Select sequencing platform
- Identify whether sequencing was performed on:
  - **NextSeq500** (`irisXXX`)
  - **AVITI** (`arisXXX`)
- Use:
  - `IRIS_scRNA_QC_pipeline_AVITI_R1_16bp.sh` (for AVITI)
  - `IRIS_scRNA_QC_pipeline_NextSeq_R1_16bp.sh` (for NextSeq500)

## Run pipeline

**Input:**
- Config file  
- Location of Fastq files  

**Run in background:**
```bash
nohup bash /home/userID/sequencing/1_scripts/iris_scRNAseq/1_1_scRNA-seq/1_pipeline_scRNAseq/IRIS_scRNA_QC_pipeline_NextSeq_R1_16bp.sh \
  -c /home/userID/sequencing/2_config_files/FolderName/config_file_example.txt &