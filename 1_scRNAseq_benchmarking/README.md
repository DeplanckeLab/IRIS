# scRNA-seq benchmarking
Here we compare the sensitivity of the IRIS scRNA-seq biochemistry was compared against DisCO (Drop-seq derivative) and 10X Genomics V3 on human HEK293T.

## [Preprocessing](./preprocessing)
We added libraries for two competing technologies: 10x and DISCO. For each of them, we've run the following steps:
- **STAR solo** mapping. The foundation to obtain all QC parameters per cell.
- Extraction of number of sequenced reads and mapped reads per cell from BAM.
- Subsampling from BAM.
- Compilation of subsampling for visualization.

## [Visualization](./visualization)
Here are the R-based scripts to generate plots from **[preprocessing](./preprocessing)** for 10x and DISCO libs, or from [0_scRNAseq_preprocessing](../0_scRNAseq_preprocessing) for IRIS libs
- **CELLLINEperCC** encapsulation of murine NIH/3T3 with one Cell Code and human HEK293T with a different second CellCode, with thereby one mouse and one human cell per well.
- **CELLSperHOUR** continuous profiling of HEK293T cells for 1h.
- **CELLSperPLATE** continuous profiling of HEK293T cells with 13 different Cell Codes yielding *droplet consortia* of 13 cells per well.
- **WellCoding** species-mixing experiments with NIH/3T3 and HEK293T mixed at equal ratios to optimize Well Coding primer exclusion during cDNA generation.
- **Downsampling** visualization of downsampled IRIS, 10X V3 and Drop-seq/DisCO data on HEK293T.
