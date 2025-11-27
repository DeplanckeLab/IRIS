# scRNA-seq benchmarking
- Identifying molecular signatures underlying the T cell architectures in healthy human naive CD8+ T cells.

## Preprocessing
- Output of **0_scRNAseq_preprocessing**.
- UMI count matrices and fastq available on GEO upon publication

## Figure_6
- Jupyter Notebook for all plots in Figure 6 and Supplement figure.

## Visualization
- **PBMC** profiling of four healthy donors.
    - **image_annotation** contains per cell image annotation for nuclear invagination for each donor.
    - "PBMC.R" for all plots of PBMCs in Figure 6 combining manually annotated nuclear morphologies with scRNA-seq of four donors.
- **Tcells_CD8/CD8_GeneOntology/CD8_GeneOntology.R** Gene ontology analysis of differentially upregulated genes in 'stripy' naive T cell morphology as compared to 'conventional' morphology. Differential gene expression analysis was performed in scanpy. See "Figure_6"


