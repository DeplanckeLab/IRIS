# Cell cycle analysis
## Introduction
This folder contains the workflows needed to recreate the results shown in Figures 4 and 5, including the publication revision:

1. `1_scRNA-seq_analysis`: FUCCI scRNA-seq analysis (Figure 4)
2. `2_DeepLearning_pipeline`: Deep-learning pipeline for cell shape, nuclear, and angular models (Figure 5).
3. `3_DeepLearning_analysis`: Analysis of the original deep-learning results (Figure 5).
4. `04_Feature_analysis`: Morphology-feature RF models and revised Figure 5 analyses.

All additional necessary files are deposited in the following folders:
- `input_files`: Pregenerated files to run scripts.
- `output_files`: Comprehensive output list of scRNA-seq analysis.
- `utils`: Utility files/modules.

## Input data
### .h5 dataset of preprocessed imaging and expression data
The imaging and count matrices are available as an hdf5 dataset. This dataset was generated from the raw sequencing data as described in the Methods section, and from the raw imaging data that was preprocessed as shown in Supplementary Figure 3. The hdf5 dataset contains each cell as group with a predefined cell id (cid). For each cell the following data is available: 

- preprocces imaging data (center cropped) of the in focus plane and +/- 1 plane are available
- the cell segmentation mask
- the nucleus segmentation mask
- the gene expression data.

The preprocessed data is deposited on Zenodo under the link: [fucci_final.h5 and fucci_final_nucmask.h5](https://zenodo.org/records/17864950?preview=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MTc2NTI4MjM1NiwiZXhwIjoxNzk4Njc1MTk5fQ.eyJpZCI6IjRkNDU0MjMyLTgzYTUtNDFhZS04OGQ0LTBhY2ZjZjcxZjQ5NyIsImRhdGEiOnt9LCJyYW5kb20iOiI5MTFiOTUyMWQ0YTliMGRkZWNkZTQwMTk5NDA4NjRiYyJ9.KeNdIvWJE2Nnd6SCAcCYRsYgMW51cRx8G8j5_K7GxfpPJ3WJHVgnyaJrILvHe7bRYv8thKsp2MWDsoLN8JxWtA) and needs to be deposited in the `.\input_files`. The feature-analysis revision requires `fucci_final_nucmask.h5`, as documented in `04_Feature_analysis/README.md`.

### Additional data
Precomputed files are available in the `.\input_files` folder so all scripts should run independently. For rerunning the whole workflow paths in the scripts where indicated need to be changed, or files in the `.\input_files` folder replaced to run files. Files that are generated from scratch are available in `.\output_files` upon execution of scripts.

## Execution instructions (order)
> :warning: The original workflow requires `fucci_final.h5` or `fucci_final_nucmask.h5`; `04_Feature_analysis` requires `fucci_final_nucmask.h5`.
### Execution order
The workflows are partially dependent on each other. Precomputed intermediate results allow them to run independently; for full regeneration, use the following order:

1. `1_scRNA-seq_analysis`
2. `2_DeepLearning_pipeline`
3. `3_DeepLearning_analysis`
4. `04_Feature_analysis`

### Dependencies
Each workflow provides a `requirements.txt` or `environment.yml` for a compatible environment.
