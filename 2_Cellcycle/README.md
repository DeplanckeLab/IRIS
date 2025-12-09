# Cell cycle analysis
## Introduction
This folder contains the scripts that are necessary to recreate the results shown in Figure 4 and 5. The scripts are separated into 3 discrete scripts:

1. `1_scRNA-seq_analysis`: FUCCI scRNA-seq analysis (Figure 4)
2. `2_DeepLearnging_pipeline`: Deep-learning pipeline for cell shape, nuclear, and angular model (Figure 5).
3. `3_DeepLearnging_analysis`: Analysis of deep-learning results (Figure5).

All additional necessary files are deposited in the following folders:
- `input_files`: Pregenerated files to run scripts.
- `output_files`: Comprehensive output list of scRNA-seq analysis.
- `utils`: Utility files/modules.

## Input data
### .h5 dataset of preprocessed imaging and expression data
The imaging and count matrices are available as an hdf5 dataset. The dataset contains each cell as group. For each cell the following data is available: 
- preprocces imaging data (center cropped) of the in focus plane and +/- 1 plane are available
- the segmentatiom mask
- the gene expression data.
The preprocessed data is deposited on `Zenodo` under the link: [fucci_final.h5](https://zenodo.org/records/17864950?preview=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MTc2NTI4MjM1NiwiZXhwIjoxNzk4Njc1MTk5fQ.eyJpZCI6IjRkNDU0MjMyLTgzYTUtNDFhZS04OGQ0LTBhY2ZjZjcxZjQ5NyIsImRhdGEiOnt9LCJyYW5kb20iOiI5MTFiOTUyMWQ0YTliMGRkZWNkZTQwMTk5NDA4NjRiYyJ9.KeNdIvWJE2Nnd6SCAcCYRsYgMW51cRx8G8j5_K7GxfpPJ3WJHVgnyaJrILvHe7bRYv8thKsp2MWDsoLN8JxWtA) and needs to be deposited in the `.\input_files`.

### Additional data
Precomputed files are available in the `.\input_files` folder so all scripts should run independently. For rerunning the whole workflow paths in the scripts where indicated need to be changed, or files in the `.\input_files` folder replaced to run files. Files that are generated from scratch are available in `.\output_files` upon execution of scripts.

## Execution instructions (order)
> :warning: In order to run the scripts the `fucci_final.h5` needs to be deposited in `.\input_files` as described above!
### Execution order
The scripts are partially dependent on eachother. We provide all temporary results necessary that make the scripts executable individually. The only file necessary for this is `fucci_final.h5` For regeneration of all files we suggest the following order of execution:

1. `1_scRNA-seq_analysis`
2. `2_DeepLearning_pipeline`
3. `3_DeepLearning_analysis` 

### Dependencies
For each notebook a requirements.txt is provided to create an environment compatible with the notebooks.
