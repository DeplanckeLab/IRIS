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
### Preprocessed imaging and expression data
The imaging and count matrices are available as an hdf5 dataset. The dataset contains each cell as group. For each cell the following data is available: 
- preprocces imaging data (center cropped) of the in focus plane and +/- 1 plane are available
- the segmentatiom mask
- the gene expression data.
The preprocessed data is deposited on `Zenodo` under the link: `...` and needs to be deposited in the `.\input_files`.

## Execution extructions
> :warning: In order to run the scripts the `fucci_final.h5` needs to be deposited in `.\input_files` as described above!
### Execution order
The scripts are partially dependent on eachother. We provide all temporary results necessary that make the scripts executable individually. The only file necessary for this is X.h5. For regeneration of all files we suggest the following order of execution:

1. FUCCI_scRNA-seq can be fully executed except for the last cell. Produces the following files:
2. Execution of the machine learning pipeline.
3. Analysis of deep-learning 

### Dependencies
For each notebook a requirements.txt is provided to create an environment compatible with the notebooks.
