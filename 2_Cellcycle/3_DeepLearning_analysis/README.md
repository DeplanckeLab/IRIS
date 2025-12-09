# Deep learning analysis
## Introduction
In the previous step '' three models were trained and tested:
- Nucleus CNN model
- Cell shape CNN model
- Angular score MLP
The produced results are analyzed in this script as in Figure 5 of the manuscript.

### Analysis script
The notebook is dividied in 2 subsections:
1. General analysis of morphology -> gene expression association: Figure 5 A-F
2. In depth analysis of expression predictions: Figure 5 G-K

## Input files:
### Necessary files
The input files for the analysis are deposited in Requires the following files:
- `fucci_final.h5`: The raw dataset containing both preprocessed imaging and transcriptomics data. Available as described in README parent folder.
- `all_channels_results_qc_clean.csv`: Compiled test results of models. From section `2_DeepLearning_pipeline`
- `test_preds_hoechst_s1.h5ad` : Prediction on test set for nucleus model. From section `2_DeepLearning_pipeline`
- `test_preds_seg_s1.h5ad` : Prediction on test set for cell shape model. From section `2_DeepLearning_pipeline`
- `DEGs_x_window.tsv` : DE-genes per angular model. From section: `1_scRNA-seq_analysis`
- `cid_whitelist.csv` : Cell ids that pass quality filters. From section: `1_scRNA-seq_analysis`
- `cid_nuc_toremove.csv`: Cell ids to filter that are debris or fragmented cells. From section: `1_scRNA-seq_analysis`

### Analysis of rerun results
In order to run the analysis script with results obtained from rerunning either '' or '' paths have to be adapted. 
- ML training results:
  - `test_preds_hoechst_s1.h5ad` has to be replaced by the h5ad from nucleus model from one seed e.g. '<TRAIN_FOLDER_NUCLEUS_MODEL>/seed_0/test_preds.h5ad'
  - `test_preds_seg_s1.h5ad` has to be replaced by the h5ad from cell shape model from one seed e.g. '<TRAIN_FOLDER_CELLSHAPE_MODEL>/seed_0/test_preds.h5ad'
  - `all_channels_results_qc_clean.csv` from compilation script.
  > ⚠️ File `all_channels_results_qc_clean.csv`, `test_preds_hoechst_s1.h5ad`, `test_preds_seg_s1.h5ad` have to be exchanged together!
- scRNA-seq analysis results:
  - `cid_whitelist.csv`, `cid_nuc_toremove.csv`, and `DEGs_x_window.tsv` have to be replaced
  > ⚠️ Changes in the scRNA-seq analysis will change files `cid_whitelist.csv` and `cid_nuc_toremove.csv` which will impact the whole analysis.


