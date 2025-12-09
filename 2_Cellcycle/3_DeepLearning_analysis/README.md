# Deep learning analysis
In the previous step '' three models were trained and tested:
- Nucleus CNN model
- Cell shape CNN model
- Angular score MLP
The produced results are analyzed in this script as in Figure 5 of the manuscript.

## Input files:
### Necessary files
The input files for the analysis are deposited in Requires the following files:
- `fucci_final.h5`: The raw dataset containing both preprocessed imaging and transcriptomics data. Available as described in README parent folder.
- `all_channels_results_qc_clean.csv`: Compiled test results of models. From section `2_DeepLearnin_pipeline`
- `test_preds_hoechst_s1.h5ad` : Prediction on test set for nucleus model. From section `2_DeepLearnin_pipeline`
- `test_preds_seg_s1.h5ad` : Prediction on test set for cell shape model. From section `2_DeepLearnin_pipeline`
- `DEGs_x_window.tsv` : DE-genes per angular model. From section: `1_scRNA-seq_analysis`
- `cid_whitelist.csv` : Cell ids that pass quality filters. From section: `1_scRNA-seq_analysis`
- `cid_nuc_toremove.csv`: Cell ids to filter that are debris or fragmented cells. From section: `1_scRNA-seq_analysis`

### Analysis of rerun results
In order to run the analysis script with results obtained from rerunning either '' or '' paths have to be adapted. 


