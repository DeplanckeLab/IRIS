# Deep learning analysis
In the previous step '' three models were trained and tested:
- Nucleus CNN model
- Cell shape CNN model
- Angular score MLP
The produced results are analyzed in this script as in Figure 5 of the manuscript.

## Input files:
### Necessary files
The input files for the analysis are deposited in Requires the following files:
- fucci_final.h5: The raw dataset containing both preprocessed imaging and transcriptomics data.
- all_channels_results_qc_clean.csv: Compiled results from 2_DeepLearnin_pipeline
- test_preds_hoechst_s1.h5ad : Prediction on test set for nucleus model
- test_preds_seg_s1.h5ad : Prediction on test set for cell shape model
- DEGs_x_window.tsv : DE-genes per angular model from 
- cid_whitelist.csv : 
- cid_nuc_toremove.csv

### Analysis of rerun results
In order to run the analysis script with results obtained from rerunning either '' or '' paths have to be adapted. 

