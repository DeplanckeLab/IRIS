# Publication feature workflow

Run these notebooks from this directory in numeric order:

1. `00_extract_morphology_features.ipynb` extracts and filters morphology features.
2. `01_train_rf_models_and_compute_shap.ipynb` trains the RF models and computes SHAP values.
3. `02_compare_cnn_and_feature_rf_models.ipynb` compares CNN and feature-RF models (Figure 5B; Supplementary Figure 5A-C).
4. `03_analyze_morphology_predicted_programs.ipynb` analyzes predicted transcriptional programs (Figure 5C-K; Supplementary Figure 5D-H).

## Required inputs

- `inputs/fucci_final_nucmask.h5`, or set `FUCCI_H5_PATH` to the HDF5 file elsewhere.
- `inputs/all_channels_results_qc_clean.csv`, containing the frozen CNN/MLP gene metrics.
- `inputs/cid_whitelist.csv` and `inputs/cid_nuc_toremove.csv`, which are bundled curation resources.
- `inputs/FeatureListAnnotated.csv`, defining membership of the four feature-based CNN comparators.

`MANIFEST.tsv` records the provenance and SHA-256 checksum of each bundled input and publication result.

## Rerunning analysis

1. Set `FUCCI_H5_PATH` if the HDF5 file is not at `inputs/fucci_final_nucmask.h5`:

```bash
export FUCCI_H5_PATH=/path/to/fucci_final_nucmask.h5
```

2. Run every cell in notebook 00. It writes `outputs/features_combined_raw.csv` and `outputs/features_combined_corrected.csv`.
3. Run every cell in notebook 01, including the final validation cell. Regenerated model files are staged in `outputs/rerun/`; the provided files in `outputs/` are not replaced. A complete successful run writes `outputs/rerun/rerun_validation_summary.csv`.
4. Inspect the staged results with notebook 03 before promotion:

```bash
PUBLICATION_RESULTS_DIR="$PWD/outputs/rerun" jupyter notebook 03_analyze_morphology_predicted_programs.ipynb
```

5. If validation and inspection succeed, promote the files listed below from `outputs/rerun/` to `outputs/`.
6. Run notebooks 02 and 03 against the promoted files to regenerate the final panels. Figures are written to `figures/outputs/`.

Representative-cell montages in notebook 03 require `FUCCI_H5_PATH`; other sections still run without the image store.

## Promoting a validated rerun

After `outputs/rerun/rerun_validation_summary.csv` has been created and reviewed, replace these files in `outputs/` with their counterparts from `outputs/rerun/`:

- `cell_metadata.csv`
- `rf_gene_metrics.csv`
- `rf_test_predictions.csv.gz`
- `Result_RF_nuc_imgfdr.csv`
- `Result_RF_nuc_shapefdr.csv`
- `Result_RF_cell_imgfdr.csv`
- `Result_RF_cell_shapefdr.csv`
- `individual_gene_fucci_residual_fold_metrics.csv`
- `individual_gene_fucci_residual_gene_mean_metrics.csv`
- `individual_gene_fucci_residual_summary.csv`
- `shap_significant_genes.csv.gz`
