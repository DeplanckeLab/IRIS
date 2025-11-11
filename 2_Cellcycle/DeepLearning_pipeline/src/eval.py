from scipy.stats import pearsonr, spearmanr
import pandas as pd
import scanpy as sc
from sklearn.metrics import r2_score

# Function to calculate results
def calculate_results(adata_subset, preds, gts):
    corrs = []
    pvals = []
    corrs_nonzero = []
    pvals_nonzero = []
    spearman = []
    spearman_pvals = []
    gene_names = []
    r2s = []
    masks = gts > 0
    gene_counts = masks.sum(axis=0)
    for idx in range(adata_subset.shape[-1]):
        if gene_counts[idx] > 0:
            corr, pval = pearsonr(gts[:, idx], preds[:, idx])
            spear, spear_pval = spearmanr(gts[:, idx], preds[:, idx])
            r2 = r2_score(gts[:, idx], preds[:, idx])
            
            m = masks[:, idx]
            if m.sum() > 1:
                corr_nonzero, pval_nonzero = pearsonr(gts[m, idx], preds[m, idx])
                gene_names.append(adata_subset.var.index[idx])
                corrs.append(corr)
                pvals.append(pval)
                spearman.append(spear)
                spearman_pvals.append(spear_pval)
                corrs_nonzero.append(corr_nonzero)
                pvals_nonzero.append(pval_nonzero)
                r2s.append(r2)

    results_df = pd.DataFrame({
        'genes': gene_names,
        'R2': r2s,
        'Pearson': corrs,
        'pvals': pvals,
        'Spearman': spearman,
        'pvals (spearman)': spearman_pvals,
        'nonzero_Pearson': corrs_nonzero,
        'nonzero_pvals': pvals_nonzero
    })
    results_df = results_df.dropna().sort_values(by=['Pearson'], ascending=False)# .set_index('genes')
    return results_df

def calculate_results_per_batch(adata_subset, preds, gts):
    labels = adata_subset.obs['exp_id']
    unique_labels = labels.unique()
    df_final = None
    for label in unique_labels:
        m = labels == label
        results_batch_df = calculate_results(adata_subset[m], preds[m], gts[m])
        columns = {k: f'{k}_{label}' for k in results_batch_df.columns if k != 'genes'}
        results_batch_df = results_batch_df.rename(columns=columns)
        
        # Merge results
        if df_final is None:
            df_final = results_batch_df
        else:
            df_final = pd.merge(df_final, results_batch_df, on='genes', suffixes=('', ''))

    # Average scores
    columns = [c for c in df_final.columns if c.startswith('Pearson')]
    df_final['avg_Pearson'] = df_final[columns].mean(axis=1)
    columns = [c for c in df_final.columns if c.startswith('Spearman')]
    df_final['avg_Spearman'] = df_final[columns].mean(axis=1)
    columns = [c for c in df_final.columns if c.startswith('R2')]
    df_final['avg_R2'] = df_final[columns].mean(axis=1)
    df_final = df_final.sort_values(by=['avg_Pearson'], ascending=False)

    return df_final

def adata_results(adata_subset, preds, gts, output_features=None):
    adata_pred = sc.AnnData(gts)
    adata_pred.var = adata_subset.var
    adata_pred.layers['pred'] = preds
    
    if output_features is not None:
        adata_pred.obsm['features'] = output_features
    
    """
    if results_df is not None:
        for c in set(results_df.columns) - set(['genes']):
            adata_pred.var[c] = results_df.set_index('genes').loc[adata_pred.var.index][c]
    """
    adata_pred.obs = adata_subset.obs

    return adata_pred
