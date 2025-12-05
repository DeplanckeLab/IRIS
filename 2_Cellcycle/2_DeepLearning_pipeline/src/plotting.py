import matplotlib.pyplot as plt

def scatter_top_genes(preds, gts, results_batch_df, adata_subset, k=10):
    labels = adata_subset.obs['exp_id']
    unique_labels = labels.unique()

    pval_cols = [c for c in results_batch_df.columns if c.startswith('pvals_')]
    for i, gene in enumerate(results_batch_df['genes'].values[:k]):
        gene_mask = adata_subset.var.index == gene
        plt.subplot((k - 1) // 5 + 1, 5, i+1)
        corr = results_batch_df.iloc[i]['avg_Pearson']
        pval = results_batch_df.iloc[i][pval_cols].mean() # pearsonr(gts[:, idx], preds[:, idx])

        for label in unique_labels:
            m = labels == label
            plt.scatter(gts[m, gene_mask], preds[m, gene_mask], s=5, label=label)

        plt.legend()
        plt.title(f'{i+1}. {gene} (corr={corr:.2f}, p={pval:.2f})')
        plt.xlabel('True expression')
        plt.ylabel('Predicted expression')


def plot_example_images(dataset_original, dataset_comparison, out_dir, key, channel_idx=0):
        # Plot some example images with and without augmentation
    fig, axs = plt.subplots(2, 5, figsize=(20, 8))
    channel_idx = 0
    for i in range(5):
        img, _ = dataset_original[i]
        axs[0, i].imshow(img.permute(1, 2, 0).numpy()[..., channel_idx, None])
        axs[0, i].set_title(f'Image {i} (no augmentation)')
        axs[0, i].axis('off')

        img_aug, _ = dataset_comparison[i]
        axs[1, i].imshow(img_aug.permute(1, 2, 0).numpy()[..., channel_idx, None])
        axs[1, i].set_title(f'Image {i} (with augmentation)')
        axs[1, i].axis('off')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/example_{key}_{channel_idx}_images.pdf', bbox_inches='tight')
