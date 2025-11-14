##libraries

#general
import os
import csv
import h5py
import pandas as pd
from pandas import HDFStore
import matplotlib.pyplot as plt
import re
import seaborn as sns
import numpy as np
import scipy.stats as stats

#scanpy
import scanpy as sc
import anndata as ad
from matplotlib import rcParams
import gseapy as gp


##functions to deal with hdf5 files 
def gen_indexdf(dataset):
    "JJ code snippets to generate an indexed table"
    df_dict = {}
    for cell_id in dataset.keys():
        if cell_id == 'index_table': continue
        cell_ds = dataset[cell_id]
        if not df_dict:
            df_dict = {attr: [] for attr in cell_ds.attrs.keys()}
            df_dict["cellid_h5"] = []
        for attr in dataset[cell_id].attrs.keys():
            df_dict[attr].append(cell_ds.attrs[attr])
        df_dict["cellid_h5"].append(cell_id)

    cell_index_df = pd.DataFrame(df_dict)
    cell_index_df.set_index('cellid_h5', inplace=True)
    return cell_index_df

def generate_pd_dataframe(ds_dir):
    index_pddf = None
    index_exists = False

    with pd.HDFStore(ds_dir) as ds: index_exists = '/index_table' in ds

    if not index_exists:
        print("Generating index. This might take a minute.")
        with h5py.File(ds_dir) as dataset:
            index_pddf = gen_indexdf(dataset)
        # Remove the 'centroid' column if it exists
        if 'centroid' in index_pddf.columns:
            index_pddf.drop(columns=['centroid'], inplace=True)
        # Save the DataFrame to HDF5
        ds_h5pd = pd.HDFStore(ds_dir)
        ds_h5pd.append('/index_table', index_pddf, format='t', append=False)
        ds_h5pd.close()
    else:
        print("Index exists.")
        with pd.HDFStore(ds_dir, mode='r') as h5store_handle:
            index_pddf = h5store_handle.get('/index_table')
            
    return index_pddf

def extract_plate_and_cc(folder_name):
    # Use regex to extract the number after 'P' and 'CC'
    plate_match = re.search(r'_P(\d+)', folder_name)
    cc_match = re.search(r'CC(\d+)', folder_name)
    
    # If found, return both plate number and cell count
    plate = plate_match.group(1) if plate_match else None
    cc = cc_match.group(1) if cc_match else None
    return plate, cc



##filtering and QC functions 
def perform_qc(adata):
    """
    Perform QC on the provided AnnData object, including:
    - Calculating mitochondrial gene metrics
    - Plotting violin plots for key metrics
    - Plotting scatter plots of total counts vs number of genes, colored by mitochondrial content.
    
    Parameters:
    adata: AnnData object
        Single-cell data to perform QC on.
        
    Returns:
    None. QC metrics are plotted as violin plots and added as columns to adata.obs. 
    """
    # Step 1: Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    
    # Step 2: Calculate QC metrics, including mitochondrial percentages
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
    
    # Step 3: Plot violin plots for key QC metrics
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        size=2,
    )
    
    # Step 4: Plot scatter plot of total counts vs number of genes, colored by mitochondrial percentage
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", size=30)
    
    # Return message to indicate the QC process is complete
    print("QC completed: Violin and scatter plots generated.")
    


##filtering and QC functions 
def perform_qc2(adata, 
               verbose=True
):
    """
    Perform QC on the provided AnnData object, including:
    - Calculating mitochondrial, ribosomal and heat-shock gene metrics
    - Plotting violin plots for key metrics
    - Plotting scatter plots of total counts vs number of genes, colored by mitochondrial content.
    Differently from perform_qc(), it computes also the percentage of 'ribo' and 'hsp' genes.
    
    Parameters:
    adata: AnnData object
        Single-cell data to perform QC on.
        
    Returns:
    None. QC metrics are plotted as violin plots and added as columns to adata.obs.
    """
    
    # Step 1: identify mitochondrial, ribosomal and heat-shock genes
    # Add .str.toupper() to perfrom case-insensitive search
    
    # mitochondrial genes, "MT-" for human
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    # Heatshock 
    adata.var["hsp"] = adata.var_names.str.upper().str.startswith(("HSP"))
    
    # Step 2: Calculate QC metrics, including mitochondrial percentages
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hsp"], inplace=True, log1p=True)

    if verbose: 
        # Step 3: Plot violin plots for key QC metrics
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            size=2,
        )
        # Some QCs: number of genes expressed in the count matrix, total counts per cell, percentage of counts in mitochondrial genes
        sc.pl.violin(
            adata,
            ["pct_counts_ribo", "pct_counts_hsp"],
            jitter=0.4,
            multi_panel=True,
        )
    
        # Step 4: Plot scatter plot of total counts vs number of genes, colored by mitochondrial percentage
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", size=30)
        
        # Return message to indicate the QC process is complete
        print("QC completed: Violin and scatter plots generated.")
    
    print("QC completed")



def preprocess_single_cell_data(
    adata, 
    mt_threshold=20,           # Mitochondrial percentage threshold
    min_genes=2000,            # Minimum number of genes per cell
    min_cells=3,               # Minimum number of cells per gene
    n_top_genes=2000,          # Number of highly variable genes to select
    max_gene_presence=0.99    # Maximum presence of a gene in cells   
):
    """
    Preprocess single-cell data, including filtering, normalization, and selection of highly variable genes.
    
    Parameters:
    - adata (AnnData): Single-cell AnnData object.
    - mt_threshold (float): Maximum mitochondrial percentage allowed per cell.
    - min_genes (int): Minimum number of genes per cell.
    - min_cells (int): Minimum number of cells per gene.
    - n_top_genes (int): Number of highly variable genes to select.
    - max_gene_presence (float): Maximum fraction of cells in which a gene is allowed.
    - batch_key (str): Batch key for highly variable gene selection.
    
    Returns:
    AnnData: Processed AnnData object.
    """
    # Step 1: Filter cells based on mitochondrial percentage
    adata = adata[adata.obs['pct_counts_mt'] < mt_threshold, :].copy()
    
    # Step 2: Filter cells with fewer than 'min_genes' genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Step 3: Filter genes expressed in fewer than 'min_cells' cells
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Step 4: Normalize and log-transform the data
    adata.layers["counts"] = adata.X.copy()  # Save raw counts
    
    sc.pp.normalize_total(adata)
    
    # Logarithmize the data to reduce the effect of large differences, making the data more normally distributed 
    sc.pp.log1p(adata) #applies log tranform: log(1+x) to prevent taking the log of zero  

    # Step 5: Filter genes expressed in more than max_gene_presence * total cells
    sc.pp.filter_genes(adata, max_cells=int(adata.n_obs * max_gene_presence))
    
    # Step 6: Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key="sample")
    
    # Step 7: Plot highly variable genes
    sc.pl.highly_variable_genes(adata)
    
    # Step 8: Subset the data to only include highly variable genes
    adata = adata[:, adata.var['highly_variable']].copy()
    
    print(f"Data preprocessed: {adata.shape[0]} cells and {adata.shape[1]} highly variable genes remaining.")
    
    return adata


def filter_adata(
    adata,           
    min_genes=1000,            
    min_cells_pct=0.02,     
    mt_pct_thresh=15  
):
    """
    Filter single-cell data based on QC metrics. 

    Parameters: 
    - adata : Single-cell AnnData object.
    - min_genes : minimum number of genes per cell 
    - min_cells_pct : minimum number of cells per gene
    - mt_pct_thresh : mitochondrial percentage threshold - maximum pct of counts in mitochrondrial genes 

    Returns:
    Filtered adata object. 
    """
    print(f'N. of "cells x genes" before filtering: {adata.shape}' + "\n")
    adata_shape_init = adata.shape

    # Step 1: Filter genes expressed in fewer than 'min_cells' cells
    min_cells = np.round(adata.obs.shape[0] * min_cells_pct, 3)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f'N. of genes after filtering for min cells: {adata.shape[1]} - (Removed {adata_shape_init[1] - adata.shape[1]} genes)')
        
    # Step 1: Filter cells with fewer than 'min_genes' genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print(f'N. of cells after filtering for min genes: {adata.shape[0]} - Removed {adata_shape_init[0] - adata.shape[0]} cells')
    adata_shape_new = adata.shape
    
    # Step 3: Filter cells based on mitochondrial percentage
    adata = adata[adata.obs['pct_counts_mt'] < mt_pct_thresh, :].copy()
    print(f'N. of cells after filtering for pct mito: {adata.shape[0]} - Removed {adata_shape_new[0] - adata.shape[0]} cells')
    
    print("\n" + f'Final N. of "cells x genes" after filtering: {adata.shape}')

    return adata
    

def normalize_adata(
    adata, 
    target_sum=None, 
    n_top_genes=2000, 
    subset_to_variable=False, 
    scale=True, 
    batch_key='sample'
):
    """
    Normalize single-cell expression data by total counts over all genes and log-transformation.
    Then find most variable genes and scale data to unit mean-variance. 

    Parameters: 
    - adata : Single-cell AnnData object.
    - target_sum : Size factor for normalization. Default: median count; for CMPs: 1e6.
    - n_top_genes : number of highly-variable genes to select. 
    - subset_to_variable : bool. Whether to subset the AnnData object to retain only the highly-variable genes. Default: False. 
    - scale : bool. Whether the data should be scaled to unit mean-variance or not.

    Returns: 
    Modified adata object with normalized data and highly variable genes. 
    If scaling is performed, the unscaled log-normlized data is stored in adata.layers["log_norm"].
    """
    
    # Step 1: Normalize and log-transform the data (similar to Seurat NormalizeData())
    print("Normalizing data with count depth scaling and log plus one (log1p) transformation")
    print(f'Scaling factor: {target_sum}' + "\n")
    
    # Save raw counts
    adata.layers["counts"] = adata.X.copy() 

    # Normalize each cell by total counts over all genes
    sc.pp.normalize_total(adata, 
                          target_sum=target_sum,         # Size factor for normalization. Default: median count; for CMPs: 1e6
                          exclude_highly_expressed=False # Esclude highly expr. genes from the computation of the size factor for each cell
                         )
    
    # Logarithmize the data to reduce the effect of large differences, making the data more normally distributed 
    sc.pp.log1p(adata) # applies log tranform: log(1+x) to prevent taking the log of zero  

    # Save log-normlized counts (for future DGE)
    adata.layers["log_norm"] = adata.X.copy()
    print('Log-nomalized unscaled data is stored in adata.layers["log_norm"]' + "\n")
    
    # Step 2: Identify highly variable genes
    print(f'Detecting {n_top_genes} most highly variable genes for dim. reduction' + "\n")
    
    # Reduce the dimensionality of the dataset and only include the most informative genes.
    sc.pp.highly_variable_genes(adata, 
                            n_top_genes=n_top_genes, 
                            batch_key=batch_key,   # If specified, highly-variable genes are selected within each batch separately and merged.
                                                  # avoids the selection of batch-specific genes (lightweight batch correction method)
                            flavor='seurat'       # This reproduces the R-implementations of Seurat, Cell Ranger, or Seurat v3.
                           )
        
    # Plot highly variable genes
    sc.pl.highly_variable_genes(adata)

    # Subset the data to only include highly variable genes
    if subset_to_variable==True:
        adata = adata[:, adata.var['highly_variable']].copy()
        print(f"Data preprocessed: {adata.shape[0]} cells and {adata.shape[1]} highly variable genes remaining.")
    else:
        print(f"Data preprocessed: {adata.shape[0]} cells and {adata.shape[1]} genes remaining.")

    
    # Scale each gene to unit variance and zero mean
    if scale:    
        print("Scaling data to unit variance." + "\n")
        sc.pp.scale(adata, max_value=10)
    

    return adata


def plot_split_umap(
    adata, 
    split_by = "sample", 
    color_by = 'Cluster', 
    size=20, 
    palette=None, 
    figsize=(5,5)
):
    '''
    Scatter plot in UMAP basis, split according to a variable.
    It substitues Seurat's DimPlot(..., split.by = 'var')
    
    Parameters:
    - adata : adata object to plot.
    - split_by : a variable in the object metadata to split the plot by, i.e. pass 'batch' to split the plot by batch. 
    - color_by : a variable in the object metadata to color cells by. I.e. pass 'cluster' to color cells in each panel by cluster.
    - size : size of points representing cells. 
    - palette : color palette to color the UMAPs based on color_by. Default is None, colors arealy mapped to cateogries in adata will be used.  

    Returns: 
    - A series of scatter plot in UMAP basis, each one representing a set of cells split according to split_by. 
      Cells are colored by color_by and the number of panels in the figure will be as many as split_by levels are.

    Examples: 
    - plot_split_umap(adata, split_by='sample', color_by='Cluster_names', size=30)
    '''

    # TODO: 
    # - add step to check if split_by can be converted t category or is of wrong type (i.e. numeric)
    # - if number of categories is > 4, change number of rows/columns in figure. 
    
    # Split_by variable should be of type 'category'
    if not isinstance(adata.obs[split_by].dtype, pd.CategoricalDtype):
        adata.obs[split_by] = adata.obs[split_by].astype('category')
    
    categories = adata.obs[split_by].unique()
    num_categories = len(categories)

    # Plot split UMAP
    fig, axes = plt.subplots(1, num_categories, figsize=(figsize[0] * num_categories, figsize[1]))

    if palette is not None:
        palette = sns.color_palette(palette, len(adata.obs[color_by].unique()))

    # Iterate through each category and plot the UMAP for that category in a subplot
    for i, category in enumerate(categories):
        subset = adata[adata.obs[split_by] == category].copy()
        sc.pl.umap(subset, 
                   color=color_by, 
                   #title=f"{split_by}: {category}",
                   title=category,
                   size=size, 
                   palette=palette, 
                   ax=axes[i], 
                   show=False)
    
        # Remove legend manually for all but the first plot
        if i < num_categories-1:
            axes[i].get_legend().remove()

    # Add a legend title to the first subplot
    last_legend = axes[num_categories-1].get_legend()
    if last_legend:
        last_legend.set_title(color_by)

    return fig 



def plot_subsets_distribution(
    adata, 
    category1, 
    category2, 
    x_axis_order=None, 
    y_axis_order=None, 
    colors=None, 
    title=None,
    figsize=(4, 5), 
    ax=None
):
    """
    Plot distribution of a category across another category.
    I.e. Plot the distribution cell-cycle Phases (annotated according to cell-cycle gene signatures)
    across the cell-cycle phases annotated via FUCCI annotation. 

    Parameters:
    - adata: AnnData object with `.obs` containing the categories of interest. 
             Or adata.obs as a Data Frame. 
    - category1 : variable to split the data by. I.e. the name of the column contaning FUCCI cell-cycle phases annotation. 
    - category2 : variable whose distribution needs to be plotted across categroy one. 
                  I.e. the name of the column contaning transcriptomic cell-cycle phases annotation. 
    - x_axis_order: List of ordered categories for the x-axis.
    - y_axis_order: List of levels of category2 in the desired order for the stacked bar plot.
    - colors: Optional dictionary specifying colors for each level of categroy2. 
    - figsize: Tuple for figure size (default is (6, 6)).
    - ax : 

    Returns:
    - None (displays the plot).
    Categroy1 will be on the x-axis. Percentages of cells belonging to the levels of category2 will be displayed on the y-axis.

    Examples: 
    # Plot the distribution of cc-phases across FUCCI cc annotations
    colors = {'G1': 'brown', 'S': 'orange', 'G2M': 'lime'}
    utils_scanpy.plot_subsets_distribution(adata=adata_tidy, 
                                           category1='cyc_im_lenient', 
                                           category2='phase', 
                                           colors=colors, 
                                           y_axis_order=['G2M', 'S', 'G1'])
    """
    
    # If no axes are passed, create a new one, otherwise use the provided ax
    if ax is None:
        #ax = plt.gca()
        fig, ax = plt.subplots(figsize=figsize)
        
    if isinstance(adata, ad.AnnData):
        df = adata.obs.copy()
    else:
        df = adata

    # Ensure the annotation column is categorical - and potentially sorted
    if x_axis_order is None: 
        df[category1] = pd.Categorical(df[category1])
    else:
        df[category1] = pd.Categorical(df[category1], categories=x_axis_order, ordered=True)

    # Group by category1 and categroy2, then count occurrences
    cc_plot = df.groupby([category1, category2]).size().reset_index(name='count')

    # Calculate total cells per category1 to normalize counts (percentage)
    cc_plot['total'] = cc_plot.groupby(category1)['count'].transform('sum')
    cc_plot['percentage'] = np.round((cc_plot['count'] / cc_plot['total']) * 100, 2)

    # Pivot data for stacked bar plot
    pivot_df = cc_plot.pivot(index=category1, columns=category2, values='percentage').fillna(0)

    
    # Sort category2 for consistent color mapping
    if y_axis_order is not None:    
        pivot_df = pivot_df[y_axis_order]

    if colors is not None:        
        # Plot using matplotlib for a stacked bar plot
        pivot_df.plot(kind='bar', stacked=True, width=0.8, 
                      color=[colors[cat2] for cat2 in pivot_df.columns], ax=ax)
    else:
        pivot_df.plot(kind='bar', stacked=True, figsize=figsize, width=0.8, ax=ax)

    if title is None:
        title = f'{category2} distribution by {category1}'
    
    # Customize plot labels and appearance
    ax.set_ylabel('Percentage (%) of cells')
    ax.set_xlabel(category1)
    ax.set_title(title)
    ax.legend(title=category2, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ax.set_yticks(np.arange(0, 101, 10))

    return ax


    

def plot_volcano_DEGs(input_DEGs,
                      sig_metric = 'pvals_adj', # 'pvals_adj' or 'pvals'
                      sig_thresh=0.05,
                      mag_metric = 'logfoldchanges', # 'scores' or 'logfoldchnages'
                      ax=None, 
                      title_key="cluster", 
                      contrasts=None
):
    '''
    Plot volcano of differentially expressed genes between two groups, given a table containing
    results from differential gene expression analysis (with wilcoxon test).

    Parameters: 
    - input_DEGs: Data frame containing DGE results.
    - sig_metric: metric used to evaluate significance. Default: pvals_adj (adjusted P-value).
    - sig_thresh: significance threshold for sig_metris. Default: 0.05.
    - mag_metric : metric used to evaluate magnitude of differential expression. Default: logFC.
    - ax: axis for plotting, if present.
    - title_key: keyword to be added to the default title to understand the result. 
                 Usually the name of the cluster which cells belong to.
    - contrasts : list containing the 2 categories that were used to perform DGE, ordered. I.e. ['group1', 'group2']

    Returns: 
    - Volcano plot showing up and down regulated genes in group1 of contrasts vs. group2. 

    Examples: 
    utils_scanpy.plot_volcano_DEGs(DEGs_stripy_clust_filt, 
                                   sig_metric='pvals_adj', 
                                   mag_metric = 'logFC', 
                                   contrasts=['stripy', 'normi'])
    '''
    
    # If no axes are passed, create a new one, otherwise use the provided ax
    if ax is None:
        ax = plt.gca()
    
    data = input_DEGs

    if contrasts is not None:
        group1 = contrasts[0]
        group2 = contrasts[1]
    else:
        # Select groups to plot one against the other:
        if len(data.group.unique()) == 2:
            group1 = DEGs_all.group.unique()[0]
            group2 = DEGs_all.group.unique()[1]

            data = data[ (data['group'] == group1) | (data['group'] == group2) ]
    
    ax.scatter(data[mag_metric], -np.log10(data[sig_metric]), color='black', edgecolor='grey', s=40)
    ax.set_xlabel(mag_metric, fontsize=12)
    ax.set_ylabel(f'-log10 ({sig_metric})', fontsize=12)
    
    ax.set_title(f'Volcano of Diff. Expressed Genes - {title_key}')
    ax.axvline(0, color='grey', linestyle='--', linewidth=0.5)
    
    # Add horizontal line for significance level (p-value < 0.05)
    significance_level = -np.log10(sig_thresh)  
    ax.axhline(sig_thresh, color='red', linestyle='--', linewidth=1, label='p < 0.05')
    
    ax.annotate('', xy=(0.5, -0.125), xycoords='axes fraction', xytext=(0, -0.125),
            arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=2))
    ax.annotate(contrasts[0], xy=(0.75, -0.165), xycoords='axes fraction', fontsize=11)
    
    ax.annotate('', xy=(0.5, -0.125), xycoords='axes fraction', xytext=(1, -0.125),
            arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=2))
    ax.annotate(contrasts[1], xy=(0.25, -0.165), xycoords='axes fraction', fontsize=11)



# TODO: optimize. Make more general and comment. 
def plot_go_enrichment(
    go_enrichment_results,  
    n_top_terms=10,
    sig_metric = 'Adjusted P-value',
    figsize=(14, 12)
):
    '''
    Plot a heatmap with GO enrichment analysis results for different clusters. 

    Parameters:
    - go_enrichment_results : data frame storing go enrichemnt analysis results for each of the clusters considered. 
    - n_top_terms : number of top GO terms to consider. Default: 10.

    Returns:
    A heatmap with top enriched GO terms for each cluster and their significance. 
    Rows represent top enriched GO terms. Columns represent clusters. The color scale is indicative of significance.     
    '''
    
    # Store DataFrames for GO terms and clusters
    go_terms_list = []

    # Loop through each cluster and extract the top n GO terms
    for cluster, result in go_enrichment_results.items():
        
        # Get the top n GO terms for the cluster
        top_terms = result.head(n_top_terms)[['Term', sig_metric]]   
        top_terms['Cluster'] = cluster
        
        # Append the DataFrame to the list (for later concatenation)
        go_terms_list.append(top_terms)
    
    go_terms_clusters = pd.concat(go_terms_list, ignore_index=True)

    # Pivot the DataFrame to create a matrix where rows are GO terms and columns are clusters
    go_matrix = go_terms_clusters.pivot(index='Term', columns='Cluster', values=sig_metric)
    # Fill missing values with a large p-value (e.g., 1.0) to indicate no enrichment
    go_matrix = go_matrix.fillna(1.0)
    
    # Log-transform the p-values for better visualization
    go_matrix = -go_matrix.applymap(lambda x: np.log10(x))
    
    # Set the minimum and maximum log-transformed p-value for better contrast
    vmin = go_matrix.min().min()  # minimum value in the matrix
    vmax = go_matrix.max().max()  # maximum value in the matrix

    # Cap P-values ?
    
    # Plot a heatmap of the GO terms vs clusters using seaborn
    plt.figure(figsize=figsize)
    
    ax = sns.heatmap(
        go_matrix,
        cmap='Blues',      
        annot=False,       # Remove p-value annotations
        linewidths=0.8,    # Add thicker lines between cells
        linecolor='black', # Color of the lines between cells
        cbar_kws={'label': f'-log10({sig_metric})'},  # Add a color bar with a label
        vmin=vmin,          # Set min value for color scaling
        vmax=vmax             # Set max value for color scaling
    )
    
    # Adjust the labels and title
    plt.xticks(rotation=45, ha='right', fontsize=12)  # Rotate x-axis labels for better readability
    plt.yticks(fontsize=10)  # Increase y-axis label size
    plt.title('GO Term Enrichment Across Clusters', fontsize=14, pad=20)  # Add a title and increase the font size
    
    # Show the plot
    plt.tight_layout()


def _plot_bubble_single(
    result, 
    cluster,
    n_top_terms=10, 
    figsize=(12, 6)
):
    '''
    Private function:
    Plot GO enrichment results as a bubble plot, for a single cluster.
    
    Parameters:
    - result : data frame containing GO enrichment results for one cluster only. 
    - cluster : the name of the cluster considered. It will be used for plot's title. 
    - n_top_terms : number of top GO terms to consider. Default: 10.
    - figsize : size of the plot figure. Default: (12, 6)

    Returns: 
    A bubble plot with top enriched GO terms for one cluster. Bubbles are colored based on significance. 
    Bubbles' size represents the fraction of genes belonging to GO term overlapping with DEGs.     
    '''
    # Calculate GeneRatio as a fraction from the 'Overlap' column (e.g., "25/90")
    result['GeneRatio'] = np.round(result['Overlap'].apply(lambda x: eval(x)), 2)  # Converts '25/90' -> 25/90 
    
    # Get the top n GO terms for the cluster
    top_terms = result.head(n_top_terms)[['Term', 'Adjusted P-value', 'GeneRatio']]
    
    # -log10 transform Adjusted P-values for better visualization
    top_terms['log_adj_pval'] = -np.log10(top_terms['Adjusted P-value'])
    top_terms['log_adj_pval'] = np.round(top_terms['log_adj_pval'], 2)

    # Plot bubble plot using seaborn
    plt.figure(figsize=figsize)
    bubble_plot = sns.scatterplot(
        data=top_terms,
        x='GeneRatio',           # X-axis: Fraction of genes in the category
        y='Term',                # Y-axis: GO terms
        size='GeneRatio',        # Bubble size: Fraction of genes
        hue='-log_adj_pval',           # Bubble color: -log10(adjusted p-value)
        palette='coolwarm',      # Color palette for significance
        sizes=(50, 300),         # Size scale for bubbles
        edgecolor='black',       # Add border to bubbles
        alpha=0.8                # Set transparency for better readability
    )

    # Adjust axis labels and title
    plt.xlabel('Gene Ratio', fontsize=11)
    plt.ylabel('GO Term', fontsize=11)
    plt.title(f'GO Term Enrichment Bubble Plot - Cluster {cluster}', fontsize=14, pad=20)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()


def plot_go_enrichment_bubble(
    go_enrichment_results,
    cluster=None,
    n_top_terms=10,
    figsize=(12, 6)
):
    """
    Plot GO enrichment results as bubble plots, for a set of clusters.
        
    Parameters:
    - go_enrichment_results : a dict of data frames containing GO enrichment results, each one for one cluster. 
                              Or a single data frame containing GO enrichment results for one cluster.
    - cluster : the name of the cluster considered, if only one cluster is present in go_enrichment_results. 
                It will be used for plot's title. 
    - n_top_terms : number of top GO terms to consider. Default: 10.
    - figsize : size of the plot figure. Default: (12, 6)

    Returns: 
    One or a series of bubble plots with top enriched GO terms for each cluster. Bubbles are colored based on significance. 
    Bubbles' size represents the fraction of genes belonging to GO term overlapping with DEGs.     
    """
    
    # If results is a dictionry of results for multiple clusters/subsets
    if type(go_enrichment_results) == dict:

        # Extract the top n GO terms for each cluster
        for cluster, result in go_enrichment_results.items():
            _plot_bubble_single(
                result=result, 
                cluster=cluster, 
                n_top_terms=n_top_terms, 
                figsize=figsize
            )
    
    else:
        _plot_bubble_single(
            result=go_enrichment_results, 
            cluster=cluster, 
            n_top_terms=n_top_terms, 
            figsize=figsize
        )

def plot_go_enrichment_bubble_sided(
    go_enrichment_results,
    clusters,
    n_top_terms=10,
    figsize=(18, 6), 
    default_gene_ratio=0.01,
    default_adj_pval=1.0
):
    """
    Plot two bubble plots for GO enrichment results side by side, showing top GO terms for two clusters.
    
    Parameters:
    - go_enrichment_results : dict of data frames containing GO enrichment results for each cluster.
    - clusters : list of two cluster names to compare.
    - n_top_terms : number of top GO terms to consider. Default: 10.
    - figsize : size of the plot figure. Default: (18, 6).

    Returns:
    Side-by-side bubble plots of top enriched GO terms for the two clusters.
    """
    if len(clusters) != 2:
        raise ValueError("Please provide exactly two clusters for comparison.")

    cluster_1, cluster_2 = clusters
    
    # Extract data for each cluster and calculate GeneRatio
    result_1 = go_enrichment_results[cluster_1].copy()
    result_1['GeneRatio'] = np.round(result_1['Overlap'].apply(lambda x: eval(x)), 2)
    
    result_2 = go_enrichment_results[cluster_2].copy()
    result_2['GeneRatio'] = np.round(result_2['Overlap'].apply(lambda x: eval(x)), 2)
    
    # Get top n terms for both clusters
    top_terms_1 = result_1.head(n_top_terms)['Term']
    top_terms_2 = result_2.head(n_top_terms)['Term']
    top_terms_all = pd.concat([top_terms_1, top_terms_2], axis=0)
    
    # Initialize empty dataframe for combining both cluster data
    combined_df = pd.DataFrame()

    # Fill the data for both clusters
    for t in top_terms_all:
        # Extract metrics for the term t from result_1 and result_2, if present
        one_term_1 = result_1.loc[result_1['Term'] == t][['Term', 'Adjusted P-value', 'GeneRatio']]
        one_term_2 = result_2.loc[result_2['Term'] == t][['Term', 'Adjusted P-value', 'GeneRatio']]
    
        # If the term is missing in one of the clusters, fill with default values
        if one_term_1.empty:
            one_term_1 = pd.DataFrame({'Term': [t], 'Adjusted P-value': [default_adj_pval], 'GeneRatio': [default_gene_ratio]})
        if one_term_2.empty:
            one_term_2 = pd.DataFrame({'Term': [t], 'Adjusted P-value': [default_adj_pval], 'GeneRatio': [default_gene_ratio]})
    
        # Add a cluster column to differentiate
        one_term_1['Cluster'] = cluster_1
        one_term_2['Cluster'] = cluster_2
    
        # Concatenate the results for each term
        combined_df = pd.concat([combined_df, one_term_1, one_term_2], axis=0)
    


    # -log10 transform Adjusted P-values
    combined_df['log_adj_pval'] = -np.log10(combined_df['Adjusted P-value'])

    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    
    # Create the bubble plots for both clusters
    for ax, cluster in zip(axes, clusters):
        cluster_data = combined_df[combined_df['Cluster'] == cluster]
        
        bubble_plot = sns.scatterplot(
            data=cluster_data,
            x='GeneRatio',
            y='Term',
            size='GeneRatio',
            hue='log_adj_pval',
            palette='coolwarm',
            sizes=(50, 300),
            edgecolor='black',
            alpha=0.8,
            ax=ax
        )
    
        # Customizing plot
        ax.set_title(f'Cells cluster: {cluster}', fontsize=14, pad=20)
        ax.set_xlabel('')  # Remove GeneRatio label from the x-axis
        if ax != axes[0]:  # Remove y-axis label for the second plot
            ax.set_ylabel('')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
    
    # Suptitle
    plt.suptitle('GO-BP Terms Enrichment', fontsize=16, y=1.02, x=0.6, ha='center')
    
    # Remove duplicate legends
    handles, labels = axes[1].get_legend_handles_labels()
    axes[0].legend().remove()
    axes[1].legend().remove()
    custom_labels = labels
    custom_labels = [l if l != 'GeneRatio' else 'Gene Ratio' for l in custom_labels]
    custom_labels = [l if l != 'log_adj_pval' else '-log10(adj. p-value)' for l in custom_labels]
        
    legend = fig.legend(
        handles=handles,
        labels=custom_labels,
        loc='upper right',
        bbox_to_anchor=(0.98, 0.86)
    )

    # Adjust layout and display
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)  # Make room for a shared legend



def get_significant_degs(
    adata, 
    p_value_thresh=0.05, 
    sort_by='pval_adj'
):
    '''
    Extract differentially expressed genes (DEGs) from an adata object on which differential gene
    expression was performed, using sc.tl.rank_genes_groups().  

    Parameters: 
    - adata: scanpy adata object.
    - p_value_thresh: p-value threshold for filtering DEGs.  
    - sort_by = variable to sort the output dataframe. Usually one between: pval, pval_adj, logFC.

    Returns: 
    - A pandas data frame of DEGs with columns: ['gene_name', 'score', 'logFC', 'pval', 'pval_adj', 'cluster']

    Examples: 
    significant_genes_df = get_significant_degs(adata)
    '''
    
    # Retrieve relevant data from the object
    gene_names = adata.uns['rank_genes_groups']['names']
    scores = adata.uns['rank_genes_groups']['scores']
    logFC = adata.uns['rank_genes_groups']['logfoldchanges']
    pvals_adj = adata.uns['rank_genes_groups']['pvals_adj']
    pvals = adata.uns['rank_genes_groups']['pvals']
    
    # Convert rec.array to a pandas DataFrame
    clusters = gene_names.dtype.names
    significant_degs = []

    for cluster in clusters:
        genes = gene_names[cluster]
        scores_i = scores[cluster]
        logFC_i = logFC[cluster]
        pvals_adj_i = pvals_adj[cluster]
        pvals_i = pvals[cluster]
        
        # Create a DataFrame for the cluster
        df = pd.DataFrame({
            'gene_name': genes,
            'score' : scores_i, 
            'logFC' : logFC_i,
            'pval': pvals_i,
            'pval_adj': pvals_adj_i
        })

        # Filter by the significant adjusted p-value
        df = df[df['pval'] <= p_value_thresh]
        df['cluster'] = cluster

        if sort_by in ['pval', 'pval_adj']:            
            # Sort by p-value - ascending order
            df = df.sort_values(by=sort_by)
            
        elif sort_by == "logFC":
            # Sort by logFC - descending order
            df = df.sort_values(by=sort_by, ascending=False)

        # Append results
        significant_degs.append(df)

    # Combine all clusters into a single DataFrame
    result = pd.concat(significant_degs, ignore_index=True)
    
    return result


def filter_degs(
    input_df, 
    pvals_adj_thresh=0.05, 
    logFC_thresh=0.2, 
    min_pct_thresh=0.2, 
    only_pos = False
):
    '''
    Filter DEGs based on significance and other metrics.  
    '''

    # Filter by the significant adjusted p-value
    output_df = input_df[
    (input_df['pvals_adj'] <= pvals_adj_thresh) &
    (input_df['logFC'].abs() >= logFC_thresh) &
    (input_df['pct_nz_group'] >= min_pct_thresh)
    ]

    if only_pos:
        output_df = output_df[output_df['logFC'] > 0]

    return output_df

    

def layered_umap_plot(
    adata,
    umap_key='X_umap',
    category_key='phase',
    phase_order=None,
    colors=None,
    title=None,
    save_path=None,
    save_name=None,
    figsize=(6, 6),
    s=20,
    alpha=0.8,
    legend_title=None
):
    """
    Creates a layered UMAP scatter plot using custom order and colors for categorical phases.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with UMAP coordinates in obsm and categorical variable in obs.
    umap_key : str
        Key in adata.obsm containing UMAP coordinates.
    category_key : str
        Key in adata.obs containing categorical phase variable.
    phase_order : list
        List specifying the order in which to layer phases.
    colors : list
        List of colors corresponding to phase_order.
    title : str
        Title of the plot.
    save_path : str or None
        If provided, saves the figure to this path.
    figsize : tuple
        Figure size.
    s : int
        Marker size.
    alpha : float
        Marker transparency.
    legend_title : str or None
        Title for the legend.
    """
    
    fig, ax = plt.subplots(figsize=figsize)
    
    for phase, color in zip(phase_order, colors):
        print(f"Plotting phase: {phase}")
        
        mask = adata.obs[category_key] == phase
        x = adata.obsm[umap_key][mask, 0]
        y = adata.obsm[umap_key][mask, 1]
        
        ax.scatter(x, y, color=color, label=phase, s=s, alpha=alpha)
    
    # Clean axes
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Legend and title
    if legend_title is None:
        legend_title = category_key
    ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    if title is not None:
        ax.set_title(title)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(os.path.join(save_path, f"{save_name}.pdf"), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_path, f"{save_name}.eps"), dpi=300, bbox_inches='tight')

    


def plot_stripy_umap(adata, group_key="stripy_manual_solid", umap_key="X_umap", 
                     colors=None, size=20, title="UMAP subset", ax=None, figsize=(4, 5)):
    """
    Custom UMAP scatter plot highlighting stripy/normi/others.

    Parameters
    ----------
    adata : AnnData
        AnnData object with UMAP coordinates in .obsm[umap_key] 
        and categorical annotations in .obs[group_key].
    group_key : str
        Column in .obs containing group labels (default: 'stripy_manual_solid').
    umap_key : str
        Key in .obsm with UMAP coordinates (default: 'X_umap').
    colors : dict
        Dictionary mapping categories to hex colors.
        Example: {'stripy_center': '#FFB347', 'normi_circle': '#9381FF', 'Undefined': '#e6e6e6'}
    size : int
        Point size for scatter plot.
    title : str
        Plot title.
    """
    if colors is None:
        colors = {
            "TØ Stripy": "#FFB347",
            "TO Conventional": "#9381FF",
            "Undefined": "#e6e6e6"
        }

    # Extract UMAP coordinates
    umap_df = adata.obsm[umap_key]
    categories = adata.obs[group_key]

    # Separate data
    stripy_data = umap_df[categories == "TØ Stripy"]
    normi_data = umap_df[categories == "TO Conventional"]
    other_data = umap_df[~categories.isin(["TØ Stripy", "TO Conventional"])]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot undefined background
    ax.scatter(other_data[:, 0], other_data[:, 1], 
               s=size, c=colors.get("Undefined", "#e6e6e6"), 
               label="Undefined", alpha=0.7)

    # Plot normies
    ax.scatter(normi_data[:, 0], normi_data[:, 1], 
               s=size, c=colors.get("TO Conventional", "#9381FF"), 
               label="TO Conventional")

    # Plot stripies
    ax.scatter(stripy_data[:, 0], stripy_data[:, 1], 
               s=size, c=colors.get("TØ Stripy", "#FFB347"), 
               label="TØ Stripy")

    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])

    # Legend with bigger markers
    legend = ax.legend(loc="upper right", frameon=False, 
                       bbox_to_anchor=(1.55, 0.6), borderaxespad=0.,
                       fontsize="large", scatterpoints=1)
    for handle in legend.legendHandles:
        handle._sizes = [40]

    return ax

