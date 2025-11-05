
## Vignette at: https://www.bioconductor.org/packages/release/bioc/vignettes/tricycle/inst/doc/tricycle.html 

library(BiocManager)
#BiocManager::install(pkgs = "tricycle", lib = '~/updepla/users/bugani/libs/')
library(tricycle, lib.loc = '~/updepla/users/bugani/libs/')
library(ggplot2)
library(scattermore)
library(scater)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(cowplot)


##


path_matrix <- '~/updepla/users/bugani/projects/from_seurat_to_scanpy/results/2023_FUCCI.3T3_irisDataManagment/2023_FUCCI.fucci_3t3_221124_filtered.scanpy.count_matrix_log.csv'


##


# Example data
data(neurosphere_example)

# Project a single cell data set to pre-learned cell cycle space
neurosphere_example <- project_cycle_space(neurosphere_example)
neurosphere_example

# Plot the projected cell cycle space 
scater::plotReducedDim(neurosphere_example, dimred = "tricycleEmbedding") +
  labs(x = "Projected PC1", y = "Projected PC2") +
  ggtitle(sprintf("Projected cell cycle space (n=%d)",
                  ncol(neurosphere_example))) +
  theme_bw(base_size = 14)


## From a matrix of counts
counts_mat <- as.matrix(neurosphere_example@assays@data$counts)
counts_mat_log <- as.matrix(neurosphere_example@assays@data$logcounts)

# Fix gene names
cleaned_gene_names <- sub("\\..*", "", rownames(counts_mat))
rownames(counts_mat) <- cleaned_gene_names
cleaned_gene_names <- sub("\\..*", "", rownames(counts_mat_log))
rownames(counts_mat_log) <- cleaned_gene_names

sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts_mat), logcounts = counts_mat_log))

sce_example <- project_cycle_space(sce)

scater::plotReducedDim(sce_example, dimred = "tricycleEmbedding") +
  labs(x = "Projected PC1", y = "Projected PC2") +
  ggtitle(sprintf("Projected cell cycle space (n=%d)",
                  ncol(sce_example))) +
  theme_bw(base_size = 14)

# Infer cell cycle position
neurosphere_example <- estimate_cycle_position(neurosphere_example)
names(colData(neurosphere_example))

#' Note that we strive to get high resolution of cell cycle state, and we think the continuous position is more appropriate 
#' when describing the cell cycle. However, to help users understand the position variable, we also note that users can approximately relate: 
#' - 0.5pi to be the start of S stage, 
#' - pi to be the start of G2M stage, 
#' - 1.5pi to be the middle of M stage, 
#' - and 1.75pi-0.25pi to be G1/G0 stage.

# Assesssing performance 

#' Plotting the projection of the data into the cell cycle embedding is shown above. Our observation is that deeper sequenced data 
#' will have a more clearly ellipsoid pattern with an empty interior. As sequencing depth decreases, the radius of the ellipsoid 
#' decreases until the empty interior disappears. So the absence of an interior does not mean the method does not work.

top2a.idx <- which(rowData(neurosphere_example)$Gene == 'Top2a')
fit.l <- fit_periodic_loess(neurosphere_example$tricyclePosition,
                            assay(neurosphere_example, 'logcounts')[top2a.idx,],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(neurosphere_example), ")"))
names(fit.l)
fit.l$fig + theme_bw(base_size = 14)

## Alternative: Infer cell cycle stages

#' - Method proposed by Schwabe et al. (2020)
#' - On average, the performance is quite similar to the original implementation in Revelio package. 
#' - In brief, we calculate the z-scores of highly expressed stage specific cell cycle marker genes, and assign the cell to the stage with the greatest z-score.

neurosphere_example <- estimate_Schwabe_stage(neurosphere_example,
                                              gname.type = 'ENSEMBL',
                                              species = 'mouse')

scater::plotReducedDim(neurosphere_example, dimred = "tricycleEmbedding",
                       colour_by = "CCStage") +
  labs(x = "Projected PC1", y = "Projected PC2",
       title = paste0("Projected cell cycle space (n=", ncol(neurosphere_example), ")")) +
  theme_bw(base_size = 14)



##


# Load FUCCI matrix from scanpy
fucci_matrix <- read.csv(path_matrix, header = T)
rownames(fucci_matrix) <- fucci_matrix$cell_id
fucci_matrix$cell_id <- NULL

# genes on the rows
fucci_matrix <- t(fucci_matrix)

# Convert gene symbols to Ensembl gene ids
gene_ids <- mapIds(org.Mm.eg.db, keys=rownames(fucci_matrix), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
sum(is.na(gene_ids)) 
length(gene_ids)
#rownames(fucci_matrix) <- gene_ids

fucci_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(fucci_matrix)))
fucci_sce_tri <- project_cycle_space(fucci_sce, gname.type = 'SYMBOL')

scater::plotReducedDim(fucci_sce_tri, dimred = "tricycleEmbedding") +
  labs(x = "Projected PC1", y = "Projected PC2") +
  ggtitle(sprintf("Projected cell cycle space (n=%d)",
                  ncol(fucci_sce_tri))) +
  theme_bw(base_size = 14)

# Infer cell cycle position
fucci_sce_tri <- estimate_cycle_position(fucci_sce_tri)
names(colData(fucci_sce_tri))
#View(data.frame(fucci_sce_tri@colData))

# Plot the expresion of 1 gene across the circular trajectory
top2a.idx <- which(rownames(rowData(fucci_sce_tri)) == 'Top2a')
fit.l <- fit_periodic_loess(fucci_sce_tri$tricyclePosition,
                            assay(fucci_sce_tri, 'logcounts')[top2a.idx,],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(fucci_sce_tri), ")"))
names(fit.l)
fit.l$fig + theme_bw(base_size = 14)

Smc2.idx <- which(rownames(rowData(fucci_sce_tri)) == 'Smc2')
fit.l <- fit_periodic_loess(fucci_sce_tri$tricyclePosition,
                            assay(fucci_sce_tri, 'logcounts')[Smc2.idx,],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Smc2 along \u03b8 (n=",
                                               ncol(fucci_sce_tri), ")"))
names(fit.l)
fit.l$fig + theme_bw(base_size = 14)

Pcna.idx <- which(rownames(rowData(fucci_sce_tri)) == 'Pcna')
fit.l <- fit_periodic_loess(fucci_sce_tri$tricyclePosition,
                            assay(fucci_sce_tri, 'logcounts')[Pcna.idx,],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                            fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                               ncol(fucci_sce_tri), ")"))
names(fit.l)
fit.l$fig + theme_bw(base_size = 14)


## Alternative to tricycle: Infer cell cycle stages
fucci_sce_tri <- estimate_Schwabe_stage(fucci_sce_tri,
                                              gname.type = 'SYMBOL',
                                              species = 'mouse')

scater::plotReducedDim(fucci_sce_tri, dimred = "tricycleEmbedding",
                       colour_by = "CCStage") +
  labs(x = "Projected PC1", y = "Projected PC2",
       title = paste0("Projected cell cycle space (n=", ncol(fucci_sce_tri), ")")) +
  theme_bw(base_size = 14)


## Plot out embedding sccater plot colored by cell cycle position
p <- plot_emb_circle_scale(fucci_sce_tri, dimred = 1,
                           point.size = 3.5, point.alpha = 0.9) +
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))


## Assign cell-cycle phases based on tricycle position scores

#' - 0.5pi to be the start of S stage, 
#' - pi to be the start of G2M stage, 
#' - 1.5pi to be the middle of M stage, 
#' - and 1.75pi-0.25pi to be G1/G0 stage.


fucci_sce_tri@colData$tricycle_CC_stage <- 'G0'

# start of M: initially pi*1.25. 

S_phase <- (fucci_sce_tri@colData$tricyclePosition >= pi*0.5 & fucci_sce_tri@colData$tricyclePosition < pi*1)
fucci_sce_tri@colData$tricycle_CC_stage[S_phase] <- 'S'

G2M_phase <- (fucci_sce_tri@colData$tricyclePosition >= pi*1 & fucci_sce_tri@colData$tricyclePosition < pi*1.25)
fucci_sce_tri@colData$tricycle_CC_stage[G2M_phase] <- 'G2M'

M_phase <- (fucci_sce_tri@colData$tricyclePosition >= pi*1.25 & fucci_sce_tri@colData$tricyclePosition < pi*1.75)
fucci_sce_tri@colData$tricycle_CC_stage[M_phase] <- 'M'

G1_phase <- (fucci_sce_tri@colData$tricyclePosition >= pi*1.75 | fucci_sce_tri@colData$tricyclePosition < pi*0.5)
fucci_sce_tri@colData$tricycle_CC_stage[G1_phase] <- 'G1'

table(fucci_sce_tri@colData$tricycle_CC_stage)

## Plot out embedding sccater plot colored by cell cycle position

# Plotting the embedding with discrete color scale
p <- plot_emb_circle_scale(fucci_sce_tri, dimred = 1,
                           point.size = 3, point.alpha = 1, 
                           color_by = "tricycle_CC_stage", 
                           hue.colors = c("#2E22EA", "#F86BE2", "#FCCE7B", "#4BBA0F"), 
                           hue.n = 4) +
  scale_color_manual(values = c("#2E22EA", "#F86BE2", "#FCCE7B", "#4BBA0F")) +
  theme_minimal() +
  labs(title = "Cell Cycle Phases")

plot_grid(p)


# Save tricycle CC state annotation
path_output <- fs::path("~/updepla/users/bugani/projects/from_seurat_to_scanpy/results/2023_FUCCI.3T3_irisDataManagment/2023_FUCCI.fucci_3t3_221124_filtered.scanpy.tricycle_cc_phases_annotation.tsv")
data.frame(fucci_sce_tri@colData) %>% 
  rownames_to_column(., "cell_id") %>% 
  write_tsv(., path_output)


##