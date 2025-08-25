#Scope | 
# 1) Perform Seurat batch correction
# 2) Batch correction experiment-based
#Script-Type | Function
#Input | count matrices
#Output | seurat object
#Author | Joern Pezoldt
#Date | 30.05.2024

#####
#Parameter examples
#####
#species = "human"
#integrate_over = "CC"

#Threshing
#min_cells_percent = 0.005
#min_gene_number = 1000
#mito_cutoff = 15
#min_nUMI = 3000
##Seurat
#dims_use = 30
#resolution = 0.3

#common variables to regress are: "n_UMI", "percent.mt", "Phase"
#vars_to_regress <- c()

#Conditionals
# 1) If-clause for species for cell-cycle genes, mt, rpl, hsp, computation

seurat_integrate <- function(
    #signature cell cycle
    signature_cell_cycle_human_mouse,
    #species
    species,
    #cond.2.rmv, #Conditions to remove
    #UMI matrix
    matrix,
    #meta data
    meta_data,
    #Variables
    # Batch over "CC" or "expID"
    integrate_over,
    
    #Seurat
    min_cells_percent,
    min_gene_number,
    mito_cutoff,
    dims_use,
    resolution,
    #common variables to regress are: "n_UMI", "percent.mt", "Phase"
    vars_to_regress,
    #Path tp save figure
    path){
    #####
    #Load libraries
    #####
    require(Matrix)
    require(Seurat)
    #####
    #Setup species intel
    #####
    if(species == "human"){
      signature_cell_cycle <- subset(signature_cell_cycle_human_mouse, species == "human")
      mito_genes = "^MT-"
      ribo_genes = "^RPL"
      heatshock_genes = "^HSP"
    }
    if(species == "mouse"){
      signature_cell_cycle <- subset(signature_cell_cycle_human_mouse, species == "mouse")
      mito_genes = "^mt-"
      ribo_genes = "^Rpl"
      heatshock_genes = "^Hsp"
    }
    print("<<<Species intel set>>>")
    
    ######
    #Run Seurat
    ######
    matrix <- Matrix(matrix, sparse = TRUE)
    #Ensures that rownames are not considered duplicated
    rownames(matrix) <- make.unique(rownames(matrix))
    
    print("<<<Count matrix>>>")
    #Create Seurat objects
    ##Compute minimum number of cells expressing a gene
    min_cells <- ncol(matrix) * min_cells_percent
    seurat.object_v5 <- CreateSeuratObject(matrix, assay = "RNA", meta.data = meta_data,
                                           min.cells = min_cells, min.features = min_gene_number)
    
    #seurat.object_v5 <- subset(seurat.object_v5, 
     #                          subset = CellType_Condition != cond.2.rmv)
    
    #Calculate QC genes
    seurat.object_v5[["percent.mt"]] <- PercentageFeatureSet(seurat.object_v5, pattern = mito_genes)
    seurat.object_v5[["percent.Rpl"]] <- PercentageFeatureSet(seurat.object_v5, pattern = ribo_genes)
    seurat.object_v5[["percent.Hsp"]] <- PercentageFeatureSet(seurat.object_v5, pattern = heatshock_genes)
    
    #Filter cells
    print(paste("All :",nrow(seurat.object_v5@meta.data)))
    seurat.object_v5 <- subset(seurat.object_v5, subset = nFeature_RNA  > min_gene_number)
    print(paste("After Feature Thresh :",nrow(seurat.object_v5@meta.data)))
    seurat.object_v5 <- subset(seurat.object_v5, subset = percent.mt < mito_cutoff)
    print(paste("After Mito Thresh :",nrow(seurat.object_v5@meta.data)))
    
    print("<<<Seurat object threshed for cutoffs>>>")
    
    #Cell cycle scoring
    s.genes <- subset(signature_cell_cycle, cycle_state == "S")$GeneSymbol
    g2m.genes <- subset(signature_cell_cycle, cycle_state == "G2/M")$GeneSymbol
    #s.genes <- toupper(s.genes)
    #g2m.genes <- toupper(g2m.genes)
    #signature_cell_cycle
    seurat.object_v5 <- NormalizeData(seurat.object_v5)
    DefaultLayer(seurat.object_v5[["RNA"]]) <- 'data'
    seurat.object_v5 <- CellCycleScoring(seurat.object_v5,s.features = s.genes,g2m.features = g2m.genes)
    print("<<<Cell cycle score computed>>>")
    
    #Extract factor over which to integrate
    integrate_by <- seurat.object_v5[[integrate_over]]
    integrate_by_FCI <- rownames(integrate_by)
    integrate_by_type <- integrate_by[,1]
    names(integrate_by_type) <- integrate_by_FCI
    
    #Split data
    seurat.object_v5[["RNA"]] <- split(seurat.object_v5[["RNA"]], f = integrate_by_type)
    
    print("<<<Seurat layers according to integration parameter set>>>")
    
    #Find variable genes, scale and PCA
    seurat.object_v5 <- FindVariableFeatures(seurat.object_v5, verbose = FALSE)
    seurat.object_v5 <- ScaleData(seurat.object_v5,
                                  assay = "RNA",
                                  features = rownames(seurat.object_v5@assays$RNA),
                                  vars.to.regress = vars_to_regress,
                                  verbose = TRUE)
    seurat.object_v5 <- RunPCA(seurat.object_v5, npcs = dims_use, verbose = FALSE)
    
    # Ensure the output directory exists
    #if (!dir.exists(path)) {
    #  dir.create(path, recursive = TRUE)
    #}
    
    #Generate elbow plot for dimension determination
    #p <- ElbowPlot(seurat.object_v5, ndims = dims_use) +
    #  theme_minimal() +
    ##  theme(
    #    plot.background = element_rect(fill = "white", color = NA),
    #    panel.background = element_rect(fill = "white", color = NA),
    #    panel.grid.major = element_line(color = "grey80"),
    #    panel.grid.minor = element_line(color = "grey90")
    #  )
    #print(p)
    #ggsave(filename = paste0(path, "/elbow_plot.png"), plot = p)
    
    #stdevs <- seurat.object_v5@reductions$pca@stdev
    #min_stdev <- min(stdevs)
    #thresh <- min_stdev * 1.05 # Make a st_dev threshold @ 2% deviation
    #Return the dimension which fulfills the above threshold criterium, assign value to parameters
    #dims_use <- (which(stdevs <= thresh)[1])
    
    #print("<<<FVB, Scale, Regression, PCA done>>>")
    
    # Adjust dims_use based on the number of cells in the smallest object
    #num_cells <- ncol(seurat.object_v5)
    #if (dims_use >= num_cells) {
    #  dims_use <- num_cells - 1
    #  warning(paste("Adjusted dims_use to", dims_use, "to be less than the number of cells:", num_cells))
    #}
    
    #Compute Integrations
    seurat.object_v5 <- IntegrateLayers(
      object = seurat.object_v5, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE)
    
    #Cluster and UMAP
    seurat.object_v5 <- FindNeighbors(seurat.object_v5, reduction = "integrated.cca", 
                                      dims = 1:dims_use, verbose = TRUE)
    seurat.object_v5 <- FindClusters(seurat.object_v5, resolution = resolution,
                                     cluster.name = "cca_clusters", verbose = TRUE)
    seurat.object_v5 <- RunUMAP(seurat.object_v5, reduction = "integrated.cca",
                                dims = 1:dims_use, reduction.name = "umap.cca", 
                                verbose = TRUE)
    
    #Join Layers
    seurat.object_v5 <- JoinLayers(seurat.object_v5)
    
    print("<<<Integration, clustering, joing layers done>>>")
    
    return(seurat.object_v5)
}