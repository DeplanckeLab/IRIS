#Author: Joern Pezoldt
#Date: 06.10.2024
#Purpose: GO analysis of differentially expressed genes

#######################
#Compute and Plot
#######################
#####
#Load
#####
stripy_scanpy_log2fc <- read.delim(".......scanpy_to_R/DEGs_CD8_Naive_stripy_scanpbyELISA.logFC_filt.txt")
normi_scanpy_log2fc <- read.delim(".......scanpy_to_R//DEGs_CD8_Naive_normi_scanpbyELISA.logFC_filt.txt")
gene_universe_all <- read.delim(".......scanpy_to_R//genes_CD8_human_donors.txt")$X0
gene_universe_forDEGs <- read.delim(".......scanpy_to_R//genes_CD8_human_donors_naive_forDEGs.txt")$names

#####
#Parameters
#####
species = "human"
#All genes expressed across all donors profiled for CD8+ naive T cells (ca. 33K)

#####
#Compute GO
#####
gene_universe <- rownames(matrix)
length(gene_universe)
#Adapt colnames to function
colnames(stripy_scanpy_log2fc) <- c("gene","scores","avg_log2FC","p_val","p_val_adj","pct_nz_group")
rownames(stripy_scanpy_log2fc) <- stripy_scanpy_log2fc$gene
colnames(normi_scanpy_log2fc) <- c("gene","scores","avg_log2FC","p_val","p_val_adj","pct_nz_group")
rownames(normi_scanpy_log2fc) <- normi_scanpy_log2fc$gene

#Make list of DEG tables
l_DEGS_naive_stripy_normi <- list(stripy_scanpy_log2fc,
                                  normi_scanpy_log2fc)
names(l_DEGS_naive_stripy_normi) <- c("stripy","normi")
lapply(l_DEGS_naive_stripy_normi,nrow)

#Compute GO terms
l_GO_naive_stripy_normi_33000 <- gene_ontology_universe(species, gene_universe, l_DEGS_naive_stripy_normi)
l_GO_naive_stripy_normi_11000 <- gene_ontology_universe(species, gene_universe_all, l_DEGS_naive_stripy_normi)
l_GO_naive_stripy_normi_800 <- gene_ontology_universe(species, gene_universe_forDEGs, l_DEGS_naive_stripy_normi)

#####
#Plot
#####
#Heatmap----------------
#Set parameters for thresholding
pval = 0.0001
max_pval = 7
gene_ontology_heatmap(l_GO_naive_scanpy_stripy_normi_log2fc, pval, max_pval)

#Bubbleplot 33K
pval = 0.00025
ratio = 7
GO_naive_stripy <- l_GO_naive_stripy_normi[[1]]
GO_naive_stripy$Ratio <- GO_naive_stripy$Significant / GO_naive_stripy$Annotated * 100
GO_naive_stripy$weight <- as.numeric(GO_naive_stripy$weight)
GO_naive_stripy_tobubble <- subset(GO_naive_stripy, weight < pval & Ratio > ratio)
GO_naive_stripy_tobubble$log10pval <- round(-log10(GO_naive_stripy_tobubble$weight), 1)
#GO_naive_stripy_tobubble <- GO_naive_stripy_tobubble[order(GO_naive_stripy_tobubble$ratio), ]

ggplot(GO_naive_stripy_tobubble, aes(x = Ratio, y = Term, size = log10pval)) +
  geom_point(color = "deepskyblue3") + 
  geom_point(shape = 1,colour = "black")+
  scale_size(name = "-log10(pValue)",
             range = c(3, 7)) +
  xlim(0, 27)


#############################
#Functions
#############################
#####
#Gene ontology analysis
#####
#RETRURNs: list of GOs per list of differentially expressed genes suplied with p-values, GO:IDs, names...
gene_ontology_geneuniverse <- function(
  # Species used
  species,
  # All genes expressed in dataset
  gene_universe,
  #named list of differentially expressed genes
  ## each list element contains one table with at least columns "gene", "p_val"
  list_DEGs){
  #Load libraries
  library(stringr)
  library(ggplot2)
  library(foreign)
  library(GO.db)
  library(topGO)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  
  #####
  #Make a gene2GO list
  #####
  if(species == "human"){
    x <- org.Hs.egGO}
  if(species == "mouse"){
    x <- org.Mm.egGO}
  # Get the entrez gene identifiers that are mapped to a GO ID
  mapped_genes <- mappedkeys(x)
  # Build a list GeneID and GOs
  GeneID2GO <- list()
  
  xx <- as.list(x[mapped_genes])
  
  for(i in 1:length(xx)){
    #Initiate vector to collect GOIDs for geneID i
    GO_vector <- c()
    #Get geneID
    Gene_ID <- ls(xx[i])
    #grab data for geneID_i
    temp <- xx[[i]]
    
    #check Ontology category
    for(i in 1:length(temp)){
      category <- as.character(temp[[i]]["Ontology"])
      #if ontology category matches collect GOIDs  (flex)
      if(category == "BP"){
        temp_GOID <- as.character(temp[[i]]["GOID"])
        GO_vector <- c(GO_vector, temp_GOID)
      }
      #print(GO_vector)
    }
    #Generate list name geneID, content
    GO_IDs <- GO_vector 
    GeneID2GO[[Gene_ID]] <- GO_IDs 
  }
  
  GO2GeneID <- inverseList(GeneID2GO)
  
  #####
  #Setup input for GO analysis
  #####
  l_input_DEGs <- list_DEGs
  new.cluster.ids <- names(l_input_DEGs)
  myfiles <- l_input_DEGs
  mygenes <- gene_universe
  #GeneSymbol to GeneID
  
  if(species == "human"){
    for(i in 1:length(myfiles)){
      #i =1
      #Translate DEGs
      idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Hs.egSYMBOL)
      SYMBOL <- toTable(org.Hs.egSYMBOL)
      head(SYMBOL)
      m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
      GENE_ID <- SYMBOL$gene_id[m]
      store_i <- cbind(GENE_ID, myfiles[[i]])
      colnames(store_i) <- c("GENE_ID",colnames(store_i)[2:ncol(store_i)])
      #store_i <- store_i[complete.cases(store_i), ]
      myfiles[[i]] <- store_i
    }
    #Translate detected genes
    idfound_mygenes <- mygenes %in% mappedRkeys(org.Hs.egSYMBOL)
    SYMBOL <- toTable(org.Hs.egSYMBOL)
    head(SYMBOL)
    m <- match(mygenes, SYMBOL$symbol)
    GENE_ID <- SYMBOL$gene_id[m]
    store <- cbind(GENE_ID, mygenes)
    colnames(store) <- c("GENE_ID",colnames(store)[2:ncol(store)])
    mygenes <- store
  }
  if(species == "mouse"){
    for(i in 1:length(myfiles)){
      #i =1
      idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Mm.egSYMBOL)
      SYMBOL <- toTable(org.Mm.egSYMBOL)
      head(SYMBOL)
      m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
      GENE_ID <- SYMBOL$gene_id[m]
      store_i <- cbind(GENE_ID, myfiles[[i]])
      colnames(store_i) <- c("GENE_ID",colnames(store_i)[2:ncol(store_i)])
      #store_i <- store_i[complete.cases(store_i), ]
      myfiles[[i]] <- store_i
    }
    #Translate detected genes
    idfound_mygenes <- mygenes %in% mappedRkeys(org.Mm.egSYMBOL)
    SYMBOL <- toTable(org.Mm.egSYMBOL)
    head(SYMBOL)
    m <- match(mygenes, SYMBOL$symbol)
    GENE_ID <- SYMBOL$gene_id[m]
    store <- cbind(GENE_ID, mygenes)
    colnames(store) <- c("GENE_ID",colnames(store)[2:ncol(store)])
    mygenes <- store
  }
  
  #Perform GO analysis single
  #determine number of clusters: split by cluster
  #link gene IDs with p-values
  #perform GO
  l_gene_pval <- list()
  l_all_gene_pval <- list()
  
  for(i in 1:length(myfiles)){
    #i= 1
    genes_of_interest_i <- as.character(myfiles[[i]]$GENE_ID)
    genes_pval_i <- myfiles[[i]]$p_val
    names(genes_pval_i) <- genes_of_interest_i
    l_gene_pval[[i]] <- genes_pval_i 
  }
  print(class(l_gene_pval))
  print(length(l_gene_pval))
  l_all_gene_pval[[length(l_all_gene_pval)+1]] <- l_gene_pval
  
  
  #get one gene list
  l_gene_pval <- l_all_gene_pval[[1]]
  
  #Gene universe
  #Genes in DEGs list
  geneNames <- unique(as.data.frame(mygenes)$GENE_ID)[2:length(mygenes)]
  #Size gene universe
  length(geneNames)
  
  #gene lists
  l_gene_List <- list()
  for(i in 1:length(l_gene_pval)){
    geneList_i <- factor(as.integer(geneNames %in% names(l_gene_pval[[i]])))
    names(geneList_i) <- geneNames
    l_gene_List[[i]] <- geneList_i 
  }
  length(l_gene_List)
  
  #####
  #GO analysis
  #####
  List_allRes <- list()
  #Do  GO statistics for all gene lists
  for(i in 1:length(l_gene_List)){
    #Access the gene lists and p-values for the differentially expressed genes
    #i = 1
    print(paste0(i, " ||| ", names(l_input_DEGs)[i]))
    geneList <- l_gene_List[[i]]
    pvalue_of_interest <- l_gene_pval[[i]]
    
    #build GOdata object
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = pvalue_of_interest,
                  annot = annFUN.gene2GO , gene2GO = GeneID2GO, nodeSize = 5)
    
    #get number of differentially expressed genes in the GOdata object
    sg <- sigGenes(GOdata)
    numSigGenes(GOdata)
    #get the number of GO_IDs that are within the applied GeneUniverse
    graph(GOdata)
    number_GOIDs <- usedGO(GOdata)
    number_nodes <- length(number_GOIDs)
    
    
    #Run statistics
    #Fisher's works with weight
    test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.stat)
    #KS works with elim but not with weight
    #test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
    #resultElim <- getSigGroups(GOdata, test.stat)
    #runTest a high level interface for testing Fisher
    resultFis <- runTest(GOdata, statistic = "fisher")
    #Kolmogorov-Smirnov includes p-values
    #test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    #resultKS <- getSigGroups(GOdata,  test.stat)
    #runTest a high level interface for testing KS
    #elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    #make table
    allRes <- GenTable(GOdata, classic = resultFis, weight = resultWeight,
                       orderBy = "weight", topNodes = number_nodes)
    #KS = resultKS, 
    #ranksOf = "KS", 
    #make list of result tables
    List_allRes[[i]] <- allRes
  }
  names(List_allRes) <- new.cluster.ids
  #return
  List_allRes
}


#####
#Plot Heatmaps of GO statistics
#####
#RETURNs: 
## 1) heatmap
## 2) table with combined statistics per GO ID and DEG list identifer
gene_ontology_heatmap <- function(
    #output gene_ontology_compute
  list_GO_per_cluster,
  #p-value cut-off to include in heatmap
  pval_cutoff,
  #Cutoff to plot on heatmap to condense color space
  max_pval){
  #make RETURN list
  list_return <- list()
  #load library
  library(pheatmap)
  #compile individual tables
  List_TopGOs <- list()
  for(i in 1:length(list_GO_per_cluster)){
    #i = 1
    table_i <- list_GO_per_cluster[[i]]
    table_i$weight <- as.numeric(table_i$weight)
    table_tophits <- subset(table_i, weight < pval_cutoff)
    #naive 5.0e-03
    #print(i)
    #print(nrow(table_tophits))
    List_TopGOs[[i]] <- table_tophits
  }
  #collect all TopGos
  TopGOs_vector <- c()
  for(i in 1:length(List_TopGOs)){
    table_i <- List_TopGOs[[i]]
    TopGOs_i <- as.character(table_i$GO.ID)
    TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
  }
  
  #Condense tables and sort by GO.ID
  l_topGO <- list()
  for(i in 1:length(list_GO_per_cluster)){
    #i = 1
    table_i <- list_GO_per_cluster[[i]]
    head(table_i)
    table_i_subset <- subset(table_i, table_i$GO.ID %in% TopGOs_vector)
    table_i_subset_GO_weight <- table_i_subset[,c("GO.ID","weight", "Term")]
    cluster_weight_name <- names(list_GO_per_cluster)[i]
    colnames(table_i_subset_GO_weight)[2] <- cluster_weight_name
    table_i_subset_GO_weight[,2] <- as.numeric(table_i_subset_GO_weight[,2])
    table_i_subset_GO_weight <- table_i_subset_GO_weight[order(table_i_subset_GO_weight$"GO.ID", decreasing=TRUE), ]
    l_topGO[[i]] <- table_i_subset_GO_weight
  }
  #merger
  data_TopGO_weight = Reduce(function(...) merge(..., all=T), l_topGO)
  
  ####
  #Make GO comparison heatmap
  ####
  data_heatmap <- data_TopGO_weight
  rownames(data_heatmap) <- paste0(data_heatmap$GO.ID, "_", data_heatmap$Term)
  data_heatmap[is.na(data_heatmap)] <- 1
  data_heatmap <- data_heatmap[,3:ncol(data_heatmap)]
  #log transfrom
  data_heatmap <- -log10(data_heatmap)
  #Common Minimum at 5
  data_heatmap[data_heatmap > max_pval] <- max_pval
  data_heatmap_matrix <- data.matrix(data_heatmap)
  #colnames(data_heatmap_matrix) <- new.cluster.ids
  title <- c("Differential GO scRNASeq")
  ggplot_GO_heatmap <- pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
                                treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
                                scale = "none", border_color = "black", cellwidth = 10,
                                cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
                                main = title)
  #return
  list_return[[1]] <- data_TopGO_weight
  list_return[[2]] <- ggplot_GO_heatmap
  names(list_return) <- c("table_GO_weigths","ggplot_heat_GO")
  list_return
}
