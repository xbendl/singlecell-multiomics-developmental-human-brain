library(Seurat)
library(Signac)
library(cicero)
library(tidyverse)
library(data.table)
library(Matrix)
library(monocle3)
library(SeuratWrappers)


######## Create Seurat object for neuronal cells ########

load('obj_all_processed_v3.rda')

## Select only neuronal cells 
# 22635 cells in total
neurons <- rownames(obj_all_processed@meta.data)[which(obj_all_processed$celltype %in% c('RG', 'IPC', 'EN-fetal-early', 'EN-fetal-late',
                                                                                         'IN-fetal', 'IN-MGE', 'IN-CGE', 'EN'))]
obj_neurons <- subset(obj_all_processed, cells = neurons)
rm(obj_all_processed); gc()


## Remove genes/peaks that are rarely detected in neurons
tmp <- Matrix::rowSums(obj_neurons[['SCT']]@counts > 0)
obj_neurons[['SCT']] <- subset(obj_neurons[['SCT']], features = names(which(tmp >= 10)))
tmp <- Matrix::rowSums(obj_neurons[['peaks']]@counts > 0)
obj_neurons[['peaks']] <- subset(obj_neurons[['peaks']], features = names(which(tmp >= 10)))


## Redo dimensional reduction on the new neuronal dataset
DefaultAssay(obj_neurons) <- 'SCT'
obj_neurons <- FindVariableFeatures(obj_neurons, nfeatures = 3000)
obj_neurons <- RunPCA(obj_neurons)
ElbowPlot(obj_neurons, ndims = 50, reduction = 'pca')

n.pc = 10; n.lsi = 6
obj_neurons <- RunUMAP(obj_neurons, reduction = 'pca', dims = 1:n.pc, assay = 'SCT',
                  reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DimPlot(obj_neurons, reduction = 'umap.rna', group.by = 'sample.id', shuffle = T)


## Remove doublets (n=530) between IN-MGE and IN-CGE
DefaultAssay(obj_neurons) <- 'SCT'
n.pc = 10
obj_neurons <- FindNeighbors(obj_neurons, reduction = 'pca', dims = 1:n.pc, assay = 'SCT')
obj_neurons <- FindClusters(obj_neurons, graph.name = 'SCT_snn', algorithm = 3, resolution = 0.2)
DimPlot(obj_neurons, reduction = 'umap.rna', group.by = 'SCT_snn_res.0.2', shuffle = T, label = T)

Idents(obj_neurons) <- obj_neurons$SCT_snn_res.0.2
obj_neurons <- FindSubCluster(obj_neurons, cluster = 6, graph.name = 'SCT_snn', algorithm = 3, resolution = 0.2)
Idents(obj_neurons) <- obj_neurons$sub.cluster
obj_neurons <- FindSubCluster(obj_neurons, cluster = 10, graph.name = 'SCT_snn', algorithm = 3, resolution = 0.2)
DimPlot(obj_neurons, reduction = 'umap.rna', group.by = 'sub.cluster', shuffle = T, label = T)
table(obj_neurons$sub.cluster)

doublets <- colnames(obj_neurons)[which(obj_neurons$sub.cluster %in% c('6_3', '10_2'))]

obj_neurons = subset(obj_all_processed, cells = setdiff(colnames(obj_neurons), doublets))



######## Apply Monocle3 for pseudotime inference #######

DefaultAssay(obj_neurons) <- 'SCT'
cds <- as.cell_data_set(obj_neurons)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(obj_neurons[["SCT"]])

## Switch UMAP for UMAP.RNA
reducedDims(cds)$UMAP <- reducedDims(cds)$UMAP.RNA

cds <- cluster_cells(cds = cds, reduction_method = "UMAP")

cds <- learn_graph(
  cds,
  learn_graph_control = list(minimal_branch_len = 20),
  use_partition = FALSE,
  close_loop = FALSE
  )

## Use RGs as the starting point
cell_ids <- colnames(cds)[which(obj_neurons$celltype == 'RG')]
closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds),])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

## Compute the trajectory
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = TRUE
  )  

save(cds, file = '../data_processed/Monocle3_cds_allneurons_0429.rda')


## Save pseudotime information into Seurat object
obj_neurons <- AddMetaData(
  object = obj_neurons,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

# 22105 cells
save(obj_neurons, file='../data_processed/obj_neurons_0429.rda')



######## Use Monocle3 results to create input for running tradeSeq ########

## Extract the pseudotimes and cell weights for tradeSeq
# Refer to tradeSeq tutorial at github

library(magrittr)

y_to_cells <- principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells = rownames(y_to_cells)
y_to_cells$Y = y_to_cells$V1

# Get the root vertices
root = cds@principal_graph_aux$UMAP$root_pr_nodes

# Get the other endpoints
mst = principal_graph(cds)$UMAP
endpoints = names(which(igraph::degree(mst) == 1))
endpoints = endpoints[!endpoints %in% root]

# For each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  
  # Find the path between the endpoint and the root
  path = igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path = as.character(path)
  
  # Find the cells that map along that path
  df = y_to_cells[y_to_cells$Y %in% path,]
  df = data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) = endpoint
  return(df)
  
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()

rownames(cellWeights) = colnames(cds)
pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights), nrow = ncol(cds), byrow = FALSE)
rownames(pseudotime) <- rownames(cellWeights)


## Exclude cells whose cell type annotations are inconsistent with lineage assignment

# Assign cells to lineages
EN_lineage = names(which(rowSums(cellWeights[,c(1:3,6)]) > 0))
IN_CGE_lineage = names(which(cellWeights[, 5] > 0))
IN_MGE_lineage = names(which(cellWeights[, 4] > 0))

# Filter out cells for each lineage
EN_lineage = setdiff(EN_lineage, colnames(obj_neurons)[which(obj_neurons$celltype %in% c('IN-MGE', 'IN-CGE', 'IN-fetal'))])
IN_MGE_lineage = setdiff(IN_MGE_lineage, colnames(obj_neurons)[which(obj_neurons$celltype %in% c('EN', 'EN-fetal-early', 'EN-fetal-late', 'IN-CGE'))])
IN_CGE_lineage = setdiff(IN_CGE_lineage, colnames(obj_neurons)[which(obj_neurons$celltype %in% c('EN', 'EN-fetal-early', 'EN-fetal-late', 'IN-MGE'))])

# Add lineage information to Seurat object
lineage = matrix(0, nrow = ncol(obj_neurons), ncol = 3,
                 dimnames = list(colnames(obj_neurons), c('EN', 'IN-MGE', 'IN-CGE')))
lineage[EN_lineage, 'EN'] = 1
lineage[IN_MGE_lineage, 'IN-MGE'] = 1
lineage[IN_CGE_lineage, 'IN-CGE'] = 1
obj_neurons$EN_lineage = lineage[,'EN']
obj_neurons$IN_MGE_lineage = lineage[,'IN-MGE']
obj_neurons$IN_CGE_lineage = lineage[,'IN-CGE']

cellWeights_simplified = lineage
obj_neurons <- subset(obj_neurons, cells = names(which(rowSums(cellWeights_simplified) > 0)))

# 21579 cells
save(obj_neurons, file = '../data_processed/obj_neurons_filtered.rda')


## Use the filtered neuronal dataset as the input for the downstream tradeSeq analysis
counts = obj_neurons[['SCT']]@counts
cellWeights_simplified = cellWeights_simplified[colnames(obj_neurons),]
pseudotime = pseudotime[colnames(obj_neurons),]

save(counts, genes, cellWeights_simplified, pseudotime, file = '../data_processed/tradeSeq_allneurons_input_0429.rda')



######## Developmental trajectory of EN lineage ########

load('../data_processed/links_filtered_allneurons.rda')
load('../data_processed/ENlineage_km_genes.rda')

## Identify genes in links which are significantly varied along the trajectory
genes = unique(links$gene)
res = graph_test(cds[, rownames(cellWeights_simplified)[which(cellWeights_simplified[,1] > 0)]],
                 neighbor_graph = 'principal_graph', cores = 8)
genes_var_EN <- intersect(row.names(subset(res, q_value < 0.01)), genes)


## Find enriched motifs in each k-means cluster 
pbmc = obj_neurons
rm(obj_neurons); gc()
Idents(pbmc) = pbmc$celltype
DefaultAssay(pbmc) <- 'peaks'

open.peaks = AccessiblePeaks(pbmc,
                             cells = colnames(pbmc)[which(pbmc$EN_lineage == 1)])
meta.feature = GetAssayData(pbmc, assay = 'peaks', slot = 'meta.features')

enriched.motifs = vector('list', length = 4)
names(enriched.motifs) = c('km4','km3','km2','km1')
for (i in 1:4) {
  
  idx = which(links$gene %in% km.genes[[i]])
  peaks <- links$peak[idx[order(links$score[idx], decreasing = T)]]
  
  peaks.matched = MatchRegionStats(
    meta.feature = meta.feature[open.peaks,],
    query.feature = meta.feature[peaks,]
  )
  
  enriched.motifs[[i]] = FindMotifs(
    pbmc,
    features = peaks,
    background = peaks.matched
  )
  
  enriched.motifs[[i]]$pvalue.adj <- p.adjust(enriched.motifs[[i]]$pvalue, method = 'fdr')
}


## GO enrichment analysis on genes in each k-means cluster
library(clusterProfiler)

ego <- enrichGO(km.genes[[1]], 
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = 'ALL',
                universe = rownames(obj_neurons[['SCT']]@data))
barplot(ego, showCategory = 20)

