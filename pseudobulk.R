library(Seurat)
library(Signac)

options(future.globals.maxSize = 5000 * 1024^2)

load('../data_processed/obj_all_processed_v3.rda')
pbmc <- obj_all_processed
rm(obj_all_processed); gc()


######## Creating aggregated pseudobulk samples ########

## Find the nearest neighbors for each cell based on WNN graph
pbmc <- FindMultiModalNeighbors(pbmc, 
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:30, 2:10),
                                k.nn = 99)

neighbors <- pbmc@neighbors$weighted.nn

celltype_info = pbmc$celltype
names(celltype_info) = colnames(pbmc)


## Initialization for the 500 pseudobulk samples
seeds_neighors <- matrix(nrow = 500, ncol = 50)
seeds <- NULL
i = 1

## Select the 49 nearest neighbors within the cells of the same cell type for seed cells
set.seed(2021)
while (i <= nrow(seeds_neighors)) {
  
  seed <- sample(colnames(pbmc), size = 1)
  if(seed %in% seeds){
    print(seed)
    next
  }

  nn <- neighbors@cell.names[neighbors@nn.idx[which(neighbors@cell.names == seed),]]
  tmp <- which(celltype_info[nn] == celltype_info[seed])
  if(length(tmp) >= 49){
    seeds_neighors[i,] <- c(seed, nn[1:49])
    seeds <- c(seeds, seed)
    i = i+1
  }else{
    next
  }

}


## Aggregating counts for pseudobulk samples each of which comprises 50 cells
cells <- unique(as.vector(seeds_neighors))
pbmc <- subset(pbmc, cells = cells)

rna.counts <- pbmc@assays$SCT@counts
atac.counts <- pbmc@assays$peaks@counts

# RNA
rna.counts.agg <- matrix(nrow = nrow(rna.counts), ncol = length(seeds))
for (i in 1:length(seeds)) {
  rna.counts.agg[,i] <- Matrix::rowSums(rna.counts[, seeds_neighors[i,]])
}
rna.counts.agg = Matrix::Matrix(rna.counts.agg, sparse = TRUE)
rownames(rna.counts.agg) <- rownames(rna.counts)
colnames(rna.counts.agg) <- paste('bulk', sep = '_', 1:length(seeds))

# ATAC
atac.counts.agg <- matrix(nrow = nrow(atac.counts), ncol = length(seeds))
for (i in 1:length(seeds)) {
  atac.counts.agg[,i] <- Matrix::rowSums(atac.counts[, seeds_neighors[i,]])
}
atac.counts.agg = Matrix::Matrix(atac.counts.agg, sparse = TRUE)
rownames(atac.counts.agg) <- rownames(atac.counts)
colnames(atac.counts.agg) <- paste('bulk', sep = '_', 1:length(seeds))

# metadata
metadata <- data.frame(seed = seeds)
tmp = as.data.frame(seeds_neighors)
metadata$nearest.neighbors = tmp

celltype <- celltype_info[seeds]
agegroup = sapply(1:length(seeds), function(x) names(which.max(table(pbmc$age.group[which(colnames(pbmc) %in% seeds_neighors[x,])]))))
sampleid = sapply(1:length(seeds), function(x) names(which.max(table(pbmc$orig.ident[which(colnames(pbmc) %in% seeds_neighors[x,])]))))
clusterid = sapply(1:length(seeds), function(x) names(which.max(table(pbmc$wsnn_res.0.2[which(colnames(pbmc) %in% seeds_neighors[x,])]))))

metadata$celltype = celltype
metadata$age.group = agegroup
metadata$sample.id = sampleid
metadata$cluster.id = clusterid
rownames(metadata) <- paste('bulk', sep = '_', 1:length(seeds))

# create Seurat object for pseudobulk samples
pseudobulk <- CreateSeuratObject(counts = rna.counts.agg, assay = 'RNA', 
                                 meta.data = metadata)
pseudobulk[['ATAC']] <- CreateChromatinAssay(counts = atac.counts.agg)
rm(rna.counts); rm(atac.counts)
gc()

pseudobulk <- SCTransform(pseudobulk)
pseudobulk <- RunPCA(pseudobulk)
ElbowPlot(pseudobulk, ndims = 50)

DefaultAssay(pseudobulk) <- "ATAC"
pseudobulk <- RunTFIDF(pseudobulk, method = 3)
pseudobulk <- FindTopFeatures(pseudobulk, min.cutoff = 'q75')
pseudobulk <- RunSVD(pseudobulk)
DepthCor(pseudobulk)

pseudobulk[['ATAC']]@seqinfo <- pbmc[['ATAC']]@seqinfo
pseudobulk[['ATAC']]@annotation <- pbmc[['ATAC']]@annotation

DefaultAssay(pseudobulk) <- 'SCT'

save(pseudobulk, file = '../data_processed/pseudobulk_500_50_0119.rda')


## UMAP visualization of pseudobulk samples
DimPlot(pseudobulk, reduction = 'wnn.umap', label = T, repel = T, pt.size = 0.5, shuffle = T) + 
  NoLegend() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank()) + 
  ggtitle('pseudobulk-celltype')



######## Get 'pseudo-age' for each cell type ########
age = rep(NA, length = ncol(pbmc))

## Use log-scale weights for different developmental stages
age[which(pbmc$age.group %in% c('early fetal'))] = 0
age[which(pbmc$age.group %in% c('late fetal'))] = log10(3)  # 0.48
age[which(pbmc$age.group %in% c('infancy'))] = log10(5)     # 0.70
age[which(pbmc$age.group %in% c('childhood'))] = log10(7)   # 0.85
age[which(pbmc$age.group %in% c('adolescence'))] = log10(9) # 0.95
age[which(pbmc$age.group %in% c('adulthood'))] = 1

age_per_celltype = rep(NA, length(levels(pbmc$celltype)))
names(age_per_celltype) = levels(pbmc$celltype)
for (c in levels(pbmc$celltype)) {
  
  samples = which(pbmc$celltype == c)
  age_per_celltype[c] = sum(age[samples])/length(samples)
  
}
