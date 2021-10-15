library(Seurat)
library(Signac)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(ComplexHeatmap)


######## Running tradeSeq ########

## Load tradeSeq input created from Monocle3 results
load('../data_processed/tradeSeq_allneurons_input_0429.rda')


## Parallel computing
BPPARAM <- BiocParallel::bpparam()
BPPARAM
BPPARAM$workers <- 12


## Select optimal number of knots
set.seed(5)
icMat <- evaluateK(counts = counts, 
                   k = 6:12,
                   pseudotime = pseudotime[,1:3], 
                   cellWeights = cellWeights_simplified, 
                   parallel = TRUE, BPPARAM = BPPARAM)

save(icMat, file = '../data_processed/tradeSeq_allneurons_icMat_0429.rda')


## Run tradeSeq by using the optimal number of knots = 9
# by changing 'counts' to compute fitted results for different modalities (RNA, ATAC, chromvar)
sce <- fitGAM(counts = counts,
              pseudotime = pseudotime[cells,1:3], cellWeights = cellWeights[cells,1:3], 
              nknots = 9, parallel = TRUE, BPPARAM = BPPARAM)

save(sce, file = 'tradeSeq_allneurons_output.rda')


## save knots information into obj_neurons
knots = rep(1, length = ncol(sce))
knots.time = sce@metadata$tradeSeq$knots
t = pseudotime[cells,1]
for (i in 2:8) {
  
  knots[which(t > knots.time[i])] <- i
  
}
names(knots) = names(cells)
table(knots)

tmp = rep(NA, length = ncol(obj_neurons))
names(tmp) = colnames(obj_neurons)
tmp[names(cells)] = knots
obj_neurons$knots = tmp
DimPlot(obj_neurons, reduction = 'umap.rna', group.by = 'knots')



######## Lineage priming analysis ########

## Function of min-max normalization
min_max_normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


## Get the DORC scores and the associated gene expression smoothed over the inferred pseudotime

# Load neuronal peak-gene links
load('../data_processed/links_filtered_allneurons.rda')
peaks.per.gene <- table(links$gene)
genes = names(which(peaks.per.gene >= 5))   

# see how to create DORC x cell matrix in links.R
sce_DORC = fitGAM(counts = DORC,
                  pseudotime = pseudotime[,1:3],
                  cellWeights = cellWeights_simplified,
                  nknots = 9, parallel = TRUE, BPPARAM = BPPARAM)

sce_exprs <- fitGAM(counts = counts,
                    pseudotime = pseudotime[,1:3],
                    cellWeights = cellWeights_simplified,
                    genes = genes,
                    nknots = 9, parallel = TRUE, BPPARAM = BPPARAM)

yhatSmooth_DORC <- predictSmooth(sce_DORC, gene = rownames(sce_DORC), tidy = FALSE)
yhatSmooth_exprs <- predictSmooth(sce_exprs, gene = rownames(sce_exprs), tidy = FALSE)


## Perform min-max normalization
exprs_norm = t(apply(yhatSmooth_exprs, 1, min_max_normalize))
DORC_norm = t(apply(yhatSmooth_DORC, 1, min_max_normalize))
chromvar_norm = t(apply(yhatSmooth_chromvar, 1, min_max_normalize))

## Get residuals between DORC and gene expression
residual <- DORC_norm - exprs_norm


## Plot residual vs. #correlated peaks
df <- data.frame(count = as.vector(peaks.per.gene[rownames(residual)]), 
                 residual = rowMeans(residual),
                 gene = rownames(residual))
df$DORC <- ifelse(df$count >= 5, "DORC", "Not DORC")
df$residual.sign = ifelse(df$residual > 0, "Positive", "Negative")
df1 <- df[which(df$count >= 5),]

p = ggplot(df1, aes(x = count, y = residual)) +
  geom_point(aes(color = residual.sign), size = 1) +
  scale_color_manual(values = c("darkgrey", "red")) +
  xlim(1,22) +
  ylim(-0.1, 0.38)+
  theme_bw(base_size = 12) + 
  NoLegend() +
  geom_text_repel(
    data = subset(df1, residual.sign == 'Positive'),
    aes(label = gene),
    size = 3,
    direction = 'both',
    max.overlaps = 10,
    box.padding = unit(0.25, "lines")
  ) +
  xlab('Number of correlated peaks') +
  ylab('Residual')

p
ggsave('residual_neurons_vs_peaknumbers.tiff')


