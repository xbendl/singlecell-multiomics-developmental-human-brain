library(Seurat)
library(Signac)
library(gaston)
library(ggplot2)
library(reshape2)
library(data.table)
library(Matrix)
library(GenomicRanges)


######## Function for variance component model ########

varcomp <- function(object, shuffle = FALSE, seed = 2021, ncores = 10){
 
  ## CollapseToLongestTranscript(), DistanceToTSS() functions can be found in utils.R
  
  gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = object[['ATAC']]))
  
  rna = object[['SCT']]@data
  if(shuffle){
    set.seed(seed)
    rna = rna[,sample(1:ncol(rna), ncol(rna), replace = F)]
  }
  
  genes = rownames(rna)[which(rowMeans(rna > 0) >= 0.1)]
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
  
  peaks <- granges(object[['ATAC']])
  atac <- object[['ATAC']]@data
  peaks.keep <- which(rowMeans(atac > 0) >= 0.1)
  peaks <- peaks[peaks.keep]
  atac <- atac[peaks.keep,]
  
  ## Identify peaks at genes' promoter regions, i.e. peaks_TSS
  gene.promoters <- promoters(gene.coords.use, upstream = 1000, downstream = 100)
  hits <- findOverlaps(query = peaks, subject = gene.promoters)
  peaks_TSS <- rownames(atac)[unique(queryHits(hits))]
  
  
  ## Identify peaks at genes' enhancer regions, i.e. peaks_DE
  peak_distance_DE <- DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = 5e+05
  )
  peaks_DE <- rownames(atac)[which(rowSums(peak_distance_DE) >= 1)]
  peaks_DE <- setdiff(peaks_DE, peaks_TSS)
  
  
  ## Generate correlation matrices for different components
  P = cor(as.matrix(atac[peaks_TSS,]))
  
  E = cor(as.matrix(atac[peaks_DE,]))
  
  I = matrix(0, nrow = nrow(P), ncol = ncol(P), dimnames = list(rownames(P), colnames(P)))
  samples <- object$sample.id
  for (i in samples) {
    cells <- which(object$sample.id == i)
    I[cells, cells] <- 1
  }
  
  A = matrix(0, nrow = nrow(P), ncol = ncol(P), dimnames = list(rownames(P), colnames(P)))
  ageGroups <- object$age.group
  for (i in ageGroups) {
    cells <- which(object$age.group == i)
    A[cells, cells] <- 1
  }
  
  
  ## Apply AIREML to infer variance parameters
  vals <- parallel::mclapply(genes, function(i){
    Y <- scale(rna[i,])
    mod <- lmm.aireml(Y = Y, K = list(P, E, I, A), verbose = FALSE)
    round(c(mod$sigma2, mod$tau), 3)
  }, mc.cores = ncores)
  
  vals <- matrix(unlist(vals), ncol = 5, byrow = TRUE)
  
  rownames(vals) = genes
  
  return(vals)
  
}


######## Perform variance component model to the original and shuffled data ########
load('../data_processed/pseudobulk_500_50_0119.rda')
pbmc <- pseudobulk
rm(pseudobulk); gc()

val_all <- varcomp(pbmc)
val_shuffle <- varcomp(pbmc, shuffle = T)

save(val_all, val_shuffle, file = '../data_processed/varcomp_vals_0624.rda')



######## Plotting ########
vals = val_all    
# vals = val_shuffle

vdf <- data.frame(vals/rowSums(vals)*100)
names(vdf) <- c('Unexplained', 'Promoter', 'Enhancer', 'Individual', 'Age')
vdf$gene <- rownames(vals)
vdf <- vdf[order((vdf$Promoter + vdf$Enhancer), decreasing = T),]
vdf2 <- rbind(vdf[(vdf$Promoter+vdf$Enhancer) > 1,],
              (vdf[(vdf$Promoter+vdf$Enhancer) < 1,])[(order(vdf[(vdf$Promoter + vdf$Enhancer) < 1, 'Age'], decreasing = T)),]
)

vdf2$rank <- 1:dim(vdf2)[1]

extLabels = sapply(colnames(vdf2)[1:5], function(x) paste0(x, " (", round(mean(vdf2[,x]),2), "%)"))

ldf <- melt(vdf2, id.var=c("gene", "rank"))
ldf$variable = factor(ldf$variable, levels = c('Unexplained', 'Individual', 'Age', 'Promoter', 'Enhancer'), ordered = T)


ldf_all = ldf         #original 
# ldf_random = ldf    #shuffled

p1 <- ggplot(ldf_all, aes(x = rank, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(fill = "Component", x = "Genes", y = "% Variance Explained") +
  scale_fill_manual(values = c("grey", "blue", "#E4288A", "#E4AB00", "#67A61A"), 
                    labels=extLabels) +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        legend.position = c(0.2, 0.3), legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "white")) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_x_continuous( expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Original')

p1


p2 <- ggplot(ldf_random, aes(x = rank, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(fill = "Component", x = "Genes") +
  scale_fill_manual(values = c("grey", "blue", "#E4288A", "#E4AB00", "#67A61A"), 
                    labels=extLabels_random) +
  theme(panel.grid = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.7),
        legend.background = element_rect(fill = "white", color = "white")) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_x_continuous( expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  ggtitle('Permuted')

p2

p1+p2

