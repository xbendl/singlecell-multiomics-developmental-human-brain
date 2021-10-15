library(Signac)
library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)

options(future.globals.maxSize = 50000 * 1024^2)

load('../data_processed/obj_all.rda')


######## Quality control to filter cells ########

## Filter cells
# RNA-seq: nCount_RNA, percent.mt
# ATAC-seq: nCount_ATAC, nucleosome_signal, TSS.enrichment
pbmc <- subset(x = obj_all, subset = nCount_RNA > 2e2 & 
                 nCount_RNA < 5e4 &
                 nCount_ATAC > 2e2 &
                 nCount_ATAC < 1e5 &
                 percent.mt < 5 &
                 nucleosome_signal < 3 &
                 TSS.enrichment > 1)



######## Precessing and Clustering ########

preprocess <- function(object, n.pc = 30, n.lsi = 10){
  
  
  ## Filter peaks/genes which are detected in < 10 cells
  
  tmp <- Matrix::rowSums(object[['RNA']]@counts > 0)
  object[['RNA']] <- subset(object[['RNA']], features = names(which(tmp >= 10)))
  tmp <- Matrix::rowSums(object[['ATAC']]@counts > 0)
  object[['ATAC']] <- subset(object[['ATAC']], features = names(which(tmp >= 10)))
  
  
  ## Normalization, dimensional reduction, and clustering on RNA-seq and ATAC-seq separately
  
  # RNA-seq
  DefaultAssay(object) <- 'RNA'
  object <- SCTransform(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, reduction = 'pca', dims = 1:n.pc, assay = 'SCT',
                    reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  object <- FindNeighbors(object, reduction = 'pca', dims = 1:n.pc, assay = 'SCT')
  object <- FindClusters(object, graph.name = 'SCT_snn', algorithm = 3, resolution = 0.2)
  
  # ATAC-seq
  DefaultAssay(object) <- "ATAC"
  object <- RunTFIDF(object, method = 3)
  object <- FindTopFeatures(object, min.cutoff = 'q75')
  object <- RunSVD(object)
  object <- RunUMAP(object, reduction = 'lsi', dims = 2:n.lsi, assay = 'ATAC',
                    reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  object <- FindNeighbors(object, reduction = 'lsi', dims = 2:n.lsi, assay = 'ATAC')
  object <- FindClusters(object, graph.name = 'ATAC_snn', algorithm = 3, resolution = 0.2)
  
  
  # Weighted nearest neighbor (WNN) analysis using both modalities
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("pca", "lsi"),
                                    dims.list = list(1:n.pc, 2:n.lsi),
                                    modality.weight.name = 'RNA.weight')
  object <- RunUMAP(object, nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  object <- FindClusters(object, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.2)
  object <- FindMultiModalNeighbors(object,
                                    reduction.list = list("harmony.pca", "harmony.lsi"),
                                    dims.list = list(1:n.pc, 2:n.lsi),
                                    modality.weight.name = 'RNA.weight.harmony',
                                    knn.graph.name = 'wknn.harmony',
                                    snn.graph.name = 'wsnn.harmony',
                                    weighted.nn.name = 'weighted.nn.harmony')
  

  DefaultAssay(object) <- 'SCT'
  
  return(object)
  
}

pbmc = preprocess(obj_all)
rm(obj_all); gc()



######## Cell types annotation ########

## WNN-derived
celltype <- rep(NA, length = ncol(pbmc))
Idents(pbmc) <- pbmc$wsnn_res.0.2

celltype[which(Idents(pbmc) %in% c(11))] <- 'RG'
celltype[which(Idents(pbmc) %in% c('5_3'))] <- 'IPC'
celltype[which(Idents(pbmc) %in% c(8))] <- 'IN-MGE'
celltype[which(Idents(pbmc) %in% c(10))] <- 'IN-CGE'
celltype[which(Idents(pbmc) %in% c(3))] <- 'IN-fetal'
celltype[which(Idents(pbmc) %in% c(0,14))] <- 'EN-fetal-late'
celltype[which(Idents(pbmc) %in% c(6,12))] <- 'EN'
celltype[which(Idents(pbmc) %in% c('5_0','5_1','5_2','5_4'))] <- 'EN-fetal-early'

celltype[which(Idents(pbmc) %in% c(4,9,15))] <- 'Astrocytes'
celltype[which(Idents(pbmc) %in% c(1,21,22))] <- 'Oligodendrocytes'
celltype[which(Idents(pbmc) %in% c(2,19))] <- 'OPC'
celltype[which(Idents(pbmc) %in% c(7,13,17,23))] <- 'Microglia'
celltype[which(Idents(pbmc) %in% c(16))] <- 'Endothelial'
celltype[which(Idents(pbmc) %in% c(18))] <- 'Pericytes'
celltype[which(Idents(pbmc) %in% c(20))] <- 'VSMC'

pbmc$celltype = celltype


## RNA-derived
celltype <- rep(NA, length = ncol(pbmc))
Idents(pbmc) <- pbmc$SCT_snn_res.0.2.subcluster

celltype[which(Idents(pbmc) %in% c('12_0','12_1','12_3'))] <- 'RG'
celltype[which(Idents(pbmc) %in% c('12_2'))] <- 'IPC'
celltype[which(Idents(pbmc) %in% c(9))] <- 'IN-MGE'
celltype[which(Idents(pbmc) %in% c(10))] <- 'IN-CGE'
celltype[which(Idents(pbmc) %in% c(6))] <- 'IN-fetal'
celltype[which(Idents(pbmc) %in% c(8))] <- 'EN-fetal-early'
celltype[which(Idents(pbmc) %in% c(0,11))] <- 'EN-fetal-late'
celltype[which(Idents(pbmc) %in% c(4))] <- 'EN'

celltype[which(Idents(pbmc) %in% c(5,7,16))] <- 'Astrocytes'
celltype[which(Idents(pbmc) %in% c(1,17))] <- 'Oligodendrocytes'
celltype[which(Idents(pbmc) %in% c(2,15))] <- 'OPC'
celltype[which(Idents(pbmc) %in% c(3,18))] <- 'Microglia'
celltype[which(Idents(pbmc) %in% c(14))] <- 'Endothelial'
celltype[which(Idents(pbmc) %in% c('13_0','13_1','13_4'))] <- 'Pericytes'
celltype[which(Idents(pbmc) %in% c('13_2', '13_3'))] <- 'VSMC'

celltype <- factor(celltype, 
                   levels = c('RG', 'IPC', 'EN-fetal-early', 'EN-fetal-late', 'EN', 
                              'IN-fetal', 'IN-MGE', 'IN-CGE',
                              'OPC', 'Astrocytes', 'Oligodendrocytes','Microglia', 
                              'Endothelial', 'Pericytes', 'VSMC'), 
                   ordered = T)

pbmc$celltype.rna <- celltype


## ATAC-derived
celltype <- rep(NA, length = ncol(pbmc))
Idents(pbmc) <- pbmc$ATAC_snn_res.0.2.subcluster
celltype[which(Idents(pbmc) %in% c('5_1'))] <- 'RG'
celltype[which(Idents(pbmc) %in% c('5_3'))] <- 'IPC'
celltype[which(Idents(pbmc) %in% c(4))] <- 'IN'
celltype[which(Idents(pbmc) %in% c(7))] <- 'IN-fetal'
celltype[which(Idents(pbmc) %in% c(0))] <- 'EN-fetal-late'
celltype[which(Idents(pbmc) %in% c(8,10))] <- 'EN'
celltype[which(Idents(pbmc) %in% c('5_0','5_2','5_4'))] <- 'EN-fetal-early'

celltype[which(Idents(pbmc) %in% c(1,11))] <- 'Astrocytes'
celltype[which(Idents(pbmc) %in% c(2,12))] <- 'Oligodendrocytes'
celltype[which(Idents(pbmc) %in% c(3))] <- 'OPC'
celltype[which(Idents(pbmc) %in% c(6))] <- 'Microglia'
celltype[which(Idents(pbmc) %in% c(9))] <- 'Endothelial'

celltype <- factor(celltype, 
                   levels = c('RG', 'IPC', 'EN-fetal-early', 'EN-fetal-late', 'EN', 
                              'IN-fetal', 'IN',
                              'OPC', 'Astrocytes', 'Oligodendrocytes','Microglia', 
                              'Endothelial'), 
                   ordered = T)

pbmc$celltype.atac <- celltype




######## Additional processing ########

## Re-call peaks for each annotated cell type using MACS2
peaks <- CallPeaks(pbmc, assay = 'ATAC', 
                   macs2.path = '/home/kaiyi/anaconda3/bin/macs2',
                   group.by = 'celltype',
                   outdir = 'MACS2_output',
                   fragment.tempdir = 'MACS2_output',
                   cleanup = FALSE)

# Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_hg38_unified, 
                          invert = TRUE)

# Quantify counts in each peak
DefaultAssay(pbmc) <- 'ATAC'
frags <- Fragments(pbmc)

macs_count <- FeatureMatrix(
  fragments = frags,
  features = peaks,
  cells = colnames(pbmc)
)

# Create a new assay using the MACS2 peak set and add it to the Seurat object
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

pbmc[['peaks']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frags,
  annotation = annotations
)

DefaultAssay(pbmc) <- 'peaks'
pbmc <- RunTFIDF(pbmc, method = 3)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q75')


## Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)
pbmc[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(pbmc, assay = 'GeneActivity')


## Add developmental stages information
ageGroup <- rep('early fetal', length = ncol(pbmc))
sampleID <- pbmc$orig.ident
ageGroup[which(sampleID %in% c('4', '8'))] <- 'late fetal'
ageGroup[which(sampleID %in% c('4413', '4422'))] <- 'infancy'
ageGroup[which(sampleID %in% c('6032', '5977'))] <- 'childhood'
ageGroup[which(sampleID %in% c('6007', '5936'))] <- 'adolescence'
ageGroup[which(sampleID %in% c('150666', '150656'))] <- 'adulthood'
ageGroup <- factor(ageGroup, 
                   levels = c('early fetal', 'late fetal', 'infancy', 'childhood', 'adolescence', 'adulthood'), 
                   ordered = T)
pbmc$age.group = ageGroup


obj_all_processed <- pbmc
save(obj_all_processed, file = '../data_processed/obj_all_processed_v3.rda')

