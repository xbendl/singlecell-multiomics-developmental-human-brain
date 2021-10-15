library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)


######## Generate a Seurat object for one single sample ########

## Select one sample and set the working directory to the folder including its raw data
sampleID = "150656"
setwd(paste("../data", sampleID, sep = "/"))

## Read the 10x hdf5 file which contains both RNA-seq and ATAC-seq data. 
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
metadata <- read.csv('per_barcode_metrics.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

## Create Seurat object
pbmc <- CreateSeuratObject(counts = inputdata.10x$`Gene Expression`, 
                           assay = 'RNA',
                           project = sampleID,
                           meta.data = metadata)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


## Now add in the ATAC-seq data
# Only use peaks in standard chromosomes
atac_counts <- inputdata.10x$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get gene annotations for hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "atac_fragments.tsv.gz"
pbmc[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

DefaultAssay(pbmc) <- 'ATAC'
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$pct_reads_in_peaks <- pbmc$atac_peak_region_fragments / pbmc$atac_fragments * 100
pbmc$blacklist_fraction <- FractionCountsInRegion(pbmc, assay = 'ATAC', regions = blacklist_hg38_unified)


## Peak calling using MACS2

peaks <- CallPeaks(pbmc, 
                   assay = 'ATAC',
                   # macs2.path = '/Library/Frameworks/Python.framework/Versions/3.9/bin/macs2'
                   macs2.path = '/home/kaiyi/anaconda3/bin/macs2')

peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_hg38_unified, 
                          invert = TRUE)

# quantify counts in each peak
macs_count <- FeatureMatrix(fragments = Fragments(pbmc),
                            features = peaks,
                            cells = colnames(pbmc))


pbmc[['ATAC']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

## Save the sample-specific Seurat object 
obj_150656 = pbmc
save(obj_150656, file = paste('obj', sampleID, 'processed.rda', sep = '_'))



######## Combine data for all the samples (n=12) to a single Seurat object ########
setwd('../../data_processed')

obj_all <- merge(obj_6007, y = c(obj_6032, obj_11, obj_150656, obj_5977, obj_150666,
                                 obj_4413, obj_4, obj_8, obj_5936, obj_4422, obj_16), 
                 add.cell.id = c("6007","6032","11","150656", "5977", "150666",
                                 "4413", "4", "8", "5936", "4422", "16"), 
                 project = "brain_all")

save(obj_all, file = 'obj_all.rda')

