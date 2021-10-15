# Joint transcriptome and regulome analysis in the developing human brain cortex at the single cell level

## Description
- ```merge_samples.R``` describes the process of creating a single Seurat object from Cell Ranger ARC outputs of different samples. 
- ```preprocessing.R``` includes quality control, pre-processing and clustering steps on the raw data.
- ```pseudobulk.R``` shows how to create the aggregated pseudobulk samples.
- ```varcomp.R``` implements the variance component model.
- ```links.R``` includes peak-gene links identification, and DORC-related quantification.
- ```trajectory.R``` describes trajectory pseudotime analysis on neuronal populations.
- ```run_tradeseq.R``` contains downstream pseudotime analysis by using tradeSeq.
