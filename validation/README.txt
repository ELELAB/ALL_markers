* INTRODUCTION *

This directory, validation, contains scripts used to validate the predicted
gene expression markers in an independent cohort of pediatric patients with
ALL from a Danish hospital. As the associated data cannot be granted, these scripts are not runnable. 


* REQUIREMENTS *

Following R and R packages are needed: 
- R 4.2.0
- R packages:
	* tidyverse (version 2.0.0)
	* readxl (version 1.4.3)
	* SummarizedExperiment (version 1.28.0)
	* GenomicRanges (version 1.50.2)
	* TCGAbiolinks (version 2.25.3)
	* sva (version 3.46.0)
	* stats (version 4.2.0)
	* patchwork (version 1.1.3)
	* factoextra (version 1.0.7)
	* FactoMineR (version 2.6)
	* viridis (version 0.6.4)
	* gplots (version 3.1.3)
	* survival (version 3.4.0)
	* survminer (version 0.4.9)
	* survMisc (version 0.5.6)
	* limma (version 3.54.2)


* OVERVIEW OF SCRIPTS *

get_iCOPE_data.R:
- loads expression data and wrangles the expression and meta data to be on correct formats for subsequent analysis
- extracts meta data about samples
- visualizes age distribution of patients 

iCOPE_TCGAbiolinks_processing.R:
- processes (preprocessing, normalization and filtering) on expression data using TCGAbiolinks

iCOPE_voom_transform.R:
- performs voom transformation of raw and filtered expression data

iCOPE_init_pca.R:
- performs MDS of raw and filtered expression data coloring for year of diagnosis and subtype

iCOPE_batch_correction.R:
- applies batch correction correcting expression data for year of diagnosis batch factor

iCOPE_candidate_genes_pca.R:
- performs and visualizes PCA and MDS of expression data of predicted gene expression markers, of all genes except these predicted markers, and of all genes in expression data

iCOPE_clustering_candidates.R:
- performs unsupervised clustering of expression data of predicted expression markers and visualizes this in heatmaps

iCOPE_survival.R:
- performs survival analysis using Cox regression on predicted gene expression markers 

