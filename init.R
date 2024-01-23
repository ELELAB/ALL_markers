### This script runs all analyses associated with the discovery of gene
### expression markers in two acute lymphoblastic leukemia subtypes. The data
### that is analyzed is the TARGET-ALL-P2 gene expression project from the
### TARGET database.


# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Download data from OSF --------------------------------------------------
library("osfr")
node <- osf_retrieve_node("kgfpv")
files <- osf_ls_files(node)
osf_download(files, 
             recurse = TRUE, 
             conflicts = "skip")

# Restore packages from lockfile ------------------------------------------
library(renv)
renv::settings$external.libraries(c("env_ALL/lib/R/library/"))
renv::restore()
renv::activate()
detach("package:renv", 
       unload = TRUE)

# Set up for Cairo rendering, for headless machines -----------------------
options(bitmapType = "cairo")


## Run scripts ------------------------------------------------------------


## GET DATA ---------------------------------------------------------------

# Get TARGET-ALL-P2 gene expression data, investigate metadata and subset data 
# to contain only bone marrow samples
source("TARGET_data/get_TARGET_data.R")

# -------------------------------------------------------------------------


## REPLICATES -------------------------------------------------------------

# Investigate replicate samples in data
source("TARGET_replicates/TARGET_ALL_P2_bmp_replicates.R")

# Adjust replicates in data
source("TARGET_replicates/TARGET_ALL_P2_bmp_replicates_adj.R")

# -------------------------------------------------------------------------


## PROCESSING -------------------------------------------------------------

# Process data by removing outliers (preprocessing), normalizing for GC content 
# and library size (normalization) and removing lowly expressed genes (filtering)
source("TARGET_processing/bone_marrow_primary_comp/TARGET_process_bmp_comp.R")

# Investigate lost genes following normalization step
source("TARGET_processing/bone_marrow_primary_comp/TARGET_norm_genes_bmp_comp.R")

# Visualize gene types of lost genes following normalization step
source("TARGET_processing/bone_marrow_primary_comp/TARGET_norm_genes_plot_bmp_comp.R")

# --------------------------------------------------------------------------


## VOOM TRANSFORMATION  ----------------------------------------------------

# Voom transform raw data
source("TARGET_transform/voom_transform_bmp_raw/TARGET_ALL_P2_bmp_raw_voom.R")

# Voom transform filtered data
source("TARGET_transform/voom_transform_bmp_comp/TARGET_ALL_P2_bmp_voom_comp.R")

# --------------------------------------------------------------------------


## MULTIDIMENSIONAL SCALING ANALYSIS BEFORE BATCH CORRECTION ---------------

# Extract information about samples such as subtype, metadata etc
source("TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_info_comp.R")

# Perform MDS on raw gene expression data
source("TARGET_pca/bone_marrow_primary_pca_raw/TARGET_ALL_P2_bmp_raw_pca.R")

# Perform MDS on processed data before batch correction
source("TARGET_pca/bone_marrow_primary_pca_comp/init_pca/TARGET_ALL_P2_bmp_pca_comp.R")

# --------------------------------------------------------------------------


## BATCH CORRECTION  -------------------------------------------------------

# Apply batch correction
source("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_batch_comp.R")

# --------------------------------------------------------------------------


## MULTIDIMENSIONAL SCALING ANALYSIS AFTER BATCH CORRECTION  ---------------

# Perform MDS on batch corrected data
source("TARGET_pca/bone_marrow_primary_pca_comp/batch_corr_pca/TARGET_ALL_P2_bmp_pca_batch_comp.R")

# --------------------------------------------------------------------------


## DIFFERENTIAL EXPRESSION ANALYSIS ----------------------------------------

# Perform DEA using DESeq2 method
source("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/TARGET_ALL_P2_bmp_dea_deseq2.R")

# Perform DEA using edgeR method
source("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/TARGET_ALL_P2_bmp_dea_edgeR.R")

# Perform DEA using limma-voom method
source("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/TARGET_ALL_P2_bmp_dea_limma.R")

# Make visualizations of gene counts for the three methods
source("TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/TARGET_ALL_P2_bmp_dea_gene_counts.R")

# Make visualizations of DEA results found using the three methods
source("TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/TARGET_ALL_P2_bmp_dea_visualization.R")

# Compare results of DEA results found using the three methods
source("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/TARGET_ALL_P2_bmp_dea_compare.R")

# Intersect identified DEGs across the three methods
source("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/TARGET_ALL_P2_bmp_dea_gene_intersection_lists.R")

# --------------------------------------------------------------------------


## HOUSEKEEPING GENES ------------------------------------------------------

# Find overlap between identified DEGs and housekeeping genes
source("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_housekeeping.R")

# Visualize identified housekeeping DEGs
source("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_housekeeping_visualization.R")

# --------------------------------------------------------------------------


## ENRICHMENT ANALYSIS -----------------------------------------------------

# Perform enrichment analysis of consensus DEGs
source("TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/TARGET_ALL_P2_bmp_DEG_enrichment.R")

# --------------------------------------------------------------------------


## REGULARIZED LOGISTIC REGRESSION -----------------------------------------

# Perform regularized logistic regression of gene expression data
source("TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso.R")

# Visualize genes found from logistic regression 
source("TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/TARGET_ALL_P2_bmp_consensus_biotypes.R")

# --------------------------------------------------------------------------


## GENE CONTRIBUTIONS FROM PCA ---------------------------------------------

# Perform PCA on gene expression data to find feature (gene) contributions to 
# principal components
source("TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/TARGET_ALL_P2_bmp_pca_gene_contrib.R")

# Visualize feature (gene) contributions to principal components
source("TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes.R")

# --------------------------------------------------------------------------

## COMPARE RESULTS ACROSS METHODS ------------------------------------------

# Compare results across methods
source("TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_compare.R")

# --------------------------------------------------------------------------


## UNSUPERVISED CLUSTERING -------------------------------------------------

# Perform unsupervised clustering using cola framework
source("TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_run_cola.R")

# Analyze and visualize results from unsupervised clustering
source("TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_analyze_cola.R")

# --------------------------------------------------------------------------


## RANDOM FOREST ON PREDICTED CLUSTERS -------------------------------------

# Perform variable selecting using random forest on predicted clusters
source("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_clusters_random_forest.R")

# Analyze selected variables from random forest
source("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R")

# --------------------------------------------------------------------------


## SURVIVAL ANALYSIS -------------------------------------------------------

# Perform survival analysis on predicted subtype-related gene expression markers
source("TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_survival_subtype_markers.R")

# Perform survival analysis on predicted cluster-related gene expression markers
source("TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_survival_cluster_markers.R")

# --------------------------------------------------------------------------
