### This script performs unsupervised clustering on TARGET ALL P2 gene expression data using the cola framework


# Load libraries
library(cola)
library(SummarizedExperiment)


### Load data ---------------------------------------------------------------------------------------------------------

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))

# Load raw gene count data (where replicates have been adjusted for)
ALL_P2_raw_rep <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Extract gene count matrix of raw data
ALL_P2_raw_rep_mat <- assay(ALL_P2_raw_rep)

# Set seed
set.seed(173)

# Run consensus partitioning on batch corrected data
rl <- run_all_consensus_partition_methods(
  ALL_P2_corr, 
  top_value_method = all_top_value_methods(),
  partition_method = all_partition_methods(),
  cores = 4)

# Generate report of whole analysis done on batch corrected data
cola_report(rl, 
            output_dir = "TARGET_clustering/bone_marrow_primary_clustering_comp", 
            title = "TARGET_ALL_P2_bmp_batch_cola_report",
            cores = 4)

# Adjust raw gene count matrix
ALL_P2_raw_rep_mat_adj <- adjust_matrix(ALL_P2_raw_rep_mat)

# Set seed
set.seed(173)

# Run consensus partitioning on raw data
rl_raw <- run_all_consensus_partition_methods(
  ALL_P2_raw_rep_mat_adj, 
  top_value_method = all_top_value_methods(),
  partition_method = all_partition_methods(),
  cores = 4)

# Generate report of whole analysis done on raw data
cola_report(rl_raw, 
            output_dir = "TARGET_clustering/bone_marrow_primary_clustering_comp/cola_on_raw_data", 
            title = "TARGET_ALL_P2_bmp_raw_cola_report",
            cores = 4)


### Save data ---------------------------------------------------------------------------------------------------------

# Save adjusted raw data
saveRDS(ALL_P2_raw_rep_mat_adj,
        file = "TARGET_clustering/bone_marrow_primary_clustering_comp/cola_on_raw_data/TARGET_ALL_P2_bmp_raw_adj_data_cola.rds")

# Save results from consensus partitioning of batch corrected data
saveRDS(rl, 
        file = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_batch_cola_results.rds")

# Save results from consensus partitioning of raw data
saveRDS(rl_raw, 
        file = "TARGET_clustering/bone_marrow_primary_clustering_comp/cola_on_raw_data/TARGET_ALL_P2_bmp_raw_cola_results.rds")




