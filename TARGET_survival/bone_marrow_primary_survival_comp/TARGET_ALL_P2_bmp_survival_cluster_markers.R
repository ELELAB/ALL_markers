### This script performs survival analysis on predicted gene expression
### markers separating two predicted clusters of T-ALL


# Load libraries
library(tidyverse)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(survMisc)

# Source functions
source("TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_survival_functions.R")



### Load data ---------------------------------------------------------------------------------------------------------

# Load clinical data
clinical <- get(load("TARGET_data/TARGET_ALL_P2_clinical_data.rda"))

# Load predicted markers between predicted T-ALL clusters
cluster_genes_info <- read_csv("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_selected_features_T_clusters_rf.csv")

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))

# Load predicted cluster/subgroup labels of all samples
clusters_pred <- read_csv("TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_samples_clusters.csv")


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Extract subtype genes from full table
cluster_genes <- cluster_genes_info$ENSEMBL_ID 

# Get samples belonging to clusters 2 and 3 (T-ALL)
cluster_samples <- clusters_pred %>% 
  dplyr::filter(class == 2 | class == 3) 

# Subset gene expression matrix to only contain samples belonging to
# clusters 2 and 3
ALL_P2_corr_clusters <- ALL_P2_corr[, cluster_samples$barcodes]

# Convert column names of gene expression matrix to go from full barcode
# level to patient level
colnames(ALL_P2_corr_clusters) <- colnames(ALL_P2_corr_clusters) %>% 
  str_extract(string = ., pattern = "^[^-]+-[^-]+-[^-]+")

# Remove duplicate "submitter_id" columns from clinical data
clinical_subset <- clinical[, -c(73, 95)]

# Subset the clinical data to contain only patients of clusters 2 and 3
# and select only variables related to survival
# Remove patients with "unknown" vital status
clinical_survival <- clinical_subset %>% 
  as_tibble() %>%
  dplyr::filter(submitter_id %in% colnames(ALL_P2_corr_clusters)) %>%
  dplyr::select(submitter_id, days_to_last_follow_up,
                days_to_death, vital_status, gender, 
                age_at_diagnosis, primary_diagnosis) %>%
  dplyr::filter(!vital_status == "Unknown") %>%
  dplyr::mutate(days_to_death = as.numeric(days_to_death), 
                time = case_when(vital_status == "Dead" ~ days_to_death,
                                 TRUE ~ days_to_last_follow_up)) %>%
  dplyr::filter(!is.na(time)) %>%
  dplyr::mutate(time_years = time/365,
                age_years = age_at_diagnosis/365,
                vital_status_binary = case_when(vital_status == "Dead" ~ 1,
                                                TRUE ~ 0),
                gender_binary = case_when(gender == "female" ~ 1,
                                          gender == "male" ~ 0),
                primary_diagnosis = case_when(primary_diagnosis == "T lymphoblastic leukemia/lymphoma" ~ "T-ALL"))


### Analyze data ---------------------------------------------------------------------------------------------------------

## Survival analysis 

## Test for Cox proportionality assumption before doing Cox regression analysis

# For each gene i
# test the Cox proportionality assumption
cox_ph_test <- map(cluster_genes, function(i) {
  
  cox_reg_ph(exp_data = ALL_P2_corr_clusters, gene = i, clinical_data = clinical_survival)
  
})

## Proceed with only genes where the proportional hazards assumption is met

# For each tibble containing proportional hazards assumption test of a gene set
sig_genes <- map(cox_ph_test, function(x) {
  
  # Keep only genes that meet the assumption
  x %>% 
    group_by(gene) %>% 
    dplyr::filter(all(p > 0.05)) %>% 
    ungroup %>% 
    pull(gene) %>% 
    unique()
  
}) %>% 
  unlist()


## Fit a univariate Cox regression model for each gene

# For each gene i
# fit a univariate Cox regression model 
cox_reg_uni_results <- map(sig_genes, function(i) {
  
  cox_reg_gene_set <- cox_reg_univariate(exp_data = ALL_P2_corr_clusters, gene = i, clinical_data = clinical_survival) 
  
}) %>% 
  bind_rows() %>% 
  # Correct p-values from univariate Cox regression for multiple testing using fdr method
  dplyr::mutate(fdr = p.adjust(`Pr(>|z|)`, method = "fdr"))


## Find genes with a significant effect on survival from univariate Cox
## regression 
genes_for_multi <- cox_reg_uni_results %>% 
  dplyr::filter(fdr < 0.05) %>% 
  pull(gene)


## Perform Kaplan-Meier survival analysis on the cluster genes 
km_clusters <- map(cluster_genes, function(i) {
  
  km_survival(exp_data = ALL_P2_corr_clusters, gene = i, clinical_data = clinical_survival)
  
})
names(km_clusters) <- cluster_genes


## Save data ------------------------------------

# Save results of proportional hazards assumption test
save(cox_ph_test, file = "TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_cluster_markers_cox_ph.rda")

# Save results of Cox univariate regression analysis
save(cox_reg_uni_results, file = "TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_cluster_markers_cox_univariate.rda")

# Save all Kaplan-Meier survival plots 
map(seq.int(length(km_clusters)), function(x) {
  
  pdf(str_c("TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_cluster_markers_km_", names(km_clusters)[x], ".pdf"),
      width = 10, height = 10)
  print(km_clusters[[x]], newpage = FALSE)
  dev.off()
  
})






