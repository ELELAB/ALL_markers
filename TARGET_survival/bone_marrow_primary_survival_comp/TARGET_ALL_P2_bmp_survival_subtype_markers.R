### This script performs survival analysis on predicted small set of 
### subtype-related gene expression markers 


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

# Load defined small gene set of subtype-related markers
subtype_genes_info <- read_csv("TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_consensus_DEGs_elastic_pca.csv")

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Extract subtype genes from full table
subtype_genes <- subtype_genes_info$ENSEMBL_ID 

# Convert column names of gene expression matrix to go from full barcode
# level to patient level
colnames(ALL_P2_corr) <- colnames(ALL_P2_corr) %>% 
  str_extract(string = ., pattern = "^[^-]+-[^-]+-[^-]+")

# Remove duplicate "submitter_id" columns from clinical data
clinical_subset <- clinical[, -c(73, 95)]

# Subset the clinical data to contain only tumor patients of interest
# and select only variables related to survival
# Remove patients with "unknown" vital status
clinical_survival <- clinical_subset %>% 
  as_tibble() %>%
  dplyr::filter(submitter_id %in% colnames(ALL_P2_corr)) %>%
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
                primary_diagnosis = case_when(primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" ~ "B-ALL",
                                              primary_diagnosis == "T lymphoblastic leukemia/lymphoma" ~ "T-ALL"))


### Analyze data ---------------------------------------------------------------------------------------------------------

## Survival analysis 

## Test for Cox proportionality assumption before doing Cox regression analysis

# For each gene i
# test the Cox proportionality assumption
cox_ph_test <- map(subtype_genes, function(i) {
  
  cox_reg_ph(exp_data = ALL_P2_corr, gene = i, clinical_data = clinical_survival)
  
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
  
  cox_reg_gene_set <- cox_reg_univariate(exp_data = ALL_P2_corr, gene = i, clinical_data = clinical_survival) 
  
}) %>% 
  bind_rows() %>% 
  # Correct p-values from univariate Cox regression for multiple testing using fdr method
  dplyr::mutate(fdr = p.adjust(`Pr(>|z|)`, method = "fdr"))


## Find genes with a significant effect on survival from univariate Cox
## regression 
genes_for_multi <- cox_reg_uni_results %>% 
  dplyr::filter(fdr < 0.05) %>% 
  pull(gene)

## Do multivariate Cox regression on genes with a significant effect on 
## survival from univariate Cox regression 
sig_genes_multivariate <- map(genes_for_multi, function(i) {
  
  cox_reg_multivariate(exp_data = ALL_P2_corr, gene = i, clinical_data = clinical_survival) 
  
}) %>% 
  bind_rows()

## Extract genes in whose expression has a significant effect
## on survival at multivariate level
sig_genes_exp_multi <- sig_genes_multivariate %>% 
  dplyr::filter(variable == "exp_value" & `Pr(>|z|)` < 0.05) %>% 
  pull(gene) 

## Perform Kaplan-Meier survival analysis on genes with a significant
## effect on survival at multivariate level
km_sig_genes_multi <- map(sig_genes_exp_multi, function(i) {
  
  km_survival(exp_data = ALL_P2_corr, gene = i, clinical_data = clinical_survival)
  
})
names(km_sig_genes_multi) <- sig_genes_exp_multi


## Save data ------------------------------------

# Save results of proportional hazards assumption test
save(cox_ph_test, file = "TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_subtype_markers_cox_ph.rda")

# Save results of Cox univariate regression analysis
save(cox_reg_uni_results, file = "TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_subtype_markers_cox_univariate.rda")

# Save results of Cox multivariate regression analysis
save(sig_genes_multivariate, file = "TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_subtype_markers_cox_multivariate.rda")

# Save all Kaplan-Meier survival plots 
map(seq.int(length(km_sig_genes_multi)), function(x) {
  
  pdf(str_c("TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_subtype_markers_km_", names(km_sig_genes_multi)[x], ".pdf"),
      width = 10, height = 10)
  print(km_sig_genes_multi[[x]], newpage = FALSE)
  dev.off()
  
})

