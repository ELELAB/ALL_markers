### This script performs survival analysis on predicted small set of 
### subtype-related gene expression markers 


# Load libraries
library(tidyverse) #version 2.0.0
library(readxl) #version 1.4.3
library(TCGAbiolinks) #version 2.25.3
library(survival) #version 3.4.0
library(survminer) #version 0.4.9
library(survMisc) #version 0.5.6

# Source functions
source("survival/iCOPE_survival_functions.R")


### Load data -----------------------------------------------------------------

# Load survival data
survival_data <- read_excel("survival/survivaldata_leukemia_dtu.xlsx")

# Load voom transformed filtered data
voom_filt_data <- get(load("transform/iCOPE_ALL_filt_data_voom.rda"))

# Load meta data
sample_meta_data_df <- get(load("get_data/iCOPE_ALL_meta_data.rda"))


### Wrangle data --------------------------------------------------------------

# Prepare survival and metadata for finding survival times
survival_data_wrangled <- survival_data %>% 
  dplyr::mutate(vital_status_binary = case_when(`alive XX.XX.XXXX` == "Y" | `alive XX.XX.XXXX` == "y" ~ 0,
                                                `alive XX.XX.XXXX` == "N" | `alive XX.XX.XXXX` == "n" ~ 1,
                                                TRUE ~ NA)) %>% 
  dplyr::rename("MRDNR" = "MRD") %>% 
  full_join(x = ., y = sample_meta_data_df, by = "MRDNR") %>% 
  dplyr::rename("age_years" = "age",
                "gender_binary" = "sex") %>% 
  dplyr::select(-c("ids", "diagnosis", "definition", "biopsy_tissue_type")) %>% 
  dplyr::filter(!is.na(vital_status_binary)) %>% 
  dplyr::mutate(diag_date = str_extract(DIAGDATO, "\\d{4}-\\d{2}-\\d{2}")) 

survival_data_wrangled[1, "diag_date"] <- "XXXX-XX-XX"

# Add survival times to wrangled survival data tibble
survival_data_wrangled <- survival_data_wrangled %>% 
  dplyr::mutate("time_years" = case_when(vital_status_binary == 0 ~ 
                                           as.numeric(difftime("XXXX-XX-XX", diag_date, units = "days"))/365,
                                         vital_status_binary == 1 ~ 
                                           as.numeric(difftime(`date of death`, diag_date, units = "days"))/365))

# Define candidate genes
subtype_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                   "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                   "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                   "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                   "ENSG00000227706",  "ENSG00000271394")


### Analyze data --------------------------------------------------------------

## Survival analysis

## Test for Cox proportionality assumption before doing Cox regression analysis

# For each gene i
# test the Cox proportionality assumption
cox_ph_test <- map(subtype_genes, function(i) {
  
  cox_reg_ph(exp_data = voom_filt_data$E, gene = i, clinical_data = survival_data_wrangled)
  
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
  
  cox_reg_gene_set <- cox_reg_univariate(exp_data = voom_filt_data$E, 
                                         gene = i, clinical_data = survival_data_wrangled) 
  
}) %>% 
  bind_rows() %>% 
  # Correct p-values from univariate Cox regression for multiple testing using fdr method
  dplyr::mutate(fdr = p.adjust(`Pr(>|z|)`, method = "fdr"))


## Find genes with a significant effect on survival from univariate Cox
## regression 
genes_for_multi <- cox_reg_uni_results %>% 
  dplyr::filter(fdr < 0.05) %>% 
  pull(gene)


### Save data --------------------------------------------------------------

# Save results of proportional hazards assumption test
save(cox_ph_test, file = "survival/iCOPE_candidate_genes_cox_ph.rda")

# Save results of Cox univariate regression analysis
write_csv(cox_reg_uni_results, file = "survival/iCOPE_candidate_genes_cox_univariate.csv")




