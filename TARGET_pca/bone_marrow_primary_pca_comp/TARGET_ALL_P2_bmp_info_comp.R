### This script saves information about TARGET-ALL-P2 expression data 
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered 
### 6) voom transformed


# Load libraries
library(tidyverse)


## Load data ------------------------------------------------------------------

# Load in gene expression data of TARGET-ALL-P2 bone marrow (primary) samples where replicates have been removed
ALL_P2_bmp_exp_data_no_rep <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))

# Load in voom transformed data
ALL_P2_bmp_voom <- get(load("TARGET_transform/voom_transform_bmp_comp/TARGET_ALL_P2_bmp_voom.rda"))

# Load clinical data of TARGET-ALL-P2 project 
clinical_ALL_P2 <- get(load("TARGET_data/TARGET_ALL_P2_clinical_data.rda"))


## Wrangle data ------------------------------------------------------------------

# Create tibble containing barcodes, subtype, year of diagnosis, tissue portion, and number isolation from aliquot
ALL_P2_info <- tibble(
  barcodes = colnames(ALL_P2_bmp_voom), 
  patient_barcodes = str_split(barcodes, "-") %>%
    map(function(x) str_c(x[1:3], collapse = "-")) %>%
    unlist,
  diagnosis_year = clinical_ALL_P2$year_of_diagnosis[match(patient_barcodes, clinical_ALL_P2$submitter_id)],
  subtype_full = ALL_P2_bmp_exp_data_no_rep$primary_diagnosis[match(colnames(ALL_P2_bmp_voom), 
                                                                    ALL_P2_bmp_exp_data_no_rep$barcode)],
  subtype_abbr = case_when(subtype_full == "Precursor B-cell lymphoblastic leukemia" ~ "B-cell ALL",
                           subtype_full == "T lymphoblastic leukemia/lymphoma" ~ "T-cell ALL"),
  tissue_portion_aliquot = str_split(barcodes, "-") %>%
    map(function(x) x[4]) %>%
    str_split("") %>%
    map(function(x) x[3]) %>%
    unlist,
  number_isolation_aliquot = str_split(barcodes, "-") %>%
    map(function(x) x[5]) %>%
    str_split("") %>%
    map(function(x) x[2]) %>%
    unlist
)
write_csv(ALL_P2_info, 
          file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")








