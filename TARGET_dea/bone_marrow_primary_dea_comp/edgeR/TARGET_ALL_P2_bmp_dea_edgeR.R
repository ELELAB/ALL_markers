### This script performs DEA of TARGET-ALL-P2 gene expression data using the edgeR pipeline
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates

# Load libraries
library(edgeR)
library(tidyverse)
library(SummarizedExperiment)
library(biomaRt)

# Source functions
source("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/TARGET_ALL_P2_bmp_dea_edgeR_functions.R")

## Load data ----------------------

# Load filtered data (limma voom)
ALL_P2_filt <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda"))

# Load expression data (used for deseq2, this version of edgeR)
exp_data <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))

# Load batch factor data
ALL_P2_batch_info <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")

## Wrangle data ---------------------------------------------------------------------------------------------------------

## Filter expression data to include only prefiltered genes
exp_data <- assays(exp_data[rownames(exp_data) %in% rownames(ALL_P2_filt),])$"HTSeq - Counts"

# Vector of batch factor designs
batch_factor_designs <- list.dirs("TARGET_dea/bone_marrow_primary_dea_comp/edgeR", 
                                  full.names = FALSE, 
                                  recursive = FALSE)

# Reorder batch info
ALL_P2_batch_info <- ALL_P2_batch_info[match(colnames(exp_data), ALL_P2_batch_info$barcodes), ]

# Define group 
group <- ALL_P2_batch_info$subtype_abbr

# Create DGEList with a grouping variable
y <- DGEList(counts = exp_data, group = group)

# Filter genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# Calculation of norm factors
y <- calcNormFactors(y)

## Create design matrices
design_matrix_list <- map(batch_factor_designs, function(batch_factor_designs) {
  
  create_design_matrix(batch_factor = batch_factor_designs, 
                       condition_types = ALL_P2_batch_info$subtype_abbr, 
                       tissue_portion_aliquot = ALL_P2_batch_info$tissue_portion_aliquot, 
                       years_of_diagnosis = as.character(ALL_P2_batch_info$diagnosis_year))
  
})

# Assign names to resulting list of design matrices
names(design_matrix_list) <- batch_factor_designs

## DEA analysis -------------------------------------------------------------------------------------------

## Do DEA using edgeR for each DEA design differing in included batch factors
DEA_edgeR_list <- map(design_matrix_list, function(design_matrix_list){
  
  edgeR_DEA(design_matrix = design_matrix_list,
            dge_list = y)
})

## Convert ENSEMBL IDs to gene names in DEA matrix

DEA_edgeR_list_gene <- map(seq(DEA_edgeR_list), function(x) {
  
  # Convert full DEA matrix to tibble
  DEA_tbl <- DEA_edgeR_list[[x]][[1]] %>%
    as_tibble(rownames = "ENSEMBL_ID")
  
  # Convert ENSEMBL IDs to gene names and add gene names to full DEA matrix
  ensembl_hugo_conversion(ensembl_id = DEA_tbl$ENSEMBL_ID) %>%
    as_tibble() %>%
    dplyr::rename("ENSEMBL_ID" = `Gene stable ID`) %>%
    full_join(x = .,
              y = DEA_tbl, 
              by = "ENSEMBL_ID")
  
})
names(DEA_edgeR_list_gene) <- batch_factor_designs



## Find significant DEGs based on FDR 
DEA_significant_DEG <- map(DEA_edgeR_list_gene, function(DEA_edgeR_list_gene){
  
  DEA_edgeR_list_gene %>%
    filter(FDR < 0.05)
})


## Split significantly DEGs into up- and downregulated DEGs

DEG_sig_up <- map(DEA_significant_DEG, function(DEA_significant_DEG) {
  
  DEA_significant_DEG %>%
    filter(logFC > 0) %>%
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})

DEG_sig_down <- map(DEA_significant_DEG, function(DEA_significant_DEG) {
  
  DEA_significant_DEG %>%
    filter(logFC < 0) %>%
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})


## Save data ---------------------------------------------------------------------------------------------------------

# Save DEA data for all DEA designs

save_data <- map(seq(batch_factor_designs), function(x) {
  
  # Save DEA results
  saveRDS(object = DEA_edgeR_list_gene[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_full_", batch_factor_designs[x], ".rds"))
  
  # Save DEA of significant DEGs
  saveRDS(object = DEA_significant_DEG[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
  # Save list of significant DEGs
  write_csv(x = DEA_significant_DEG[[x]] %>%
              dplyr::select(ENSEMBL_ID, `Gene name`),
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_", batch_factor_designs[x], ".csv"))
  
  # Save list of significantly upregulated DEGs
  write_csv(x = DEG_sig_up[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
  # Save list of significantly downregulated DEGs
  write_csv(x = DEG_sig_down[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
})




