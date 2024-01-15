### This script performs DEA of TARGET-ALL-P2 gene expression data using the DESeq2 pipeline
### The expression data is raw, but filtered to contain genes that are left from the following workflow:
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered

# Load libraries
library(DESeq2)
library(tidyverse)
library(biomaRt)

# Source functions
source("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/TARGET_ALL_P2_bmp_dea_deseq2_functions.R")

## Load data ----------------------

# Load filtered data (used for edgeR/limma voom)
ALL_P2_filt <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda"))

# Load expression data (used for deseq2)
exp_data <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))

# Load batch factor data
batch_factor_data <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")


## Wrangle data ----------------------

## Filter expression data to include only prefiltered genes
exp_data <- exp_data[rownames(exp_data) %in% rownames(ALL_P2_filt),]

# Wrangle batch factor data
batch_factor_data <- batch_factor_data %>%
  dplyr::select(barcodes, diagnosis_year, tissue_portion_aliquot) %>%
  mutate(diagnosis_year = as.factor(diagnosis_year)) %>%
  mutate(tissue_portion_aliquot = as.factor(tissue_portion_aliquot)) %>%
  dplyr::rename(barcode = barcodes)

# Add batch factor data to expression data
colData(exp_data) <- merge(colData(exp_data), batch_factor_data, by = "barcode", sort = FALSE)

# Get batch factor combinations
batch_factor_designs <- list.dirs("TARGET_dea/bone_marrow_primary_dea_comp/deseq2", 
                                  full.names = FALSE, 
                                  recursive = FALSE)

# Set dds dataset designs
dds_design_list <- map(batch_factor_designs, create_DESeqDataSet_design)

# Assign names to  list of designs
names(dds_design_list) <- batch_factor_designs

## Run analysis -------------------

## Run differential expression analysis
DEA_deseq2_list <- map(dds_design_list, function(dds_design_list){
  
  deseq2_DEA(expression_data = exp_data,
             design = dds_design_list,
             lfc_threshold = 1,
             alpha = 0.05)
})

## Convert ENSEMBL IDs to gene names in DEA matrix

DEA_deseq2_list_gene <- map(seq(DEA_deseq2_list), function(x) {
  
  # Convert full DEA matrix to tibble
  DEA_tbl <- DEA_deseq2_list[[x]] %>%
    as_tibble(rownames = "ENSEMBL_ID")
  
  # Convert ENSEMBL IDs to gene names and add gene names to full DEA matrix
  ensembl_hugo_conversion(ensembl_id = DEA_tbl$ENSEMBL_ID) %>%
    as_tibble() %>%
    dplyr::rename("ENSEMBL_ID" = `Gene stable ID`) %>%
    full_join(x = .,
              y = DEA_tbl, 
              by = "ENSEMBL_ID")
  
})
names(DEA_deseq2_list_gene) <- batch_factor_designs


# Filter results: Find significant DEGs based on fdr
DEA_deseq2_filtered_list_gene <- map(DEA_deseq2_list_gene, function(DEA_deseq2_list_gene){
  
  filter_dea_results(dea_result = DEA_deseq2_list_gene,
                     fdr.cut = 0.05)
})

# Split significantly DEGs into up- and downregulated DEGs

DEG_sig_up <- map(DEA_deseq2_filtered_list_gene, function(DEA_filtered){
  index_up <- which(DEA_filtered$log2FoldChange > 0)
  DEA_filtered <- DEA_filtered[index_up,]
})

DEG_sig_down <- map(DEA_deseq2_filtered_list_gene, function(DEA_filtered){
  index_down <- which(DEA_filtered$log2FoldChange < 0)
  DEA_filtered <- DEA_filtered[index_down,]
})



## Save data -------------------

# Save DEA data for all DEA designs

save_data <- map(seq(batch_factor_designs), function(x) {
  
  # Save DEA results
  saveRDS(object = DEA_deseq2_list_gene[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_full_", batch_factor_designs[x], ".rds"))
  
  # Save DEA of significant DEGs
  saveRDS(object = DEA_deseq2_filtered_list_gene[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
  # Save list of significantly DEGs
  write.csv(x = DEA_deseq2_filtered_list_gene[[x]] %>%
              dplyr::select(ENSEMBL_ID, `Gene name`), 
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_", batch_factor_designs[x], ".csv"),
            quote = FALSE,
            row.names = FALSE)
  
  # Save list of significantly upregulated DEGs
  write.csv(x = DEG_sig_up[[x]] %>%
              dplyr::select(ENSEMBL_ID, `Gene name`),
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"),
            quote = FALSE,
            row.names = FALSE)
  
  # Save list of significantly downregulated DEGs
  write.csv(x = DEG_sig_down[[x]] %>%
              dplyr::select(ENSEMBL_ID, `Gene name`),
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"),
            quote = FALSE,
            row.names = FALSE)
  
})

