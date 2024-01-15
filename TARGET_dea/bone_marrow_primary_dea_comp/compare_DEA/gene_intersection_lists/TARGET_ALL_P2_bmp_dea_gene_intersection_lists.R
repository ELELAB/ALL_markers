### This script creates gene name intersection lists DEA results of TARGET-ALL-P2 gene expression data from three different DEA pipelines
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered


# Load libraries
library(tidyverse)

## Load data ---------------------------------------------------------------------------------------------------------

# Vector of batch factor designs
batch_factor_designs <- list.dirs("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom", 
                                  full.names = FALSE, 
                                  recursive = FALSE)


# Load DEA of limma-voom pipeline of all batch designs
limma_voom_DEA <- map(seq(batch_factor_designs), function(x) {
  
  # Read DEA of significant DEGs
  readRDS(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
})
names(limma_voom_DEA) <- batch_factor_designs

# Load DEA of edgeR pipeline of all batch designs
edger_DEA <- map(seq(batch_factor_designs), function(x) {
  
  # Read DEA of significant DEGs
  readRDS(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
})
names(edger_DEA) <- batch_factor_designs

# Load DEA of DESeq2 pipeline of all batch designs
deseq2_DEA <- map(seq(batch_factor_designs), function(x) {
  
  # Read DEA of significant DEGs
  readRDS(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
})
names(deseq2_DEA) <- batch_factor_designs

# Load list of significantly upregulated DEGs for all DEA designs performed using limma-voom
DEGs_up_limma_voom <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly upregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_up_limma_voom) <- str_c(batch_factor_designs, "_up")

# Load list of significantly downregulated DEGs for all DEA designs performed using limma-voom
DEGs_down_limma_voom <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly downregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_down_limma_voom) <- str_c(batch_factor_designs, "_down")

# Load list of significantly upregulated DEGs for all DEA designs performed using edgeR
DEGs_up_edger <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly upregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
})
names(DEGs_up_edger) <- str_c(batch_factor_designs, "_up")

# Load list of significantly downregulated DEGs for all DEA designs performed using edgeR
DEGs_down_edger <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly downregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
})
names(DEGs_down_edger) <- str_c(batch_factor_designs, "_down")

# Load list of significantly upregulated DEGs for all DEA designs performed using DESeq2
DEGs_up_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly upregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
})
names(DEGs_up_deseq2) <- str_c(batch_factor_designs, "_up")

# Load list of significantly downregulated DEGs for all DEA designs performed using DESeq2
DEGs_down_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly downregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
})
names(DEGs_down_deseq2) <- str_c(batch_factor_designs, "_down")


## Wrangle data ---------------------------------------------------------------------------------------------------------

## Select only gene names for significant gene lists

# Select only ENSEMBL_ID and `Gene name` for limma voom significant genes
limma_voom_DEA <- map(limma_voom_DEA, function(limma_voom_DEA) {
  
  limma_voom_DEA %>% 
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})

# Select only ENSEMBL_ID and `Gene name` for edgeR significant genes
edger_DEA <- map(edger_DEA, function(edger_DEA) {
  
  edger_DEA %>% 
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})

# Select only ENSEMBL_ID and `Gene name` for limma voom significant genes
deseq2_DEA <- map(deseq2_DEA, function(deseq2_DEA) {
  
  deseq2_DEA %>% 
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})

## Intersections for all significant genes -------------------------

# Intersection of all three pipelines
DEGs_all_pipelines <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- limma_voom_DEA[[x]] %>%
    inner_join(edger_DEA[[x]]) %>%
    inner_join(deseq2_DEA[[x]])
})
names(DEGs_all_pipelines) <- batch_factor_designs

# Intersection for limma voom and edgeR
DEGs_limma_edger <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- limma_voom_DEA[[x]] %>%
    inner_join(edger_DEA[[x]])
})
names(DEGs_limma_edger) <- batch_factor_designs

# Intersection for edgeR and deseq2
DEGs_edger_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- edger_DEA[[x]] %>%
    inner_join(deseq2_DEA[[x]])
})
names(DEGs_edger_deseq2) <- batch_factor_designs

# Intersection for deseq2 and limma voom
DEGs_deseq2_limma <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- deseq2_DEA[[x]] %>%
    inner_join(limma_voom_DEA[[x]])
})
names(DEGs_deseq2_limma) <- batch_factor_designs


## Intersections for upregulated significant genes -------------------------
DEGs_up_all_pipelines <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_up_limma_voom[[x]] %>%
    inner_join(DEGs_up_edger[[x]]) %>%
    inner_join(DEGs_up_deseq2[[x]])
})
names(DEGs_up_all_pipelines) <- batch_factor_designs

# Intersection for limma voom and edgeR
DEGs_up_limma_edger <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_up_limma_voom[[x]] %>%
    inner_join(DEGs_up_edger[[x]])
})
names(DEGs_up_limma_edger) <- batch_factor_designs

# Intersection for edgeR and deseq2
DEGs_up_edger_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_up_edger[[x]] %>%
    inner_join(DEGs_up_deseq2[[x]])
})
names(DEGs_up_edger_deseq2) <- batch_factor_designs

# Intersection for deseq2 and limma voom
DEGs_up_deseq2_limma <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_up_deseq2[[x]] %>%
    inner_join(DEGs_up_limma_voom[[x]])
})
names(DEGs_up_deseq2_limma) <- batch_factor_designs


## Intersections for upregulated significant genes -------------------------
DEGs_down_all_pipelines <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_down_limma_voom[[x]] %>%
    inner_join(DEGs_down_edger[[x]]) %>%
    inner_join(DEGs_down_deseq2[[x]])
})
names(DEGs_down_all_pipelines) <- batch_factor_designs

# Intersection for limma voom and edgeR
DEGs_down_limma_edger <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_down_limma_voom[[x]] %>%
    inner_join(DEGs_down_edger[[x]])
})
names(DEGs_down_limma_edger) <- batch_factor_designs

# Intersection for edgeR and deseq2
DEGs_down_edger_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_down_edger[[x]] %>%
    inner_join(DEGs_down_deseq2[[x]])
})
names(DEGs_down_edger_deseq2) <- batch_factor_designs

# Intersection for deseq2 and limma voom
DEGs_down_deseq2_limma <- map(seq(batch_factor_designs), function(x) {
  
  DEG_tbl <- DEGs_down_deseq2[[x]] %>%
    inner_join(DEGs_down_limma_voom[[x]])
})
names(DEGs_down_deseq2_limma) <- batch_factor_designs


## Countings

# Counts of all significant genes
num_DEG_sig <- tibble("Design" = batch_factor_designs,
                      "DEG_all_pipelines" = map(DEGs_all_pipelines, nrow) %>% 
                        unlist,
                      "DEG_limma_edgeR" = map(DEGs_limma_edger, nrow) %>%
                        unlist,
                      "DEG_edgeR_DESeq2" = map(DEGs_edger_deseq2, nrow) %>%
                        unlist,
                      "DEG_DESeq2_limma" = map(DEGs_deseq2_limma, nrow) %>%
                        unlist)

# Counts of all upregulated significant genes
num_DEG_sig_up <- tibble("Design" = batch_factor_designs,
                      "DEG_up_all_pipelines" = map(DEGs_up_all_pipelines, nrow) %>% 
                        unlist,
                      "DEG_up_limma_edgeR" = map(DEGs_up_limma_edger, nrow) %>%
                        unlist,
                      "DEG_up_edgeR_DESeq2" = map(DEGs_up_edger_deseq2, nrow) %>%
                        unlist,
                      "DEG_up_DESeq2_limma" = map(DEGs_up_deseq2_limma, nrow) %>%
                        unlist)

# Counts of all downregulated significant genes
num_DEG_sig_down <- tibble("Design" = batch_factor_designs,
                      "DEG_down_all_pipelines" = map(DEGs_down_all_pipelines, nrow) %>% 
                        unlist,
                      "DEG_down_limma_edgeR" = map(DEGs_down_limma_edger, nrow) %>%
                        unlist,
                      "DEG_down_edgeR_DESeq2" = map(DEGs_down_edger_deseq2, nrow) %>%
                        unlist,
                      "DEG_down_DESeq2_limma" = map(DEGs_down_deseq2_limma, nrow) %>%
                        unlist)


## Save results -------------------------------------------

# Save gene intersection lists
save_data <- map(seq(batch_factor_designs), function(x) {
  
  # Save list of all common significant DEGs
  write_csv(x = DEGs_all_pipelines[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_all_dea_pipelines_", batch_factor_designs[x], ".csv"))
  
  # Save list of common significant DEGs for limma voom and edgeR
  write_csv(x = DEGs_limma_edger[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_limma_edger_", batch_factor_designs[x], ".csv"))
  
  # Save list of common significant DEGs for edgeR and deseq2
  write_csv(x = DEGs_edger_deseq2[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_edger_deseq2_", batch_factor_designs[x], ".csv"))
  
  # Save list of common significant DEGs for deseq2 and limma voom
  write_csv(x = DEGs_deseq2_limma[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_deseq2_limma_", batch_factor_designs[x], ".csv"))
  
  # -------
  
  # Save list of all common upregulated significant DEGs
  write_csv(x = DEGs_up_all_pipelines[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_all_dea_pipelines_", batch_factor_designs[x], ".csv"))
  
  # Save list of common upregulated significant DEGs for limma voom and edgeR
  write_csv(x = DEGs_up_limma_edger[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_limma_edger_", batch_factor_designs[x], ".csv"))
  
  # Save list of common upregulated significant DEGs for edgeR and deseq2
  write_csv(x = DEGs_up_edger_deseq2[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_edger_deseq2_", batch_factor_designs[x], ".csv"))
  
  # Save list of common upregulated ignificant DEGs for deseq2 and limma voom
  write_csv(x = DEGs_up_deseq2_limma[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_deseq2_limma_", batch_factor_designs[x], ".csv"))
  
  # -------
  
  # Save list of all common downregulated significant DEGs
  write_csv(x = DEGs_down_all_pipelines[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_all_dea_pipelines_", batch_factor_designs[x], ".csv"))
  
  # Save list of common downregulated significant DEGs for limma voom and edgeR
  write_csv(x = DEGs_down_limma_edger[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_limma_edger_", batch_factor_designs[x], ".csv"))
  
  # Save list of common downregulated significant DEGs for edgeR and deseq2
  write_csv(x = DEGs_down_edger_deseq2[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_edger_deseq2_", batch_factor_designs[x], ".csv"))
  
  # Save list of common downregulated ignificant DEGs for deseq2 and limma voom
  write_csv(x = DEGs_down_deseq2_limma[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_deseq2_limma_", batch_factor_designs[x], ".csv"))
  
})


# Save count matrices
write_csv(num_DEG_sig, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/ALL_P2_bmp_numbers_sig.csv")
write_csv(num_DEG_sig_up, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/ALL_P2_bmp_numbers_sig_up.csv")
write_csv(num_DEG_sig_down, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/gene_intersection_lists/ALL_P2_bmp_numbers_sig_down.csv")






