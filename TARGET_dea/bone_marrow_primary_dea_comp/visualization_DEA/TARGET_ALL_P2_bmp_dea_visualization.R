### This script visualizes DEA results of TARGET-ALL-P2 gene expression data 
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered


# Load libraries
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(TCGAbiolinks)
library(VennDiagram)
library(viridis)
library(biomaRt)

# Source functions
source("TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/TARGET_ALL_P2_bmp_dea_visualization_functions.R")


## Load data ---------------------------------------------------------------------------------------------------------

# Vector of batch factor designs
batch_factor_designs <- list.dirs("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom", 
                                  full.names = FALSE, 
                                  recursive = FALSE)

# Load list of significantly DEGs for all DEA designs performed using limma-voom
DEGs_limma_voom <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_limma_voom) <- batch_factor_designs

# Load list of significantly DEGs for all DEA designs performed using edgeR
DEGs_edger <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_", batch_factor_designs[x], ".csv"))
  
})
names(DEGs_edger) <- batch_factor_designs

# Load list of significantly DEGs for all DEA designs performed using deseq2
DEGs_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_deseq2) <- batch_factor_designs

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


## Wrangle data - for UpSet plots ---------------------------------------------------------------------------------------------------------

# Combine up- and downregulated DEGs identified using limma-voom into one list
DEGs_up_down_limma_voom_comb <- c(DEGs_up_limma_voom %>%
                                    map(., function(x) { x[[1]] }),
                                  DEGs_down_limma_voom %>%
                                    map(., function(x) { x[[1]] }))

# Combine up- and downregulated DEGs identified using edgeR into one list
DEGs_up_down_edger_comb <- c(DEGs_up_edger %>%
                               map(., function(x) { x[[1]] }),
                             DEGs_down_edger %>%
                               map(., function(x) { x[[1]] }))

# Combine up- and downregulated DEGs identified using DESeq2 into one list
DEGs_up_down_deseq2_comb <- c(DEGs_up_deseq2 %>%
                                map(., function(x) { x[[1]] }),
                              DEGs_down_deseq2 %>%
                                map(., function(x) { x[[1]] }))

# Combine DEGs identified using limma-voom, edgeR and DESeq2 using no batch design into one list
DEGs_all_no_batch_comb <- list("limma_voom_no_batch" = DEGs_limma_voom$no_batch$ENSEMBL_ID, 
                               "edger_no_batch" = DEGs_edger$no_batch$ENSEMBL_ID, 
                               "deseq2_no_batch" = DEGs_deseq2$no_batch$ENSEMBL_ID)

# Combine DEGs identified using limma-voom, edgeR and DESeq2 using portion as batch factor into one list
DEGs_all_portion_batch_comb <- list("limma_voom_portion_batch" = DEGs_limma_voom$portion_batch$ENSEMBL_ID, 
                                    "edger_portion_batch" = DEGs_edger$portion_batch$ENSEMBL_ID, 
                                    "deseq2_portion_batch" = DEGs_deseq2$portion_batch$ENSEMBL_ID)

# Combine DEGs identified using limma-voom, edgeR and DESeq2 using portion and years of diagnosis as batch factors into one list
DEGs_all_portion_years_batch_comb <- list("limma_voom_portion_years_batch" = DEGs_limma_voom$portion_years_batch$ENSEMBL_ID,
                                          "edger_portion_years_batch" = DEGs_edger$portion_years_batch$ENSEMBL_ID,
                                          "deseq2_portion_years_batch" = DEGs_deseq2$portion_years_batch$ENSEMBL_ID)

# Combine DEGs identified using limma-voom, edgeR and DESeq2 using years of diagnosis as batch factors into one list
DEGs_all_years_batch_comb <- list("limma_voom_years_batch" = DEGs_limma_voom$years_batch$ENSEMBL_ID,
                                  "edger_years_batch" = DEGs_edger$years_batch$ENSEMBL_ID,
                                  "deseq2_years_batch" = DEGs_deseq2$years_batch$ENSEMBL_ID)

# Combine up- and downregulated DEGs identified using limma-voom, edgeR and DESeq2 and no batch factor design into one list
DEGs_up_down_all_no_batch_comb <- list("limma_voom_up_no_batch" = DEGs_up_limma_voom$no_batch_up$ENSEMBL_ID,
                                       "edger_up_no_batch" = DEGs_up_edger$no_batch_up$ENSEMBL_ID,
                                       "deseq2_up_no_batch" = DEGs_up_deseq2$no_batch_up$ENSEMBL_ID,
                                       "limma_voom_down_no_batch" = DEGs_down_limma_voom$no_batch_down$ENSEMBL_ID,
                                       "edger_down_no_batch" = DEGs_down_edger$no_batch_down$ENSEMBL_ID,
                                       "deseq2_down_no_batch" = DEGs_down_deseq2$no_batch_down$ENSEMBL_ID)

# Combine up- and downregulated DEGs identified using limma-voom, edgeR and DESeq2 and portion as batch factor into one list
DEGs_up_down_all_portion_batch_comb <- list("limma_voom_up_portion_batch" = DEGs_up_limma_voom$portion_batch_up$ENSEMBL_ID,
                                            "edger_up_portion_batch" = DEGs_up_edger$portion_batch_up$ENSEMBL_ID,
                                            "deseq2_up_portion_batch" = DEGs_up_deseq2$portion_batch_up$ENSEMBL_ID,
                                            "limma_voom_down_portion_batch" = DEGs_down_limma_voom$portion_batch_down$ENSEMBL_ID,
                                            "edger_down_portion_batch" = DEGs_down_edger$portion_batch_down$ENSEMBL_ID,
                                            "deseq2_down_portion_batch" = DEGs_down_deseq2$portion_batch_down$ENSEMBL_ID)

# Combine up- and downregulated DEGs identified using limma-voom, edgeR and DESeq2 and portion and years as batch factor into one list
DEGs_up_down_all_portion_years_batch_comb <- list("limma_voom_up_portion_years_batch" = DEGs_up_limma_voom$portion_years_batch_up$ENSEMBL_ID,
                                                  "edger_up_portion_years_batch" = DEGs_up_edger$portion_years_batch_up$ENSEMBL_ID,
                                                  "deseq2_up_portion_years_batch" = DEGs_up_deseq2$portion_years_batch_up$ENSEMBL_ID,
                                                  "limma_voom_down_portion_years_batch" = DEGs_down_limma_voom$portion_years_batch_down$ENSEMBL_ID,
                                                  "edger_down_portion_years_batch" = DEGs_down_edger$portion_years_batch_down$ENSEMBL_ID,
                                                  "deseq2_down_portion_years_batch" = DEGs_down_deseq2$portion_years_batch_down$ENSEMBL_ID)

# Combine up- and downregulated DEGs identified using limma-voom, edgeR and DESeq2 and years as batch factor into one list
DEGs_up_down_all_years_batch_comb <- list("limma_voom_up_years_batch" = DEGs_up_limma_voom$years_batch_up$ENSEMBL_ID,
                                          "edger_up_years_batch" = DEGs_up_edger$years_batch_up$ENSEMBL_ID,
                                          "deseq2_up_years_batch" = DEGs_up_deseq2$years_batch_up$ENSEMBL_ID,
                                          "limma_voom_down_years_batch" = DEGs_down_limma_voom$years_batch_down$ENSEMBL_ID,
                                          "edger_down_years_batch" = DEGs_down_edger$years_batch_down$ENSEMBL_ID,
                                          "deseq2_down_years_batch" = DEGs_down_deseq2$years_batch_down$ENSEMBL_ID)


## Visualize data - UpSet plots ---------------------------------------------------------------------------------------------------------

# UpSet plot of DEGs found using limma-voom pipeline and different batch designs
upset_DEGs_limma_voom <- upset_plot(file_name = "ALL_P2_bmp_DEG_limma_voom_batch_designs_upset.pdf", 
                                    dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/limma_voom", 
                                    sets_list = DEGs_limma_voom %>%
                                      map(., function(x) { x[[1]] }), 
                                    names_sets = names(DEGs_limma_voom), 
                                    title_plot = "Intersections of identified DEGs using different designs in the limma-voom pipeline")

# UpSet plot of DEGs found using edgeR pipeline and different batch designs
upset_DEGs_edger <- upset_plot(file_name = "ALL_P2_bmp_DEG_edger_batch_designs_upset.pdf", 
                               dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/edgeR", 
                               sets_list = DEGs_edger %>%
                                 map(., function(x) { x[[1]] }), 
                               names_sets = names(DEGs_edger), 
                               title_plot = "Intersections of identified DEGs using different designs in the edgeR pipeline")

# UpSet plot of DEGs found using DESeq2 pipeline and different batch designs
upset_DEGs_deseq2 <- upset_plot(file_name = "ALL_P2_bmp_DEG_deseq2_batch_designs_upset.pdf", 
                                dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/deseq2", 
                                sets_list = DEGs_deseq2 %>%
                                  map(., function(x) { x[[1]] }), 
                                names_sets = names(DEGs_deseq2), 
                                title_plot = "Intersections of identified DEGs using different designs in the DESeq2 pipeline")

# UpSet plot of up and down DEGs found using limma-voom pipeline and different batch designs
upset_DEGs_up_down_limma_voom <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_limma_voom_batch_designs_upset.pdf",
                                            dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/limma_voom",
                                            sets_list = DEGs_up_down_limma_voom_comb,
                                            names_sets = names(DEGs_up_down_limma_voom_comb),
                                            title_plot = "Intersections of identified up- and downregulated DEGs using different designs in the limma-voom pipeline")

# UpSet plot of up and down DEGs found using edgeR pipeline and different batch designs
upset_DEGs_up_down_edger <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_edger_batch_designs_upset.pdf",
                                       dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/edgeR",
                                       sets_list = DEGs_up_down_edger_comb,
                                       names_sets = names(DEGs_up_down_edger_comb),
                                       title_plot = "Intersections of identified up- and downregulated DEGs using different designs in the edgeR pipeline")

# UpSet plot of up and down DEGs found using DESeq2 pipeline and different batch designs
upset_DEGs_up_down_deseq2 <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_deseq2_batch_designs_upset.pdf",
                                       dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/deseq2",
                                       sets_list = DEGs_up_down_deseq2_comb,
                                       names_sets = names(DEGs_up_down_deseq2_comb),
                                       title_plot = "Intersections of identified up- and downregulated DEGs using different designs in the DESeq2 pipeline")

# UpSet plot of DEGs found using limma-voom, edgeR and DESeq2 pipelines using no batch design 
upset_DEGs_all_no_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_all_pipelines_no_batch_upset.pdf",
                                      dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                      sets_list = DEGs_all_no_batch_comb,
                                      names_sets = names(DEGs_all_no_batch_comb),
                                      title_plot = "Intersections of identified DEGs using limma-voom, edgeR and DESeq2 and no batch factor design")

# UpSet plot of DEGs found using limma-voom, edgeR and DESeq2 pipelines using portion as batch factor
upset_DEGs_all_portion_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_all_pipelines_portion_batch_upset.pdf",
                                           dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                           sets_list = DEGs_all_portion_batch_comb,
                                           names_sets = names(DEGs_all_portion_batch_comb),
                                           title_plot = "Intersections of identified DEGs using limma-voom, edgeR and DESeq2 and portion as batch factor")

# UpSet plot of DEGs found using limma-voom, edgeR and DESeq2 pipelines using portion and year of diagnosis as batch factor
upset_DEGs_all_portion_years_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_all_pipelines_portion_years_batch_upset.pdf",
                                                 dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                                 sets_list = DEGs_all_portion_years_batch_comb,
                                                 names_sets = names(DEGs_all_portion_years_batch_comb),
                                                 title_plot = "Intersections of identified DEGs using limma-voom, edgeR and DESeq2 and portion and years as batch factor")

# UpSet plot of DEGs found using limma-voom, edgeR and DESeq2 pipelines using year of diagnosis as batch factor
upset_DEGs_all_years_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_all_pipelines_years_batch_upset.pdf",
                                         dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                         sets_list = DEGs_all_years_batch_comb,
                                         names_sets = names(DEGs_all_years_batch_comb),
                                         title_plot = "Intersections of identified DEGs using limma-voom, edgeR and DESeq2 and years as batch factor")

# UpSet plot of up- and downregulated DEGs found using limma-voom, edgeR and DESeq2 pipelines using no batch factor design
upset_DEGs_up_down_all_no_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_all_pipelines_no_batch_upset.pdf",
                                              dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                              sets_list = DEGs_up_down_all_no_batch_comb,
                                              names_sets = names(DEGs_up_down_all_no_batch_comb),
                                              title_plot = "Intersections of identified up- and downregulated DEGs using limma-voom, edgeR and DESeq2 and no batch factor")

# UpSet plot of up- and downregulated DEGs found using limma-voom, edgeR and DESeq2 pipelines using portion as batch factor
upset_DEGs_up_down_all_portion_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_all_pipelines_portion_batch_upset.pdf",
                                                   dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                                   sets_list = DEGs_up_down_all_portion_batch_comb,
                                                   names_sets = names(DEGs_up_down_all_portion_batch_comb),
                                                   title_plot = "Intersections of identified up- and downregulated DEGs using limma-voom, edgeR and DESeq2 and portion as batch factor")

# UpSet plot of up- and downregulated DEGs found using limma-voom, edgeR and DESeq2 pipelines using portion and years as batch factor
upset_DEGs_up_down_all_portion_years_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_all_pipelines_portion_years_batch_upset.pdf",
                                                         dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                                         sets_list = DEGs_up_down_all_portion_years_batch_comb,
                                                         names_sets = names(DEGs_up_down_all_portion_years_batch_comb),
                                                         title_plot = "Intersections of identified up- and downregulated DEGs using limma-voom, edgeR and DESeq2 and portion and years as batch factor")

# UpSet plot of up- and downregulated DEGs found using limma-voom, edgeR and DESeq2 pipelines using years as batch factor
upset_DEGs_up_down_all_years_batch <- upset_plot(file_name = "ALL_P2_bmp_DEG_up_down_all_pipelines_years_batch_upset.pdf",
                                                 dir_output = "TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/all_dea_pipelines",
                                                 sets_list = DEGs_up_down_all_years_batch_comb,
                                                 names_sets = names(DEGs_up_down_all_years_batch_comb),
                                                 title_plot = "Intersections of identified up- and downregulated DEGs using limma-voom, edgeR and DESeq2 and years as batch factor")








