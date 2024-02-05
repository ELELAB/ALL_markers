### This script performs enrichment analysis on consensus DEGs
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered
### 6) subject to DEA


# Load libraries
library(tidyverse)
library(enrichR)
library(patchwork)

# Source functions
source("TARGET_enrichment/bone_marrow_primary_enrichment_comp/TARGET_ALL_P2_bmp_enrichment_functions.R")


## Load data ---------------------------------------------------------------------------------------------------------

# Load consensus DEA table
DEA_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))


## Wrangle data ---------------------------------------------------------------------------------------------------------

# Get vector of consensus DEGs
DEG_consensus <- DEA_consensus$Gene_name

# Get vector of upregulated consensus DEGs
DEG_consensus_up <- DEA_consensus %>% 
  dplyr::filter(logFC_mean > 0) %>%
  dplyr::select(Gene_name) %>%
  pull

# Get vector of downregulated consensus DEGs
DEG_consensus_down <- DEA_consensus %>% 
  dplyr::filter(logFC_mean < 0) %>%
  dplyr::select(Gene_name) %>%
  pull

# Create vector of databases to use for enrichment analysis 
dbs <- c("GO_Molecular_Function_2021",
         "GO_Biological_Process_2021",
         "MSigDB_Hallmark_2020")

# Perform enrichment analysis on consensus DEGs
enrich_analysis <- enrichr(genes = DEG_consensus, 
                           databases = dbs)

# Perform enrichment analysis on upregulated consensus DEGs
enrich_analysis_up <- enrichr(genes = DEG_consensus_up, 
                              databases = dbs)

# Perform enrichment analysis on downregulated consensus DEGs
enrich_analysis_down <- enrichr(genes = DEG_consensus_down, 
                                databases = dbs)


## Visualize data ---------------------------------------------------------------------------------------------------------

## Visualize enrichment analysis of consensus DEGs

# Visualize enrichment analysis results of GO molecular function terms for consensus DEGs
DEG_plot_mf <- goplot(data = enrich_analysis[[1]], 
                      title = paste0(str_replace_all(names(enrich_analysis)[[1]],
                                                     "_",
                                                     " "),
                                     " consensus DEGs"), 
                      top = 10)

# Visualize enrichment analysis results of GO biological processes terms for consensus DEGs
DEG_plot_bp <- goplot(data = enrich_analysis[[2]], 
                      title = paste0(str_replace_all(names(enrich_analysis)[[2]],
                                                     "_",
                                                     " "),
                                     " consensus DEGs"), 
                      top = 10)

# Visualize enrichment analysis results of MSigDB hallmark terms for consensus DEGs
DEG_plot_msig <- goplot(data = enrich_analysis[[3]], 
                        title = paste0(str_replace_all(names(enrich_analysis)[[3]],
                                                       "_",
                                                       " "),
                                       " consensus DEGs"), 
                        top = 10)

# Plot the above three plots in one plot
DEG_enrich_all <- DEG_plot_mf + DEG_plot_bp + DEG_plot_msig


## Visualize enrichment analysis of upregulated consensus DEGs

# Visualize enrichment analysis results of GO molecular function terms for upregulated consensus DEGs
DEG_up_plot_mf <- goplot(data = enrich_analysis_up[[1]], 
                         title = paste0(str_replace_all(names(enrich_analysis_up)[[1]],
                                                        "_",
                                                        " "),
                                        " upregulated consensus DEGs"), 
                         top = 10)

# Visualize enrichment analysis results of GO biological processes terms for consensus DEGs
DEG_up_plot_bp <- goplot(data = enrich_analysis_up[[2]], 
                         title = paste0(str_replace_all(names(enrich_analysis_up)[[2]],
                                                        "_",
                                                        " "),
                                        " upregulated consensus DEGs"), 
                         top = 10)

# Visualize enrichment analysis results of MSigDB hallmark terms for consensus DEGs
DEG_up_plot_msig <- goplot(data = enrich_analysis_up[[3]], 
                           title = paste0(str_replace_all(names(enrich_analysis_up)[[3]],
                                                          "_",
                                                          " "),
                                          " upregulated consensus DEGs"), 
                           top = 10)

# Plot the above three plots in one plot
DEG_up_enrich_all <- DEG_up_plot_mf + DEG_up_plot_bp + DEG_up_plot_msig


## Visualize enrichment analysis of downregulated consensus DEGs

# Visualize enrichment analysis results of GO molecular function terms for downregulated consensus DEGs


# Visualize enrichment analysis results of GO biological processes terms for downregulated consensus DEGs
DEG_down_plot_bp <- goplot(data = enrich_analysis_down[[2]], 
                           title = paste0(str_replace_all(names(enrich_analysis_down)[[2]],
                                                          "_",
                                                          " "),
                                          " downregulated consensus DEGs"), 
                           top = 10)

# Visualize all enrichment plots of up- and downregulated consensus DEGs together
DEG_up_down_enrich_plot <- (DEG_up_plot_mf + DEG_up_plot_bp) + (DEG_up_plot_msig / DEG_down_plot_bp)


## Save data ---------------------------------------------------------------------------------------------------------

# Save list containing enrichment analysis results of consensus DEGs for all specified databases 
saveRDS(enrich_analysis, file = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/TARGET_ALL_P2_bmp_enrichment_analyses_DEG.rds")

# Save list containing enrichment analysis results of upregulated consensus DEGs for all specified databases
saveRDS(enrich_analysis_up, file = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/TARGET_ALL_P2_bmp_enrichment_analyses_DEG_up.rds")

# Save list containing enrichment analysis results of downregulated consensus DEGs for all specified databases
saveRDS(enrich_analysis_down, file = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/TARGET_ALL_P2_bmp_enrichment_analyses_DEG_down.rds")

# Save plot of enrichment analysis results of consensus DEGs
ggsave(filename = "TARGET_ALL_P2_bmp_consensus_DEG_enrich.pdf",  
       plot = DEG_enrich_all, 
       path = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/", 
       width = 30, 
       height = 15)

# Save plot of enrichment analysis results of upregulated consensus DEGs
ggsave(filename = "TARGET_ALL_P2_bmp_consensus_DEG_up_enrich.pdf",  
       plot = DEG_up_enrich_all, 
       path = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/", 
       width = 30, 
       height = 15)

# Save plot of enrichment analysis results of downregulated consensus DEGs
ggsave(filename = "TARGET_ALL_P2_bmp_consensus_DEG_down_enrich.pdf",  
       plot = DEG_down_plot_bp, 
       path = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/", 
       width = 8, 
       height = 10)

# Save plot of enrichment analysis results of up- and downregulated consensus DEGs
ggsave(filename = "TARGET_ALL_P2_bmp_consensus_DEG_up_down_enrich.pdf",  
       plot = DEG_up_down_enrich_plot, 
       path = "TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/", 
       width = 30, 
       height = 15)




