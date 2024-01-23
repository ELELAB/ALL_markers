### This script compares results, which areselected ENSEMBL gene IDs, from various gene discovery analyses:
### DEA, dimensionality reduction (PCA) and elastic net logistic regression which have been performed on 
### TARGET ALL P2 gene expression data 


# Load libraries
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(viridis)


# Source functions
source("TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_compare_functions.R")


### Load data ---------------------------------------------------------------------------------------------------------

# Load consensus DEA table
DEA_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load list of housekeeping genes
DEG_eisenberg_hkg <- read_csv("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/DEG_vs_eisenberg.csv")

# Load selected ENSEMBL gene IDs from elastic net logistic regression
elastic_net <- read_csv("TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/TARGET_ALL_P2_bmp_elastic_consensus_biotype.csv")

# Load top 40 ENSEMBL gene IDs found from PCA
pca_results <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes_info.csv")

# Load cancer genes from NCG
NCG_annotation <- read_tsv("TARGET_compare_genes/bone_marrow_primary_compare_comp/NCG_cancerdrivers_annotation_supporting_evidence.tsv")
NCG_properties <- read_tsv("TARGET_compare_genes/bone_marrow_primary_compare_comp/NCG_cancerdrivers_systemslevelproperties.tsv")


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Extract ENSEMBL gene IDs from all loaded tables 
DEA_consensus_genes <- DEA_consensus %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

DEA_HKG_genes <- DEG_eisenberg_hkg %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

elastic_net_genes <- elastic_net %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

pca_genes <- pca_results %>%
  dplyr::filter(PCA_dim == "PC_1") %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

# Extract gene symbols from all loaded tables 
# This is to compare results of analyses with NCG drivers as NCG drivers
# are on the symbol format
DEA_consensus_genes_symbols <- DEA_consensus %>%
  dplyr::select(Gene_name) %>%
  dplyr::filter(Gene_name != "") %>%
  distinct() %>%
  pull

DEA_HKG_genes_symbols <- DEG_eisenberg_hkg %>%
  dplyr::select(Gene_name.x) %>%
  distinct() %>%
  pull

elastic_net_genes_symbols <- elastic_net %>%
  dplyr::select(gene_name) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  distinct() %>%
  pull

pca_genes_symbols <- pca_results %>%
  dplyr::filter(PCA_dim == "PC_1") %>%
  dplyr::select(gene_name) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  distinct() %>%
  pull

# Extract NCG drivers which are divided into canonical and candidate drivers
NCG_genes_canonical <- NCG_annotation %>%
  dplyr::filter(type == "Canonical Cancer Driver") %>%
  dplyr::select(symbol) %>%
  distinct() %>%
  pull

NCG_genes_candidate <- NCG_annotation %>%
  dplyr::filter(!symbol %in% NCG_genes_canonical) %>%
  dplyr::select(symbol) %>%
  distinct() %>%
  pull


## Combine above gene lists into one list

# For upset plot quantified using ENSEMBL gene IDs (i.e. without NCG drivers)
genes_list <- list("consensus_DEGs" = DEA_consensus_genes,
                   "HKG_DEGs" = DEA_HKG_genes,
                   "elastic_net" = elastic_net_genes,
                   "PCA" = pca_genes)

# For upset plot quantified using gene symbols (i.e. with NCG drivers)
genes_list_symbols <- list("consensus_DEGs" = DEA_consensus_genes_symbols,
                           "HKG_DEGs" = DEA_HKG_genes_symbols,
                           "elastic_net" = elastic_net_genes_symbols,
                           "PCA" = pca_genes_symbols,
                           "NCG_canonical" = NCG_genes_canonical,
                           "NCG_candidate" = NCG_genes_candidate)


## Create lists of ENSEMBL gene IDs that intersect between different combinations of gene lists

# Intersect consensus DEGs with elastic net logistic regression genes
consensus_DEGs_elastic <- inner_join(x = DEA_consensus, 
                                     y = elastic_net,
                                     by = "ENSEMBL_ID")

# Intersect consensus DEGs with PCA genes
consensus_DEGs_pca <- inner_join(x = DEA_consensus, 
                                 y = pca_results %>%
                                   filter(PCA_dim == "PC_1"),
                                 by = "ENSEMBL_ID")

# Intersect elastic net logistic regression genes with PCA genes
elastic_pca <- inner_join(x = elastic_net,
                          y = pca_results %>%
                            filter(PCA_dim == "PC_1"),
                          by = "ENSEMBL_ID")

# Intersect elastic net logistic regression genes with HKG DEGs
elastic_HKG <- inner_join(x = elastic_net,
                          y = DEG_eisenberg_hkg,
                          by = "ENSEMBL_ID")

# Intersect PCA genes with HKG DEGs
pca_HKG <- inner_join(x = pca_results %>%
                        filter(PCA_dim == "PC_1"),
                      y = DEG_eisenberg_hkg,
                      by = "ENSEMBL_ID")

# Intersect consensus DEGs, elastic net genes and PCA genes
consensus_DEGs_elastic_pca <- inner_join(x = DEA_consensus,
                                         y = elastic_net,
                                         by = "ENSEMBL_ID") %>%
  inner_join(x = .,
             y = pca_results %>%
               filter(PCA_dim == "PC_1"),
             by = "ENSEMBL_ID")


## Create lists of gene symbols that intersect between different combinations of gene lists

# Create binary combination matrix
comb_mat_symbols <- make_comb_mat(genes_list_symbols, 
                                  mode = "intersect")

# Extract elements in each combination of intersections 
intersections <- map(comb_name(comb_mat_symbols), 
                     function(x) extract_comb(comb_mat_symbols, x))

# Set name to each intersection
names(intersections) <- comb_name(comb_mat_symbols)
names(intersections) <- c("conDEG_elastic_PCA_NCGcanon",
                          "conDEG_HKG_NCGcanon",
                          "conDEG_HKG_NCGcand",
                          "conDEG_elastic_PCA",
                          "conDEG_elastic_NCGcanon",
                          "conDEG_PCA_NCGcanon",
                          "conDEG_PCA_NCGcand",
                          "elastic_PCA_NCGcanon",
                          "conDEG_HKG",
                          "conDEG_elastic",
                          "conDEG_PCA",
                          "conDEG_NCGcanon",
                          "conDEG_NCGcand",
                          "HKG_NCGcanon",
                          "HKG_NCGcand",
                          "elastic_PCA",
                          "elastic_NCGcanon",
                          "PCA_NCGcanon",
                          "PCA_NCGcand",
                          "conDEG",
                          "HKG",
                          "elastic",
                          "PCA",
                          "NCGcanon",
                          "NCGcand")


### Visualize data ---------------------------------------------------------------------------------------------------------

# Visualize intersections between consensus DEGs, HKG DEGs, elastic net genes and PCA genes in UpSet plot
# Quantified using ENSEMBL gene IDs
upset_genes <- upset_plot(file_name = "TARGET_ALL_P2_bmp_compare_ensembl_ids_upset.pdf", 
                          dir_output = "TARGET_compare_genes/bone_marrow_primary_compare_comp", 
                          sets_list = genes_list, 
                          names_sets = names(genes_list), 
                          title_plot = "Intersections of ENSEMBL gene IDs from various analyses")

# Visualize intersections between consensus DEGs, HKG DEGs, elastic net genes, PCA genes and NCG drivers in UpSet plot
# Quantified using gene symbols
upset_genes_symbols <- upset_plot(file_name = "TARGET_ALL_P2_bmp_compare_symbols_upset.pdf", 
                                  dir_output = "TARGET_compare_genes/bone_marrow_primary_compare_comp", 
                                  sets_list = genes_list_symbols, 
                                  names_sets = names(genes_list_symbols), 
                                  title_plot = "Intersections of gene symbols from various analyses")


### Save data ---------------------------------------------------------------------------------------------------------

## Save intersections of ENSEMBL gene IDs

# Save intersection between consensus DEGs and elastic net logistic regression genes
write_csv(consensus_DEGs_elastic, file = "TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_consensus_DEGs_elastic.csv") 

# Save intersection between consensus DEGs and PCA genes
write_csv(consensus_DEGs_pca, file = "TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_consensus_DEGs_pca.csv")

# Save intersection between elastic net logistic regression genes and PCA genes
write_csv(elastic_pca, file = "TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_elastic_pca.csv")

# Save intersection between elastic net logistic regression genes and HKG DEGs
write_csv(elastic_HKG, file = "TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_elastic_HKG.csv") 

# Save intersection between PCA genes and HKG DEGs
write_csv(pca_HKG, file = "TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_pca_HKG.csv") 

# Save intersection between consensus DEGs, elastic net genes and PCA genes
write_csv(consensus_DEGs_elastic_pca, file = "TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_consensus_DEGs_elastic_pca.csv")

## Save intersections of gene symbols
map(1:length(intersections), function(x) write.csv(intersections[[x]], file = str_c("TARGET_compare_genes/bone_marrow_primary_compare_comp/",names(intersections[x]),"_symbols.csv")))




