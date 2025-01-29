# Load libraries
library(GeoDE)
library(tidyverse)
library(viridis)
library(UpSetR)
library(ComplexHeatmap)

# Source UpSet plotting function
source("TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_compare_functions.R")

# Load results
chdir_results <- get(load("TARGET_linear_discriminant_analysis/TARGET_ALL_P2_bmp_lda_results.rda"))

# Load consensus DEA table
DEA_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load selected ENSEMBL gene IDs from elastic net logistic regression
elastic_net <- read_csv("TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/TARGET_ALL_P2_bmp_elastic_consensus_biotype.csv")

# Load top 40 ENSEMBL gene IDs found from PCA
pca_results <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes_info.csv")


# Examine the results 
# Get the first 10 most important genes
lapply(chdir_results$results, function(x) x[1:10])

# We can also extract the results of the \code{chdirSig} function
# for example chdir_analysis_example$chdirprops[[1]] gives the whole
# characteristic direction vector for each value of gamma:

# Extract the characteristic direction vector
chdir_vector <- chdir_results$chdirprops[[1]]

# Get the number of significant genes 
chdir_num_sig_genes <- chdir_results$chdirprops$number_sig_genes

# Get the significant genes
chdir_sig_genes <- chdir_results[["results"]][[1]] %>% 
  names()

# Compare significant genes with DEA consensus results, PCA and logistic regression

# Extract ENSEMBL gene IDs from all loaded tables 
DEA_consensus_genes <- DEA_consensus %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

elastic_net_genes <- elastic_net %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

pca_genes <- pca_results %>%
  dplyr::filter(PCA_dim == "PC_1") %>%
  dplyr::select(ENSEMBL_ID) %>%
  pull

# Define predicted markers
cand_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                "ENSG00000227706",  "ENSG00000271394")

# Combine above gene lists into one list
genes_list <- list("chdir" = chdir_sig_genes,
                   "consensus_DEGs" = DEA_consensus_genes,
                   "elastic_net" = elastic_net_genes,
                   "PCA" = pca_genes,
                   "pred_markers" = cand_genes)

# Upset plot of all gene lists
upset_genes <- upset_plot(file_name = "TARGET_ALL_P2_bmp_compare_genes_chdir_upset.pdf", 
                          dir_output = "TARGET_linear_discriminant_analysis", 
                          sets_list = genes_list, 
                          names_sets = names(genes_list), 
                          title_plot = "Intersections of ENSEMBL gene IDs from various analyses")

# Find the intersected genes among all 5 gene lists 
intersect_all <- Reduce(intersect, genes_list)
identical(sort(intersect_all), sort(cand_genes)) # TRUE
cand_genes %in% chdir_sig_genes # All TRUE

