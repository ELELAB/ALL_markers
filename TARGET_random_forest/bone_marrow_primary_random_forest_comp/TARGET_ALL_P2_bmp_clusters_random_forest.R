### This script performs variable selection using random forest on TARGET ALL P2 gene expression data
### The random forest variable selection is performed within the predicted subgroups as found using
### unsupervised clustering and the Cola framework


# Load libraries
library(tidyverse)
library(randomForest)
library(varSelRF)
library(caret)


# Source functions script
source("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_random_forest_functions.R")


### Load data ---------------------------------------------------------------------------------------------------------

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))

# Load table containing batch information about samples
ALL_P2_info <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")

# Load predicted cluster/subgroup labels of all samples
subgroups_pred <- read_csv("TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_samples_clusters.csv")


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Binarize the two predicted subgroups within the B-cell ALL subgroup
group_B <- subgroups_pred %>%
  dplyr::filter(subtype_abbr == "B-cell ALL") %>%
  dplyr::select(class) %>%
  dplyr::mutate(class = case_when(class == "1" ~ 1,
                                  class == "4" ~ 2)) %>%
  pull %>%
  as.factor()

# Binarize the two predicted subgroups within the T-cell ALL subgroup
group_T <- subgroups_pred %>%
  dplyr::filter(subtype_abbr == "T-cell ALL") %>%
  dplyr::select(class) %>%
  dplyr::mutate(class = case_when(class == "2" ~ 1,
                                  class == "3" ~ 2)) %>%
  pull %>%
  as.factor()

# Subset batch corrected data to contain only the predicted clusters 
# within the actual B-ALL subgroup
B_samples <- subgroups_pred %>%
  dplyr::filter(subtype_abbr == "B-cell ALL") %>%
  pull(barcodes)
ALL_P2_corr_B <- ALL_P2_corr[,B_samples]

# Subset batch corrected data to contain only the predicted clusters 
# within the actual T-ALL subgroup
T_samples <- subgroups_pred %>%
  dplyr::filter(subtype_abbr == "T-cell ALL") %>%
  pull(barcodes)
ALL_P2_corr_T <- ALL_P2_corr[,T_samples]


# Run random forest feature selection on the batch corrected data containing only B-ALL samples
rf_ALL_B <- RunRF(data = ALL_P2_corr_B, 
                  group = group_B, 
                  split.size = 20, 
                  test.train.ratio = 0.25, 
                  num.trees.init = 5000, 
                  num.trees.iterat = 2000)

# Run random forest feature selection on the batch corrected data containing only T-ALL samples
rf_ALL_T <- RunRF(data = ALL_P2_corr_T, 
                  group = group_T, 
                  split.size = 20, 
                  test.train.ratio = 0.25, 
                  num.trees.init = 5000, 
                  num.trees.iterat = 2000)



### Save data ---------------------------------------------------------------------------------------------------------

# Save random forest performed on B-cell ALL subgroups
save(rf_ALL_B, 
     file = "TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_B_clusters_random_forest_results.rda")

# Save random forest performed on T-cell ALL subgroups
save(rf_ALL_T, 
     file = "TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_T_clusters_random_forest_results.rda")


