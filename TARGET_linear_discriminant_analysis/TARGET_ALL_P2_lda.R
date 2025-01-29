# Load libraries
library(GeoDE)
library(tidyverse)

# Load the data
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))
ALL_P2_info <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")

# Wrangle data
sample_class <- ALL_P2_info %>%
  dplyr::mutate(sample_class_enc = as.integer(factor(subtype_abbr, levels = c("B-cell ALL", "T-cell ALL")))) %>% 
  dplyr::select(sample_class_enc) %>% 
  pull %>% 
  as.factor()

ALL_P2_corr_wrl <- ALL_P2_corr %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL_ID")


# Run the analysis
chdir_analysis_results <- chdirAnalysis(datain = ALL_P2_corr_wrl,
                                        sampleclass = sample_class,
                                        CalculateSig = TRUE,
                                        nnull = 10)

# Save results
save(chdir_analysis_results, file = "TARGET_linear_discriminant_analysis/TARGET_ALL_P2_bmp_lda_results.rda")
