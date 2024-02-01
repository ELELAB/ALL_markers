### This script performs batch correction of gene expression data 



## Load libraries -------------------------------------------------------------------
library(TCGAbiolinks) #2.25.3
library(tidyverse) #2.0.0
library(sva) #3.46.0

## Source functions
source("batch_correction/iCOPE_batch_correction_functions.R")


## Load data ------------------------------------------------------------------------

# Load voom transformed filtered data
voom_filt_data <- get(load("transform/iCOPE_ALL_filt_data_voom.rda"))

# Load meta data
sample_meta_data_df <- get(load("get_data/iCOPE_ALL_meta_data.rda"))


## Wrangle data ---------------------------------------------------------------------

# Add year of diagnosis to meta data table
sample_meta_data_year <- sample_meta_data_df %>% 
  mutate(year = str_extract(sample_meta_data_df$date, "\\d{4}"))

# Extract information about batches
diag_year_subtype_count <- sample_meta_data_year %>%
  dplyr::count(year, subtype)

# Apply batch correction correcting expression data for year of diagnosis batch factor
# corrected_counts_year <- batch_correction(batchTable = sample_meta_data_year, 
#                                          barcodesCol = barcode, batchFactor = year, 
#                                          subtype_group = as.factor(sample_meta_data_year$subtype), 
#                                          filePath = "batch_correction/iCOPE_combat_year_plots.pdf", 
#                                          plotWidth = 12, plotHeight = 8, expData = voom_filt_data$E)


## Save data ---------------------------------------------------------------------

# Save overview of number of samples in year of diagnosis batches
write.csv(diag_year_subtype_count, file = "batch_correction/iCOPE_samples_year_subtype.csv")




