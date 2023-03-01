### This script corrects for batch factors on the TARGET-ALL-P2 expression data 
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered
### 6) voom transformed


# Load libraries
library(TCGAbiolinks)
library(tidyverse)
library(sva)

# Source functions
source("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_batch_functions.R")


## Load data ---------------------------------------------------------------------------------------------------------

# Load in voom transformed data
ALL_P2_bmp_voom <- get(load("TARGET_transform/voom_transform_bmp_comp/TARGET_ALL_P2_bmp_voom.rda")) 

# Load table containing batch information about samples
ALL_P2_batch <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")


## Wrangle data --------------------------------------------------------------------------------------------------------

## Extract information about batches

# Number of samples in each year 
total_samples_year <- sample_number(batchTable = ALL_P2_batch, 
                                    groupNum = 1, 
                                    batchFactor = diagnosis_year)

# Number of samples of each subtype (condition) in each year
total_samples_subtype_year <- sample_number(batchTable = ALL_P2_batch, 
                                            groupNum = 2, 
                                            batchFactor = diagnosis_year, 
                                            condition = subtype_full)

# Number of samples in each tissue portion/aliquot
total_samples_portion <- sample_number(batchTable = ALL_P2_batch, 
                                       groupNum = 1, 
                                       batchFactor = tissue_portion_aliquot)

# Number of samples of each subtype (condition) in each portion
total_samples_subtype_portion <- sample_number(batchTable = ALL_P2_batch, 
                                               groupNum = 2, 
                                               batchFactor = tissue_portion_aliquot, 
                                               condition = subtype_full)

## Apply batch correction 

# Correct expression data for year of diagnosis batch factor
corrected_counts_year <- batch_correction(batchTable = ALL_P2_batch, 
                                          barcodesCol = barcodes, 
                                          batchFactor = diagnosis_year, 
                                          subtype_group = as.factor(ALL_P2_batch$subtype_abbr), 
                                          filePath = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_combat_year_plots.pdf", 
                                          plotWidth = 12, 
                                          plotHeight = 8, 
                                          expData = ALL_P2_bmp_voom$E)

# Correct expression data for tissue portion/aliquot batch factor 
# corrected_counts_portion <- batch_correction(batchTable = ALL_P2_batch, 
#                                             barcodesCol = barcodes, 
#                                             batchFactor = tissue_portion_aliquot, 
#                                             filePath = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_combat_portion_plots.pdf", 
#                                             plotWidth = 12, 
#                                             plotHeight = 8, 
#                                             expData = ALL_P2_bmp_voom$E)


## Save data ---------------------------------------------------------------------------------------------------------

# Save tables containing number of samples for condition (subtype) and batch factor (year and portion)
write_csv(total_samples_year, 
          file = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_samples_year.csv")
write_csv(total_samples_subtype_year, 
          file = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_samples_year_subtype.csv")
write_csv(total_samples_portion, 
          file = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_samples_portion.csv")
write_csv(total_samples_subtype_portion, 
          file = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_samples_portion_subtype.csv")

# Save batch effect (year) corrected counts
save(corrected_counts_year, file = "TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda")








  

