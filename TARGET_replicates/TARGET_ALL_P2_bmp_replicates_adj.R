### This script adjusts for technical replicates found in the TARGET-ALL-P2 gene expression data 
### The gene expression data has been subsetted to contain only bone marrow (primary blood derived cancer) samples 
### The adjustment of replicates is done by selecting the sample for a patient with the highest total read count 


# Load libraries
library(tidyverse)


## Load data ---------------------------------------------------------------------

# Load in gene expression data of TARGET-ALL-P2 bone marrow (primary) samples
ALL_P2_bmp_exp_data <- get(load("TARGET_data/ALL-P2_bone_marrow_primary_exp_data.rda"))

# Load in tibble of patient barcode, ful barcode and total read count of replicates
rep_info <- read_csv(file = "TARGET_replicates/TARGET_ALL_P2_bmp_rep_read_count.csv")


## Wrangle data ---------------------------------------------------------------------

# Group replicates by patient and keep only sample among each patient with highest total read count
keep_rep <- rep_info %>%
  group_by(patient) %>%
  slice_max(order_by = sum_read_counts) %>%
  select(-c(portion,
            isolation,
            aliquot))

# Find replicates that are discarded
discard_rep <- anti_join(x = rep_info, 
                         y = keep_rep, 
                         by = "barcode") %>%
  select(-c(portion,
            isolation,
            aliquot))

# Subset list of barcodes in ALL P2 bmp data to filter out discarded replicate samples 
barcodes <- tibble(barcode = ALL_P2_bmp_exp_data$barcode) %>%
  anti_join(x = ., 
            y = discard_rep, 
            by = "barcode") %>%
  pull(barcode)

# Subset ALL P2 bmp gene expression data to filter out discarded replicate samples
ALL_P2_bmp_exp_data_no_rep <- subset(ALL_P2_bmp_exp_data,
                                     select = (barcode %in% barcodes))


## Save data ---------------------------------------------------------------------

# Save list of kept replicate patients
write_csv(x = keep_rep, 
          file = "TARGET_replicates/TARGET_ALL_P2_bmp_kept_rep.csv")

# Save list of discarded replicate patients
write_csv(x = discard_rep, 
          file = "TARGET_replicates/TARGET_ALL_P2_bmp_discarded_rep.csv")

# Save ALL P2 bone marrow primary gene expression data where replicate samples have been filtered
save(ALL_P2_bmp_exp_data_no_rep, file = "TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda")



