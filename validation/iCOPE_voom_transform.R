### This script performs a voom transformation of raw and filtered 
### gene expression data



## Load libraries -------------------------------------------------------------------
library(limma) #3.54.2


## Load data ------------------------------------------------------------------------

# Load filtered expression data 
filt_data <- get(load("processing/iCOPE_ALL_filtered_data.rda"))

# Load raw expression data
raw_data <- get(load("get_data/iCOPE_ALL_gene_counts.rda"))


## Wrangle data ---------------------------------------------------------------------

# Voom transform filtered data
voom_filt_data <- voom(filt_data)

# Voom transform raw data
voom_raw_data <- voom(raw_data)


## Save data ------------------------------------------------------------------

# Save voom transformed filtered data
save(voom_filt_data, file = "transform/iCOPE_ALL_filt_data_voom.rda")

# Save voom transformed raw data
save(voom_raw_data, file = "transform/iCOPE_ALL_raw_data_voom.rda")


