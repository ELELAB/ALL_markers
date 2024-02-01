### This script performs processing on gene count matrix that is output from
### snakemake gencode based pipeline


# Source functions
source("processing/TCGAanalyze_Normalization_revised.R")

## Load libraries -------------------------------------------------------------------

library(TCGAbiolinks) #2.25.3
library(SummarizedExperiment) #1.28.0


## Load data ------------------------------------------------------------------------

# Load in SummarizedExperiment object of iCOPE ALL gene expression data
iCOPE_ALL_gene_exp_se <- get(load("get_data/iCOPE_ALL_gene_exp_SE.rda"))


## Process data ------------------------------------------------------------------------

# Preprocessing step using TCGAbiolinks
prep_data <- TCGAanalyze_Preprocessing(object = iCOPE_ALL_gene_exp_se, 
                                       cor.cut = 0.6)

# Normalization step using TCGAbiolinks
norm_data <- TCGAanalyze_Normalization_revised(tabDF = prep_data, 
                                               geneInfo = geneInfoHT, 
                                               method = "gcContent")

# Filtering step using TCGAbiolinks
filt_data <- TCGAanalyze_Filtering(tabDF = norm_data, 
                                   method = "quantile",
                                   qnt.cut = 0.25)


## Save data ------------------------------------------------------------------------

# Save preprocessed data
save(prep_data, file = "processing/iCOPE_ALL_preprocessed_data.rda")

# Save normalized data
save(norm_data, file = "processing/iCOPE_ALL_normalized_data.rda")

# Save filtered data 
save(filt_data, file = "processing/iCOPE_ALL_filtered_data.rda")




