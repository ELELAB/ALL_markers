### This script runs all analyses associated with the discovery of gene
### expression markers in two acute lymphoblastic leukemia subtypes. The data
### that is analyzed is the TARGET-ALL-P2 gene expression project from the
### TARGET database.


# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Download data from OSF --------------------------------------------------
library("osfr")
node <- osf_retrieve_node("kgfpv")
files <- osf_ls_files(node)
osf_download(files, 
             recurse = TRUE, 
             conflicts = "skip")

# Restore packages from lockfile ------------------------------------------
library(renv)
renv::settings$external.libraries(c("env_ALL/lib/R/library/"))
renv::restore()
renv::activate()
detach("package:renv", 
       unload = TRUE)

# Set up for Cairo rendering, for headless machines -----------------------
options(bitmapType = "cairo")


## Run scripts ------------------------------------------------------------


## GET DATA ---------------------------------------------------------------

# Get TARGET-ALL-P2 gene expression data, investigate metadata and subset data 
# to contain only bone marrow samples
source("TARGET_data/get_TARGET_data.R")

# -------------------------------------------------------------------------


## REPLICATES -------------------------------------------------------------

# Investigate replicate samples in data
source("TARGET_replicates/TARGET_ALL_P2_bmp_replicates.R")

# Adjust replicates in data
source("TARGET_replicates/TARGET_ALL_P2_bmp_replicates_adj.R")

# -------------------------------------------------------------------------


## PROCESSING -------------------------------------------------------------

# Process data by removing outliers (preprocessing), normalizing for GC content 
# and library size (normalization) and removing lowly expressed genes (filtering)
source("TARGET_processing/bone_marrow_primary_comp/TARGET_process_bmp_comp.R")

# Investigate lost genes following normalization step
source("TARGET_processing/bone_marrow_primary_comp/TARGET_norm_genes_bmp_comp.R")

# Visualize gene types of lost genes following normalization step
source("TARGET_processing/bone_marrow_primary_comp/TARGET_norm_genes_plot_bmp_comp.R")

# --------------------------------------------------------------------------


## VOOM TRANSFORMATION  ----------------------------------------------------

# Voom transform raw data
source("TARGET_transform/voom_transform_bmp_raw/TARGET_ALL_P2_bmp_raw_voom.R")

# Voom transform filtered data
source("TARGET_transform/voom_transform_bmp_comp/TARGET_ALL_P2_bmp_voom_comp.R")

# --------------------------------------------------------------------------










