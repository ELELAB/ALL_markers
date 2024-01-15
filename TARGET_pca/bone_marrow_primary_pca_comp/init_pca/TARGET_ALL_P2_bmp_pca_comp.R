### This script performs PCA on TARGET-ALL-P2 expression data 
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered 
### 6) voom transformed
### The MDSPlot function from the CAMPP pipeline is used for this


# Load libraries
library(tidyverse)
library(patchwork)

# Source functions script
source("TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_pca_functions.R")


## Load data ------------------------------------------------------------------

# Load in voom transformed data 
ALL_P2_bmp_voom <- get(load("TARGET_transform/voom_transform_bmp_comp/TARGET_ALL_P2_bmp_voom.rda"))

# Load table containing information about samples
ALL_P2_info <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv") %>%
  mutate(diagnosis_year = as.character(diagnosis_year),
         number_isolation_aliquot = as.character(number_isolation_aliquot))


## Visualize data ------------------------------------------------------------------

## MDS plots

# MDSPlot requires a vector of group IDs that is used for coloring and labeling. 
# The length of this vector should be equal to the number of columns in the exp dataframe.


# Create MDS plot using MDSPlot function from CAMPP pipeline coloring for subtypes
MDS_subtype <- MDSPlot(my.data = ALL_P2_bmp_voom$E, 
                       my.group = ALL_P2_info$subtype_abbr, 
                       my.labels = colnames(ALL_P2_bmp_voom), 
                       my.cols = c("orange", "darkgreen"))
ggsave(filename = "ALL_P2_bmp_subtype_MDS_init.pdf", 
       plot = MDS_subtype, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for year of diagnosis of samples
MDS_year <- MDSPlot(my.data = ALL_P2_bmp_voom$E, 
                    my.group = ALL_P2_info$diagnosis_year, 
                    my.labels = colnames(ALL_P2_bmp_voom), 
                    my.cols = c("#FFEDA0", "#FED975", "#FEC44F", "#FD8D3C", "#FC4E2A", "#E3201B", "#BD1926", "#7F0E26"))
ggsave(filename = "ALL_P2_bmp_year_MDS_init.pdf", 
       plot = MDS_year, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for tissue portion / aliquot
MDS_portion <- MDSPlot(my.data = ALL_P2_bmp_voom$E, 
                       my.group = ALL_P2_info$tissue_portion_aliquot, 
                       my.labels = colnames(ALL_P2_bmp_voom), 
                       my.cols = c("orange", "darkgreen"))
ggsave(filename = "ALL_P2_bmp_portion_MDS_init.pdf", 
       plot = MDS_portion, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for number isolation from aliquot
MDS_isolation <- MDSPlot(my.data = ALL_P2_bmp_voom$E, 
                         my.group = ALL_P2_info$number_isolation_aliquot, 
                         my.labels = colnames(ALL_P2_bmp_voom), 
                         my.cols = c("orange", "darkgreen"))
ggsave(filename = "ALL_P2_bmp_isolation_MDS_init.pdf", 
       plot = MDS_isolation, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/", 
       width = 12, 
       height = 8)

## Combine subtype and batch effect plots side by side 
subtype_year <- MDS_subtype + MDS_year
ggsave(filename = "ALL_P2_bmp_subtype_year_MDS_init.pdf",
       plot = subtype_year, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/",
       width = 12, 
       height = 8)

subtype_portion <- MDS_subtype + MDS_portion
ggsave(filename = "ALL_P2_bmp_subtype_portion_MDS_init.pdf",
       plot = subtype_portion, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/",
       width = 12, 
       height = 8)

subtype_isolation <- MDS_subtype + MDS_isolation
ggsave(filename = "ALL_P2_bmp_subtype_isolation_MDS_init.pdf",
       plot = subtype_isolation, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/init_pca/",
       width = 12, 
       height = 8)








