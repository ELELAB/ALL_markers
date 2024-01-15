### This script performs PCA on TARGET-ALL-P2 expression data 
### This data has been 
### subsetted to contain only bone marrow (primary blood derived cancer) samples
### voom transformed
### The MDSPlot function from the CAMPP pipeline is used for this


# Load libraries
library(tidyverse)
library(SummarizedExperiment)
library(patchwork)
library(ggrepel)

# Source functions script
source("TARGET_pca/bone_marrow_primary_pca_raw/TARGET_ALL_P2_bmp_raw_pca_functions.R")


## Load data ------------------------------------------------------------------

# Load in gene expression data of TARGET-ALL-P2 bone marrow (primary) samples
ALL_P2_bone_marrow_primary_exp_data <- get(load("TARGET_data/ALL-P2_bone_marrow_primary_exp_data.rda"))

# Load in voom transformed raw read counts of TARGET-ALL-P2 bone marrow (primary) samples
ALL_P2_bone_marrow_primary_raw_voom <- get(load("TARGET_transform/voom_transform_bmp_raw/TARGET_ALL_P2_bmp_raw_voom.rda"))

# Load clinical data of TARGET-ALL-P2 project 
clinical_ALL_P2 <- get(load("TARGET_data/TARGET_ALL_P2_clinical_data.rda"))

# Load information about replicates 
ALL_P2_bmp_rep_info <- read_csv("TARGET_replicates/TARGET_ALL_P2_bmp_rep_read_count.csv")


## Wrangle data ------------------------------------------------------------------

# Create tibble containing 
# barcodes 
# subtype 
# year of diagnosis 
# tissue portion
# number isolation from aliquot
# replicate info
ALL_P2_info <- tibble(
  barcodes = colnames(ALL_P2_bone_marrow_primary_raw_voom), 
  patient_barcodes = str_split(barcodes, "-") %>%
    map(function(x) str_c(x[1:3], collapse = "-")) %>%
    unlist,
  diagnosis_year = clinical_ALL_P2$year_of_diagnosis[match(patient_barcodes, clinical_ALL_P2$submitter_id)],
  subtype_full = ALL_P2_bone_marrow_primary_exp_data$primary_diagnosis[match(colnames(ALL_P2_bone_marrow_primary_raw_voom), 
                                                                             ALL_P2_bone_marrow_primary_exp_data$barcode)],
  subtype_abbr = case_when(subtype_full == "Precursor B-cell lymphoblastic leukemia" ~ "B-cell ALL",
                           subtype_full == "T lymphoblastic leukemia/lymphoma" ~ "T-cell ALL"),
  tissue_portion_aliquot = str_split(barcodes, "-") %>%
    map(function(x) x[4]) %>%
    str_split("") %>%
    map(function(x) x[3]) %>%
    unlist,
  number_isolation_aliquot = str_split(barcodes, "-") %>%
    map(function(x) x[5]) %>%
    str_split("") %>%
    map(function(x) x[2]) %>%
    unlist,
  aliquot = str_split(barcodes, "-") %>% 
    map(function(x) x[4:5]) %>% 
    map(function(x) str_c(x, collapse = "-")) %>% 
    unlist,
  replicate_groups = case_when(barcodes %in% ALL_P2_bmp_rep_info$barcode ~ patient_barcodes,
                               TRUE ~ "Non replicate")
)
write_csv(ALL_P2_info, 
          file = "TARGET_pca/bone_marrow_primary_pca_raw/TARGET_ALL_P2_bmp_batch_info.csv")

# Subset raw read counts to contain only replicates
ALL_P2_bmp_raw_rep <- ALL_P2_bone_marrow_primary_raw_voom %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL_names") %>% 
  select("ENSEMBL_names", ALL_P2_bmp_rep_info$barcode) %>% 
  column_to_rownames(var = "ENSEMBL_names")

# Subset tibble containing information about samples to include only replicates
ALL_P2_info_rep <- ALL_P2_info %>%
  filter(barcodes %in% ALL_P2_bmp_rep_info$barcode) %>%
  arrange(patient_barcodes)


## Visualize data ------------------------------------------------------------------

## MDS plots

# MDSPlot requires a vector of group IDs that is used for coloring and labeling. 
# The length of this vector should be equal to the number of columns in the exp dataframe.


# Create MDS plot using MDSPlot function from CAMPP pipeline coloring for subtypes
MDS_subtype <- MDSPlot(my.data = ALL_P2_bone_marrow_primary_raw_voom$E, 
                       my.group = ALL_P2_info$subtype_abbr, 
                       my.labels = colnames(ALL_P2_bone_marrow_primary_raw_voom), 
                       my.cols = c("orange", "darkgreen"))
ggsave(filename = "TARGET_ALL_P2_bmp_raw_subtype_MDS.pdf", 
       plot = MDS_subtype, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for year of diagnosis of samples
MDS_year <- MDSPlot(my.data = ALL_P2_bone_marrow_primary_raw_voom$E, 
                    my.group = ALL_P2_info$diagnosis_year, 
                    my.labels = colnames(ALL_P2_bone_marrow_primary_raw_voom), 
                    my.cols = c("#FFEDA0", "#FED975", "#FEC44F", "#FD8D3C", "#FC4E2A", "#E3201B", "#BD1926", "#7F0E26"))
ggsave(filename = "TARGET_ALL_P2_bmp_raw_year_MDS.pdf", 
       plot = MDS_year, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for tissue portion / aliquot
MDS_portion <- MDSPlot(my.data = ALL_P2_bone_marrow_primary_raw_voom$E, 
                       my.group = ALL_P2_info$tissue_portion_aliquot, 
                       my.labels = colnames(ALL_P2_bone_marrow_primary_raw_voom), 
                       my.cols = c("orange", "darkgreen"))
ggsave(filename = "TARGET_ALL_P2_bmp_raw_portion_MDS.pdf", 
       plot = MDS_portion, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for number isolation from aliquot
MDS_isolation <- MDSPlot(my.data = ALL_P2_bone_marrow_primary_raw_voom$E, 
                         my.group = ALL_P2_info$number_isolation_aliquot, 
                         my.labels = colnames(ALL_P2_bone_marrow_primary_raw_voom), 
                         my.cols = c("orange", "darkgreen"))
ggsave(filename = "TARGET_ALL_P2_bmp_raw_isolation_MDS.pdf", 
       plot = MDS_isolation, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP pipeline) coloring for replicate samples
MDS_raw <- MDSPlot_transparent(my.data = ALL_P2_bone_marrow_primary_raw_voom$E,
                               my.group = ALL_P2_info$replicate_groups,
                               my.labels = colnames(ALL_P2_bone_marrow_primary_raw_voom),
                               my.cols = c("gray60", 
                                           "#E97D72", 
                                           "#C99532",
                                           "#97A934", 
                                           "#53B74C", 
                                           "#57BEA1",
                                           "#53B6DF", 
                                           "#6E9BF8", 
                                           "#CD78F4", 
                                           "#ED6CBF")) +
  theme(legend.position = "right")
ggsave(filename = "TARGET_ALL_P2_bmp_raw_MDS.pdf", 
       plot = MDS_raw, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function (from CAMPP) pipeline of only replicate samples 
MDS_raw_rep <- MDSPlot(my.data = ALL_P2_bmp_raw_rep,
                       my.group = ALL_P2_info_rep$replicate_groups,
                       my.labels = colnames(ALL_P2_bmp_raw_rep),
                       my.cols = c("#E97D72", 
                                   "#C99532",
                                   "#97A934", 
                                   "#53B74C", 
                                   "#57BEA1",
                                   "#53B6DF", 
                                   "#6E9BF8", 
                                   "#CD78F4", 
                                   "#ED6CBF")) +
  theme(legend.position = "right") +
  geom_text_repel(mapping = aes(label = ALL_P2_info_rep$aliquot))
ggsave(filename = "TARGET_ALL_P2_bmp_raw_rep_MDS.pdf", 
       plot = MDS_raw_rep, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/", 
       width = 12, 
       height = 8)


## Combine subtype and batch effect plots side by side 

subtype_year <- MDS_subtype + MDS_year
ggsave(filename = "TARGET_ALL_P2_bmp_raw_subtype_year_MDS.pdf",
       plot = subtype_year, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/",
       width = 12, 
       height = 8)

subtype_portion <- MDS_subtype + MDS_portion
ggsave(filename = "TARGET_ALL_P2_bmp_raw_subtype_portion_MDS.pdf",
       plot = subtype_portion, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/",
       width = 12, 
       height = 8)

subtype_isolation <- MDS_subtype + MDS_isolation
ggsave(filename = "TARGET_ALL_P2_bmp_raw_subtype_isolation_MDS.pdf",
       plot = subtype_isolation, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/",
       width = 12, 
       height = 8)

## Combine replicate plots side by side
rep_MDS <- (MDS_raw_rep + 
              theme(legend.position = "none")) + MDS_raw
ggsave(filename = "TARGET_ALL_P2_bmp_raw_rep_side_MDS.pdf",
       plot = rep_MDS, 
       path = "TARGET_pca/bone_marrow_primary_pca_raw/",
       width = 18, 
       height = 8)




