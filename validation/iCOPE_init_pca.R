### This script performs PCA of raw and filtered 
### gene expression data



## Load libraries -------------------------------------------------------------------
library(stats) #4.2.0
library(patchwork) #1.1.3

## Source functions
source("pca/iCOPE_pca_functions.R")


## Load data ------------------------------------------------------------------------

# Load voom transformed raw data
voom_raw_data <- get(load("transform/iCOPE_ALL_raw_data_voom.rda"))

# Load voom transformed filtered data
voom_filt_data <- get(load("transform/iCOPE_ALL_filt_data_voom.rda"))

# Load meta data
sample_meta_data_df <- get(load("get_data/iCOPE_ALL_meta_data.rda"))


## Wrangle data ---------------------------------------------------------------------

# MDSPlot requires a vector of group IDs that is used for coloring and labeling. 
# The length of this vector should be equal to the number of columns in the exp dataframe.

# Get year of samples by extracting the year from the date column
year <- str_extract(sample_meta_data_df$date, "\\d{4}")

# Get vector of subtypes
subtypes <- sample_meta_data_df$subtype


## Visualize data ---------------------------------------------------------------------

## For voom transformed filtered data

# Create MDS plot using MDSPlot function from CAMPP pipeline coloring for subtypes
MDS_subtype_filt <- MDSPlot(my.data = voom_filt_data$E,
                            my.group = subtypes, 
                            my.labels = sample_meta_data_df$ids, 
                            my.cols = c("orange", "darkgreen"))
ggsave(filename = "iCOPE_ALL_filt_subtype_MDS_init.pdf", 
       plot = MDS_subtype_filt, 
       path = "pca/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function from CAMPP pipeline coloring for year
MDS_year_filt <- MDSPlot(my.data = voom_filt_data$E,
                         my.group = year, 
                         my.labels = sample_meta_data_df$ids, 
                         my.cols = c("#FFEDA0", "#FED975", "#FEC44F", "#FD8D3C", 
                                     "#FC4E2A", "#E3201B", "#BD1926", "#7F0E26"))
ggsave(filename = "iCOPE_ALL_filt_year_MDS_init.pdf", 
       plot = MDS_year_filt, 
       path = "pca/", 
       width = 12, 
       height = 8)


## For voom transformed raw data

# Create MDS plot using MDSPlot function from CAMPP pipeline coloring for subtypes
MDS_subtype_raw <- MDSPlot(my.data = voom_raw_data$E,
                           my.group = subtypes, 
                           my.labels = sample_meta_data_df$ids, 
                           my.cols = c("orange", "darkgreen"))
ggsave(filename = "iCOPE_ALL_raw_subtype_MDS_init.pdf", 
       plot = MDS_subtype_raw, 
       path = "pca/", 
       width = 12, 
       height = 8)

# Create MDS plot using MDSPlot function from CAMPP pipeline coloring for year
MDS_year_raw <- MDSPlot(my.data = voom_raw_data$E,
                        my.group = year, 
                        my.labels = sample_meta_data_df$ids, 
                        my.cols = c("#FFEDA0", "#FED975", "#FEC44F", "#FD8D3C", 
                                    "#FC4E2A", "#E3201B", "#BD1926", "#7F0E26"))
ggsave(filename = "iCOPE_ALL_raw_year_MDS_init.pdf", 
       plot = MDS_year_raw, 
       path = "pca/", 
       width = 12, 
       height = 8)

# Plot subtype and year MDS plots of voom transformed filtered data side by side
MDS_subtype_year_filt <- MDS_subtype_filt + MDS_year_filt
ggsave(filename = "iCOPE_ALL_filt_subtype_year_MDS_init.pdf", 
       plot = MDS_subtype_year_filt, 
       path = "pca/", 
       width = 12, 
       height = 8)

# Plot subtype and year MDS plots of voom transformed raw data side by side
MDS_subtype_year_raw <- MDS_subtype_raw + MDS_year_raw
ggsave(filename = "iCOPE_ALL_raw_subtype_year_MDS_init.pdf", 
       plot = MDS_subtype_year_raw, 
       path = "pca/", 
       width = 12, 
       height = 8)



