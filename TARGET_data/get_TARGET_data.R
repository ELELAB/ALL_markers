### This script loads in the TARGET data and saves information about the data ###

# Load libraries
library(tidyverse)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(patchwork)

# Source functions script
source("TARGET_data/get_TARGET_data_functions.R")


## Load data ------------------------------------------------------------

# Load gene expression data of project TARGET-ALL-P2
ALL_P2_exp_data <- get(load("TARGET_data/ALL-P2_exp_data.rda"))


## Wrangle data of TARGET-ALL-P2 project -------------------------------------

# Get clinical data of TARGET-ALL-P2 project
TARGET_ALL_P2_clinical <- GDCquery_clinic(project = "TARGET-ALL-P2",
                                          type = "clinical")

# Extract information about gene expression data of TARGET-ALL-P2 project
info_exp_data <- get_meta_data(SE_exp_data = ALL_P2_exp_data, 
                               var_vector = c("gender", 
                                              "primary_diagnosis", 
                                              "vital_status"))

# Get age of patients and convert age from days to years
ALL_P2_age <- tibble("age_years_since_birth" = ALL_P2_exp_data$age_at_diagnosis / 365)

# Convert column data from summarized experiment object into tibble
ALL_P2_col_data <- as_tibble(ALL_P2_exp_data@colData@listData)

# Count number of samples belonging to combinations of subtype (primary_diagnosis) and type of tissue source (definition)
# and combine the two columns (primary_diagnosis and definition) into one column
count_subtype_tissue <- dplyr::count(ALL_P2_col_data, primary_diagnosis, definition) %>%
  unite(combination, primary_diagnosis, definition, sep = " \n ") %>%
  dplyr::rename(count = n)


## Wrangle data of TARGET-ALL-P2 project (create summaries) ------------------

## 1. Create tables of meta data/variables available in clinical data of TARGET-ALL-P2 project

# Convert clinical data to tibble while repairing duplicate column names
tbl_clin_ALL_P2 <- as_tibble(TARGET_ALL_P2_clinical, 
                             .name_repair = "unique")

# Find number of NA values for each variable (i.e. column name)
number_NA_values <- tbl_clin_ALL_P2 %>% 
  map_int(~ sum(is.na(.))) %>% 
  as_tibble() %>%
  dplyr::rename(Number_of_NA_values = value)

# Find number of unique values for each variable (i.e. column name)
number_unique_values <- tbl_clin_ALL_P2 %>% 
  map_int(n_distinct) %>% 
  as_tibble() %>%
  dplyr::rename(Number_of_unique_values = value)

# Create tibble to contain summary of NA and unique values
summary_clin_ALL_P2 <- tibble(
  "Variable" = colnames(tbl_clin_ALL_P2),
  number_NA_values,
  number_unique_values
)

# Add if a variable/column name contains only NA values 
# Add unique values for variables
summary_clin_ALL_P2 <- summary_clin_ALL_P2 %>%
  mutate("All_NA_values" = case_when(Number_of_NA_values == nrow(tbl_clin_ALL_P2) ~ "All values in variable are NA",
                                     Number_of_NA_values != nrow(tbl_clin_ALL_P2) ~ "All values in variable are not NA"),
         map_df(tbl_clin_ALL_P2, ~tibble("Unique_values" = toString(unique(.)))))

## 2. Create tables of meta data/variables available in SummarizedExperiment object of TARGET-ALL-P2 data

# Find number of NA values for each variable (i.e. column name)
number_NA_values_exp <- ALL_P2_col_data %>% 
  map_int(~ sum(is.na(.))) %>% 
  as_tibble() %>%
  dplyr::rename(Number_of_NA_values = value)

# Find number of unique values for each variable (i.e. column name)
number_unique_values_exp <- ALL_P2_col_data %>% 
  map_int(n_distinct) %>% 
  as_tibble() %>%
  dplyr::rename(Number_of_unique_values = value)

# Create tibble to contain summary of NA and unique values
summary_exp_ALL_P2 <- tibble(
  "Variable" = colnames(ALL_P2_col_data),
  number_NA_values_exp,
  number_unique_values_exp
)

# Add if a variable/column name contains only NA values 
# Add unique values for variables
summary_exp_ALL_P2 <- summary_exp_ALL_P2 %>%
  mutate("All_NA_values" = case_when(Number_of_NA_values == nrow(ALL_P2_col_data) ~ "All values in variable are NA",
                                     Number_of_NA_values != nrow(ALL_P2_col_data) ~ "All values in variable are not NA"),
         map_df(ALL_P2_col_data, ~tibble("Unique_values" = toString(unique(.)))))


## Visualize data of TARGET-ALL-P2 project ------------------------------------

# Plot the age distribution of TARGET-ALL-P2 patients
# Age_at_diagnosis refers to age at the time of diagnosis expressed in number of days since birth
ALL_P2_age_p <- ggplot(ALL_P2_age, 
                       aes(x = age_years_since_birth)) + 
  geom_histogram(binwidth = 1,
                 boundary = 0,
                 color = "black", 
                 fill = "darkgreen") +
  scale_x_continuous(breaks = seq(from = 0, 
                                  to = 30, 
                                  by = 2)) + 
  scale_y_continuous(breaks = seq(from = 0, 
                                  to = 60, 
                                  by = 5)) + 
  theme_light() + 
  xlab("Age in years since birth") + 
  ylab("Number of patients") +
  ggtitle("Age distribution of TARGET-ALL-P2 patients expressed in years since birth") + 
  theme(panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(size = 16,
                                   color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16,
                                   color = "black"),
        plot.title = element_text(size = 16, 
                                  hjust = 0.5)) 

# Plot number of samples belonging to different combinations of subtype (primary_diagnosis) and tissue source (definition)
ALL_P2_subtype_tissue_p <- ggplot(data = count_subtype_tissue, 
                                  mapping = aes(x = combination, 
                                                y = count)) +
  geom_bar(stat = "identity", 
           fill = "darkgreen") +
  geom_text(mapping = aes(label = count), 
            hjust = 1,
            size = 6,
            color = "black") +
  theme_light() + 
  xlab("") +
  ylab("Number of samples") +
  ggtitle("Number of samples in TARGET-ALL-P2 project \n of combinations of subtype and tissue source") + 
  theme(panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(size = 14,
                                   color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14,
                                   color = "black"),
        plot.title = element_text(size = 14, 
                                  hjust = 0.5)) +
  coord_flip()


## Wrangle data subsetted to bone marrow primary samples -------------------

# Subset gene expression data to contain only bone marrow primary samples 
bone_marrow_primary_exp_data <- subset(ALL_P2_exp_data,
                                       select = (definition == "Primary Blood Derived Cancer - Bone Marrow"))

# Obtain raw read counts of gene expression data subsetted to contain only bone marrow (primary) samples
bmp_raw_read_counts <- assay(bone_marrow_primary_exp_data)

# Extract information about bone marrow primary data 
info_bmp_exp_data <- get_meta_data(SE_exp_data = bone_marrow_primary_exp_data, 
                                   var_vector = c("gender", 
                                                  "primary_diagnosis", 
                                                  "vital_status"))

# Get age of patients and convert age from days to years
ALL_P2_bmp_age <- tibble("age_years_since_birth" = bone_marrow_primary_exp_data$age_at_diagnosis / 365)


## Visualize data subsetted to bone marrow primary samples -------------------

# Plot age distribution of patients subsetted to contain only bone marrow primary data
ALL_P2_bmp_age_p <- ggplot(ALL_P2_bmp_age, 
                           aes(x = age_years_since_birth)) + 
  geom_histogram(binwidth = 1,
                 boundary = 0,
                 color = "black", 
                 fill = "darkgreen") +
  scale_x_continuous(breaks = seq(from = 0, 
                                  to = 30, 
                                  by = 2)) + 
  scale_y_continuous(breaks = seq(from = 0, 
                                  to = 60, 
                                  by = 5)) + 
  theme_light() + 
  xlab("Age in years since birth") + 
  ylab("Number of patients") +
  ggtitle("Age distribution of TARGET-ALL-P2 patients expressed in years since birth") + 
  theme(panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_line(size = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(size = 16,
                                   color = "black"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16,
                                   color = "black"),
        plot.title = element_text(size = 16, 
                                  hjust = 0.5)) 

# Plot side by side distribution of tissue source/subtype and age of bone marrow primary data
ALL_P2_age_subtype_tissue <- ALL_P2_subtype_tissue_p + ALL_P2_bmp_age_p 


## Save data ------------------------------------------------------------

# Save clinical data of TARGET-ALL-P2 project
save(TARGET_ALL_P2_clinical, 
     file = "TARGET_data/TARGET_ALL_P2_clinical_data.rda")

# Save headers of clinical data of TARGET-ALL-P2 project
write_csv(data.frame(colnames(TARGET_ALL_P2_clinical)), 
          file = "TARGET_data/TARGET_ALL_P2_clinical_data_headers.csv")

# Save summary table of meta data/variables available in clinical data of TARGET-ALL-P2 project
write_csv(summary_clin_ALL_P2, "TARGET_data/TARGET_ALL_P2_clinical_data_summary.csv", col_names = TRUE)

# Save headers of colData in SummarizedExperiment of downloaded expression data of TARGET-ALL-P2 project
write_csv(data.frame(names(ALL_P2_exp_data@colData@listData)), 
          file = "TARGET_data/TARGET_ALL_P2_exp_data_headers.csv")

# Save summary table tables of meta data/variables available in SummarizedExperiment object of TARGET-ALL-P2 data
write_csv(summary_exp_ALL_P2, "TARGET_data/TARGET_ALL_P2_exp_data_summary.csv", col_names = TRUE)

# Save tibble containing information about TARGET-ALL-P2 project
write_csv(x = info_exp_data, 
          file = "TARGET_data/TARGET_ALL_P2_exp_data_info.csv") 

# Save tibble containing ages of TARGET-ALL-P2 patients
write_csv(ALL_P2_age, 
          file = "TARGET_data/TARGET_ALL_P2_exp_data_age_table.csv") 

# Save table of number of tissue source and subtypes as csv file
write_csv(x = count_subtype_tissue, 
          file = "TARGET_data/TARGET_ALL_P2_exp_data_subtype_tissue.csv")

# Save gene expression data of bone marrow samples (primary blood derived cancer)
save(bone_marrow_primary_exp_data, 
     file = "TARGET_data/ALL-P2_bone_marrow_primary_exp_data.rda")

# Save raw read counts of gene expression data subsetted to contain only bone marrow (primary) samples
save(bmp_raw_read_counts, file = "TARGET_data/ALL_P2_bmp_raw_read_counts.rda")

# Save table containing information about number of samples
write_csv(x = info_bmp_exp_data, 
          file = "TARGET_data/TARGET_ALL_P2_bmp_exp_data_info.csv")

# Save tibble containing ages of TARGET-ALL-P2 patients
write_csv(ALL_P2_bmp_age, 
          file = "TARGET_data/TARGET_ALL_P2_bmp_exp_data_age_table.csv")

# Save plot of histogram of ages of TARGET-ALL-P2 patients
ggsave(filename =  "TARGET_ALL_P2_exp_data_age_plot.pdf",
       plot = ALL_P2_age_p, 
       path = "TARGET_data/", 
       width = 12, 
       height = 8)

# Save plot of combinations of samples by subtype and tissue source
ggsave(filename = "TARGET_ALL_P2_exp_data_subtype_tissue_plot.pdf",
       plot = ALL_P2_subtype_tissue_p, 
       path = "TARGET_data/", 
       width = 12, 
       height = 8)

# Save plot of histogram of ages of TARGET-ALL-P2 patients in bone marrow
# primary subsetted data
ggsave(filename = "TARGET_ALL_P2_bmp_exp_data_age_plot.pdf", 
       plot = ALL_P2_bmp_age_p, 
       path = "TARGET_data/", 
       width = 12, 
       height = 8)

# Save plot containing the above two plots side by side
ggsave(filename = "TARGET_ALL_P2_bmp_exp_data_age_subtype_tissue_plot.pdf", 
       plot = ALL_P2_age_subtype_tissue, 
       path = "TARGET_data/", 
       width = 20, 
       height = 8)

