### This script investigates technical replicates found in the TARGET-ALL-P2 gene expression data 
### The gene expression data has been subsetted to contain only bone marrow (primary blood derived cancer) samples 


# Load libraries
library(tidyverse)
library(ggrepel)
library(SummarizedExperiment)


## Load data ---------------------------------------------------------------------

# Load in gene expression data of TARGET-ALL-P2 bone marrow (primary) samples
ALL_P2_bone_marrow_primary_exp_data <- get(load("TARGET_data/ALL-P2_bone_marrow_primary_exp_data.rda"))


## Wrangle data ---------------------------------------------------------------------

# Extract patient barcodes and full barcodes of all samples 
barcodes_all_samples <- ALL_P2_bone_marrow_primary_exp_data@colData %>% 
  as_tibble() %>% 
  select(patient, barcode) 

# Calculate total read count of each sample 
# Join with above tibble containing barcode information
# Add aliquot labels 
read_counts_sum <- assay(ALL_P2_bone_marrow_primary_exp_data) %>% 
  as_tibble() %>% 
  summarise_all(sum) %>% 
  pivot_longer(cols = everything(), 
               names_to = "barcode", 
               values_to = "sum_read_counts") %>% 
  inner_join(x = ., 
             y = barcodes_all_samples, 
             by = "barcode") %>% 
  mutate(portion = str_split(barcode, "-") %>% 
           map(function(x) x[4]) %>% 
           str_split("") %>% 
           map(function(x) x[3]) %>% 
           unlist,
         isolation = str_split(barcode, "-") %>% 
           map(function(x) x[5]) %>% 
           str_split("") %>% 
           map(function(x) x[2]) %>% 
           unlist,
         aliquot = str_split(barcode, "-") %>% 
           map(function(x) x[4:5]) %>% 
           map(function(x) str_c(x, collapse = "-")) %>% 
           unlist)

# Find duplicated patients 
duplicate_patients <- barcodes_all_samples %>% 
  group_by(patient) %>% 
  dplyr::count() %>% 
  filter(n > 1) %>% 
  select(patient) %>% 
  ungroup()

# Subset read_counts_sum for duplicate patients 
dup_patients_sum <- inner_join(x = duplicate_patients, 
                               y = read_counts_sum, 
                               by = "patient") 


## Visualize data ---------------------------------------------------------------------

# Plot total sum of read counts in duplicate samples
dup_read_counts <- dup_patients_sum %>% 
  group_by(patient) %>% 
  ggplot(data = ., 
         mapping = aes(x = barcode, 
                       y = sum_read_counts,
                       fill = patient)) +
  geom_col(width = 0.5) +
  scale_y_continuous(breaks = seq(from = 0, 
                                  to = max(dup_patients_sum$sum_read_counts), 
                                  by = 10e+06)) +
  labs(x = "Sample barcode", 
       y = "Total read count", 
       title = "Total read count of TARGET ALL replicates") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 65, 
                                   hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(color = "black"))
dup_read_counts

# Plot total sum of read counts of all samples 
read_counts_plot <- read_counts_sum %>% 
  ggplot(data = ., 
         mapping = aes(x = reorder(barcode, 
                                   -sum_read_counts), 
                       y = sum_read_counts)) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = seq(from = 0, 
                                  to = max(read_counts_sum$sum_read_counts), 
                                  by = 10e+06)) +
  labs(x = "Samples", 
       y = "Total read count", 
       title = "Total read count of TARGET ALL samples") +
  theme_minimal(base_size = 24) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
read_counts_plot

# Highlight duplicate samples in plot
read_counts_duplicates <- read_counts_plot +
  geom_point(data = group_by(dup_patients_sum, 
                             patient),
             mapping = aes(x = reorder(barcode, 
                                       -sum_read_counts),
                           y = sum_read_counts,
                           color = patient),
             size = 5) +
  geom_text_repel(data = group_by(dup_patients_sum,
                                  patient),
                  mapping = aes(label = aliquot),
                  size = 5,
                  vjust = 1,
                  hjust = 0.5)
read_counts_duplicates

# Highlight portion B samples in plot
read_counts_portion_B <- read_counts_plot +
  geom_point(data = read_counts_sum %>% 
               filter(portion == "B"),
             mapping = aes(x = reorder(barcode, 
                                       -sum_read_counts),
                           y = sum_read_counts,
                           color = portion),
             size = 5) +
  scale_color_manual("Portion", 
                     values = c("darkorange"))
read_counts_portion_B

# Highlight in read_counts_portion_B which of portion B samples are replicates
portion_B_rep <- read_counts_sum %>% 
  filter(portion == "B") %>% 
  mutate(is_replicate = case_when(barcode %in% dup_patients_sum$barcode ~ "Is replicate",
                                  TRUE ~ "Non replicate"))
portion_B_rep_plot <- read_counts_plot +
  geom_point(data = portion_B_rep,
             mapping = aes(x = reorder(barcode,
                                       -sum_read_counts),
                           y = sum_read_counts,
                           color = is_replicate),
             size = 3) +
  scale_color_manual("Portion B", 
                     values = c("deepskyblue3", 
                                "darkorange"))
portion_B_rep_plot

# Highlight number isolation from aliquot in plot
read_counts_isolation_2 <- read_counts_plot +
  geom_point(data = read_counts_sum %>% 
               filter(isolation == 2),
             mapping = aes(x = reorder(barcode, 
                                       -sum_read_counts),
                           y = sum_read_counts,
                           color = isolation),
             size = 5) +
  scale_color_manual("Isolation",
                     values = c("darkorange"))
read_counts_isolation_2

# Highlight in read_counts_isolation_2 plot which of the isolation 2 samples are replicates
isolation_2_rep <- read_counts_sum %>% 
  filter(isolation == 2) %>% 
  mutate(is_replicate = case_when(barcode %in% dup_patients_sum$barcode ~ "Is replicate",
                                  TRUE ~ "Non replicate"))
isolation_2_rep_plot <- read_counts_plot +
  geom_point(data = isolation_2_rep,
             mapping = aes(x = reorder(barcode,
                                       -sum_read_counts),
                           y = sum_read_counts,
                           color = is_replicate),
             size = 5) +
  scale_color_manual("Isolation 2", 
                     values = c("deepskyblue3", 
                                "darkorange"))
isolation_2_rep_plot



## Save data ---------------------------------------------------------------------

# Save tibble of total read counts and barcode information of all samples
write_csv(x = read_counts_sum,
          file = "TARGET_replicates/TARGET_ALL_P2_bmp_all_read_count.csv")

# Save list of patients that are replicates
write_csv(x = duplicate_patients, 
          file = "TARGET_replicates/TARGET_ALL_P2_bmp_replicates.csv")

# Save tibble of total read counts and barcode information of replicates
write_csv(x = dup_patients_sum, 
          file = "TARGET_replicates/TARGET_ALL_P2_bmp_rep_read_count.csv")

# Save barplot of total read count of technical replicates 
ggsave(filename = "TARGET_ALL_P2_bmp_rep_read_count.pdf", 
       plot = dup_read_counts, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)

# Save scatter plot of total read count of all samples 
ggsave(filename = "TARGET_ALL_P2_bmp_all_read_count.pdf", 
       plot = read_counts_plot, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)

# Save scatter plot of total read count of all samples, highlighting replicates 
ggsave(filename = "TARGET_ALL_P2_bmp_all_rep_read_count.pdf", 
       plot = read_counts_duplicates, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)

# Save scatter plot of total read count of all samples, highlighting portion B samples 
ggsave(filename = "TARGET_ALL_P2_bmp_portion_B_read_count.pdf", 
       plot = read_counts_portion_B, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)

# Save scatter plot of total read count of all samples, highlighting portion B samples and replicates
ggsave(filename = "TARGET_ALL_P2_bmp_portion_B_rep_read_count.pdf", 
       plot = portion_B_rep_plot, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)

# Save scatter plot of total read count of all samples, highlighting isolation 2 samples
ggsave(filename = "TARGET_ALL_P2_bmp_iso_2_read_count.pdf", 
       plot = read_counts_isolation_2, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)

# Save scatter plot of total read count of all samples, highlighting isolation 2 samples and replicates
ggsave(filename = "TARGET_ALL_P2_bmp_iso_2_rep_read_count.pdf", 
       plot = isolation_2_rep_plot, 
       path = "TARGET_replicates/", 
       width = 12, 
       height = 8)













