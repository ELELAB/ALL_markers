### This script loads in and makes ready the gene count matrix that is output from
### snakemake gencode based pipeline



## Load libraries -------------------------------------------------------------------

library(tidyverse) #2.0.0
library(readxl) #1.4.3
library(SummarizedExperiment) #1.28.0
library(GenomicRanges) #1.50.2


## Load data ------------------------------------------------------------------------

# Load in feature counts matrix from GENCODE (ensembl IDs as rows and samples as columns)
gene_count_data <- read_delim("gene_counts/gencode/campp_input_featureCounts_gencode.txt", col_names = TRUE)

# Load in sample metadata from GENCODE (ensembl IDs as rows)
sample_meta_data <- read_delim("gene_counts/gencode/sample_metadata_2.txt", col_names = TRUE)

# Load larger metadata table
large_meta_data <- read_excel("get_data/pheno_for_adrian_566.xlsx")


## Wrangle data ---------------------------------------------------------------------

# Set gene names (in first column) to rownames and convert to matrix
gene_count_df <- gene_count_data %>% 
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

# Truncate ENSEMBL IDs in rownames
rownames(gene_count_df) <- str_replace(string = rownames(gene_count_df), 
                                       pattern = "\\..*", 
                                       replacement = "")

# Replace , in Age variable in sample_meta_data to .
# Copy ids into rownames and into a barcode column
# Create a biopsy/tissue type column using the barcode nomenclature indicating
#   the biopsy/tissue type of patient (3rd component in barcode)
# Create a definition column corresponding to the biopsy/tissue type 
# Add MRD numbers to table
sample_meta_data_df <- sample_meta_data %>%
  mutate(age = str_replace_all(age, ",", "."),
         age = as.numeric(age), 
         ids_rownames = ids,
         barcode = ids,
         biopsy_tissue_type = map_chr(str_split(barcode, "[:punct:]"), 3),
         definition = case_when(biopsy_tissue_type == "G" ~ "Germline, i.e. blood leukocyte DNA",
                                biopsy_tissue_type == "T" ~ "Tumor, i.e. any biopsy of tumor cells",
                                biopsy_tissue_type == "H" ~ "Heel tendon biopsy (post-mortem)",
                                biopsy_tissue_type == "F" ~ "Fibroblast culture from skin",
                                biopsy_tissue_type == "P" ~ "Paraffin-embedded tissue"),
         MRDNR = as.numeric(str_extract(barcode, "^[^_]+"))) %>%
  column_to_rownames(var = "ids_rownames") %>% 
  as_tibble()

# Join two metadata tables
meta_data_join <- large_meta_data %>%
  dplyr::select(FamilyID, MRDNR, DIAGDATO, gender)
meta_data_join <- left_join(x = sample_meta_data_df, y = meta_data_join, by = "MRDNR") 

meta_data_join[1, "DIAGDATO"] <- meta_data_join[1, "date"]$date
meta_data_join[1, "gender"] <- "male"

# Get information about data

# Get number of gender samples
gender_counts <- meta_data_join %>%
  group_by(gender) %>%
  dplyr::count() %>%
  mutate(gender = as.character(gender),
         gender = paste0("sex_", gender)) %>%
  group_split()
gender_1 <- gender_counts[[1]]
gender_2 <- gender_counts[[2]]

# Get number of subtype samples
subtype_counts <- meta_data_join %>%
  group_by(subtype) %>%
  dplyr::count() %>%
  mutate(subtype = paste0("subtype_", subtype)) %>%
  group_split()
subtype_1 <- subtype_counts[[1]] 
subtype_2 <- subtype_counts[[2]]

# Get number of diagnosis samples
diagnosis_counts <- meta_data_join %>%
  group_by(diagnosis) %>%
  dplyr::count() %>%
  mutate(diagnosis = paste0("diagnosis_", diagnosis))

# Get number of biopsy tissue type
biopsy_tissue_counts <- meta_data_join %>%
  group_by(biopsy_tissue_type) %>%
  dplyr::count() %>%
  mutate(biopsy_tissue_type = paste0("biopsy_tissue_type_", biopsy_tissue_type))

# Create tibble containing above information:
# Number of genes, number of samples, number of gender samples, number of subtype samples, number of diagnosis samples, number of biopsy tissue type samples
info_exp_data <- tibble(
  "number_of_genes" = nrow(gene_count_df),
  "number_of_samples" = ncol(gene_count_df),
  !!gender_1$gender := gender_1$n,
  !!gender_2$gender := gender_2$n,
  !!subtype_1$subtype := subtype_1$n,
  !!subtype_2$subtype := subtype_2$n,
  !!diagnosis_counts$diagnosis := diagnosis_counts$n,
  !!biopsy_tissue_counts$biopsy_tissue_type := biopsy_tissue_counts$n
)

# Create SummarizedExperiment object

# Create a rowRanges object needed in SummarizedExperiment
row_data <- GRanges(rep(c("NA", "NA"), c(0, nrow(gene_count_df))),
                    IRanges(c(0)))

# Create a single SummarizedExperiment object from the meta data and gene counts
se_counts_meta <- SummarizedExperiment(gene_count_df, colData = DataFrame(meta_data_join),
                                       rowRanges = row_data)

# Add name to counts data in SE
names(se_counts_meta@assays@data@listData) <- "Counts"


## Visualize data ------------------------------------------------------------------

# Create plot of age distribution of patients 

# File path for saving plot of ages
file_path_ALL_age <- c("get_data/iCOPE_ALL_exp_data_age_plot.pdf")

# Save plot as pdf
pdf(file = file_path_ALL_age, 
    width = 8,
    height = 8)

# Plot age of iCOPE ALL patients in histogram
ggplot(as.data.frame(meta_data_join$age), 
       aes(x = meta_data_join$age)) + 
  geom_histogram(binwidth = 1,
                 boundary = 0,
                 color = "black", 
                 fill = "darkgreen") +
  scale_x_continuous(breaks = seq(0, 30, 2)) + 
  scale_y_continuous(breaks = seq(0, 60, 5)) + 
  theme_light() + 
  xlab("Age in years since birth") + 
  ylab("Number of patients") +
  ggtitle("Age distribution of iCOPE ALL patients expressed in years since birth") + 
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

# End with dev.off()
dev.off()


## Save data ------------------------------------------------------------------

# Save list of MRD numbers
write.csv(meta_data_join$MRDNR, file = "get_data/iCOPE_ALL_MRD.csv")

# Save wrangled gene count dataframe 
save(gene_count_df, file = "get_data/iCOPE_ALL_gene_counts.rda")

# Save wrangled meta data dataframe
save(meta_data_join, file = "get_data/iCOPE_ALL_meta_data.rda")

# Save tibble containing information about gene expression data
write_csv(info_exp_data, file = "get_data/iCOPE_ALL_exp_data_info.csv")

# Save SummarizedExperiment object of gene expression and meta data 
save(se_counts_meta, file = "get_data/iCOPE_ALL_gene_exp_SE.rda")




