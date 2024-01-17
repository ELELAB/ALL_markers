### Functions for housekeeping analysis

### ---------------------------------------------------------------------------------------------------------

### This function, wrangle_exp_data, wrangles gene expression data to prepare it
### for visualization in a boxplot stratified on subtype. 
### This function takes five inputs:
### exp_data:           gene expression data with genes in rows and samples in columns
### HK_data:            table of housekeeping genes including ENSEMBL gene ID and external gene name
### gene_row_num:       row number of gene whose expression values are to be wrangled
### B_ALL_barcodes:     vector of barcodes defining B-ALL subtype samples
### T_ALL_barcodes:     vector of barcodes defining T-ALL subtype samples
### This function returns a tibble containing barcodes of all samples, expression values,
### subtype, ENSEMBL gene ID and external gene name of respective gene. 
wrangle_exp_data <- function(exp_data, 
                             HK_data,
                             gene_row_num,
                             B_ALL_barcodes,
                             T_ALL_barcodes) {
  
  # Get gene expression values of gene of T subtype samples
  gene_exp_T <- exp_data[gene_row_num, 
                         T_ALL_barcodes] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "barcodes") %>% 
    as_tibble() %>% 
    dplyr::mutate(exp_value = as.numeric(.),
           subtype = "T-cell ALL",
           ENSEMBL_ID = rownames(exp_data)[gene_row_num]) %>% 
    dplyr::select(-.)
  
  # Get gene expression values of gene of B subtype samples
  # Combine into one tibble with gene_exp_T 
  gene_exp_T_B <- exp_data[gene_row_num,
                           B_ALL_barcodes] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "barcodes") %>% 
    as_tibble() %>% 
    dplyr::mutate(exp_value = as.numeric(.),
           subtype = "B-cell ALL",
           ENSEMBL_ID = rownames(exp_data)[gene_row_num]) %>% 
    dplyr::select(-.) %>% 
    bind_rows(., 
              gene_exp_T) %>%
    dplyr::mutate(HK_data %>% 
             filter(ENSEMBL_ID == rownames(exp_data)[gene_row_num]) %>% 
             select(Gene_name)) 
  
  # Return tibble of gene expression values for respective gene   
  return(gene_exp_T_B)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, boxplot_exp, visualizes gene expression data in a boxplot
### stratified on subtype. 
### This function takes seven inputs:
### wrangled_exp_data:      wrangled gene expression returned by function wrangle_exp_data
### subtype_col:            column indicating subtype annotations of samples 
### exp_col:                column containing expression values to be plotted
### ensembl_col:            column containing ENSEMBL gene ID of respective gene
### gene_name_col:          column containing external gene name of respective gene
### ylab:                   string defining label for y axis
### plot_title:             string defining title of plot
### This function returns a boxplot of expression values of the respective gene stratified
### on subtype. 
boxplot_exp <- function(wrangled_exp_data,
                        subtype_col,
                        exp_col,
                        ensembl_col,
                        gene_name_col,
                        ylab,
                        plot_title) {
  
  # Create boxplot of gene stratified on subtype
  boxplot_gene <- ggplot(data = wrangled_exp_data, 
                         mapping = aes(x = subtype_col, 
                                       y = exp_col,
                                       fill = subtype_col)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#21908CFF", 
                                 "#FDE725FF")) +
    labs(x = "",
         y = ylab,
         title = plot_title,
         fill = "Subtype") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 16,
                                     color = "black"),
          axis.text.y = element_text(size = 16,
                                     color = "black"),
          axis.title.y = element_text(size = 16,
                                      color = "black"),
          panel.border = element_rect(size = 1.5))
  
  return(boxplot_gene)
  
}

### ---------------------------------------------------------------------------------------------------------





