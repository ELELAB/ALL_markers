### This script contains functions used for DEA of TARGET-ALL-P2 gene expression data using the limma-voom pipeline
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered



### ---------------------------------------------------------------------------------------------------------

### This function, create_design_matrix, creates a design matrix for DEA based on 
### different designs of batch effects. This function takes four 
### inputs:
### batch_factor:             A string character indicating which batch factor should be included
###                           in the design matrix. Options are "no_batch", "portion_batch", 
###                           "years_batch", or "portion_years_batch".
### condition_types:          A vector indicating the condition type of each sample. 
### tissue_portion_aliquot:   A vector indicating the tissue portion of each sample.
### years_of_diagnosis:       A vector indicating the year of diagnosis of each sample. 
### The function outputs a design matrix depending on the chosen batch factor. 
create_design_matrix <- function(batch_factor,
                                 condition_types,
                                 tissue_portion_aliquot,
                                 years_of_diagnosis) {
  
  # If no batch correction is desired
  if (batch_factor == "no_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types)
    colnames(design_matrix)[1:2] <- c("B_cell", "T_cell")
    
  }
  
  # If tissue portion as batch factor is desired
  else if (batch_factor == "portion_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types+tissue_portion_aliquot)
    colnames(design_matrix)[1:2] <- c("B_cell", "T_cell")
    
    # If years of diagnosis as batch factor is desired
  } else if (batch_factor == "years_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types+years_of_diagnosis)
    colnames(design_matrix)[1:2] <- c("B_cell", "T_cell")
    
    # If both tissue portion and years of diagnosis as batch factors are desired 
  } else if (batch_factor == "portion_years_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types+tissue_portion_aliquot+years_of_diagnosis)
    colnames(design_matrix)[1:2] <- c("B_cell", "T_cell")
    
    # Print message if user enters invalid batch factor option 
  } else {
    
    print("Batch factor not valid. Select one of the following batch correction options: no_batch, portion_batch, years_batch or portion_years_batch")
    
  }
  
  return(design_matrix)
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, voom_transform, voom transforms processed gene expression data. 
### This function takes four inputs:
### file_name:          File name of plot to be saved displaying the mean-variance
###                     trend (voom transformation)
### dir_output:         Directory where the plot should be saved
### processed_data:     A gene expression data matrix with genes in rows and samples
###                     in columns. 
### design_matrix:      A design matrix 
### The function outputs a voom object. 
voom_transform <- function(file_name,
                           dir_output,
                           processed_data,
                           design_matrix) {
  
  # Save a plot of the mean-variance trend
  png(filename = paste0(dir_output, 
                        "/", 
                        file_name))
  
  # Do voom transformation
  voom_data <- voom(counts = processed_data,
                    design = design_matrix,
                    plot = TRUE, 
                    save.plot = TRUE)
  
  # End with dev.off() to save plot
  dev.off()
  
  return(voom_data)
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, limma_DEA, performs DEA using the limma-voom pipeline. This 
### function takes two inputs:
### design_matrix:        A design matrix to be used in the DEA
### voom_object:          A voom object containing the voom transformed gene expression data   
### The function outputs a table containing significantly DEGs. 
limma_DEA <- function(design_matrix,
                      voom_object) {
  
  # Make group contrasts
  subtype_contrast <- makeContrasts("B_cell-T_cell", 
                                    levels = design_matrix)
  
  ### Do DEA with limma 
  
  # Fit linear model to processed data incorporating design matrix
  lm_fit <- lmFit(object = voom_object, 
                  design = design_matrix)
  
  # Compute estimated coefficients for contrasts
  contr_fit <- contrasts.fit(fit = lm_fit, 
                             contrasts = subtype_contrast)
  
  # Fit empirical Bayes model for differential expression using treat
  treat_fit <- treat(fit = contr_fit, 
                     lfc = 1)
  
  # Extract table of genes from fit
  DEA_table <- topTreat(fit = treat_fit, 
                        coef = 1, 
                        n = Inf)
  
  return(DEA_table)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, DEG_significant, finds significantly DEGs according to a cutoff
### of FDR values. This function takes two inputs:
### DEA_table:            A matrix containing results from DEA among which 
###                       significant DEGs should be found
### fdr.cut:              A cutoff for adjusted p values to find significantly DEGs 
### The function outputs the DEA table containing only significantly DEGs. 
DEG_significant <- function(DEA_table,
                            fdr.cut) {
  
  # Filter DEA_table according to FDR cutoff to extract significant results 
  DEA_table_extract <- DEA_table[which(DEA_table$adj.P.Val < fdr.cut),]
  
  return(DEA_table_extract)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, ensembl_hugo_conversion, converts ENSEMBL IDs to gene names.
### The function takes one input:
### ensembl_id:     A list of ENSEMBL IDs to be converted into gene names.
### The function returns the ENSEMBL IDs and their corresponding gene name.
### The function is based on code from the CAMPP pipeline.
ensembl_hugo_conversion <- function(ensembl_id) {
  
  # Connect to biomart database and dataset hosted by ENSEMBL (release 106)
  ensembl_106 <- useEnsembl(biomart = 'genes', 
                            dataset = 'hsapiens_gene_ensembl',
                            version = 106)
  
  # Convert ENSEMBL IDs to HUGO gene names
  gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      mart = ensembl_106,
                      values = ensembl_id,
                      filters = 'ensembl_gene_id',
                      uniqueRows = TRUE,
                      bmHeader = T)
  
  return(gene_names)
  
}

### ---------------------------------------------------------------------------------------------------------



