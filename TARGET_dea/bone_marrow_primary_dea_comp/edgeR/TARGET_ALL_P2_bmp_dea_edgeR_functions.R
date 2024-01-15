### This script contains functions used for DEA of TARGET-ALL-P2 gene expression data using the edgeR pipeline
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates




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
    
  }
  
  # If tissue portion as batch factor is desired
  else if (batch_factor == "portion_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types+tissue_portion_aliquot)
    
    # If years of diagnosis as batch factor is desired
  } else if (batch_factor == "years_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types+years_of_diagnosis)
    
    # If both tissue portion and years of diagnosis as batch factors are desired 
  } else if (batch_factor == "portion_years_batch") {
    
    # Create design matrix
    design_matrix <- model.matrix(~0+condition_types+tissue_portion_aliquot+years_of_diagnosis)
    
    # Print message if user enters invalid batch factor option 
  } else {
    
    print("Batch factor not valid. Select one of the following batch correction options: no_batch, portion_batch, years_batch or portion_years_batch")
    
  }
  
  colnames(design_matrix)[1:2] <- c("B_cell", "T_cell")
  
  return(design_matrix)
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, edgeR_DEA, performs DEA using the edgeR pipeline. This 
### function takes two inputs:
### design_matrix:        A design matrix to be used in the DEA
### dge_list:             A DGEList data object containing the gene expression data counts   
### The function outputs a table with DEA results for all genes. 
edgeR_DEA <- function(design_matrix,
                      dge_list){
  
  # Create contrast
  subtype_contrast <- makeContrasts("B_cell-T_cell", levels = design_matrix)
  
  # Estimate dispersion factors 
  y <- estimateDisp(y, design_matrix)
  
  # Fit model to count data
  fit <- glmQLFit(y, design_matrix)
  
  # Conduct threshold-specific statistical testing
  tr <- glmTreat(fit, contrast=subtype_contrast, lfc=1)
  
  # Get results
  DEA_table <- topTags(tr, n = Inf, sort.by = "p.value", adjust.method = "fdr")
  
  
  return(DEA_table)
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




