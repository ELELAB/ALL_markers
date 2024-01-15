#### This script performs DEA of TARGET-ALL-P2 gene expression data using the DESeq2 pipeline
### The expression data is raw, but filtered to contain genes that are left from the following workflow:
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered

### ---------------------------------------------------------------------------------------------------------

### This function, create_DESeqDataSet_design, creates a design formula for DEA based on 
### different designs of batch effects. This function takes one 
### input:
### batch_factor:             A string character indicating which batch factor should be included
###                           in the design matrix. Options are "no_batch", "portion_batch", 
###                           "years_batch", or "portion_years_batch".
### The function outputs a design formula depending on the chosen batch factor. 
create_DESeqDataSet_design <- function(batch_factor){
  
  # No batch correction is desired
  if (batch_factor == "no_batch"){
    dds_design <- "~ primary_diagnosis"
  }
  
  # If tissue portion as batch factor is desired
  else if (batch_factor == "portion_batch"){
    dds_design <- "~ tissue_portion_aliquot + primary_diagnosis"
  }
  
  # If year of diagnosis as batch factor is desired
  else if (batch_factor == "years_batch"){
    dds_design <- "~ diagnosis_year + primary_diagnosis"
  }
  
  # If both tissue portion and years of diagnosis as batch factors are desired 
  else if (batch_factor == "portion_years_batch"){
    dds_design <- "~ diagnosis_year + tissue_portion_aliquot + primary_diagnosis"
  }
  
  return(dds_design)
}


### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, deseq2_DEA, performs DEA using the deseq2 pipeline. This
### function takes four inputs:
### expression_data:      A RangedSummarizedExperiment with counts data to be used in the DEA
### design:               A string character with the design formula for the design of the DEA
### lfc_threshold:        A floating point number of which we want to test for log2 fold changes greater
###                       in absolute value than this value
### alpha:                A floating point number indicating the desired significance level for independent
###                       filtering
### The function outputs a table containing results for all tested genes.
deseq2_DEA <- function(expression_data, 
                       design, 
                       lfc_threshold, 
                       alpha){
  
  # Create deseq data set
  dds <- DESeqDataSet(expression_data, design = as.formula(design))
  
  # Rewrite level names
  levels(dds$primary_diagnosis) <- sub("-", "_", levels(dds$primary_diagnosis))
  levels(dds$primary_diagnosis) <- sub("/.*", "", levels(dds$primary_diagnosis))
  levels(dds$primary_diagnosis) <- gsub(" ", "_", levels(dds$primary_diagnosis))
  
  # Factor the reference condition and set the order                                         
  dds$primary_diagnosis <- factor(dds$primary_diagnosis, levels = c("Precursor_B_cell_lymphoblastic_leukemia", "T_lymphoblastic_leukemia"))
  
  # Run differential expression analysis
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds, maxit = 500)
  
  # Get results and set contrasts (logFC: B-cell/T-cell) 
  res = results(dds, 
                contrast = c("primary_diagnosis", "Precursor_B_cell_lymphoblastic_leukemia", "T_lymphoblastic_leukemia"), 
                alpha = alpha, 
                lfcThreshold = lfc_threshold)
  
  # Order by ascending adjusted p-value
  resOrdered <- res[order(res$padj),]
  
  return(resOrdered)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, filter_dea_results, finds significantly DEGs according to a cutoff
### of FDR values. This function takes two inputs:
### dea_results:          A matrix containing results from DEA among which 
###                       significant DEGs should be found
### fdr.cut:              A cutoff for adjusted p values to find significantly DEGs 
### The function outputs the DEA table containing only significantly DEGs.
filter_dea_results <- function(dea_results, 
                               fdr.cut){
  
  # Find and collect index for significant log2FoldChanges
  index_significant <- which(dea_results$padj < fdr.cut)
  
  # Filter the DEA results
  filtered_dea_results <- dea_results[index_significant,]
  
  # Sort filtered results by adjusted p-value and then lfc
  filtered_dea_results <- filtered_dea_results[order(filtered_dea_results$padj),]
  
  return(filtered_dea_results)
  
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
