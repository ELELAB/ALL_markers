* INTRODUCTION *

This directory (TARGET_pca) performs MDS analysis of TARGET-ALL-P2 
gene expression data. This directory contains the following subdirectories
where MDS analysis is performed for different purposes and in
different workflows:

1. bone_marrow_primary_pca_comp
-  Here, MDS is performed on TARGET-ALL-P2 gene expression data
   where data has been subject to two workflows. These two workflows
   share the same six first steps which are the following and
   which define the first workflow:
	a. subset to bone marrow primary samples
	b. adjust for replicates
	c. preprocess
	d. normalize using updated geneInfoHT table
	e. filter
	f. voom transform
   In the second workflow, a subsequent batch correction is done:
	g. batch corrected for year of diagnosis

2. bone_marrow_primary_pca_raw
-  Here, MDS is performed on raw TARGET-ALL-P2 gene expression 
   data where data has been subject to the following workflow:
	a. subset to bone marrow primary samples
