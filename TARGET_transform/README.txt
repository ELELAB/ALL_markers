* INTRODUCTION *

This directory (TARGET_transform) contains transformed TARGET-ALL-P2
data. This directory contains two subdirectories used for 
transforming data of different purposes and workflows:

1. voom_transform_bmp_comp
-  Here, TARGET-ALL-P2 gene expression data is voom transformed
   where data has been subject to workflow:
	* subset to bone marrow primary samples
	* adjust for replicates
	* preprocess
	* normalize using updated geneInfoHT table
	* filter

2. voom_transform_bmp_raw
-  Here, raw TARGET-ALL-P2 gene expression data is voom transformed
   where data has been subsetted to bone marrow primary samples

