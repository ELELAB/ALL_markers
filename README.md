Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen

# ALL_markers

## Introduction

This repository contains scripts related to the discovery of gene expression markers 
that distinguish two acute lymphoblastic leukemia subtypes: B-cell and T-cell ALL. 

## Installation requirements and guideline for reproducing the analyses

All analyses have been performed on a GNU/Linux server. 

### Computing environment

To reproduce the data and results, you will need to set up a conda environment which
contains the expected R version and the required R packages. This requires being able
to run Anaconda by means of the `conda` executable. 

If you don't have access to `conda`, please see the [Miniconda installer page](https://docs.conda.io/en/latest/miniconda.html)
for instructions on how to install Miniconda. 

Once you have access to `conda`, follow the below guidelines to reproduce the results
and data:

1. Clone the GitHub repository into a local directory on your local machine:

```
git clone https://github.com/ELELAB/ALL_markers.git
cd ALL_markers
```

2. Create a virtual environment using `conda` and activate it afterwards. 
The environment should be placed in the ALL_markers folder:

```
conda create --prefix ./env_ALL -c conda-forge r-base=4.2 r-pacman=0.5.1 r-curl=4.3.3 r-ragg=1.2.5 r-renv=0.16.0 r-osfr=0.2.9 r-cairo=1.6.0
conda activate ./env_ALL
```

3. Run the analyses:

```
bash ./run_all.sh
```

**WARNING**: our scripts use the [renv](https://rstudio.github.io/renv/articles/renv.html) 
R package to handle automatic dependency installation. `Renv` writes packages in 
its own cache folder, which is by default in the user's home directory. This might not be 
desirable if free space in the home directory is limited. You can change the location of 
the `Renv` root folder by setting a system environment variable - please see comments 
in the `run_all.sh` script.

The `run_all.sh` will perform the following steps to reproduce all results and data:

1. Download data from the corresponding [OSF repository](https://osf.io/kgfpv/) which
contains the required data to run the analyses and all results associated with the analyses. 

2. Install in the environment all necessary packages to run the analyses.

3. Perform all analyses. 

## Structure and content of GitHub and OSF repositories

The GitHub and OSF repositories contain scripts and data/results associated
with this publication, respectively. Both repositories are structured in the
same way with a main folder named after the main analyses which then contains 
all scripts and data/results associated with the main analysis. Below is an
overview of these main folders. See README files in each main folder for more
details.

TARGET_data:
- This directory contains gene expression data of the TARGET-ALL-P2 project and
  information about its metadata such as age of patients, number of samples,
  number of genes, subtype information etc. Moreover, here the data is
  subsetted to contain only Primary Blood Derived Cancer - Bone Marrow samples 

TARGET_replicates:
- This directory investigates replicate samples found in the TARGET-ALL-P2
  data and adjusts the data for these replicates
  
TARGET_processing:
- This directory processes gene expression data of the TARGET-ALL-P2
  project (preprocessing, normalization, and filtering using TCGAbiolinks) 
  
TARGET_transform:
- This directory voom transforms TARGET-ALL-P2 gene expression data
- In here, the directory voom_transform_bmp_raw voom transforms the raw
  gene expression data and the directory voom_transform_bmp_comp voom
  transforms processed data

TARGET_pca:
- This directory performs dimensionality reduction of TARGET-ALL-P2 gene expression
  data
- In here, the directory bone_marrow_primary_pca_raw performs MDS of the raw gene
  expression data and the directory bone_marrow_primary_pca_comp performs MDS of 
  the processed data (both before and after batch correction). Moreover, in 
  bone_marrow_primary_pca_comp, PCA is performed to find the contributions 
  of genes to principal components 1 and 2 

TARGET_batch:
- This directory performs batch correction of TARGET-ALL-P2 gene expression data

TARGET_dea:
- This directory performs differential expression analysis of TARGET-ALL-P2
  gene expression data

TARGET_housekeeping
- This directory investigates overlaps of discovered differentially expressed 
  genes with lists of housekeeping genes

TARGET_enrichment
- This directory performs enrichment analyses of consensus differentially
  expressed genes
