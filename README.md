# fusion-organoid-single-cell

## This repository contains code to process single-cell RNAseq data presented in: 
Samarasinghe, R.A., Miranda, O.A., Buth, J.E. et al. Identification of neural oscillations and epileptiform changes in human brain organoids. Nat Neurosci 24, 1488–1500 (2021). https://doi.org/10.1038/s41593-021-00906-5

## Data Availability:
- The raw dataset is available at GEO under the accession number GSE165577: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165577
- A Seurat object containing the processed single-cell RNA sequencing data is available at zenodo: https://doi.org/10.5281/zenodo.5732813

## Contents:
- 01: Bash scripts to run cellranger count on the fastq files
- 02: R code for filtering, normalization, and differential expression using Seurat
- 03: R code to run single cell differential composition analysis (scDC; see https://sydneybiox.github.io/scDC/index.html for additional information)
