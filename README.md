# fusion-organoid-single-cell

## This repository contains code to process single-cell RNAseq data presented in: 
Samarasinghe, R.A., Miranda, O.A., Buth, J.E. et al. Identification of neural oscillations and epileptiform changes in human brain organoids. Nat Neurosci 24, 1488–1500 (2021). https://doi.org/10.1038/s41593-021-00906-5

## Dataset Description:
- Cortex and Medial Ganglionic Eminence (MGE) organoids were generated from human induced pluripotent stem cell lines generated from individuals with Rett syndrome (MECP2-deficient) and control (WT) human pluripotent stem cell lines (Mut and iCTRL) and collected:
    - day 56 (prior to fusion)
    - day 70 (fusion)
    - day 100 (fusion)
- Annotation: refdata-gex-GRCh38-2020-A (10X Genomics)

## Data Availability:
- The raw dataset is available at GEO under the accession number GSE165577: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165577
- A Seurat object containing the processed single-cell RNA sequencing data is available at zenodo: https://doi.org/10.5281/zenodo.5732813

## Contents:
- 01: Bash scripts to run the Cell Ranger 4.0.0 pipeline (fastq -> counts)
- 02: R code for filtering, normalization, and differential expression using Seurat (version 3.2.0)
- 03: R code to run single cell differential composition analysis to compare Mut and iCTRL fusion organoids  (scDC; see https://sydneybiox.github.io/scDC/index.html for additional information)
