# fusion-organoid-single-cell

## This repository contains code to process single-cell RNAseq data presented in: 
Samarasinghe, R.A., Miranda, O.A., Buth, J.E. et al. Identification of neural oscillations and epileptiform changes in human brain organoids. Nat Neurosci 24, 1488–1500 (2021). https://doi.org/10.1038/s41593-021-00906-5

## Dataset Description:
- Cortex and Medial Ganglionic Eminence (MGE) organoids were generated from human induced pluripotent stem cell lines generated from individuals with Rett syndrome (MECP2-deficient; "Mut") and control ("iCTRL") and collected:
    - day 56 (prior to fusion)
    - day 70 (fusion)
    - day 100 (fusion)
- Annotation: refdata-gex-GRCh38-2020-A (10X Genomics)

## Data Availability:
- The raw dataset is available at GEO under the accession number GSE165577: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165577
- A Seurat object containing the processed single-cell RNA sequencing data is available at zenodo: https://doi.org/10.5281/zenodo.5732813

## Code Overview:
- 01: Bash scripts to process fastq with the Cell Ranger 4.0.0 pipeline (FASTQ -> count matrix)
- 02: R code for filtering, batch correction (LIGER), normalization, and differential expression using Seurat (version 3.2.0)
- 03: R code to calculate bootstrapped confidence intervals for cell-type proportions using single cell differential composition analysis (scDC; see https://sydneybiox.github.io/scDC/index.html for additional information)

## Directory Structure:
- fusion-organoid-single-cell/
    - code/
      - log
      - 01_cellranger_batch_qsub.sh
      - 01_cellranger_count.sh
    - bin/
      - cellranger-4.0.0
    - annotation/
      - refdata-gex-GRCh38-2020-A
    - fastq/
    - cellranger/
    - R/
      - 02_filtering_normalization_differentialexpression_GO
      - 03_scDC


