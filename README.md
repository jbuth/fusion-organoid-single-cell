# fusion-organoid-single-cell

This repository contains code to process single-cell RNA-seq data from:  

Samarasinghe, R.A., Miranda, O.A., Buth, J.E. et al. Identification of neural oscillations and epileptiform changes in human brain organoids. Nat Neurosci 24, 1488–1500 (2021). https://doi.org/10.1038/s41593-021-00906-5

## Dataset Description:
- Cortex and Medial Ganglionic Eminence (MGE) organoids were generated from human induced pluripotent stem cell lines generated from individuals with Rett syndrome (MECP2-deficient; "Mut") and control ("iCTRL") and collected:
    - day 56 (prior to fusion)
    - day 70 (fusion)
    - day 100 (fusion)
- Annotation: refdata-gex-GRCh38-2020-A (10X Genomics)

## Data Availability:
- The raw dataset is available at GEO under the accession number GSE165577: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165577
- A Seurat object containing the processed single-cell RNA-seq data is available at Zenodo: https://doi.org/10.5281/zenodo.5732813

## Pipeline Overview:
FASTQ → Cell Ranger (count matrix generation) → Seurat (filtering, batch correction, normalization, clustering, differential expression, pathway analysis) → scDC (cell-type proportion analysis)

## Code Overview:
- 01: Bash scripts to run Cell Ranger 4.0.0 (FASTQ → gene expression count matrix)
- 02: R code using Seurat (v3.2.0) for filtering, batch correction (LIGER), normalization, clustering, differential expression, pathway analysis, and scDC (v0.1.0) for calculation of bootstrapped confidence intervals for cell-type proportions (see https://sydneybiox.github.io/scDC/index.html)

## Directory Structure:
The bash scripts (01) were run on an HPC cluster (UCLA Hoffman2). The R script (02) was run locally, with the Cell Ranger output copied to the local project folder.
```
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
    - cellranger/       # Cell Ranger output directories (one per sample)
    - R/
      - 02_Seurat.R
```

