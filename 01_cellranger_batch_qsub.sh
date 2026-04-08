#!/bin/bash

## Save this script as:
## /path/to/fusion-organoid-single-cell/code/01_cellranger_batch_qsub.sh

## To execute, cd into the code folder, then type:
## ./01_cellranger_batch_qsub.sh

## Description:
## Loops through sample folders in /fastq and submits one qsub job per sample
## qsub runs cellranger count via the script 01_cellranger_count.sh

## ------ Setup directories ------- ##

## Base directory: /path/to/
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

## Subdirectories
CODE_DIR="${BASE_DIR}/fusion-organoid-single-cell/code"
DATA_DIR="${BASE_DIR}/fusion-organoid-single-cell/fastq"

## --- Data directory structure --- ##

## Example:
## ${DATA_DIR}/D100_Ctrl_docked_1
##      /__D100_Ctrl_docked_1.mro
##      /D100_Ctrl_docked_1_S5_L001_I1_001.fastq.gz
##      /D100_Ctrl_docked_1_S5_L001_R1_001.fastq.gz
##      /D100_Ctrl_docked_1_S5_L001_R2_001.fastq.gz
##      /D100_Ctrl_docked_1_S5_L002_I1_001.fastq.gz
##      /D100_Ctrl_docked_1_S5_L002_R1_001.fastq.gz
##      /D100_Ctrl_docked_1_S5_L002_R2_001.fastq.gz

## --- Change to data directory --- ##

cd "${DATA_DIR}" || exit

## ---- Batch submit qsub jobs ---- ##

for file in */; do
  
  name="${file%/}"  ## save the sample folder name
  echo "${name}"
  
  ## directory for log stdout (-o) and stderr (-e)
  ## runtime, memory per core, highp runs on lab purchased compute nodes
  ## shared memory request for 8 cores with 16G per core (~128G total)
  qsub \
    -o "${CODE_DIR}/log" \
    -e "${CODE_DIR}/log" \
    -l h_rt=8:00:00,h_data=16G,highp \
    -pe shared 8 \
    "${CODE_DIR}/01_cellranger_count.sh" "${name}" "${BASE_DIR}"
  
done
