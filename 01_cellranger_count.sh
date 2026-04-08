#!/bin/bash

name="$1" # sample folder name
BASE_DIR="$2" # base directory
DATA_DIR="${BASE_DIR}/fusion-organoid-single-cell/fastq"

REF_DATA="${BASE_DIR}/fusion-organoid-single-cell/annotation/refdata-gex-GRCh38-2020-A"
CELLRANGER_DIR="${BASE_DIR}/fusion-organoid-single-cell/bin/cellranger-4.0.0"
FASTQ_PATH="${DATA_DIR}/${name}"

"${CELLRANGER_DIR}/cellranger" count \
  --id="${name}" \
  --transcriptome="${REF_DATA}" \
  --fastqs="${FASTQ_PATH}" \
  --localcores=8 \
  --localmem=120
