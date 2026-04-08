#!/bin/bash

name="$1" # sample folder name
BASE_DIR="$2" # base directory
PROJECT_DIR="${BASE_DIR}/fusion-organoid-single-cell"

CELLRANGER_DIR="${PROJECT_DIR}/bin/cellranger-4.0.0"
REF_DATA="${PROJECT_DIR}/annotation/refdata-gex-GRCh38-2020-A"
DATA_DIR="${PROJECT_DIR}/fastq"
FASTQ_PATH="${DATA_DIR}/${name}"

mkdir -p "${PROJECT_DIR}/cellranger"
OUT_DIR="${PROJECT_DIR}/cellranger"
cd "${OUT_DIR}" || exit

echo "Running cellranger count for sample: ${name}"

"${CELLRANGER_DIR}/cellranger" count \
  --id="${name}" \
  --transcriptome="${REF_DATA}" \
  --fastqs="${FASTQ_PATH}" \
  --localcores=8 \
  --localmem=120
