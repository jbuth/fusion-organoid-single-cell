#!/bin/bash

name=$1 # sample folder name from previous script
BASE_DIR=$2 # saved base directory from previous script
DATA_DIR=${BASE_DIR}/fusion-organoid-single-cell/fastq

ref_data=${BASE_DIR}/fusion-organoid-single-cell/annotation/refdata-gex-GRCh38-2020-A
cellranger_dir=${BASE_DIR}/fusion-organoid-single-cell/bin/cellranger-4.0.0
fastq_path=${DATA_DIR}/${name}

cd "${DATA_DIR}"/"${name} || exit" || exit

"${cellranger_dir}"/cellranger count --id="${name}" --transcriptome="${ref_data}" --fastqs="${fastq_path}"
