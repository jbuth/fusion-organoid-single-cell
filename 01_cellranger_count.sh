#!/bin/bash

DATA_DIR=/path/to/FASTQ
name=$1 # sample name
cellranger_ref_data=/bin/refdata-gex-GRCh38-2020-A
fastq_path=${DATA_DIR}/${name}
cellranger_dir=/path/to/cellranger-4.0.0

cd ${DATADIR}/${name}

${cellranger_dir}/cellranger count --id=${name} --transcriptome=${cellranger_ref_data} --fastqs=${fastq_path}
