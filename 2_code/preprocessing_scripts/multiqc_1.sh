#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH -J multiqc_1
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/multiqc/multiqc_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/multiqc/multiqc_%j.err

set -euo pipefail

module load MultiQC

BASE_DIR=/home/yuch3531/genome_analysis_project/3_preprocessing
FASTQC_DIR=${BASE_DIR}/fastqc
MULTIQC_DIR=${BASE_DIR}/multiqc

mkdir -p ${MULTIQC_DIR}

for subdir in illumina pacbio rna_bh_raw rna_bh_trimmed rna_serum_raw rna_serum_trimmed tnseq_bhi tnseq_hserum tnseq_serum
do
    echo "======================================"
    echo "Running MultiQC for: ${subdir}"

    mkdir -p ${MULTIQC_DIR}/${subdir}

    multiqc ${FASTQC_DIR}/${subdir} \
        --force \
        --outdir ${MULTIQC_DIR}/${subdir} \
        --dirs

    echo "Finished: ${subdir}"
done

echo "======================================"
echo "All MultiQC reports finished."
