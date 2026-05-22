#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH -J re_fastqc
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc/re_fastqc_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc/re_fastqc_%j.err

set -euo pipefail

module load FastQC

# -----------------------------
# Input directories
# -----------------------------
TRIM_BASE="/home/yuch3531/genome_analysis_project/3_preprocessing/trimming"
ILLUMINA_DIR="${TRIM_BASE}/illumina"
RNA_BH_DIR="${TRIM_BASE}/rna_bh"
RNA_SERUM_DIR="${TRIM_BASE}/rna_serum"

# -----------------------------
# Output directories
# -----------------------------
OUT_BASE="/home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc"
ILLUMINA_OUT="${OUT_BASE}/illumina"
RNA_BH_OUT="${OUT_BASE}/rna_bh"
RNA_SERUM_OUT="${OUT_BASE}/rna_serum"

mkdir -p "${OUT_BASE}" "${ILLUMINA_OUT}" "${RNA_BH_OUT}" "${RNA_SERUM_OUT}"

echo "Starting FastQC on trimmed data at $(date)"

# -----------------------------
# Illumina trimmed reads
# -----------------------------
if compgen -G "${ILLUMINA_DIR}/*.fastq.gz" > /dev/null; then
    echo "Running FastQC for trimmed Illumina data"
    fastqc -t 2 -o "${ILLUMINA_OUT}" "${ILLUMINA_DIR}"/*.fastq.gz
else
    echo "No FASTQ files found in ${ILLUMINA_DIR}"
fi

# -----------------------------
# RNA_BH trimmed reads
# -----------------------------
if compgen -G "${RNA_BH_DIR}/*.fastq.gz" > /dev/null; then
    echo "Running FastQC for trimmed RNA_BH data"
    fastqc -t 2 -o "${RNA_BH_OUT}" "${RNA_BH_DIR}"/*.fastq.gz
else
    echo "No FASTQ files found in ${RNA_BH_DIR}"
fi

# -----------------------------
# RNA_Serum trimmed reads
# -----------------------------
if compgen -G "${RNA_SERUM_DIR}/*.fastq.gz" > /dev/null; then
    echo "Running FastQC for trimmed RNA_Serum data"
    fastqc -t 2 -o "${RNA_SERUM_OUT}" "${RNA_SERUM_DIR}"/*.fastq.gz
else
    echo "No FASTQ files found in ${RNA_SERUM_DIR}"
fi

echo "FastQC on trimmed data finished at $(date)"
