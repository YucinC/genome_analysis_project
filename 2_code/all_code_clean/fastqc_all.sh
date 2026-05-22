#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -p haswell
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J fastqc_all
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/fastqc/fastqc_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/fastqc/fastqc_%j.err

set -euo pipefail

module load FastQC

# -----------------------------
# Input directories
# -----------------------------
ILLUMINA_DIR="/home/yuch3531/genome_analysis_project/1_data/genomics_data/Illumina"
RNABH_RAW_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_BH/raw"
RNABH_TRIM_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_BH/trimmed"
RNASERUM_RAW_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_Serum/raw"
RNASERUM_TRIM_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_Serum/trimmed"
TNSEQ_BHI_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/Tn-Seq_BHI"
TNSEQ_HSERUM_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/Tn-Seq_HSerum"
TNSEQ_SERUM_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/Tn-Seq_Serum"
PACBIO_DIR="/home/yuch3531/genome_analysis_project/1_data/genomics_data/PacBio"

# -----------------------------
# Output base directory
# -----------------------------
OUT_BASE="/home/yuch3531/genome_analysis_project/3_preprocessing/fastqc"

# Create output subdirectories
mkdir -p "${OUT_BASE}/illumina"
mkdir -p "${OUT_BASE}/rna_bh_raw"
mkdir -p "${OUT_BASE}/rna_bh_trimmed"
mkdir -p "${OUT_BASE}/rna_serum_raw"
mkdir -p "${OUT_BASE}/rna_serum_trimmed"
mkdir -p "${OUT_BASE}/tnseq_bhi"
mkdir -p "${OUT_BASE}/tnseq_hserum"
mkdir -p "${OUT_BASE}/tnseq_serum"
mkdir -p "${OUT_BASE}/pacbio"

echo "Starting FastQC at $(date)"

# -----------------------------
# Illumina genomic reads
# -----------------------------
echo "Running FastQC: Illumina genomic reads"
fastqc -o "${OUT_BASE}/illumina" "${ILLUMINA_DIR}"/*.fq.gz

# -----------------------------
# RNA-seq raw
# -----------------------------
echo "Running FastQC: RNA-Seq BH raw"
fastqc -o "${OUT_BASE}/rna_bh_raw" "${RNABH_RAW_DIR}"/*.fastq.gz

echo "Running FastQC: RNA-Seq Serum raw"
fastqc -o "${OUT_BASE}/rna_serum_raw" "${RNASERUM_RAW_DIR}"/*.fastq.gz

# -----------------------------
# RNA-seq trimmed
# -----------------------------
echo "Running FastQC: RNA-Seq BH trimmed"
fastqc -o "${OUT_BASE}/rna_bh_trimmed" "${RNABH_TRIM_DIR}"/*.fastq.gz

echo "Running FastQC: RNA-Seq Serum trimmed"
fastqc -o "${OUT_BASE}/rna_serum_trimmed" "${RNASERUM_TRIM_DIR}"/*.fastq.gz

# -----------------------------
# Tn-seq
# -----------------------------
echo "Running FastQC: Tn-Seq BHI"
fastqc -o "${OUT_BASE}/tnseq_bhi" "${TNSEQ_BHI_DIR}"/*.fastq.gz

echo "Running FastQC: Tn-Seq HSerum"
fastqc -o "${OUT_BASE}/tnseq_hserum" "${TNSEQ_HSERUM_DIR}"/*.fastq.gz

echo "Running FastQC: Tn-Seq Serum"
fastqc -o "${OUT_BASE}/tnseq_serum" "${TNSEQ_SERUM_DIR}"/*.fastq.gz

# -----------------------------
# PacBio subreads
# -----------------------------
echo "Running FastQC: PacBio subreads"
fastqc -o "${OUT_BASE}/pacbio" "${PACBIO_DIR}"/*.fastq.gz

echo "FastQC finished at $(date)"
