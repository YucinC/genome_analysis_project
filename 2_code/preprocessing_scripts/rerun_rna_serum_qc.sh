#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH -J rerun_rna_serum_qc
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/rerun_rna_serum_qc_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/rerun_rna_serum_qc_%j.err

set -euo pipefail

module load FastQC
module load MultiQC

# -----------------------------
# Input directory
# -----------------------------
RNA_SERUM_DIR="/home/yuch3531/genome_analysis_project/3_preprocessing/trimming/rna_serum"

# -----------------------------
# FastQC output directory
# -----------------------------
FASTQC_OUT="/home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc/rna_serum"

# -----------------------------
# MultiQC output directory
# -----------------------------
MULTIQC_OUT="/home/yuch3531/genome_analysis_project/3_preprocessing/re_multiqc/rna_serum"

# Create output directories if they do not already exist
mkdir -p "$FASTQC_OUT"
mkdir -p "$MULTIQC_OUT"

echo "======================================"
echo "Start RNA serum FastQC + MultiQC rerun"
date
echo "Input: $RNA_SERUM_DIR"
echo "FastQC output: $FASTQC_OUT"
echo "MultiQC output: $MULTIQC_OUT"
echo "======================================"

# -----------------------------
# Check input FASTQ files
# -----------------------------
FASTQ_FILES=("${RNA_SERUM_DIR}"/*.fastq.gz)

if [ ! -e "${FASTQ_FILES[0]}" ]; then
    echo "ERROR: No FASTQ files found in ${RNA_SERUM_DIR}"
    exit 1
fi

echo "Detected ${#FASTQ_FILES[@]} FASTQ files:"
printf '%s\n' "${FASTQ_FILES[@]}"

# -----------------------------
# Run FastQC on all RNA serum trimmed FASTQ files
# -----------------------------
echo "Running FastQC for all RNA serum trimmed data..."
fastqc -t "${SLURM_CPUS_PER_TASK:-2}" -o "$FASTQC_OUT" "${FASTQ_FILES[@]}"

echo "FastQC finished at $(date)"

# -----------------------------
# Run MultiQC for RNA serum FastQC results
# -----------------------------
echo "Running MultiQC for RNA serum..."
multiqc "$FASTQC_OUT" \
    --outdir "$MULTIQC_OUT" \
    --filename rna_serum_multiqc_report.html \
    --force

echo "MultiQC finished at $(date)"
echo "Done."
