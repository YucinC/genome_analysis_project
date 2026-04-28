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

mkdir -p "$FASTQC_OUT"
mkdir -p "$MULTIQC_OUT"

echo "======================================"
echo "Start rerun RNA serum FastQC + MultiQC"
date
echo "Input: $RNA_SERUM_DIR"
echo "FastQC output: $FASTQC_OUT"
echo "MultiQC output: $MULTIQC_OUT"
echo "======================================"

# -----------------------------
# Clean old RNA serum FastQC results
# -----------------------------
echo "Removing old RNA serum FastQC results..."
rm -f "$FASTQC_OUT"/*

# -----------------------------
# Run FastQC on all RNA serum trimmed FASTQ files
# -----------------------------
if compgen -G "${RNA_SERUM_DIR}/*.fastq.gz" > /dev/null; then
    echo "Running FastQC for RNA serum trimmed data..."
    fastqc -t 2 -o "$FASTQC_OUT" "${RNA_SERUM_DIR}"/*.fastq.gz
else
    echo "ERROR: No FASTQ files found in ${RNA_SERUM_DIR}"
    exit 1
fi

echo "FastQC finished at $(date)"

# -----------------------------
# Clean old RNA serum MultiQC results
# -----------------------------
echo "Removing old RNA serum MultiQC results..."
rm -rf "${MULTIQC_OUT:?}"/*

# -----------------------------
# Run MultiQC only for RNA serum
# -----------------------------
echo "Running MultiQC for RNA serum..."
multiqc "$FASTQC_OUT" \
    --outdir "$MULTIQC_OUT" \
    --filename rna_serum_multiqc_report.html \
    --force

echo "MultiQC finished at $(date)"
echo "Done."
