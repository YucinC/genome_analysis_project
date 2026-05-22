#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -J trim_illumina
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/trimming/trim_illumina_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/trimming/trim_illumina_%j.err

set -euo pipefail

# Load Trimmomatic.
# If this exact module name fails on your account, run: module spider Trimmomatic
module load Trimmomatic

# Input files
IN1="/home/yuch3531/genome_analysis_project/1_data/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz"
IN2="/home/yuch3531/genome_analysis_project/1_data/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz"

# Output directory
OUTDIR="/home/yuch3531/genome_analysis_project/3_preprocessing/trimming/illumina"
mkdir -p "${OUTDIR}"

# Output files
P1="${OUTDIR}/trim_paired_E745_1.fastq.gz"
U1="${OUTDIR}/trim_unpaired_E745_1.fastq.gz"
P2="${OUTDIR}/trim_paired_E745_2.fastq.gz"
U2="${OUTDIR}/trim_unpaired_E745_2.fastq.gz"

echo "Starting Trimmomatic at $(date)"
echo "Input R1: ${IN1}"
echo "Input R2: ${IN2}"
echo "Output dir: ${OUTDIR}"

trimmomatic PE \
  -threads 2 \
  -phred33 \
  "${IN1}" "${IN2}" \
  "${P1}" "${U1}" \
  "${P2}" "${U2}" \
  TRAILING:20 \
  SLIDINGWINDOW:4:20 \
  MINLEN:36

echo "Finished Trimmomatic at $(date)"
