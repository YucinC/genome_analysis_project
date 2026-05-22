#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 6:00:00
#SBATCH -J spades_asm
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/4_genome_assembly/spades_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/4_genome_assembly/spades_%j.err

set -euo pipefail

module load SPAdes

ILLUMINA_DIR="/home/yuch3531/genome_analysis_project/1_data/genomics_data/Illumina"
NANOPORE_READS="/home/yuch3531/genome_analysis_project/1_data/genomics_data/Nanopore/E745_all.fasta.gz"
ASSEMBLY_DIR="/home/yuch3531/genome_analysis_project/4_genome_assembly"

ILLUMINA_R1="${ILLUMINA_DIR}/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz"
ILLUMINA_R2="${ILLUMINA_DIR}/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz"

OUTDIR="${ASSEMBLY_DIR}/spades_${SLURM_JOB_ID}"

mkdir -p "${ASSEMBLY_DIR}"
mkdir -p "${OUTDIR}"

echo "Job started at: $(date)"
echo "Running on: $(hostname)"

echo "Running SPAdes with Illumina + Nanopore..."
spades.py \
  --isolate \
  -1 "${ILLUMINA_R1}" \
  -2 "${ILLUMINA_R2}" \
  --nanopore "${NANOPORE_READS}" \
  -t 2 \
  -o "${OUTDIR}"

echo "Job finished at: $(date)"
