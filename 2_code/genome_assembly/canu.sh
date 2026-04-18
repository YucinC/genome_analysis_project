#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J canu_asm
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/4_genome_assembly/canu_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/4_genome_assembly/canu_%j.err

set -euo pipefail

module load canu
module load SAMtools

PACBIO_DIR="/home/yuch3531/genome_analysis_project/1_data/genomics_data/PacBio"
ASSEMBLY_BASE="/home/yuch3531/genome_analysis_project/4_genome_assembly"
OUTDIR="${ASSEMBLY_BASE}/canu_${SLURM_JOB_ID}"
PREFIX="E745_canu"
GENOME_SIZE="2.8m"

PACBIO_MERGED="${ASSEMBLY_BASE}/pacbio_all_${SLURM_JOB_ID}.fastq.gz"

mkdir -p "${ASSEMBLY_BASE}"
mkdir -p "${OUTDIR}"

echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "samtools path: $(which samtools)"

echo "Merging PacBio reads..."
cat "${PACBIO_DIR}"/*.fastq.gz > "${PACBIO_MERGED}"

echo "Running Canu with PacBio reads..."
canu \
  -p "${PREFIX}" \
  -d "${OUTDIR}" \
  genomeSize="${GENOME_SIZE}" \
  useGrid=false \
  maxThreads=2 \
  corThreads=2 \
  corConcurrency=1 \
  stopOnLowCoverage=2 \
  minInputCoverage=1 \
  samtools="$(which samtools)" \
  -pacbio "${PACBIO_MERGED}"

echo "Job finished at: $(date)"
