#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J prokka_spades
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/5_genome_annotation/Prokka/prokka_spades_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/5_genome_annotation/Prokka/prokka_spades_%j.err

set -euo pipefail

module load prokka

THREADS=${SLURM_NTASKS:-2}

ASSEMBLY="/home/yuch3531/genome_analysis_project/4_genome_assembly/spades_4947505/spades_contigs.fasta"

OUTDIR="/home/yuch3531/genome_analysis_project/5_genome_annotation/Prokka/E745_spades"
PREFIX="E745_spades"

mkdir -p "$OUTDIR"

if [[ ! -s "$ASSEMBLY" ]]; then
    echo "ERROR: Assembly file not found or empty:"
    echo "$ASSEMBLY"
    exit 1
fi

echo "Running Prokka annotation"
echo "Input assembly: $ASSEMBLY"
echo "Output directory: $OUTDIR"
echo "Threads: $THREADS"

prokka \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX" \
  --kingdom Bacteria \
  --genus Enterococcus \
  --species faecium \
  --strain E745 \
  --locustag E745S \
  --gcode 11 \
  --cpus "$THREADS" \
  --addgenes \
  --force \
  "$ASSEMBLY"

echo "Prokka annotation finished"
echo "Main output files:"
ls -lh "$OUTDIR"
