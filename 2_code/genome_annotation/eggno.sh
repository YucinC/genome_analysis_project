#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J eggnog_spades_hmmer
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/5_genome_annotation/eggNOG/eggnog_spades_hmmer_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/5_genome_annotation/eggNOG/eggnog_spades_hmmer_%j.err

set -euo pipefail

module load eggnog-mapper

THREADS=${SLURM_NTASKS:-2}

PROKKA_DIR="/home/yuch3531/genome_analysis_project/5_genome_annotation/Prokka/E745_spades"
PROTEIN_FILE="${PROKKA_DIR}/E745_spades.faa"

OUTDIR="/home/yuch3531/genome_analysis_project/5_genome_annotation/eggNOG/E745_spades"
PREFIX="E745_spades_eggnog_hmmer"

mkdir -p "$OUTDIR"

echo "Checking input file..."
echo "Protein file: $PROTEIN_FILE"

if [[ ! -s "$PROTEIN_FILE" ]]; then
    echo "ERROR: Protein file not found or empty:"
    echo "$PROTEIN_FILE"
    echo "Please check the Prokka output directory."
    exit 1
fi

echo "Input protein file found."
echo "Running eggNOG-mapper functional annotation"
echo "Search mode: HMMER"
echo "Target database: bact"
echo "Output directory: $OUTDIR"
echo "Output prefix: $PREFIX"
echo "Threads: $THREADS"

emapper.py \
    -i "$PROTEIN_FILE" \
    --itype proteins \
    -m hmmer \
    -d bact \
    --cpu "$THREADS" \
    --output "$PREFIX" \
    --output_dir "$OUTDIR" \
    --override

echo "eggNOG-mapper HMMER annotation finished."
echo "Output files:"
ls -lh "$OUTDIR"
