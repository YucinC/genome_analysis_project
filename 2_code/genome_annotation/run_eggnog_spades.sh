#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 15:00:00
#SBATCH -J eggnog_spades_hmmer
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/5_genome_annotation/eggNOG/eggnog_spades_hmmer_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/5_genome_annotation/eggNOG/eggnog_spades_hmmer_%j.err

set -euo pipefail

module load eggnog-mapper

THREADS=${SLURM_NTASKS:-2}

PROKKA_DIR="/home/yuch3531/genome_analysis_project/5_genome_annotation/Prokka/E745_spades"
PROTEIN_FILE="${PROKKA_DIR}/E745_spades.faa"

# Shared eggNOG database on UPPMAX.
# This directory must directly contain eggnog.db and eggnog.taxa.db.
EGGNOG_DATA_BASE="/sw/data/eggNOG/5.0"

# The public eggNOG data are stored in cluster-named subdirectories.
# Pelle may not have its own subdirectory here, so use an existing one with eggnog.db.
EGGNOG_DATA=""

for candidate in pelle rackham snowy miarka bianca; do
    if [[ -s "${EGGNOG_DATA_BASE}/${candidate}/eggnog.db" ]]; then
        EGGNOG_DATA="${EGGNOG_DATA_BASE}/${candidate}"
        break
    fi
done

if [[ -z "$EGGNOG_DATA" ]]; then
    echo "ERROR: Could not find eggnog.db under:"
    echo "$EGGNOG_DATA_BASE"
    echo "Checked candidate directories: pelle rackham snowy miarka bianca"
    exit 1
fi

OUTDIR="/home/yuch3531/genome_analysis_project/5_genome_annotation/eggNOG/E745_spades"
PREFIX="E745_spades_eggnog_hmmer"

mkdir -p "$OUTDIR"

echo "Checking input file..."
echo "Protein file: $PROTEIN_FILE"

if [[ ! -s "$PROTEIN_FILE" ]]; then
    echo "ERROR: Protein file not found or empty:"
    echo "$PROTEIN_FILE"
    exit 1
fi

echo "Checking eggNOG database directory..."
echo "eggNOG data directory: $EGGNOG_DATA"

if [[ ! -s "${EGGNOG_DATA}/eggnog.db" ]]; then
    echo "ERROR: eggnog.db not found:"
    echo "${EGGNOG_DATA}/eggnog.db"
    exit 1
fi

if [[ ! -s "${EGGNOG_DATA}/eggnog.taxa.db" ]]; then
    echo "ERROR: eggnog.taxa.db not found:"
    echo "${EGGNOG_DATA}/eggnog.taxa.db"
    exit 1
fi

echo "Input protein file and eggNOG databases found."
echo "Running eggNOG-mapper functional annotation"
echo "Search mode: HMMER"
echo "Target database taxID: 2, Bacteria"
echo "Output directory: $OUTDIR"
echo "Output prefix: $PREFIX"
echo "Threads: $THREADS"

emapper.py \
    -i "$PROTEIN_FILE" \
    --itype proteins \
    -m hmmer \
    -d 2 \
    --data_dir "$EGGNOG_DATA" \
    --cpu "$THREADS" \
    --output "$PREFIX" \
    --output_dir "$OUTDIR" \
    --override

echo "eggNOG-mapper HMMER annotation finished."
echo "Output files:"
ls -lh "$OUTDIR"
