#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J htseq_rna_canu
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH -o htseq_rna_canu_%j.out
#SBATCH -e htseq_rna_canu_%j.err

set -euo pipefail

echo "============================================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "============================================================"

# -----------------------------
# 1. Load modules
# -----------------------------

module purge

# Try this first. If it fails, run: module spider HTSeq
module load HTSeq || module load htseq || module load Python

echo "Loaded modules:"
module list

# -----------------------------
# 2. Basic settings
# -----------------------------

PROJECT_DIR="/home/yuch3531/genome_analysis_project"

RNA_BAM_BASE="${PROJECT_DIR}/6_differential_gene_analysis/rna_mapping/canu"
COUNT_OUT_BASE="${PROJECT_DIR}/6_differential_gene_analysis/rna_counting/canu"

GFF="${PROJECT_DIR}/5_genome_annotation/Prokka/E745_canu/E745_canu.no_fasta.gff"

mkdir -p "$COUNT_OUT_BASE"

echo "RNA BAM base: $RNA_BAM_BASE"
echo "Count output base: $COUNT_OUT_BASE"
echo "GFF annotation: $GFF"

# -----------------------------
# 3. Check input
# -----------------------------

if [[ ! -d "$RNA_BAM_BASE" ]]; then
    echo "ERROR: RNA BAM directory not found: $RNA_BAM_BASE" >&2
    exit 1
fi

if [[ ! -s "$GFF" ]]; then
    echo "ERROR: GFF file not found or empty: $GFF" >&2
    echo "Please check your Canu Prokka GFF path." >&2
    exit 1
fi

if ! command -v htseq-count >/dev/null 2>&1; then
    echo "ERROR: htseq-count command not found after module loading." >&2
    echo "Try checking available modules with: module spider HTSeq" >&2
    exit 1
fi

echo "HTSeq version:"
htseq-count --version || true

# -----------------------------
# 4. Run HTSeq-count
# -----------------------------

find "$RNA_BAM_BASE" -name "*.canu.sorted.bam" | sort | while read -r BAM; do
    SAMPLE="$(basename "$BAM" | sed -E 's/\.canu\.sorted\.bam$//')"
    DATASET="$(basename "$(dirname "$(dirname "$BAM")")")"

    SAMPLE_OUTDIR="${COUNT_OUT_BASE}/${DATASET}"
    mkdir -p "$SAMPLE_OUTDIR"

    OUTCOUNT="${SAMPLE_OUTDIR}/${SAMPLE}.canu.htseq.counts.txt"
    LOG="${SAMPLE_OUTDIR}/${SAMPLE}.canu.htseq.log"

    echo "------------------------------------------------------------"
    echo "Sample: $SAMPLE"
    echo "Dataset: $DATASET"
    echo "BAM: $BAM"
    echo "Output count: $OUTCOUNT"
    echo "Log: $LOG"
    echo "------------------------------------------------------------"

    if [[ -s "$OUTCOUNT" ]]; then
        echo "Count file already exists, skipping: $OUTCOUNT"
        continue
    fi

    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -t CDS \
        -i ID \
        "$BAM" \
        "$GFF" \
        > "$OUTCOUNT" \
        2> "$LOG"

    echo "Finished HTSeq-count for sample: $SAMPLE"
done

echo "============================================================"
echo "HTSeq-count for RNA-seq finished at: $(date)"
echo "============================================================"
