#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J htseq_tnseq_canu_bowtie2
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH -o htseq_tnseq_canu_bowtie2_%j.out
#SBATCH -e htseq_tnseq_canu_bowtie2_%j.err

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

module load HTSeq || module load htseq || module load Python

echo "Loaded modules:"
module list

# -----------------------------
# 2. Basic settings
# -----------------------------

PROJECT_DIR="/home/yuch3531/genome_analysis_project"

TNSEQ_BAM_BASE="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/bowtie2_canu_25bp_window"
COUNT_OUT_BASE="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/htseq_canu_bowtie2"

GFF="${PROJECT_DIR}/5_genome_annotation/Prokka/E745_canu/E745_canu.no_fasta.gff"

mkdir -p "$COUNT_OUT_BASE"

echo "Tn-seq BAM base: $TNSEQ_BAM_BASE"
echo "Count output base: $COUNT_OUT_BASE"
echo "GFF annotation: $GFF"

# -----------------------------
# 3. Check input
# -----------------------------

if [[ ! -d "$TNSEQ_BAM_BASE" ]]; then
    echo "ERROR: Tn-seq BAM directory not found: $TNSEQ_BAM_BASE" >&2
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

echo "Checking BAM files found:"
find "$TNSEQ_BAM_BASE" -type f -name "*.canu.bowtie2.sorted.bam" | sort | head -20

BAM_COUNT=$(find "$TNSEQ_BAM_BASE" -type f -name "*.canu.bowtie2.sorted.bam" | wc -l)

echo "Number of Bowtie2 BAM files found: $BAM_COUNT"

if [[ "$BAM_COUNT" -eq 0 ]]; then
    echo "ERROR: No Bowtie2 BAM files found." >&2
    echo "Expected pattern: *.canu.bowtie2.sorted.bam" >&2
    echo "Try checking files with:" >&2
    echo "find $TNSEQ_BAM_BASE -type f -name '*.bam'" >&2
    exit 1
fi

# -----------------------------
# 4. Optional sanity check: contig names
# -----------------------------

echo "Preview BAM contig names:"
FIRST_BAM=$(find "$TNSEQ_BAM_BASE" -type f -name "*.canu.bowtie2.sorted.bam" | sort | head -1)
samtools idxstats "$FIRST_BAM" | head || true

echo "Preview GFF contig names:"
grep -v "^#" "$GFF" | head || true

# -----------------------------
# 5. Run HTSeq-count
# -----------------------------

find "$TNSEQ_BAM_BASE" -type f -name "*.canu.bowtie2.sorted.bam" | sort | while read -r BAM
do
    SAMPLE="$(basename "$BAM" | sed -E 's/\.canu\.bowtie2\.sorted\.bam$//')"
    DATASET="$(basename "$(dirname "$(dirname "$BAM")")")"

    SAMPLE_OUTDIR="${COUNT_OUT_BASE}/${DATASET}"
    mkdir -p "$SAMPLE_OUTDIR"

    OUTCOUNT="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bowtie2.htseq.counts.txt"
    LOG="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bowtie2.htseq.log"

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
        -i locus_tag \
        "$BAM" \
        "$GFF" \
        > "$OUTCOUNT" \
        2> "$LOG"

    echo "Finished HTSeq-count for sample: $SAMPLE"

    echo "Special HTSeq summary rows:"
    grep "^__" "$OUTCOUNT" || true
done

echo "============================================================"
echo "HTSeq-count for Tn-seq finished at: $(date)"
echo "Output directory:"
echo "$COUNT_OUT_BASE"
echo "============================================================"
