#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J bwa_canu_tnseq_single
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH -o bwa_canu_tnseq_single_%j.out
#SBATCH -e bwa_canu_tnseq_single_%j.err

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
module load BWA/0.7.19-GCCcore-13.3.0
module load SAMtools/1.22.1-GCC-13.3.0

echo "Loaded modules:"
module list

# -----------------------------
# 2. Basic settings
# -----------------------------

THREADS="${SLURM_CPUS_PER_TASK:-2}"

PROJECT_DIR="/home/yuch3531/genome_analysis_project"

REF="${PROJECT_DIR}/4_genome_assembly/canu_4935797/E745_canu.contigs.fasta"

BASE_INPUT_DIR="$/6_differential_gene_analysis/tnseq/trimmed_magellan6_barcode"
BASE_OUTPUT_DIR="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/bwa_mapping"

DATASETS=(
    "Tn-Seq_BHI"
    "Tn-Seq_HSerum"
    "Tn-Seq_Serum"
)

echo "Using $THREADS threads"
echo "Starting single-end BWA mapping of Tn-seq reads to Canu assembly"
echo "Reference fasta: $REF"
echo "Base input directory: $BASE_INPUT_DIR"
echo "Base output directory: $BASE_OUTPUT_DIR"

mkdir -p "$BASE_OUTPUT_DIR"

# -----------------------------
# 3. Check reference and BWA index
# -----------------------------

if [[ ! -s "$REF" ]]; then
    echo "ERROR: Reference fasta does not exist or is empty: $REF" >&2
    exit 1
fi

if [[ ! -s "${REF}.bwt" || ! -s "${REF}.sa" ]]; then
    echo "BWA index not found. Building BWA index for: $REF"
    bwa index "$REF"
else
    echo "BWA index already exists for: $REF"
fi

# -----------------------------
# 4. Helper function: sample name
# -----------------------------

sample_name_from_fastq () {
    local FASTQ="$1"
    local SAMPLE

    SAMPLE="$(basename "$FASTQ")"

    SAMPLE="${SAMPLE%.fastq.gz}"
    SAMPLE="${SAMPLE%.fq.gz}"
    SAMPLE="${SAMPLE%.fastq}"
    SAMPLE="${SAMPLE%.fq}"

    SAMPLE="$(echo "$SAMPLE" \
        | sed -E 's/^trim_//' \
        | sed -E 's/^trimmed_//' \
        | sed -E 's/[^A-Za-z0-9_.-]+/_/g')"

    echo "$SAMPLE"
}

# -----------------------------
# 5. Function: map one Tn-seq dataset
# -----------------------------

map_dataset () {
    local DATASET="$1"

    # Your Tn-seq files are directly inside the dataset folder,
    # not inside a trimmed/ subfolder.
    local INPUT_DIR="${BASE_INPUT_DIR}/${DATASET}"
    local DATASET_OUTDIR="${BASE_OUTPUT_DIR}/${DATASET}"

    echo "============================================================"
    echo "Processing dataset: $DATASET"
    echo "Input directory: $INPUT_DIR"
    echo "Output directory: $DATASET_OUTDIR"
    echo "Reference: canu"
    echo "Read type: single-end"
    echo "============================================================"

    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERROR: Input directory does not exist: $INPUT_DIR" >&2
        exit 1
    fi

    mkdir -p "$DATASET_OUTDIR"

    mapfile -t FASTQ_FILES < <(
        find "$INPUT_DIR" -maxdepth 1 -type f \
            \( -name "trim_*.fastq.gz" -o -name "trim_*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) \
            ! -name "*paired*" \
            ! -name "*unpaired*" \
            ! -name "*single*" \
            | sort
    )

    echo "Number of FASTQ files found: ${#FASTQ_FILES[@]}"

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in: $INPUT_DIR" >&2
        echo "Files in this directory:" >&2
        ls -lh "$INPUT_DIR" >&2
        exit 1
    fi

    for FASTQ in "${FASTQ_FILES[@]}"; do
        local SAMPLE
        local SAMPLE_OUTDIR
        local OUTBAM
        local BWA_LOG
        local FLAGSTAT
        local IDXSTATS
        local SORT_TMP_BASE

        SAMPLE="$(sample_name_from_fastq "$FASTQ")"

        SAMPLE_OUTDIR="${DATASET_OUTDIR}/${SAMPLE}"
        mkdir -p "$SAMPLE_OUTDIR"

        OUTBAM="${SAMPLE_OUTDIR}/${SAMPLE}.canu.sorted.bam"
        BWA_LOG="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bwa.log"
        FLAGSTAT="${SAMPLE_OUTDIR}/${SAMPLE}.canu.flagstat.txt"
        IDXSTATS="${SAMPLE_OUTDIR}/${SAMPLE}.canu.idxstats.txt"

        echo "------------------------------------------------------------"
        echo "Sample: $SAMPLE"
        echo "FASTQ: $FASTQ"
        echo "Output BAM: $OUTBAM"
        echo "BWA log: $BWA_LOG"
        echo "------------------------------------------------------------"

        # Skip valid existing BAM
        if [[ -s "$OUTBAM" ]] && samtools quickcheck "$OUTBAM" 2>/dev/null; then
            echo "Valid BAM already exists. Skipping mapping for sample: $SAMPLE"

            if [[ ! -s "${OUTBAM}.bai" ]]; then
                echo "BAM index missing. Creating index for: $OUTBAM"
                samtools index "$OUTBAM"
            fi

            if [[ ! -s "$FLAGSTAT" ]]; then
                samtools flagstat "$OUTBAM" > "$FLAGSTAT"
            fi

            if [[ ! -s "$IDXSTATS" ]]; then
                samtools idxstats "$OUTBAM" > "$IDXSTATS"
            fi

            echo "Finished existing valid sample: $SAMPLE"
            continue
        fi

        # Remove corrupted/incomplete BAM
        if [[ -e "$OUTBAM" ]]; then
            echo "Existing BAM is incomplete or corrupted. Removing:"
            echo "$OUTBAM"
            rm -f "$OUTBAM" "${OUTBAM}.bai"
        fi

        rm -f "${SAMPLE_OUTDIR}/${SAMPLE}.sorttmp"*

        if [[ -n "${SNIC_TMP:-}" && -d "${SNIC_TMP:-}" ]]; then
            SORT_TMP_BASE="${SNIC_TMP}/${SAMPLE}.sorttmp"
        else
            mkdir -p "${SAMPLE_OUTDIR}/tmp"
            SORT_TMP_BASE="${SAMPLE_OUTDIR}/tmp/${SAMPLE}.sorttmp"
        fi

        echo "Sort temporary prefix: $SORT_TMP_BASE"
        echo "Running single-end BWA mem and SAMtools sort for sample: $SAMPLE"

        bwa mem \
            -t "$THREADS" \
            "$REF" \
            "$FASTQ" \
            2> "$BWA_LOG" | \
        samtools sort \
            -@ "$THREADS" \
            -m 1G \
            -T "$SORT_TMP_BASE" \
            -o "$OUTBAM" \
            -

        echo "Checking BAM integrity for sample: $SAMPLE"

        if ! samtools quickcheck -v "$OUTBAM"; then
            echo "ERROR: BAM file is incomplete or corrupted after mapping: $OUTBAM" >&2
            rm -f "$OUTBAM" "${OUTBAM}.bai"
            exit 1
        fi

        echo "Indexing BAM for sample: $SAMPLE"
        samtools index "$OUTBAM"

        echo "Generating mapping statistics for sample: $SAMPLE"
        samtools flagstat "$OUTBAM" > "$FLAGSTAT"
        samtools idxstats "$OUTBAM" > "$IDXSTATS"

        echo "Finished sample: $SAMPLE"
    done
}

# -----------------------------
# 6. Run all Tn-seq datasets
# -----------------------------

for DATASET in "${DATASETS[@]}"; do
    map_dataset "$DATASET"
done

# -----------------------------
# 7. Final global check
# -----------------------------

echo "============================================================"
echo "Running final BAM integrity check for all BAM files"
echo "============================================================"

find "$BASE_OUTPUT_DIR" -name "*.sorted.bam" -print0 | xargs -0 -r samtools quickcheck -v

echo "============================================================"
echo "Checking missing BAM index files"
echo "============================================================"

find "$BASE_OUTPUT_DIR" -name "*.sorted.bam" | while read -r BAM; do
    if [[ ! -s "${BAM}.bai" ]]; then
        echo "Missing BAM index: ${BAM}.bai"
        exit 1
    fi
done

echo "============================================================"
echo "All available BAM files passed quickcheck and have index files."
echo "Job finished successfully at: $(date)"
echo "============================================================"
