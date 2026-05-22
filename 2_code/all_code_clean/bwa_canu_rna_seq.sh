#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J bwa_canu_mapping_all_resume
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH --mail-type=ALL
#SBATCH -o bwa_canu_mapping_all_resume_%j.out
#SBATCH -e bwa_canu_mapping_all_resume_%j.err

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

THREADS=2

REF="/home/yuch3531/genome_analysis_project/4_genome_assembly/canu_4935797/E745_canu.contigs.fasta"

BASE_INPUT_DIR="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data"
BASE_OUTPUT_DIR="/home/yuch3531/genome_analysis_project/6_differential_gene_analysis/rna_mapping/canu"

echo "Using $THREADS threads"
echo "Starting BWA mapping to Canu assembly"
echo "Reference fasta: $REF"
echo "Base input directory: $BASE_INPUT_DIR"
echo "Base output directory: $BASE_OUTPUT_DIR"

# -----------------------------
# 3. Check reference and BWA index
# -----------------------------

if [[ ! -s "$REF" ]]; then
    echo "ERROR: Reference fasta does not exist or is empty: $REF" >&2
    exit 1
fi

if [[ ! -s "${REF}.bwt" ]]; then
    echo "BWA index not found. Building BWA index for: $REF"
    bwa index "$REF"
else
    echo "BWA index already exists for: $REF"
fi

# -----------------------------
# 4. Function: map one dataset
# -----------------------------

map_dataset () {
    local DATASET="$1"
    local INPUT_DIR="${BASE_INPUT_DIR}/${DATASET}/trimmed"
    local DATASET_OUTDIR="${BASE_OUTPUT_DIR}/${DATASET}"

    echo "============================================================"
    echo "Processing dataset: $DATASET"
    echo "Input directory: $INPUT_DIR"
    echo "Output directory: $DATASET_OUTDIR"
    echo "Reference: canu"
    echo "============================================================"

    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "WARNING: Input directory does not exist, skipping dataset: $INPUT_DIR" >&2
        return 0
    fi

    mkdir -p "$DATASET_OUTDIR"

    shopt -s nullglob

    local R1_FILES=("${INPUT_DIR}"/trim_paired_*_1.fastq.gz)

    if [[ ${#R1_FILES[@]} -eq 0 ]]; then
        echo "WARNING: No R1 files found in: $INPUT_DIR" >&2
        return 0
    fi

    for R1 in "${R1_FILES[@]}"; do
        local R2="${R1/_1.fastq.gz/_2.fastq.gz}"

        if [[ ! -s "$R2" ]]; then
            echo "ERROR: R2 file not found or empty for R1: $R1" >&2
            echo "Expected R2: $R2" >&2
            exit 1
        fi

        local R1_BASENAME
        R1_BASENAME=$(basename "$R1")

        local SAMPLE
        SAMPLE="${R1_BASENAME#trim_paired_}"
        SAMPLE="${SAMPLE%_1.fastq.gz}"

        local SAMPLE_OUTDIR="${DATASET_OUTDIR}/${SAMPLE}"
        mkdir -p "$SAMPLE_OUTDIR"

        local OUTBAM="${SAMPLE_OUTDIR}/${SAMPLE}.canu.sorted.bam"
        local BWA_LOG="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bwa.log"

        echo "------------------------------------------------------------"
        echo "Sample: $SAMPLE"
        echo "R1: $R1"
        echo "R2: $R2"
        echo "Output BAM: $OUTBAM"
        echo "BWA log: $BWA_LOG"
        echo "------------------------------------------------------------"

        # -----------------------------
        # 4.1 Skip valid existing BAM
        # -----------------------------

        if [[ -s "$OUTBAM" ]] && samtools quickcheck "$OUTBAM" 2>/dev/null; then
            echo "Valid BAM already exists. Skipping mapping for sample: $SAMPLE"

            if [[ ! -s "${OUTBAM}.bai" ]]; then
                echo "BAM index missing. Creating index for: $OUTBAM"
                samtools index "$OUTBAM"
            else
                echo "BAM index already exists: ${OUTBAM}.bai"
            fi

            echo "Finished existing valid sample: $SAMPLE"
            continue
        fi

        # -----------------------------
        # 4.2 Remove corrupted/incomplete BAM
        # -----------------------------

        if [[ -e "$OUTBAM" ]]; then
            echo "Existing BAM is incomplete or corrupted. Removing:"
            echo "$OUTBAM"
            rm -f "$OUTBAM" "${OUTBAM}.bai"
        fi

        # Remove possible old temporary files
        rm -f "${SAMPLE_OUTDIR}/${SAMPLE}.sorttmp"*

        # -----------------------------
        # 4.3 Temporary directory
        # -----------------------------

        local SORT_TMP_BASE

        if [[ -n "${SNIC_TMP:-}" && -d "${SNIC_TMP:-}" ]]; then
            SORT_TMP_BASE="${SNIC_TMP}/${SAMPLE}.sorttmp"
        else
            mkdir -p "${SAMPLE_OUTDIR}/tmp"
            SORT_TMP_BASE="${SAMPLE_OUTDIR}/tmp/${SAMPLE}.sorttmp"
        fi

        echo "Sort temporary prefix: $SORT_TMP_BASE"

        # -----------------------------
        # 4.4 BWA mapping + SAMtools sort
        # -----------------------------

        echo "Running BWA mem and SAMtools sort for sample: $SAMPLE"

        bwa mem \
            -t "$THREADS" \
            "$REF" \
            "$R1" \
            "$R2" \
            2> "$BWA_LOG" | \
        samtools sort \
            -@ "$THREADS" \
            -m 1G \
            -T "$SORT_TMP_BASE" \
            -o "$OUTBAM" \
            -

        # -----------------------------
        # 4.5 Check BAM integrity
        # -----------------------------

        echo "Checking BAM integrity for sample: $SAMPLE"

        if ! samtools quickcheck -v "$OUTBAM"; then
            echo "ERROR: BAM file is incomplete or corrupted after mapping: $OUTBAM" >&2
            rm -f "$OUTBAM" "${OUTBAM}.bai"
            exit 1
        fi

        # -----------------------------
        # 4.6 Index BAM
        # -----------------------------

        echo "Indexing BAM for sample: $SAMPLE"
        samtools index "$OUTBAM"

        echo "Finished sample: $SAMPLE"
    done
}

# -----------------------------
# 5. Run all datasets
# -----------------------------

map_dataset "RNA-Seq_BH"
map_dataset "RNA-Seq_Serum"
# map_dataset "Tn-Seq_BHI"
# map_dataset "Tn-Seq_HSerum"
# map_dataset "Tn-Seq_Serum"

# -----------------------------
# 6. Final global check
# -----------------------------

echo "============================================================"
echo "Running final BAM integrity check for all BAM files"
echo "============================================================"

find "$BASE_OUTPUT_DIR" -name "*.sorted.bam" -print0 | xargs -0 samtools quickcheck -v

echo "============================================================"
echo "Checking missing BAM index files"
echo "============================================================"

find "$BASE_OUTPUT_DIR" -name "*.sorted.bam" | while read -r bam; do
    if [[ ! -s "${bam}.bai" ]]; then
        echo "Missing BAM index: ${bam}.bai"
        exit 1
    fi
done

echo "============================================================"
echo "All available BAM files passed quickcheck and have index files."
echo "Job finished successfully at: $(date)"
echo "============================================================"
