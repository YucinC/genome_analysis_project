#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J bowtie2_tnseq_25bp
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH -o bowtie2_tnseq_25bp_%j.out
#SBATCH -e bowtie2_tnseq_25bp_%j.err

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
module load Bowtie2
module load SAMtools/1.22.1-GCC-13.3.0

echo "Loaded modules:"
module list

# -----------------------------
# 2. Basic settings
# -----------------------------

THREADS="${SLURM_CPUS_PER_TASK:-2}"

PROJECT_DIR="/home/yuch3531/genome_analysis_project"

REF="${PROJECT_DIR}/4_genome_assembly/canu_4935797/E745_canu.contigs.fasta"

BASE_INPUT_DIR="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/trimmed_magellan6_barcode"

BASE_OUTPUT_DIR="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/bowtie2_canu_25bp_window"

INDEX_DIR="${BASE_OUTPUT_DIR}/bowtie2_index"
INDEX_PREFIX="${INDEX_DIR}/E745_canu"

WINDOW_SIZE=25

DATASETS=(
    "Tn-Seq_BHI"
    "Tn-Seq_HSerum"
    "Tn-Seq_Serum"
)

echo "Using $THREADS threads"
echo "Reference fasta: $REF"
echo "Input base directory: $BASE_INPUT_DIR"
echo "Output base directory: $BASE_OUTPUT_DIR"
echo "Window size: ${WINDOW_SIZE} bp"

mkdir -p "$BASE_OUTPUT_DIR"
mkdir -p "$INDEX_DIR"

# -----------------------------
# 3. Check reference and build index
# -----------------------------

if [[ ! -s "$REF" ]]; then
    echo "ERROR: Reference fasta does not exist or is empty: $REF" >&2
    exit 1
fi

if [[ ! -s "${REF}.fai" ]]; then
    echo "Reference fasta index not found. Creating .fai:"
    samtools faidx "$REF"
else
    echo "Reference fasta index already exists: ${REF}.fai"
fi

if [[ ! -s "${INDEX_PREFIX}.1.bt2" && ! -s "${INDEX_PREFIX}.1.bt2l" ]]; then
    echo "Bowtie2 index not found. Building Bowtie2 index for: $REF"
    bowtie2-build "$REF" "$INDEX_PREFIX"
else
    echo "Bowtie2 index already exists."
fi

# -----------------------------
# 4. Create 25 bp genome windows
# -----------------------------

WINDOW_BED="${BASE_OUTPUT_DIR}/E745_canu.${WINDOW_SIZE}bp_windows.bed"

echo "Creating ${WINDOW_SIZE} bp genome windows:"
echo "$WINDOW_BED"

awk -v W="$WINDOW_SIZE" '
{
    contig=$1
    len=$2
    for (start=0; start<len; start+=W) {
        end=start+W
        if (end>len) end=len
        print contig "\t" start "\t" end "\t" contig ":" start+1 "-" end
    }
}
' "${REF}.fai" > "$WINDOW_BED"

# -----------------------------
# 5. Helper function: clean sample name
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
        | sed -E 's/\.cut_5bp6_magellan6$//' \
        | sed -E 's/\.cut_barcode_magellan6$//' \
        | sed -E 's/_pass$//' \
        | sed -E 's/[^A-Za-z0-9_.-]+/_/g')"

    echo "$SAMPLE"
}

# -----------------------------
# 6. Map and count one dataset
# -----------------------------

map_dataset () {
    local DATASET="$1"

    local INPUT_DIR="${BASE_INPUT_DIR}/${DATASET}"
    local DATASET_OUTDIR="${BASE_OUTPUT_DIR}/${DATASET}"

    echo "============================================================"
    echo "Processing dataset: $DATASET"
    echo "Input directory: $INPUT_DIR"
    echo "Output directory: $DATASET_OUTDIR"
    echo "============================================================"

    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERROR: Input directory does not exist: $INPUT_DIR" >&2
        exit 1
    fi

    mkdir -p "$DATASET_OUTDIR"

    mapfile -t FASTQ_FILES < <(
        find "$INPUT_DIR" -maxdepth 1 -type f \
            \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
            | sort
    )

    echo "Number of FASTQ files found: ${#FASTQ_FILES[@]}"

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No FASTQ files found in: $INPUT_DIR" >&2
        ls -lh "$INPUT_DIR" >&2
        exit 1
    fi

    for FASTQ in "${FASTQ_FILES[@]}"
    do
        local SAMPLE
        local SAMPLE_OUTDIR
        local OUTBAM
        local BOWTIE2_LOG
        local FLAGSTAT
        local IDXSTATS
        local SORT_TMP_BASE
        local RAW_WINDOW_COUNTS
        local WINDOW_COUNTS

        SAMPLE="$(sample_name_from_fastq "$FASTQ")"

        SAMPLE_OUTDIR="${DATASET_OUTDIR}/${SAMPLE}"
        mkdir -p "$SAMPLE_OUTDIR"

        OUTBAM="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bowtie2.sorted.bam"
        BOWTIE2_LOG="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bowtie2.log"
        FLAGSTAT="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bowtie2.flagstat.txt"
        IDXSTATS="${SAMPLE_OUTDIR}/${SAMPLE}.canu.bowtie2.idxstats.txt"

        RAW_WINDOW_COUNTS="${SAMPLE_OUTDIR}/${SAMPLE}.${WINDOW_SIZE}bp_window.raw_counts.tsv"
        WINDOW_COUNTS="${SAMPLE_OUTDIR}/${SAMPLE}.${WINDOW_SIZE}bp_window.counts.tsv"

        echo "------------------------------------------------------------"
        echo "Sample: $SAMPLE"
        echo "FASTQ: $FASTQ"
        echo "Output BAM: $OUTBAM"
        echo "Window count output: $WINDOW_COUNTS"
        echo "------------------------------------------------------------"

        if [[ -n "${SNIC_TMP:-}" && -d "${SNIC_TMP:-}" ]]; then
            SORT_TMP_BASE="${SNIC_TMP}/${SAMPLE}.sorttmp"
        else
            mkdir -p "${SAMPLE_OUTDIR}/tmp"
            SORT_TMP_BASE="${SAMPLE_OUTDIR}/tmp/${SAMPLE}.sorttmp"
        fi

        echo "Running Bowtie2 mapping for sample: $SAMPLE"

        bowtie2 \
            --end-to-end \
            -N 0 \
            -L 8 \
            -x "$INDEX_PREFIX" \
            -U "$FASTQ" \
            -p "$THREADS" \
            2> "$BOWTIE2_LOG" | \
        samtools sort \
            -@ "$THREADS" \
            -m 1G \
            -T "$SORT_TMP_BASE" \
            -o "$OUTBAM" \
            -

        echo "Checking BAM integrity for sample: $SAMPLE"

        if ! samtools quickcheck -v "$OUTBAM"; then
            echo "ERROR: BAM file is incomplete or corrupted: $OUTBAM" >&2
            rm -f "$OUTBAM" "${OUTBAM}.bai"
            exit 1
        fi

        echo "Indexing BAM for sample: $SAMPLE"
        samtools index "$OUTBAM"

        echo "Generating mapping statistics for sample: $SAMPLE"
        samtools flagstat "$OUTBAM" > "$FLAGSTAT"
        samtools idxstats "$OUTBAM" > "$IDXSTATS"

        echo "Counting mapped reads in ${WINDOW_SIZE} bp windows for sample: $SAMPLE"

        samtools view -F 4 "$OUTBAM" \
            | awk -v W="$WINDOW_SIZE" '
                {
                    contig=$3
                    pos=$4
                    start=int((pos-1)/W)*W
                    end=start+W
                    key=contig "\t" start "\t" end
                    count[key]++
                }
                END {
                    for (key in count) {
                        print key "\t" count[key]
                    }
                }
            ' > "$RAW_WINDOW_COUNTS"

        echo -e "condition\tsample\tcontig\twindow_start_0based\twindow_end_0based\twindow_name\tcount" > "$WINDOW_COUNTS"

        awk -v condition="$DATASET" -v sample="$SAMPLE" '
            NR==FNR {
                key=$1 "\t" $2 "\t" $3
                count[key]=$4
                next
            }
            {
                key=$1 "\t" $2 "\t" $3
                c=(key in count ? count[key] : 0)
                print condition "\t" sample "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" c
            }
        ' "$RAW_WINDOW_COUNTS" "$WINDOW_BED" >> "$WINDOW_COUNTS"

        echo "Finished sample: $SAMPLE"
    done
}

# -----------------------------
# 7. Run all datasets
# -----------------------------

for DATASET in "${DATASETS[@]}"
do
    map_dataset "$DATASET"
done

# -----------------------------
# 8. Merge all 25 bp window count files
# -----------------------------

MERGED_COUNTS="${BASE_OUTPUT_DIR}/tnseq_canu_bowtie2_${WINDOW_SIZE}bp_window_counts.all_samples.tsv"

echo "Merging all ${WINDOW_SIZE} bp window count files:"
echo "$MERGED_COUNTS"

FIRST_FILE=1

find "$BASE_OUTPUT_DIR" -type f -name "*.${WINDOW_SIZE}bp_window.counts.tsv" | sort | while read -r COUNT_FILE
do
    if [[ "$FIRST_FILE" -eq 1 ]]; then
        cat "$COUNT_FILE" > "$MERGED_COUNTS"
        FIRST_FILE=0
    else
        tail -n +2 "$COUNT_FILE" >> "$MERGED_COUNTS"
    fi
done

# -----------------------------
# 9. Final check
# -----------------------------

echo "============================================================"
echo "Final BAM integrity check"
echo "============================================================"

find "$BASE_OUTPUT_DIR" -name "*.sorted.bam" -print0 | xargs -0 -r samtools quickcheck -v

echo "============================================================"
echo "Job finished successfully at: $(date)"
echo "Output directory:"
echo "$BASE_OUTPUT_DIR"
echo "Merged window counts:"
echo "$MERGED_COUNTS"
echo "============================================================"
