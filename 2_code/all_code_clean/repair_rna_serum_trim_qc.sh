#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH -J repair_rna_serum_trim_qc
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/repair_rna_serum_trim_qc_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/repair_rna_serum_trim_qc_%j.err

set -euo pipefail
shopt -s nullglob

module load Trimmomatic
module load FastQC
module load MultiQC

threads="${SLURM_CPUS_PER_TASK:-2}"

# -----------------------------
# Input directories
# -----------------------------
SERUM_RAW="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_Serum/raw"

# -----------------------------
# Output directories
# -----------------------------
OUT_BASE="/home/yuch3531/genome_analysis_project/3_preprocessing/trimming"
SERUM_OUT="${OUT_BASE}/rna_serum"
FASTQC_OUT="/home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc/rna_serum"
MULTIQC_OUT="/home/yuch3531/genome_analysis_project/3_preprocessing/re_multiqc/rna_serum"

# Temporary rebuild directory
TMP_DIR="${OUT_BASE}/tmp_rna_serum_repair_${SLURM_JOB_ID}"

mkdir -p "${SERUM_OUT}" "${FASTQC_OUT}" "${MULTIQC_OUT}" "${TMP_DIR}"

trap 'rm -rf "${TMP_DIR}"' EXIT

# -----------------------------
# Try to find adapter file
# -----------------------------
ADAPTERS=""
for p in \
    "${EBROOTTRIMMOMATIC:-}/adapters/TruSeq3-PE.fa" \
    "${TRIMMOMATIC_HOME:-}/adapters/TruSeq3-PE.fa" \
    "${TRIMMOMATIC_ROOT:-}/adapters/TruSeq3-PE.fa"
do
    if [[ -n "$p" && -f "$p" ]]; then
        ADAPTERS="$p"
        break
    fi
done

if [[ -z "${ADAPTERS}" ]]; then
    echo "ERROR: Could not find TruSeq3-PE.fa adapter file automatically."
    echo "Please run: module spider Trimmomatic"
    echo "Then locate the adapters/TruSeq3-PE.fa file and set ADAPTERS manually in this script."
    exit 1
fi

echo "Using adapter file: ${ADAPTERS}"
echo "Starting repair workflow at $(date)"

# -----------------------------
# Same trimming parameters as your original script
# -----------------------------
TRIM_PARAMS=(
    "ILLUMINACLIP:${ADAPTERS}:2:30:10"
    "TRAILING:20"
    "SLIDINGWINDOW:4:20"
    "MINLEN:36"
)

# Only rebuild the broken samples
SAMPLES=("ERR1797969" "ERR1797970")

check_gzip_ok() {
    local f="$1"
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing file: $f"
        exit 1
    fi
    echo "Checking gzip integrity: $f"
    gzip -t "$f"
}

validate_fastq_gz() {
    local f="$1"

    if [[ ! -s "$f" ]]; then
        echo "ERROR: Output file missing or empty: $f"
        exit 1
    fi

    gzip -t "$f"

    local mod4
    mod4=$(zcat "$f" 2>/dev/null | awk 'END{print NR % 4}')

    if [[ "$mod4" -ne 0 ]]; then
        echo "ERROR: FASTQ structure broken for $f (line count mod4 = $mod4)"
        exit 1
    fi

    echo "Validated: $f"
}

echo "======================================"
echo "Raw input:     ${SERUM_RAW}"
echo "Trim output:   ${SERUM_OUT}"
echo "FastQC output: ${FASTQC_OUT}"
echo "MultiQC out:   ${MULTIQC_OUT}"
echo "Temp dir:      ${TMP_DIR}"
echo "======================================"

# -----------------------------
# Step 1: check raw files
# -----------------------------
echo
echo "===== Step 1: check raw files ====="

for SAMPLE in "${SAMPLES[@]}"
do
    IN1="${SERUM_RAW}/${SAMPLE}_1.fastq.gz"
    IN2="${SERUM_RAW}/${SAMPLE}_2.fastq.gz"

    echo "Checking raw sample: ${SAMPLE}"
    check_gzip_ok "${IN1}"
    check_gzip_ok "${IN2}"
done

# -----------------------------
# Step 2: retrim to temp directory
# -----------------------------
echo
echo "===== Step 2: retrim broken samples to temp directory ====="

for SAMPLE in "${SAMPLES[@]}"
do
    IN1="${SERUM_RAW}/${SAMPLE}_1.fastq.gz"
    IN2="${SERUM_RAW}/${SAMPLE}_2.fastq.gz"

    P1="${TMP_DIR}/trim_paired_${SAMPLE}_pass_1.fastq.gz"
    U1="${TMP_DIR}/trim_single_${SAMPLE}_pass_1.fastq.gz"
    P2="${TMP_DIR}/trim_paired_${SAMPLE}_pass_2.fastq.gz"
    U2="${TMP_DIR}/trim_single_${SAMPLE}_pass_2.fastq.gz"
    SUMMARY="${TMP_DIR}/${SAMPLE}_trimmomatic_summary.txt"

    echo "Trimming Serum sample: ${SAMPLE}"
    trimmomatic PE \
        -threads "${threads}" \
        -phred33 \
        "${IN1}" "${IN2}" \
        "${P1}" "${U1}" \
        "${P2}" "${U2}" \
        "${TRIM_PARAMS[@]}" \
        2>&1 | tee "${SUMMARY}"

    validate_fastq_gz "${P1}"
    validate_fastq_gz "${U1}"
    validate_fastq_gz "${P2}"
    validate_fastq_gz "${U2}"
done

# -----------------------------
# Step 3: remove old broken files and replace them
# -----------------------------
echo
echo "===== Step 3: replace broken trimmed files ====="

rm -f \
    "${SERUM_OUT}/trim_paired_ERR1797969_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_paired_ERR1797969_pass_2.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797969_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797969_pass_2.fastq.gz" \
    "${SERUM_OUT}/trim_paired_ERR1797970_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_paired_ERR1797970_pass_2.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797970_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797970_pass_2.fastq.gz"

mv -f "${TMP_DIR}"/trim_* "${SERUM_OUT}/"
mv -f "${TMP_DIR}"/*_trimmomatic_summary.txt "${SERUM_OUT}/"

# Final check after moving
for f in \
    "${SERUM_OUT}/trim_paired_ERR1797969_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_paired_ERR1797969_pass_2.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797969_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797969_pass_2.fastq.gz" \
    "${SERUM_OUT}/trim_paired_ERR1797970_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_paired_ERR1797970_pass_2.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797970_pass_1.fastq.gz" \
    "${SERUM_OUT}/trim_single_ERR1797970_pass_2.fastq.gz"
do
    validate_fastq_gz "$f"
done

# -----------------------------
# Step 4: clean old RNA serum FastQC / MultiQC results
# -----------------------------
echo
echo "===== Step 4: clean old RNA serum FastQC / MultiQC results ====="

rm -f "${FASTQC_OUT}"/*_fastqc.html "${FASTQC_OUT}"/*_fastqc.zip
rm -rf "${MULTIQC_OUT:?}"/*

# -----------------------------
# Step 5: rerun FastQC for all RNA serum trimmed files
# -----------------------------
echo
echo "===== Step 5: rerun FastQC for all RNA serum trimmed files ====="

files=("${SERUM_OUT}"/*.fastq.gz)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in ${SERUM_OUT}"
    exit 1
fi

fastqc \
    -t "${threads}" \
    -o "${FASTQC_OUT}" \
    "${files[@]}"

# -----------------------------
# Step 6: rerun MultiQC for RNA serum
# -----------------------------
echo
echo "===== Step 6: rerun MultiQC for RNA serum ====="

multiqc "${FASTQC_OUT}" \
    --outdir "${MULTIQC_OUT}" \
    --filename rna_serum_multiqc_report.html \
    --force

echo
echo "Finished repair workflow at $(date)"
