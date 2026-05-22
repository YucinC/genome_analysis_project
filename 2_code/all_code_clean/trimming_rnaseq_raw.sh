#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 03:00:00
#SBATCH -J trim_rnaseq
#SBATCH -o /home/yuch3531/genome_analysis_project/3_preprocessing/trimming/trim_rnaseq_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/3_preprocessing/trimming/trim_rnaseq_%j.err

set -euo pipefail

module load Trimmomatic

# -----------------------------
# Input directories
# -----------------------------
BH_RAW="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_BH/raw"
SERUM_RAW="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_Serum/raw"

# -----------------------------
# Output directories
# -----------------------------
OUT_BASE="/home/yuch3531/genome_analysis_project/3_preprocessing/trimming"
BH_OUT="${OUT_BASE}/rna_bh"
SERUM_OUT="${OUT_BASE}/rna_serum"

mkdir -p "${OUT_BASE}" "${BH_OUT}" "${SERUM_OUT}"

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
echo "Starting RNA-seq trimming at $(date)"

# -----------------------------
# Common trimming parameters
# -----------------------------
TRIM_PARAMS="ILLUMINACLIP:${ADAPTERS}:2:30:10 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"

# -----------------------------
# Trim RNA-Seq_BH raw samples
# -----------------------------
for SAMPLE in ERR1797972 ERR1797973 ERR1797974
do
    IN1="${BH_RAW}/${SAMPLE}_1.fastq.gz"
    IN2="${BH_RAW}/${SAMPLE}_2.fastq.gz"

    P1="${BH_OUT}/trim_paired_${SAMPLE}_pass_1.fastq.gz"
    U1="${BH_OUT}/trim_single_${SAMPLE}_pass_1.fastq.gz"
    P2="${BH_OUT}/trim_paired_${SAMPLE}_pass_2.fastq.gz"
    U2="${BH_OUT}/trim_single_${SAMPLE}_pass_2.fastq.gz"

    echo "Trimming BH sample: ${SAMPLE}"
    trimmomatic PE \
        -threads 2 \
        -phred33 \
        "${IN1}" "${IN2}" \
        "${P1}" "${U1}" \
        "${P2}" "${U2}" \
        ${TRIM_PARAMS}
done

# -----------------------------
# Trim RNA-Seq_Serum raw samples
# -----------------------------
for SAMPLE in ERR1797969 ERR1797970 ERR1797971
do
    IN1="${SERUM_RAW}/${SAMPLE}_1.fastq.gz"
    IN2="${SERUM_RAW}/${SAMPLE}_2.fastq.gz"

    P1="${SERUM_OUT}/trim_paired_${SAMPLE}_pass_1.fastq.gz"
    U1="${SERUM_OUT}/trim_single_${SAMPLE}_pass_1.fastq.gz"
    P2="${SERUM_OUT}/trim_paired_${SAMPLE}_pass_2.fastq.gz"
    U2="${SERUM_OUT}/trim_single_${SAMPLE}_pass_2.fastq.gz"

    echo "Trimming Serum sample: ${SAMPLE}"
    trimmomatic PE \
        -threads 2 \
        -phred33 \
        "${IN1}" "${IN2}" \
        "${P1}" "${U1}" \
        "${P2}" "${U2}" \
        ${TRIM_PARAMS}
done

echo "Finished RNA-seq trimming at $(date)"
