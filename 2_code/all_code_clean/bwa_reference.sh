#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J bwa_reference_mapping
#SBATCH -t 02:00:00
#SBATCH -n 2
#SBATCH --mail-type=ALL
#SBATCH -o bwa_reference_mapping_%j.out
#SBATCH -e bwa_reference_mapping_%j.err

set -euo pipefail

############################################
# 1. Load modules
############################################

module purge
module load bwa || module load BWA
module load samtools || module load SAMtools

THREADS="${SLURM_CPUS_PER_TASK:-8}"

############################################
# 2. Project paths
############################################

PROJECT_DIR="/home/yuch3531/genome_analysis_project"
DATA_BASE="${PROJECT_DIR}/1_data/transcriptomics_data"

REFERENCE_FASTA="${PROJECT_DIR}/1_data/reference/E745_reference.fasta"
REFERENCE_NAME="reference_genome"

RNA_OUT_BASE="${PROJECT_DIR}/6_differential_gene_analysis/rna_mapping/${REFERENCE_NAME}"
TN_OUT_BASE="${PROJECT_DIR}/6_differential_gene_analysis/tnseq_mapping/${REFERENCE_NAME}"

mkdir -p "${RNA_OUT_BASE}" "${TN_OUT_BASE}"

############################################
# 3. Dataset groups
############################################

RNA_DATASETS=(
    "RNA-Seq_BH"
    "RNA-Seq_Serum"
)

TN_DATASETS=(
    "Tn-Seq_BHI"
    "Tn-Seq_HSerum"
    "Tn-Seq_Serum"
)

############################################
# 4. Helper functions
############################################

index_reference_if_needed() {
    local ref="$1"

    if [[ ! -s "${ref}" ]]; then
        echo "ERROR: reference genome not found or empty: ${ref}" >&2
        exit 1
    fi

    if [[ ! -s "${ref}.bwt" || ! -s "${ref}.sa" ]]; then
        echo "BWA index not found for ${ref}"
        echo "Running bwa index..."
        bwa index "${ref}"
    else
        echo "BWA index already exists for ${ref}"
    fi
}

find_mate_r2() {
    local r1="$1"
    local candidate

    local candidates=(
        "${r1/_R1_/_R2_}"
        "${r1/_R1./_R2.}"
        "${r1/_R1/_R2}"
        "${r1/%_1.fastq.gz/_2.fastq.gz}"
        "${r1/%_1.fq.gz/_2.fq.gz}"
        "${r1/%_1.fastq/_2.fastq}"
        "${r1/%_1.fq/_2.fq}"
    )

    for candidate in "${candidates[@]}"; do
        if [[ "${candidate}" != "${r1}" && -s "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

sample_name_from_r1() {
    local file="$1"
    local sample

    sample="$(basename "${file}")"

    sample="${sample%.fastq.gz}"
    sample="${sample%.fq.gz}"
    sample="${sample%.fastq}"
    sample="${sample%.fq}"

    sample="$(echo "${sample}" \
        | sed -E 's/^trim_paired_//' \
        | sed -E 's/^paired_//' \
        | sed -E 's/_R?1(_[0-9]+)?$//' \
        | sed -E 's/_1$//' \
        | sed -E 's/[^A-Za-z0-9_.-]+/_/g')"

    echo "${sample}"
}

run_mapping_for_group() {
    local data_type="$1"
    local out_base="$2"
    shift 2

    local datasets=("$@")

    for dataset in "${datasets[@]}"; do
        local trimmed_dir="${DATA_BASE}/${dataset}/trimmed"

        if [[ ! -d "${trimmed_dir}" ]]; then
            echo "WARNING: trimmed directory not found, skipping: ${trimmed_dir}" >&2
            continue
        fi

        echo "============================================================"
        echo "Processing ${data_type}: ${dataset}"
        echo "Input directory: ${trimmed_dir}"
        echo "Reference: ${REFERENCE_NAME}"
        echo "============================================================"

        mapfile -t r1_files < <(
            find "${trimmed_dir}" -maxdepth 1 -type f \
                \( -name "*_1.fastq.gz" -o -name "*_1.fq.gz" -o -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" \) \
                ! -name "*single*" \
                ! -name "*unpaired*" \
                | sort
        )

        if [[ "${#r1_files[@]}" -eq 0 ]]; then
            echo "WARNING: no paired-end R1 files found in ${trimmed_dir}" >&2
            continue
        fi

        for r1 in "${r1_files[@]}"; do
            local r2
            local sample
            local sample_out_dir
            local bam
            local flagstat
            local idxstats
            local log_file
            local rg

            if ! r2="$(find_mate_r2 "${r1}")"; then
                echo "WARNING: could not find R2 mate for ${r1}, skipping." >&2
                continue
            fi

            sample="$(sample_name_from_r1 "${r1}")"
            sample_out_dir="${out_base}/${dataset}/${sample}"

            mkdir -p "${sample_out_dir}"

            bam="${sample_out_dir}/${sample}.${REFERENCE_NAME}.sorted.bam"
            flagstat="${sample_out_dir}/${sample}.${REFERENCE_NAME}.flagstat.txt"
            idxstats="${sample_out_dir}/${sample}.${REFERENCE_NAME}.idxstats.txt"
            log_file="${sample_out_dir}/${sample}.${REFERENCE_NAME}.bwa.log"

            rg="@RG\tID:${sample}_${REFERENCE_NAME}\tSM:${sample}\tPL:ILLUMINA\tLB:${data_type}_${dataset}"

            echo "Mapping sample: ${sample}"
            echo "R1: ${r1}"
            echo "R2: ${r2}"
            echo "Output: ${bam}"

            if [[ -s "${bam}" && -s "${bam}.bai" ]]; then
                echo "Output already exists, skipping: ${bam}"
                continue
            fi

            bwa mem \
                -M \
                -t "${THREADS}" \
                -R "${rg}" \
                "${REFERENCE_FASTA}" \
                "${r1}" \
                "${r2}" \
                2> "${log_file}" \
                | samtools sort \
                    -@ "${THREADS}" \
                    -o "${bam}" \
                    -

            samtools index "${bam}"
            samtools flagstat "${bam}" > "${flagstat}"
            samtools idxstats "${bam}" > "${idxstats}"
        done
    done
}

############################################
# 5. Main workflow
############################################

echo "Starting BWA mapping to downloaded reference genome"
echo "Reference fasta: ${REFERENCE_FASTA}"
echo "Threads: ${THREADS}"

index_reference_if_needed "${REFERENCE_FASTA}"

run_mapping_for_group "RNA-seq" "${RNA_OUT_BASE}" "${RNA_DATASETS[@]}"
run_mapping_for_group "Tn-seq" "${TN_OUT_BASE}" "${TN_DATASETS[@]}"

echo "BWA mapping to downloaded reference genome finished."
