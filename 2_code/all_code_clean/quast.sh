#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J quast_eval
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/4_genome_assembly/quast_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/4_genome_assembly/quast_%j.err

set -euo pipefail

module load QUAST

ASSEMBLY_BASE="/home/yuch3531/genome_analysis_project/4_genome_assembly"
QUAST_OUT="${ASSEMBLY_BASE}/quast_${SLURM_JOB_ID}"

# ===== assembly files =====
SPADES_CONTIG="${ASSEMBLY_BASE}/spades_4947505/contigs.fasta"
CANU_CONTIG="${ASSEMBLY_BASE}/canu_4935797/E745_canu.contigs.fasta"

# ===== optional reference =====
REF="/home/yuch3531/genome_analysis_project/1_data/reference/E745_reference.fasta"

mkdir -p "${QUAST_OUT}"

echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Output dir: ${QUAST_OUT}"

echo "SPAdes assembly: ${SPADES_CONTIG}"
echo "Canu assembly:   ${CANU_CONTIG}"

if [[ ! -f "${SPADES_CONTIG}" ]]; then
    echo "ERROR: SPAdes contigs file not found: ${SPADES_CONTIG}"
    exit 1
fi

if [[ ! -f "${CANU_CONTIG}" ]]; then
    echo "ERROR: Canu contigs file not found: ${CANU_CONTIG}"
    exit 1
fi

if [[ -n "${REF}" ]]; then
    if [[ ! -f "${REF}" ]]; then
        echo "ERROR: Reference file not found: ${REF}"
        exit 1
    fi

    echo "Reference genome: ${REF}"
    echo "Running QUAST with reference..."
    quast.py \
        "${SPADES_CONTIG}" \
        "${CANU_CONTIG}" \
        -o "${QUAST_OUT}" \
        -r "${REF}" \
        -l "SPAdes,Canu" \
        --min-contig 500 \
        --threads 2
else
    echo "No reference provided. Running QUAST in reference-free mode..."
    quast.py \
        "${SPADES_CONTIG}" \
        "${CANU_CONTIG}" \
        -o "${QUAST_OUT}" \
        -l "SPAdes,Canu" \
        --min-contig 500 \
        --threads 2
fi

echo "Job finished at: $(date)"
