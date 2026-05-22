#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J busco_check
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/4_genome_assembly/assembly_check/busco_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/4_genome_assembly/assembly_check/busco_%j.err

set -euo pipefail

module load BUSCO/5.8.2

BASE_DIR="/home/yuch3531/genome_analysis_project/4_genome_assembly"
OUT_BASE="${BASE_DIR}/assembly_check"

SPADES_FASTA="${BASE_DIR}/spades_4947505/contigs.fasta"
CANU_FASTA="${BASE_DIR}/canu_4935797/E745_canu.contigs.fasta"

mkdir -p "${OUT_BASE}"

busco \
  -i "${SPADES_FASTA}" \
  -m genome \
  -l bacteria_odb12 \
  -c 2 \
  -o spades_busco \
  --out_path "${OUT_BASE}" \
  -f

busco \
  -i "${CANU_FASTA}" \
  -m genome \
  -l bacteria_odb12 \
  -c 2 \
  -o canu_busco \
  --out_path "${OUT_BASE}" \
  -f

