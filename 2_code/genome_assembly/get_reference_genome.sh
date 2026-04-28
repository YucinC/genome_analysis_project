#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J get_ref
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/1_data/reference/get_ref_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/1_data/reference/get_ref_%j.err

set -euo pipefail

REF_DIR="/home/yuch3531/genome_analysis_project/1_data/reference"
REF_FASTA="${REF_DIR}/E745_reference_CP014529_CP014535.fasta"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

module load EDirect

efetch -db nuccore \
  -id CP014529,CP014530,CP014531,CP014532,CP014533,CP014534,CP014535 \
  -format fasta > "${REF_FASTA}"

echo "Reference saved to: ${REF_FASTA}"
