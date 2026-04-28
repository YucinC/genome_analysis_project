#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -J re_multiqc_by_type
#SBATCH --mail-type=ALL
#SBATCH --output=/home/yuch3531/genome_analysis_project/3_preprocessing/re_multiqc_by_type_%j.out
#SBATCH --error=/home/yuch3531/genome_analysis_project/3_preprocessing/re_multiqc_by_type_%j.err

set -euo pipefail

module load MultiQC

input_base="/home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc"
output_base="/home/yuch3531/genome_analysis_project/3_preprocessing/re_multiqc"

mkdir -p "$output_base"

echo "======================================"
echo "Start MultiQC by type"
date
echo "Input base: $input_base"
echo "Output base: $output_base"
echo "======================================"

for type in illumina rna_bh rna_serum
do
    input_dir="$input_base/$type"
    output_dir="$output_base/$type"

    echo "--------------------------------------"
    echo "Processing type: $type"
    echo "Input dir: $input_dir"
    echo "Output dir: $output_dir"
    date

    if [ ! -d "$input_dir" ]; then
        echo "Warning: $input_dir does not exist, skipping."
        continue
    fi

    mkdir -p "$output_dir"

    multiqc "$input_dir" \
        --outdir "$output_dir" \
        --filename "${type}_multiqc_report.html" \
        --force

    echo "Finished type: $type"
    date
done

echo "======================================"
echo "All MultiQC jobs finished"
date
echo "Done."
