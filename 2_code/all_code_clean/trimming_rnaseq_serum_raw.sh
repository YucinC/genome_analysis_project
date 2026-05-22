#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -J re_fastqc_multiqc
#SBATCH --mail-type=ALL
#SBATCH --output=/home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc_multiqc_%j.out
#SBATCH --error=/home/yuch3531/genome_analysis_project/3_preprocessing/re_fastqc_multiqc_%j.err

set -euo pipefail

module load FastQC
module load MultiQC

threads=$SLURM_CPUS_PER_TASK

base_dir="/home/yuch3531/genome_analysis_project/3_preprocessing"
input_dir="$base_dir/trimming"
output_dir="$base_dir/re_fastqc"

mkdir -p "$output_dir/illumina"
mkdir -p "$output_dir/rna_bh"
mkdir -p "$output_dir/rna_serum"
mkdir -p "$output_dir/multiqc"

shopt -s nullglob

echo "===== Start FastQC ====="
date

for data_type in illumina rna_bh rna_serum
do
    echo "Processing ${data_type} ..."
    
    files=("$input_dir/$data_type"/*.fastq.gz)
    
    if [ ${#files[@]} -eq 0 ]; then
        echo "No fastq.gz files found in $input_dir/$data_type"
        continue
    fi
    
    fastqc \
        -t "$threads" \
        -o "$output_dir/$data_type" \
        "${files[@]}"
done

echo "===== FastQC finished ====="
date

echo "===== Start MultiQC ====="
multiqc "$output_dir" \
    -o "$output_dir/multiqc" \
    -n multiqc_report.html \
    --force

echo "===== MultiQC finished ====="
date

echo "All done."
