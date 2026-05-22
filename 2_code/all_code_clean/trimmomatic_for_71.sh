#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -J trim_ERR1797971
#SBATCH --mail-type=ALL
#SBATCH --output=/home/yuch3531/genome_analysis_project/3_preprocessing/trim_ERR1797971_%j.out
#SBATCH --error=/home/yuch3531/genome_analysis_project/3_preprocessing/trim_ERR1797971_%j.err

set -euo pipefail

module load Trimmomatic

threads=$SLURM_CPUS_PER_TASK

in_dir="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data/RNA-Seq_Serum/raw"
out_dir="/home/yuch3531/genome_analysis_project/3_preprocessing/trimming/rna_serum"
tmp_dir="$SNIC_TMP/trim_ERR1797971"

sample="ERR1797971"
r1="$in_dir/${sample}_1.fastq.gz"
r2="$in_dir/${sample}_2.fastq.gz"

adapter="${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa"

mkdir -p "$out_dir"
mkdir -p "$tmp_dir"

echo "======================================"
echo "Start trimming sample: $sample"
date
echo "Input R1: $r1"
echo "Input R2: $r2"
echo "Output dir: $out_dir"
echo "Temp dir: $tmp_dir"
echo "Threads: $threads"
echo "Adapter file: $adapter"
echo "======================================"

if [ ! -f "$r1" ]; then
    echo "ERROR: File not found: $r1"
    exit 1
fi

if [ ! -f "$r2" ]; then
    echo "ERROR: File not found: $r2"
    exit 1
fi

if [ ! -f "$adapter" ]; then
    echo "ERROR: Adapter file not found: $adapter"
    echo "EBROOTTRIMMOMATIC=$EBROOTTRIMMOMATIC"
    exit 1
fi

trimmomatic PE \
  -threads "$threads" \
  -summary "$tmp_dir/${sample}_trimmomatic_summary.txt" \
  "$r1" "$r2" \
  "$tmp_dir/trim_paired_${sample}_1.fastq.gz" \
  "$tmp_dir/trim_single_${sample}_1.fastq.gz" \
  "$tmp_dir/trim_paired_${sample}_2.fastq.gz" \
  "$tmp_dir/trim_single_${sample}_2.fastq.gz" \
  ILLUMINACLIP:"$adapter":2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36

echo "Trimmomatic finished."
date

echo "Checking gzip integrity..."
gzip -t "$tmp_dir"/trim_*_"${sample}"_*.fastq.gz

echo "Checking FASTQ structure..."
for f in "$tmp_dir"/trim_*_"${sample}"_*.fastq.gz
do
    echo "Checking $f"
    zcat "$f" | awk 'END {print "lines=" NR, " mod4=" NR%4; if (NR%4!=0) exit 1}'
done

mv "$tmp_dir"/trim_*_"${sample}"_*.fastq.gz "$out_dir"/
mv "$tmp_dir/${sample}_trimmomatic_summary.txt" "$out_dir"/

echo "All files moved to $out_dir"
date
echo "Done."
