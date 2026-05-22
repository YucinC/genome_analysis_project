#!/bin/bash

set -euo pipefail

PROJECT_DIR="/home/yuch3531/genome_analysis_project"

TN_COUNT_DIR="${PROJECT_DIR}/6_differential_gene_analysis/tnseq_counting/canu"
OUTDIR="${PROJECT_DIR}/6_differential_gene_analysis/tnseq_counting_summary/canu"

mkdir -p "$OUTDIR"

SUMMARY="${OUTDIR}/tnseq_htseq_count_summary.tsv"

echo -e "dataset\tsample\ttotal_reads\tassigned_to_genes\tnonzero_genes\tno_feature\tambiguous\ttoo_low_aQual\tnot_aligned\talignment_not_unique\tassigned_rate_percent\tnot_aligned_rate_percent\tcount_file" > "$SUMMARY"

find "$TN_COUNT_DIR" -name "*.counts.txt" | sort | while read -r COUNT_FILE; do
    DATASET="$(basename "$(dirname "$COUNT_FILE")")"

    SAMPLE="$(basename "$COUNT_FILE")"
    SAMPLE="${SAMPLE%.canu.htseq.counts.txt}"
    SAMPLE="${SAMPLE%.htseq.counts.txt}"
    SAMPLE="${SAMPLE%.counts.txt}"

    TOTAL="$(awk '{sum += $2} END {print sum+0}' "$COUNT_FILE")"

    ASSIGNED="$(awk '$1 !~ /^__/ {sum += $2} END {print sum+0}' "$COUNT_FILE")"

    NONZERO_GENES="$(awk '$1 !~ /^__/ && $2 > 0 {n++} END {print n+0}' "$COUNT_FILE")"

    NO_FEATURE="$(awk '$1=="__no_feature" {print $2+0}' "$COUNT_FILE")"
    AMBIGUOUS="$(awk '$1=="__ambiguous" {print $2+0}' "$COUNT_FILE")"
    TOO_LOW_AQUAL="$(awk '$1=="__too_low_aQual" {print $2+0}' "$COUNT_FILE")"
    NOT_ALIGNED="$(awk '$1=="__not_aligned" {print $2+0}' "$COUNT_FILE")"
    ALIGNMENT_NOT_UNIQUE="$(awk '$1=="__alignment_not_unique" {print $2+0}' "$COUNT_FILE")"

    ASSIGNED_RATE="$(awk -v a="$ASSIGNED" -v t="$TOTAL" 'BEGIN {if (t>0) printf "%.2f", a/t*100; else print "NA"}')"
    NOT_ALIGNED_RATE="$(awk -v n="$NOT_ALIGNED" -v t="$TOTAL" 'BEGIN {if (t>0) printf "%.2f", n/t*100; else print "NA"}')"

    echo -e "${DATASET}\t${SAMPLE}\t${TOTAL}\t${ASSIGNED}\t${NONZERO_GENES}\t${NO_FEATURE}\t${AMBIGUOUS}\t${TOO_LOW_AQUAL}\t${NOT_ALIGNED}\t${ALIGNMENT_NOT_UNIQUE}\t${ASSIGNED_RATE}\t${NOT_ALIGNED_RATE}\t${COUNT_FILE}" >> "$SUMMARY"
done

echo "Summary written to:"
echo "$SUMMARY"

echo
echo "Preview:"
column -t -s $'\t' "$SUMMARY"
