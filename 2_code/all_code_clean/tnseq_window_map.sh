#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J tnseq_window_to_gene
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH -o tnseq_window_to_gene_%j.out
#SBATCH -e tnseq_window_to_gene_%j.err

set -euo pipefail

PROJECT_DIR="/home/yuch3531/genome_analysis_project"

COUNT_FILE="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/bowtie2_canu_25bp_window/tnseq_canu_bowtie2_25bp_window_counts.all_samples.tsv"

# IMPORTANT:
# This GFF must match the same Canu assembly used for Bowtie2 mapping.
# Change this path if your actual Canu Prokka GFF is elsewhere.
GFF="${PROJECT_DIR}/5_genome_annotation/Prokka/E745_canu/E745_canu.gff"

OUT_DIR="${PROJECT_DIR}/6_differential_gene_analysis/tnseq/bowtie2_canu_25bp_window/window_to_gene"

WINDOW_BED="${OUT_DIR}/tnseq_25bp_windows_with_counts.bed"
GENE_BED="${OUT_DIR}/canu_CDS_features.bed"
OVERLAP_TSV="${OUT_DIR}/tnseq_25bp_windows_overlapping_CDS.tsv"
GENE_COUNTS_LONG="${OUT_DIR}/tnseq_gene_counts_from_25bp_windows.long.tsv"
GENE_COUNTS_MATRIX="${OUT_DIR}/tnseq_gene_counts_from_25bp_windows.matrix.tsv"
SAMPLE_METADATA="${OUT_DIR}/tnseq_sample_metadata.tsv"

mkdir -p "$OUT_DIR"

echo "Count file: $COUNT_FILE"
echo "GFF file:   $GFF"
echo "Output dir: $OUT_DIR"

if [[ ! -s "$COUNT_FILE" ]]; then
    echo "ERROR: Count file does not exist or is empty: $COUNT_FILE" >&2
    exit 1
fi

if [[ ! -s "$GFF" ]]; then
    echo "ERROR: GFF file does not exist or is empty: $GFF" >&2
    echo "Please check the Canu Prokka annotation path." >&2
    exit 1
fi

# Load bedtools if needed.
# If module name differs on UPPMAX, run: module spider BEDTools
if ! command -v bedtools >/dev/null 2>&1; then
    module load BEDTools || module load BEDTools/2.31.1-GCC-13.3.0 || true
fi

if ! command -v bedtools >/dev/null 2>&1; then
    echo "ERROR: bedtools not found. Try: module spider BEDTools" >&2
    exit 1
fi

echo "Using bedtools:"
bedtools --version

echo "------------------------------------------------------------"
echo "1. Convert 25 bp window count table to BED-like format"
echo "------------------------------------------------------------"

# Input columns:
# condition sample contig window_start_0based window_end_0based window_name count
#
# Output BED columns:
# contig start end condition sample window_name count

awk -F'\t' 'BEGIN{OFS="\t"}
NR > 1 {
    print $3, $4, $5, $1, $2, $6, $7
}
' "$COUNT_FILE" > "$WINDOW_BED"

echo "Window BED written to:"
echo "$WINDOW_BED"

echo "------------------------------------------------------------"
echo "2. Convert GFF CDS features to BED"
echo "------------------------------------------------------------"

# GFF is 1-based inclusive.
# BED is 0-based half-open, so start = GFF_start - 1, end = GFF_end.
#
# For Prokka GFF, locus_tag is usually the best gene ID.
# If locus_tag is absent, fallback to ID or Name.

awk -F'\t' 'BEGIN{OFS="\t"}
function get_attr(attr, key,    n, a, i, kv) {
    n = split(attr, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) {
            return kv[2]
        }
    }
    return ""
}
!/^#/ && $3 == "CDS" {
    gene_id = get_attr($9, "locus_tag")
    if (gene_id == "") gene_id = get_attr($9, "ID")
    if (gene_id == "") gene_id = get_attr($9, "Name")
    if (gene_id == "") gene_id = "unknown_CDS_" NR

    product = get_attr($9, "product")
    if (product == "") product = "NA"

    print $1, $4 - 1, $5, gene_id, $3, $7, product
}
' "$GFF" > "$GENE_BED"

echo "CDS BED written to:"
echo "$GENE_BED"

echo "Number of CDS features:"
wc -l "$GENE_BED"

echo "------------------------------------------------------------"
echo "3. Intersect 25 bp windows with CDS features"
echo "------------------------------------------------------------"

# Output columns:
# window_contig window_start window_end condition sample window_name count
# gene_contig gene_start gene_end gene_id feature_type strand product

bedtools intersect \
    -a "$WINDOW_BED" \
    -b "$GENE_BED" \
    -wa -wb \
    > "$OVERLAP_TSV"

echo "Overlap table written to:"
echo "$OVERLAP_TSV"

echo "Number of window-CDS overlaps:"
wc -l "$OVERLAP_TSV"

echo "------------------------------------------------------------"
echo "4. Aggregate 25 bp window counts to gene-level counts"
echo "------------------------------------------------------------"

echo -e "condition\tsample\tgene_id\tfeature_type\tstrand\tproduct\tgene_count" > "$GENE_COUNTS_LONG"

awk -F'\t' 'BEGIN{OFS="\t"}
{
    condition=$4
    sample=$5
    count=$7

    gene_id=$11
    feature_type=$12
    strand=$13
    product=$14

    key=condition SUBSEP sample SUBSEP gene_id
    gene_count[key] += count

    meta[key] = condition "\t" sample "\t" gene_id "\t" feature_type "\t" strand "\t" product
}
END {
    for (key in gene_count) {
        print meta[key], gene_count[key]
    }
}
' "$OVERLAP_TSV" | sort -k1,1 -k2,2 -k3,3 >> "$GENE_COUNTS_LONG"

echo "Long gene count table written to:"
echo "$GENE_COUNTS_LONG"

echo "------------------------------------------------------------"
echo "5. Create gene x sample count matrix and sample metadata"
echo "------------------------------------------------------------"

python3 - <<PY
import csv
from collections import defaultdict

long_file = "${GENE_COUNTS_LONG}"
matrix_file = "${GENE_COUNTS_MATRIX}"
metadata_file = "${SAMPLE_METADATA}"

counts = defaultdict(dict)
products = {}
samples = []
sample_conditions = {}

with open(long_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        condition = row["condition"]
        sample = row["sample"]
        gene_id = row["gene_id"]
        product = row["product"]
        count = int(float(row["gene_count"]))

        if sample not in samples:
            samples.append(sample)
        sample_conditions[sample] = condition

        counts[gene_id][sample] = count
        products[gene_id] = product

samples = sorted(samples)

with open(matrix_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["gene_id", "product"] + samples)

    for gene_id in sorted(counts):
        row = [gene_id, products.get(gene_id, "NA")]
        for sample in samples:
            row.append(counts[gene_id].get(sample, 0))
        writer.writerow(row)

with open(metadata_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["sample", "condition"])
    for sample in samples:
        writer.writerow([sample, sample_conditions.get(sample, "NA")])
PY

echo "Gene count matrix written to:"
echo "$GENE_COUNTS_MATRIX"

echo "Sample metadata written to:"
echo "$SAMPLE_METADATA"

echo "------------------------------------------------------------"
echo "Preview: gene count matrix"
echo "------------------------------------------------------------"
head -10 "$GENE_COUNTS_MATRIX" | column -t -s $'\t'

echo "------------------------------------------------------------"
echo "Preview: sample metadata"
echo "------------------------------------------------------------"
cat "$SAMPLE_METADATA" | column -t -s $'\t'

echo "Done."
