#!/bin/bash -l
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -J make_window_to_gene
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH -o make_window_to_gene_%j.out
#SBATCH -e make_window_to_gene_%j.err

set -euo pipefail

module purge
module load BEDTools || module load bedtools

# -----------------------------
# 1. Set paths
# -----------------------------

WINDOW_DIR="/home/yuch3531/genome_analysis_project/6_differential_gene_analysis/tnseq/bowtie2_canu_25bp_window"

# Change this only if your actual Canu Prokka GFF path is different.
GFF="/home/yuch3531/genome_analysis_project/5_genome_annotation/Prokka/E745_canu/E745_canu.gff"

OUTDIR="/home/yuch3531/genome_analysis_project/6_differential_gene_analysis/tnseq/window_to_gene"
mkdir -p "${OUTDIR}"

WINDOW_BED="${OUTDIR}/windows.bed"
GENE_BED="${OUTDIR}/genes.bed"
RAW_INTERSECT="${OUTDIR}/window_gene_overlap_raw.tsv"
WINDOW_TO_GENE="${OUTDIR}/window_to_gene.tsv"

# -----------------------------
# 2. Check input files
# -----------------------------

if [ ! -d "${WINDOW_DIR}" ]; then
  echo "ERROR: WINDOW_DIR does not exist:"
  echo "${WINDOW_DIR}"
  exit 1
fi

if [ ! -f "${GFF}" ]; then
  echo "ERROR: GFF file does not exist:"
  echo "${GFF}"
  exit 1
fi

WINDOW_FILE=$(find "${WINDOW_DIR}" -type f -name "*.25bp_window.counts.tsv" | head -n 1)

if [ -z "${WINDOW_FILE}" ]; then
  echo "ERROR: No .25bp_window.counts.tsv file found in:"
  echo "${WINDOW_DIR}"
  exit 1
fi

echo "Using window file to define genomic windows:"
echo "${WINDOW_FILE}"

echo "Using GFF annotation:"
echo "${GFF}"

# -----------------------------
# 3. Convert window count file to BED
# -----------------------------

# Expected columns in window count file:
# contig  window_start_0based  window_end_0based  window_name  count
#
# BED format:
# contig  start  end  window_id

awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {
  for (i=1; i<=NF; i++) header[$i]=i
}
NR>1 {
  print $header["contig"], $header["window_start_0based"], $header["window_end_0based"], $header["window_name"]
}' "${WINDOW_FILE}" > "${WINDOW_BED}"

echo "Window BED created:"
echo "${WINDOW_BED}"
echo "Number of windows:"
wc -l "${WINDOW_BED}"

# -----------------------------
# 4. Convert Prokka GFF CDS features to BED
# -----------------------------

# GFF is 1-based inclusive.
# BED is 0-based half-open.
# Therefore:
# BED_start = GFF_start - 1
# BED_end   = GFF_end
#
# Output BED columns:
# contig  start  end  locus_tag  gene_name  product_name  strand

awk -F'\t' 'BEGIN{OFS="\t"}
  $0 !~ /^#/ && $3 == "CDS" {
    contig=$1
    start=$4-1
    end=$5
    strand=$7
    attr=$9

    locus_tag="NA"
    gene="NA"
    product="NA"

    n=split(attr, fields, ";")
    for (i=1; i<=n; i++) {
      split(fields[i], kv, "=")
      if (kv[1] == "locus_tag") locus_tag=kv[2]
      if (kv[1] == "gene") gene=kv[2]
      if (kv[1] == "product") product=kv[2]
    }

    gsub(/%20/, " ", product)
    gsub(/%2C/, ",", product)
    gsub(/%3B/, ";", product)

    print contig, start, end, locus_tag, gene, product, strand
  }
' "${GFF}" > "${GENE_BED}"

echo "Gene BED created:"
echo "${GENE_BED}"
echo "Number of CDS features:"
wc -l "${GENE_BED}"

# -----------------------------
# 5. Intersect windows with genes
# -----------------------------

# bedtools intersect -wo appends overlap length as the last column.
#
# A columns:
# 1 window_contig
# 2 window_start
# 3 window_end
# 4 window_id
#
# B columns:
# 5 gene_contig
# 6 gene_start
# 7 gene_end
# 8 locus_tag
# 9 gene_name
# 10 product_name
# 11 strand
#
# 12 overlap_bp

bedtools intersect \
  -a "${WINDOW_BED}" \
  -b "${GENE_BED}" \
  -wo > "${RAW_INTERSECT}"

echo "Raw window-gene overlaps:"
echo "${RAW_INTERSECT}"
echo "Number of overlapping records:"
wc -l "${RAW_INTERSECT}"

# -----------------------------
# 6. Keep the gene with the largest overlap for each window
# -----------------------------

awk -F'\t' 'BEGIN{OFS="\t"}
{
  window_id=$4
  overlap=$12

  if (!(window_id in best_overlap) || overlap > best_overlap[window_id]) {
    best_overlap[window_id]=overlap
    best_line[window_id]=$0
  }
}
END {
  print "window_id","contig","window_start_0based","window_end_0based","locus_tag","gene_name","product_name","strand","overlap_bp"

  for (w in best_line) {
    split(best_line[w], x, "\t")

    window_contig=x[1]
    window_start=x[2]
    window_end=x[3]
    window_id=x[4]

    locus_tag=x[8]
    gene_name=x[9]
    product_name=x[10]
    strand=x[11]
    overlap_bp=x[12]

    print window_id, window_contig, window_start, window_end, locus_tag, gene_name, product_name, strand, overlap_bp
  }
}' "${RAW_INTERSECT}" | sort -k2,2 -k3,3n > "${WINDOW_TO_GENE}"

echo "Final window-to-gene annotation table:"
echo "${WINDOW_TO_GENE}"
echo "Number of annotated windows:"
tail -n +2 "${WINDOW_TO_GENE}" | wc -l

echo "Done."
