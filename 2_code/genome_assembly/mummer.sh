#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 2
#SBATCH -t 00:10:00
#SBATCH -J mummer_pretty
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/4_genome_assembly/mummer/mummer_pretty_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/4_genome_assembly/mummer/mummer_pretty_%j.err

set -euo pipefail

module load MUMmer

REF_FASTA="/home/yuch3531/genome_analysis_project/1_data/reference/E745_reference.fasta"
SPADES_FASTA="/home/yuch3531/genome_analysis_project/4_genome_assembly/spades_4947505/spades_contigs.fasta"
CANU_FASTA="/home/yuch3531/genome_analysis_project/4_genome_assembly/canu_4935797/E745_canu.contigs.fasta"

OUTDIR="/home/yuch3531/genome_analysis_project/4_genome_assembly/mummer_2"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "Starting MUMmer comparison..."

# ----------------------------------------
# 1. Create simplified FASTA headers
# ----------------------------------------

SPADES_SHORT="${OUTDIR}/spades_contigs.shortid.fasta"
CANU_SHORT="${OUTDIR}/canu_contigs.shortid.fasta"

awk '
/^>/ {
    i++
    print ">spades_ctg_" i
    next
}
{ print }
' "${SPADES_FASTA}" > "${SPADES_SHORT}"

awk '
/^>/ {
    i++
    print ">canu_ctg_" i
    next
}
{ print }
' "${CANU_FASTA}" > "${CANU_SHORT}"

# ----------------------------------------
# 2. Function to run nucmer + filter + plot
# ----------------------------------------

run_plot () {
    local QUERY_FASTA="$1"
    local PREFIX="$2"
    local XLAB="$3"
    local YLAB="$4"

    echo "Running ${PREFIX}..."

    nucmer --prefix="${PREFIX}" "${REF_FASTA}" "${QUERY_FASTA}"

    delta-filter -q -r "${PREFIX}.delta" > "${PREFIX}.filtered.delta"

    show-coords -rcl "${PREFIX}.filtered.delta" > "${PREFIX}.filtered.coords"

    mummerplot \
        "${PREFIX}.filtered.delta" \
        -R "${REF_FASTA}" \
        -Q "${QUERY_FASTA}" \
        --png \
        --layout \
        --prefix="${PREFIX}.filtered"

    # ----------------------------------------
    # 3. Create a cleaner custom gnuplot script
    # ----------------------------------------

    GP_ORIG="${PREFIX}.filtered.gp"
    GP_CUSTOM="${PREFIX}.filtered.custom.gp"
    PNG_CUSTOM="${PREFIX}.filtered.custom.png"

    grep -vE '^(set terminal|set output|set xlabel|set ylabel|set xtics|set ytics|set key)' "${GP_ORIG}" > "${PREFIX}.filtered.body.gp"

    cat > "${GP_CUSTOM}" <<EOF
set terminal pngcairo size 2600,2600 enhanced font "Arial,18"
set output "${PNG_CUSTOM}"

set xlabel "${XLAB}" font ",24"
set ylabel "${YLAB}" font ",24"

set xtics rotate by 45 right font ",12"
set ytics font ",12"

set lmargin 35
set bmargin 12
set rmargin 5
set tmargin 5

unset key
EOF

    cat "${PREFIX}.filtered.body.gp" >> "${GP_CUSTOM}"

    gnuplot "${GP_CUSTOM}"

    echo "Finished ${PREFIX}"
    echo "Main plot: ${OUTDIR}/${PNG_CUSTOM}"
}

# ----------------------------------------
# 4. Run two comparisons
# ----------------------------------------

run_plot "${SPADES_SHORT}" "ref_vs_spades" "E745 reference genome" "SPAdes assembly"
run_plot "${CANU_SHORT}"   "ref_vs_canu"   "E745 reference genome" "Canu assembly"

echo "All MUMmer analyses finished."
