#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J replot_mummer
#SBATCH --mail-type=ALL
#SBATCH -o /gorilla/home/yuch3531/genome_analysis_project/4_genome_assembly/replot_mummer/replot_mummer_%j.out
#SBATCH -e /gorilla/home/yuch3531/genome_analysis_project/4_genome_assembly/replot_mummer/replot_mummer_%j.err

set -euo pipefail

module load gnuplot

INDIR="/gorilla/home/yuch3531/genome_analysis_project/4_genome_assembly/mummer"
OUTDIR="/gorilla/home/yuch3531/genome_analysis_project/4_genome_assembly/replot_mummer"

mkdir -p "${OUTDIR}"

cd "${INDIR}"

replot_gp() {
    local GP_ORIG="$1"
    local PREFIX="$2"
    local XLAB="$3"
    local YLAB="$4"

    local IMG_WIDTH=2600
    local IMG_HEIGHT=2600
    local BASE_FONT="Arial,30"

    local XLAB_FONT=32
    local YLAB_FONT=32
    local XTIC_FONT=24
    local YTIC_FONT=24
    local XTIC_ROTATE=45

    local LMARGIN=15
    local BMARGIN=15
    local RMARGIN=5
    local TMARGIN=5

    local GP_BODY="${OUTDIR}/${PREFIX}.replot.body.gp"
    local GP_CUSTOM="${OUTDIR}/${PREFIX}.replot.gp"
    local PNG_CUSTOM="${OUTDIR}/${PREFIX}.replot.png"

    echo "Replotting ${GP_ORIG}"
    echo "Output: ${PNG_CUSTOM}"

    awk '
    BEGIN { skip=0 }

    /^set xtics/ {
        if ($0 !~ /\)$/) skip=1
        next
    }

    /^set ytics/ {
        if ($0 !~ /\)$/) skip=1
        next
    }

    skip == 1 {
        if ($0 ~ /\)$/) skip=0
        next
    }

    /^set terminal/ { next }
    /^set output/ { next }
    /^set xlabel/ { next }
    /^set ylabel/ { next }
    /^set key/ { next }
    /^set size/ { next }

    { print }
    ' "${GP_ORIG}" > "${GP_BODY}"

    cat > "${GP_CUSTOM}" <<EOF
set terminal pngcairo size ${IMG_WIDTH},${IMG_HEIGHT} enhanced font "${BASE_FONT}"
set output "${PNG_CUSTOM}"

set xlabel "${XLAB}" font ",${XLAB_FONT}"
set ylabel "${YLAB}" font ",${YLAB_FONT}"

set xtics rotate by ${XTIC_ROTATE} right font ",${XTIC_FONT}"
set ytics font ",${YTIC_FONT}"

set lmargin ${LMARGIN}
set bmargin ${BMARGIN}
set rmargin ${RMARGIN}
set tmargin ${TMARGIN}

unset key
EOF

    cat "${GP_BODY}" >> "${GP_CUSTOM}"

    gnuplot "${GP_CUSTOM}"

    echo "Finished: ${PNG_CUSTOM}"
    echo
}

replot_gp "ref_vs_spades.filtered.gp" "ref_vs_spades.filtered" "E745 reference genome" "SPAdes assembly"
replot_gp "ref_vs_canu.filtered.gp"   "ref_vs_canu.filtered"   "E745 reference genome" "Canu assembly"

echo "All replots finished."
echo "Results are in: ${OUTDIR}"
