#!/bin/bash
#SBATCH -A uppmax2026-1-94
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -J check_tnseq
#SBATCH --mail-type=ALL
#SBATCH -o /home/yuch3531/genome_analysis_project/6_differential_gene_analysis/tnseq/check_tnseq_v3_%j.out
#SBATCH -e /home/yuch3531/genome_analysis_project/6_differential_gene_analysis/tnseq/check_tnseq_v3_%j.err

set -euo pipefail

INPUT_BASE="/home/yuch3531/genome_analysis_project/1_data/transcriptomics_data"
OUTPUT_DIR="/home/yuch3531/genome_analysis_project/6_differential_gene_analysis/tnseq"

TAG="GCCAAT"
MOTIF="ACAGGTTGGATGATAAGTCCCCGGTCT"

SUMMARY_FILE="${OUTPUT_DIR}/tnseq_structure_summary_v3.tsv"
EXAMPLE_FILE="${OUTPUT_DIR}/tnseq_structure_examples_v3.tsv"
VARLEN_FILE="${OUTPUT_DIR}/tnseq_variable_length_distribution_v3.tsv"
TAILLEN_FILE="${OUTPUT_DIR}/tnseq_tail_length_distribution_v3.tsv"

mkdir -p "$OUTPUT_DIR"

echo -e "sample\tcondition\tfile\ttotal_reads\tstart_with_TAG\tcontains_Magellan6_RC\tTAG_before_Magellan6_RC\tTAG_before_Magellan6_RC_percent" > "$SUMMARY_FILE"
echo -e "sample\tcondition\tfile\tparsed_read" > "$EXAMPLE_FILE"
echo -e "sample\tcondition\tfile\tvariable_length\tcount" > "$VARLEN_FILE"
echo -e "sample\tcondition\tfile\ttail_length_after_motif\tcount" > "$TAILLEN_FILE"

for CONDITION_DIR in "${INPUT_BASE}"/Tn-Seq_*
do
    CONDITION=$(basename "$CONDITION_DIR")

    echo "========================================"
    echo "Checking condition: $CONDITION"
    echo "Directory: $CONDITION_DIR"

    while IFS= read -r FQ
    do
        SAMPLE=$(basename "$FQ")
        TMP="${OUTPUT_DIR}/${SAMPLE}.tmp_tnseq_check.tsv"

        echo "----------------------------------------"
        echo "Checking file: $SAMPLE"
        echo "Path: $FQ"

        zcat "$FQ" \
            | awk 'NR%4==2' \
            | awk -v tag="$TAG" -v motif="$MOTIF" '
                {
                    read = $0
                    total++

                    starts = (index(read, tag) == 1)
                    motif_pos = index(read, motif)

                    if (starts) {
                        start_tag++
                    }

                    if (motif_pos > 0) {
                        motif_count++
                    }

                    if (starts && motif_pos > 0) {
                        tag_before_motif++

                        var_len = motif_pos - length(tag) - 1
                        tail_len = length(read) - motif_pos - length(motif) + 1

                        var_count[var_len]++
                        tail_count[tail_len]++

                        if (example_count < 10) {
                            var_seq = substr(read, length(tag)+1, var_len)
                            tail_seq = substr(read, motif_pos + length(motif))
                            print "EXAMPLE\t" tag " | " var_seq " | " motif " | " tail_seq
                            example_count++
                        }
                    }
                }

                END {
                    print "SUMMARY\ttotal_reads\t" total
                    print "SUMMARY\tstart_with_TAG\t" start_tag+0
                    print "SUMMARY\tcontains_Magellan6_RC\t" motif_count+0
                    print "SUMMARY\tTAG_before_Magellan6_RC\t" tag_before_motif+0

                    for (l in var_count) {
                        print "VARLEN\t" l "\t" var_count[l]
                    }

                    for (t in tail_count) {
                        print "TAILLEN\t" t "\t" tail_count[t]
                    }
                }
            ' > "$TMP"

        TOTAL=$(awk '$1=="SUMMARY" && $2=="total_reads" {print $3}' "$TMP")
        START_TAG=$(awk '$1=="SUMMARY" && $2=="start_with_TAG" {print $3}' "$TMP")
        CONTAINS_MOTIF=$(awk '$1=="SUMMARY" && $2=="contains_Magellan6_RC" {print $3}' "$TMP")
        TAG_MOTIF=$(awk '$1=="SUMMARY" && $2=="TAG_before_Magellan6_RC" {print $3}' "$TMP")

        PERCENT=$(awk -v a="$TAG_MOTIF" -v b="$TOTAL" 'BEGIN {if (b>0) printf "%.2f", a/b*100; else print "NA"}')

        echo -e "${SAMPLE}\t${CONDITION}\t${FQ}\t${TOTAL}\t${START_TAG}\t${CONTAINS_MOTIF}\t${TAG_MOTIF}\t${PERCENT}" >> "$SUMMARY_FILE"

        awk -v sample="$SAMPLE" -v condition="$CONDITION" -v file="$FQ" '
            $1=="EXAMPLE" {
                sub(/^EXAMPLE\t/, "")
                print sample "\t" condition "\t" file "\t" $0
            }
        ' "$TMP" >> "$EXAMPLE_FILE"

        awk -v sample="$SAMPLE" -v condition="$CONDITION" -v file="$FQ" '
            $1=="VARLEN" {
                print sample "\t" condition "\t" file "\t" $2 "\t" $3
            }
        ' "$TMP" | sort -k4,4n >> "$VARLEN_FILE"

        awk -v sample="$SAMPLE" -v condition="$CONDITION" -v file="$FQ" '
            $1=="TAILLEN" {
                print sample "\t" condition "\t" file "\t" $2 "\t" $3
            }
        ' "$TMP" | sort -k4,4n >> "$TAILLEN_FILE"

        rm -f "$TMP"

    done < <(find "$CONDITION_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)

done

echo "========================================"
echo "Done."
echo "Summary file: $SUMMARY_FILE"
echo "Example file: $EXAMPLE_FILE"
echo "Variable length file: $VARLEN_FILE"
echo "Tail length file: $TAILLEN_FILE"
