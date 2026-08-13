#!/usr/bin/env bash
set -euo pipefail

sample=240
REF="resources/mouse_chrM.fa"
GENES_TSV="resources/mito_genes_mouse.tsv"
OUTDIR="outputs"

mkdir -p "$OUTDIR"

samtools merge -u - /data/hestia/rkabiljo/ONT_mice/${sample}/*.bam \
    | samtools fastq -T MM,ML - \
    | seqkit seq -m 500 \
    | minimap2 -ax map-ont -y -t 8 "$REF" - \
    | samtools view -bh -o "${OUTDIR}/${sample}_aligned.bam" -
samtools sort -o "${OUTDIR}/${sample}_aligned_sorted.bam" "${OUTDIR}/${sample}_aligned.bam"
samtools index "${OUTDIR}/${sample}_aligned_sorted.bam"

python3 deletion_analysis.py \
    "${OUTDIR}/${sample}_aligned_sorted.bam" \
    "${OUTDIR}/${sample}_deletions.tsv"

Rscript plot_dels.R \
    "${OUTDIR}/${sample}_deletions.tsv" \
    "$GENES_TSV" \
    "${OUTDIR}/${sample}_deletions_plot.png"
