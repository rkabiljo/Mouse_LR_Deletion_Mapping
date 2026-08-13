#!/usr/bin/env bash
set -euo pipefail

sample=<SAMPLE_NAME> 
REF="resources/hg38_chrM_mouse.fa"
GENES_TSV="resources/mito_genes_mouse.tsv"
OUTDIR="outputs"

mkdir -p "$OUTDIR"

#merge all bam files, align, sort
samtools merge -u - <PATH_TO_UNALIGNED_BAM_FILES>/${sample}/*.bam \
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
