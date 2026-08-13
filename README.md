# Mouse mtDNA Deletion Mapping

Code used to generate the mitochondrial DNA deletion plots for:

> Diego Perez-Rodriguez, Ilaria Dalla Rosa, Anastasia Magoulopoulou, Sara Ricciardi, Samantha Davidson, Rebecca Lasalandra, Renata Kabiljo, Aine Moylett, Tamara Hill ,Tengfei Wan, Aleck W.E. Jones, Chloe F. Moss, Charlotte Zierz, Melissa L. Salazar, Radha Desai, Shar-yin N. Huang, Yves Pommier, Monika Hofer, Robert D.S. Pitceathly, Martin A.M. Reijns, Joanna Poulton, Christos Proukakis, J. Paul Simons, Mats Nilsson, Robert W. Taylor, Ian J. Holt* & Antonella Spinazzola*. 
> A stress-adapted fibre state promotes clonal expansion of deleted
> mitochondrial DNA. *Cell Reports*. 2026.

## Overview

This code aligns ONT reads (unaligned bam files) to the mouse mitochondrial genome (chrM), extracts
large deletions directly from read CIGAR strings, and plots each deletion
as an arc on a circular map of mouse chrM.

## Files

- `mice_deletions.sh` – main script: alignment, deletion extraction, plotting
- `deletion_analysis.py` – extracts deletions ≥500 bp from an aligned BAM
- `plot_dels.R` – draws the circular deletion plot

## Resources

- `mito_genes_mouse.tsv` – mouse mitochondrial gene annotation
- `hg38_chrM_mouse.fa`  mouse reference chrM and indices

## Requirements

- `samtools`, `minimap2`, `seqkit`
- Python 3 with `pysam`
- R with `circlize`, `dplyr`, `readr`

Create a conda environment for dependencies 
```bash
conda create -n map_mouse -c conda-forge -c bioconda \
    python=3.11 r-base=4.3 minimap2 samtools seqkit \
    pysam r-dplyr r-readr r-circlize
conda activate map_mouse
```

## Usage

Edit the `sample`, and input BAM path at the top of
`mice_deletions.sh`, then run:

```bash
bash mice_deletions.sh
```

`REF` should point to a mouse chrM reference FASTA (`resources/hg38_chrM_mouse.fa`).

Outputs are written to `outputs/`:
- `<sample>_aligned_sorted.bam` – aligned reads
- `<sample>_deletions.tsv` – individual deletion events (read name, breakpoints, size)
- `<sample>_deletions_plot.png` – circular deletion plot

## Notes

- Deletion detection uses a minimum size threshold of 500 bp.
- Each row in `<sample>_deletions.tsv` is one deletion event from one read;
  no clustering or filtering by read support is applied.

