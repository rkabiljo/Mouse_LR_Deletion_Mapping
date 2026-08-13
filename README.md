# Mouse mtDNA Deletion Mapping

Code used to generate the mitochondrial DNA deletion plots for:

> Diego Perez-Rodriguez1¶, Ilaria Dalla Rosa1¶, Anastasia Magoulopoulou2, Sara Ricciardi1, Samantha Davidson1, Rebecca Lasalandra1, Renata Kabiljo3, Aine Moylett3, Tamara Hill3 ,Tengfei Wan1, Aleck W.E. Jones1, Chloe F. Moss4, Charlotte Zierz5, Melissa L. Salazar1, Radha Desai4, Shar-yin N. Huang6, Yves Pommier6, Monika Hofer7, Robert D.S. Pitceathly3, Martin A.M. Reijns8, Joanna Poulton9, Christos Proukakis1, J. Paul Simons10, Mats Nilsson2, Robert W. Taylor5,11, Ian J. Holt1,12* & Antonella Spinazzola1*. 
> A stress-adapted fibre state promotes clonal expansion of deleted
> mitochondrial DNA. *Cell Reports*. 2026.

## Overview

This code aligns reads to the mouse mitochondrial genome (chrM), extracts
large deletions directly from read CIGAR strings, and plots each deletion
as an arc on a circular map of chrM.

## Files

- `mice_deletions.sh` – main script: alignment, deletion extraction, plotting
- `deletion_analysis.py` – extracts deletions ≥500 bp from an aligned BAM
- `plot_dels.R` – draws the circular deletion plot
- `resources: mito_genes_mouse.tsv` – mouse mitochondrial gene annotation  and mouse reference chrM

## Requirements

- `samtools`, `minimap2`, `seqkit`
- Python 3 with `pysam`
- R with `circlize`, `dplyr`, `readr`

```bash
conda create -n map_mouse -c conda-forge -c bioconda \
    python=3.11 r-base=4.3 minimap2 samtools seqkit \
    pysam r-dplyr r-readr r-circlize
conda activate map_mouse
```

## Usage

Edit the `sample`, `REF`, and input BAM path at the top of
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

