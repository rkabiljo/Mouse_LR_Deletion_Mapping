#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
deletions_file <- args[1]
gene_file      <- args[2]
output_png     <- args[3]

genome_length <- 16299   # mouse chrM length
soft_black <- "#595959"
default_gene_color <- "#c5d5db"
del_color <- "#8b0000"

dels <- read_tsv(deletions_file, col_types = cols())

if (nrow(dels) == 0) {
  message("No deletions to plot.")
  quit(save = "no", status = 0)
}

genes <- read_tsv(gene_file, col_types = cols())

dels <- dels %>% arrange(size)
n_del <- nrow(dels)

# Packs deletion lines tightly together, close to the outer edge of the band
pack_near_outer <- function(n, inner, outer, ring_spacing = 0.025) {
  if (n <= 1) return(outer)
  needed_range <- (n - 1) * ring_spacing
  available_range <- outer - inner
  if (needed_range <= available_range) {
    offsets <- rev(seq_len(n) - 1) * ring_spacing
    return(outer - offsets)
  } else {
    return(seq(inner, outer, length.out = n))
  }
}
y_positions <- pack_near_outer(n_del, 0.55, 0.95)


png(output_png, width = 4000, height = 4000, res = 400)
par(mar = c(4, 4, 4, 4))
circos.clear()
circos.par(clock.wise = FALSE, xaxis.clock.wise = FALSE, start.degree = 90, gap.after = 0.2)
circos.initialize(factors = "chrM", xlim = c(0, genome_length))

# Track 1: tick marks
circos.trackPlotRegion(factors = "chrM", ylim = c(0, 1), track.height = 0.15,
                        bg.border = NA,
                        panel.fun = function(x, y) {
  tick_interval <- 500
  major_ticks <- seq(0, genome_length, by = tick_interval)
  labels <- ifelse(major_ticks %% 1000 == 0, as.character(major_ticks), "")
  circos.axis(h = "top", major.at = major_ticks, labels = labels,
              labels.cex = 0.5, major.tick.length = 0.05,
              col = soft_black, labels.col = soft_black)
})

# Track 2: mitochondrial genes
for (i in 1:nrow(genes)) {
  gene_start <- genes$start[i]
  gene_end   <- genes$end[i]
  gene_name  <- genes$gene[i]

  if (gene_name %in% c("OL", "OH")) {
    border_col <- "#2C3E50"
    border_lwd <- 2.5
  } else {
    border_col <- soft_black
    border_lwd <- 1
  }

  circos.rect(gene_start, 0, gene_end, 1,
              col = default_gene_color, border = border_col, lwd = border_lwd)

  gene_length <- gene_end - gene_start
  angle_span <- (gene_length / genome_length) * 360
  label_width <- nchar(gene_name) * 6 * 0.4
  label_facing <- if (label_width < angle_span) "inside" else "reverse.clockwise"
  label_pos <- (gene_start + gene_end) / 2

  if (gene_name == "OL") {
    circos.lines(c(label_pos, label_pos), c(0.8, 1.5), col = "navy", lwd = 0.6)
    circos.text(label_pos, 1.6, gene_name, facing = "downward", niceFacing = TRUE,
                cex = 0.7, font = 2, col = "navy")
  } else if (gene_name == "OH") {
    circos.text(label_pos, 0.5, gene_name, facing = label_facing, niceFacing = TRUE,
                cex = 0.7, font = 2, col = "navy")
  } else {
    circos.text(label_pos, 0.5, gene_name, facing = label_facing, niceFacing = TRUE,
                cex = 0.5, font = 2, col = "#333333")
  }
}

# Track 3: every individual deletion event, plain red
circos.trackPlotRegion(factors = "chrM", ylim = c(0, 1),
                        track.height = 0.65, bg.border = NA,
                        panel.fun = function(x, y) {
  for (i in seq_len(n_del)) {
    start_pos <- dels$left_bp[i]
    end_pos   <- dels$right_bp[i]

    if (start_pos > end_pos) {
      bigger  <- max(start_pos, end_pos)
      smaller <- min(start_pos, end_pos)
      circos.segments(bigger, y_positions[i], smaller + genome_length, y_positions[i],
                       col = del_color, lwd = 0.6, lty = 1)
    } else {
      circos.segments(start_pos, y_positions[i], end_pos, y_positions[i],
                       col = del_color, lwd = 0.6, lty = 1)
    }
  }
})

title("Mitochondrial DNA deletions")
invisible(dev.off())

message("Wrote plot to ", output_png)
