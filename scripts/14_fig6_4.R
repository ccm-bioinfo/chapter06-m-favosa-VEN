#!/usr/bin/env Rscript
# Chapter Fig. 6.4 redrawn: (A) ITS1 and (B) 16S rRNA genus composition of Melipona favosa honey (10MFH,
# B-VEN-D9), side by side with the same height and fonts. Same data and rule as the original panels
# (amplicon/scripts/microbiome_plots.R, plot_taxa_bar): top 15 genus-level taxa, the rest "Other", taxa
# without genus "Unclassified"; RColorBrewer Paired in alphabetical order (interpolated if > 12 genera).
# Usage (Mac): Rscript scripts/14_fig6_4.R   ->  results/14_fig6_4/
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork); library(RColorBrewer) })
out <- "results/14_fig6_4"; dir.create(out, showWarnings = FALSE)
TOP <- 15

panel_data <- function(file) {
  d <- read.delim(file, stringsAsFactors = FALSE)
  d$ab <- d[[ncol(d)]]
  top <- d$GenusLabel[order(-d$ab)][seq_len(TOP)]
  d$taxon <- ifelse(!d$GenusLabel %in% top, "Other", ifelse(is.na(d$Genus) | d$Genus == "NA", "Unclassified", d$Genus))
  d |> group_by(taxon) |> summarise(ab = sum(ab), .groups = "drop")
}

panel <- function(d, tag, title) {
  special <- c("Other", "Unclassified")
  genera  <- sort(setdiff(d$taxon, special))
  n   <- length(genera)
  pal <- c(if (n <= 12) brewer.pal(max(3, n), "Paired")[seq_len(n)] else colorRampPalette(brewer.pal(12, "Paired"))(n),
           "grey65", "grey88")
  names(pal) <- c(genera, special)
  lev <- c(genera, intersect(special, d$taxon))
  d$taxon <- factor(d$taxon, levels = rev(lev))
  d <- d |> arrange(taxon) |> mutate(mid = 1 - (cumsum(ab) - ab / 2))
  ggplot(d, aes("honey", ab, fill = taxon)) +
    geom_col(width = 0.7) +
    geom_text(data = subset(d, ab >= 0.05), aes(y = mid, label = sprintf("%.0f%%", 100 * ab)), size = 4) +
    scale_fill_manual(values = pal, breaks = lev, name = "Genus",   # genera in italics, Other/Unclassified upright
                      labels = lapply(lev, function(l) if (l %in% special) l else bquote(italic(.(l))))) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 1.0001),
                       breaks = seq(0, 1, 0.25)) +
    labs(tag = tag, title = title, x = NULL, y = "Relative Abundance") +
    theme_bw(base_size = 12) +
    theme(legend.text = element_text(size = 11), legend.text.align = 0, legend.title = element_text(size = 12),
          legend.key.size = unit(0.45, "cm"), axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 13), plot.tag = element_text(size = 16, face = "bold"),
          panel.grid.major.x = element_blank())
}

a <- panel(panel_data("data/vit/ITS_genus_relabund_Mfav.tsv"), "A", "ITS1 genus composition")
b <- panel(panel_data("data/vit/16S_genus_relabund_Mfav.tsv"), "B", "16S rRNA genus composition")
p <- a + b + plot_layout(widths = c(1, 1)) &
  plot_annotation(caption = "Melipona favosa pot-honey (10MFH). Genera outside the 15 most abundant taxa are grouped as Other.",
                  theme = theme(plot.caption = element_text(size = 10, hjust = 0)))
ggsave(file.path(out, "fig6_4_ITS_16S.pdf"), p, width = 11, height = 6.5)
ggsave(file.path(out, "fig6_4_ITS_16S.png"), p, width = 11, height = 6.5, dpi = 300, bg = "white")
message("written ", out)
