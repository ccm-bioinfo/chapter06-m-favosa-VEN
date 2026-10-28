#!/usr/bin/env Rscript
# Chapter Fig. 6.4 redrawn for Melipona favosa honey (10MFH, B-VEN-D9), panels with the same height and fonts:
#  (A) ITS1 genus composition, same data and rule as the original panel (amplicon/scripts/microbiome_plots.R,
#      plot_taxa_bar): top 15 genus-level taxa, the rest "Other", no genus "Unclassified"; RColorBrewer Paired
#      in alphabetical order (interpolated if > 12 genera)
#  (B) 16S genus composition per classification reference (fig1 of 12_figures.R, replaces the original 16S bar)
# Also writes table6_4_fungi_genus.tsv: fungal genera as % of all ITS reads, of genus-assigned reads and of the
# 15 listed genera (the denominator of the submitted Table 6.4).
# Usage (Mac, after 12_figures.R): Rscript exploration/scripts/14_fig6_4.R   ->  results/14_fig6_4/
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

a <- panel(panel_data("data/vit/ITS_genus_relabund_Mfav.tsv"), "A", "Fungi (ITS1)")
saveRDS(a, file.path(out, "fig6_4A.rds"))   # panel A of the full Fig. 6.4 (17_fig6_4_full.R)
b <- readRDS("results/12_figures/fig1_genus_by_method.rds") +
  labs(tag = "B", title = "Bacteria (16S rRNA) by classification reference", subtitle = NULL) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE, title.position = "top"))

# one look for both panels: same base size, legend text and key size, legends below the bars
common <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.text.align = 0, legend.key.size = unit(0.45, "cm"),
        legend.justification = "left", axis.text.x = element_text(size = 10.5),
        plot.title = element_text(size = 13), plot.tag = element_text(size = 16, face = "bold"),
        panel.grid.major.x = element_blank())
a <- a + common + guides(fill = guide_legend(ncol = 2, title.position = "top"))
b <- b + common
p <- a + b + plot_layout(widths = c(1, 4.2)) &
  plot_annotation(caption = paste("Melipona favosa pot-honey (10MFH). (A) The 15 most abundant fungal taxa; the rest as Other.",
                                  "(B) The same 16S ASVs classified with RDP (submitted), SILVA 138.2, GTDB r226, GTDB r226 + our 16S",
                                  "from stingless bee MAGs and genomes, and EPA-ng tree placement; genera < 0.3% in every bar as Other.",
                                  sep = "\n"),
                  theme = theme(plot.caption = element_text(size = 10, hjust = 0)))
ggsave(file.path(out, "fig6_4_ITS_16S.pdf"), p, width = 15, height = 9)
ggsave(file.path(out, "fig6_4_ITS_16S.png"), p, width = 15, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out, "fig6_4_ITS_16S.svg"), p, width = 15, height = 9, device = svglite::svglite)

# fungal genus table (counts behind Table 6.4)
cnt <- read.delim("data/vit/ITS_genus_counts_Mfav.tsv", stringsAsFactors = FALSE)
cnt$n <- cnt[[ncol(cnt)]]
total <- sum(cnt$n); g <- cnt[!is.na(cnt$Genus) & cnt$Genus != "NA", ]
g <- g |> group_by(Genus) |> summarise(counts = sum(n), .groups = "drop") |> arrange(-counts)
top15 <- sum(head(g$counts, 15))
tab <- g |> mutate(pct_all_reads = round(100 * counts / total, 1),
                   pct_genus_assigned = round(100 * counts / sum(counts), 1),
                   pct_of_top15 = round(100 * counts / top15, 1))
write.table(tab, file.path(out, "table6_4_fungi_genus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("written ", out, "; ITS reads: ", total, " (genus-assigned ", sum(g$counts), ", top 15 genera ", top15, ")")
