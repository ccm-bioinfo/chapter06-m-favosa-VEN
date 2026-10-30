#!/usr/bin/env Rscript
# 10MFH bacterial community without Acetilactobacillus (final classification of 11_compare.py: organelles removed
# with SILVA, Lactobacillaceae with GTDB r226 + our 16S, the rest with SILVA 138.2): the Clade 3-4 and the
# Acetilactobacillus reads (88.8 + 0.5% of all 16S reads) and the organelle 16S are removed, the rest rescaled to 100%.
# One stacked bar; genera >= 2% of the remaining reads shown, the rest Other, no genus Unclassified.
# Colours as in fig1_genus_by_method (12_figures.R), so each genus keeps its Fig. 6.4B colour.
# Usage (Mac, after 11_compare.py): Rscript exploration/scripts/18_non_acetilactobacillus_bar.R
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
out <- "results/18_non_acetilactobacillus"; dir.create(out, showWarnings = FALSE)
CUT <- 0.02
FINAL <- "Final (SILVA + GTDB + our 16S)"
ORG <- "Organelles (mitochondria, chloroplast)"
DROP <- c("Clade 3-4 Acetilactobacillus-like", "Acetilactobacillus", ORG)

s <- read.delim("results/11_compare/summary_genus_all.tsv", stringsAsFactors = FALSE)
# palette of 12_figures.R fig1 (same genus set, same order), reused here
keep_all <- unique(s$genus[s$pct_reads >= 0.3])
special <- c("Other", "Lactobacillaceae (no genus)", "Unclassified", ORG)
CLADE <- c("Clade 3-4 Acetilactobacillus-like", "Clade 2 Nicoliella-related", "Acetilactobacillus")
fig64 <- c(Bacteroides = "#A6CEE3", Comamonas = "#1F78B4", Dysgonomonas = "#B2DF8A", Fusobacterium = "#33A02C",
           Haemophilus = "#FB9A99", Lactobacillus = "#E31A1C", Neisseria = "#FDBF6F", Rothia = "#FF7F00",
           Sphingobacterium = "#CAB2D6", Streptococcus = "#6A3D9A", Veillonella = "#FFFF99")
extra <- c("#B15928", "#8DD3C7", "#BC80BD", "#80B1D3", "#FCCDE5")
new_g <- setdiff(sort(setdiff(keep_all, c(special, CLADE))), names(fig64))
pal <- c(fig64, setNames(extra[seq_along(new_g)], new_g), Other = "grey65", Unclassified = "grey88")

g <- s |> filter(method == FINAL, !genus %in% DROP)
n_rest <- sum(g$reads)
g <- g |> mutate(pct = reads / n_rest,
                 genus = ifelse(genus == "Unclassified" | pct >= CUT, genus, "Other")) |>
  group_by(genus) |> summarise(reads = sum(reads), pct = sum(pct), .groups = "drop")
named <- g$genus[!g$genus %in% special][order(-g$pct[!g$genus %in% special])]
stopifnot(all(named %in% names(pal)))
lev <- c(named, "Other", "Unclassified")
g$genus <- factor(g$genus, levels = rev(lev))
g <- g |> arrange(desc(genus)) |> mutate(mid = cumsum(pct) - pct / 2)   # ggplot stacks the last level at the bottom
write.table(g |> arrange(desc(genus)) |> select(genus, reads, pct), file.path(out, "non_acetilactobacillus_genus.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(g, aes("10MFH", pct, fill = genus)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.6) +
  geom_text(data = subset(g, pct >= 0.03), aes(y = mid, label = sprintf("%.0f%%", 100 * pct)),
            size = 4.5, colour = "grey15") +
  scale_fill_manual(values = pal[lev], breaks = lev, name = "Genus",
                    labels = lapply(lev, function(l) if (l %in% special) l else bquote(italic(.(l))))) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 1.0001),
                     breaks = seq(0, 1, 0.25)) +
  theme_bw(base_size = 13) +
  theme(legend.text.align = 0, legend.text = element_text(size = 12), legend.title = element_text(size = 13, face = "bold"),
        panel.grid.major.x = element_blank()) +
  labs(title = "16S genera without Acetilactobacillus",
       subtitle = sprintf(paste0("Final classification; Clade 3-4, Acetilactobacillus\nand organelle 16S removed ",
                                 "(%s of %s reads\nremain, rescaled to 100%%); genera < %d%% as Other"),
                          format(n_rest, big.mark = ","), format(sum(s$reads[s$method == FINAL]), big.mark = ","),
                          100 * CUT),
       x = NULL, y = "Relative Abundance")
ggsave(file.path(out, "non_acetilactobacillus_bar.pdf"), p, width = 6, height = 7)
ggsave(file.path(out, "non_acetilactobacillus_bar.png"), p, width = 6, height = 7, dpi = 300, bg = "white")
message("written ", out)
