#!/usr/bin/env Rscript
# Step 5 (Mac): figures for the chapter.
#  fig1_classification_by_method  % of B-VEN-D9 Lactobacillaceae reads by label, per classification method
#  fig2_placement_tree            circular 16S tree (GTDB + our genomes + MG16S) with the favosa ASVs grafted
# Usage: Rscript scripts/12_figures.R   (needs results/11_compare and results/08_tree)
suppressPackageStartupMessages({ library(ape); library(ggplot2); library(ggtree); library(dplyr) })
out <- "results/12_figures"; dir.create(out, showWarnings = FALSE)

# categorical slots in fixed order (reference palette), neutrals for "not resolved"
COL <- c("Clado 3 Acetilactobacillus-like" = "#2a78d6", "Clado 4 Acetilactobacillus-like" = "#eb6834",
         "Clado 2 Nicoliella-related" = "#1baf7a", "Acetilactobacillus" = "#eda100",
         "Other genera" = "#9a9890", "Lactobacillaceae (family only)" = "#d6d4cd",
         "Not resolved to Lactobacillaceae" = "#52514e")
INK <- "#0b0b0b"; INK2 <- "#52514e"

# ---- Figure 1 ---------------------------------------------------------------------------------
s <- read.delim("results/11_compare/summary_by_method.tsv", stringsAsFactors = FALSE)
s$label[!s$label %in% names(COL)] <- "Other genera"
s <- s |> group_by(method, label) |> summarise(pct = sum(pct_LAB_reads), .groups = "drop")
methods <- c("RDP (AMPtk)", "SILVA 138.2", "GTDB r226", "GTDB + our 16S", "Tree placement (EPA-ng)")
s$method <- factor(s$method, levels = rev(methods))
s$label  <- factor(s$label, levels = names(COL))
s <- s |> arrange(method, label) |> group_by(method) |>
  mutate(mid = 100 - (cumsum(pct) - pct / 2)) |> ungroup()
tot <- sum(read.delim("results/04_queries/favosa_counts.tsv")$reads[read.delim("results/04_queries/favosa_counts.tsv")$is_LAB == "True"])

p1 <- ggplot(s, aes(pct, method, fill = label)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.6) +
  geom_text(data = subset(s, pct >= 5), aes(x = mid, label = sprintf("%.0f%%", pct)),
            colour = "white", size = 3.2, fontface = "bold") +
  scale_fill_manual(values = COL, drop = FALSE, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), labels = function(x) paste0(x, "%")) +
  labs(x = "Lactobacillaceae reads", y = NULL,
       title = "Our 16S resolve the Melipona favosa honey Lactobacillaceae to paper clades",
       subtitle = sprintf("Sample 10MFH (B-VEN-D9), 69 ASVs, %s reads; SINTAX cutoff 0.8; tree = EPA-ng best placement",
                          format(tot, big.mark = ","))) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8, colour = INK2),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = INK2), plot.title = element_text(colour = INK, face = "bold", size = 11),
        plot.subtitle = element_text(colour = INK2, size = 8.5)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))
ggsave(file.path(out, "fig1_classification_by_method.pdf"), p1, width = 8, height = 4.2)
ggsave(file.path(out, "fig1_classification_by_method.png"), p1, width = 8, height = 4.2, dpi = 300, bg = "white")

# ---- Figure 2 ---------------------------------------------------------------------------------
tr   <- read.tree("results/08_tree/graft/epa_result.newick")
info <- read.delim("results/08_tree/tip_info.tsv", stringsAsFactors = FALSE)
cls  <- read.delim("results/11_compare/favosa_LAB_classification.tsv", stringsAsFactors = FALSE)
out_tips <- info$tip[info$family %in% c("Streptococcaceae", "Enterococcaceae")]
tr <- root(tr, outgroup = intersect(out_tips, tr$tip.label), resolve.root = TRUE)

d <- data.frame(label = tr$tip.label, stringsAsFactors = FALSE) |>
  left_join(info, by = c("label" = "tip")) |>
  left_join(cls[, c("ASV", "reads")], by = c("label" = "ASV")) |>
  mutate(type = ifelse(is.na(type), "favosa ASV", type),
         clade_lab = gsub("_", " ", sub("^(Clado_\\d)_", "\\1 ", clade)),
         shown = case_when(type == "favosa ASV" ~ paste0(label, " (", format(reads, big.mark = ",", trim = TRUE), ")"),
                           type == "GTDB" ~ label.y,
                           type == "genome" ~ paste0(label.y, if_else(clade != "", paste0(" [", clade_lab, "]"), "")),
                           TRUE ~ paste0("MG16S ", label.y, if_else(clade != "", paste0(" [", clade_lab, "]"), ""))),
         who = recode(type, GTDB = "GTDB r226 reference", genome = "Our genome 16S",
                      metagenome = "16S recovered from our metagenomes", `favosa ASV` = "M. favosa ASV"),
         is_asv = type == "favosa ASV", is_ours = type %in% c("genome", "metagenome"))
# clade colour: references by their paper clade, ASVs by their EPA-ng clade (clades are not monophyletic in
# the 16S tree, so clade shading by MRCA would be misleading)
asv_clade <- setNames(sub("^(Clado_\\d)_", "\\1 ", sub(".*;(Clado_[^;]+)$", "\\1", cls$EPA_taxopath)), cls$ASV)
d$cl <- ifelse(d$is_asv, gsub("_", " ", asv_clade[d$label]), d$clade_lab)
d$cl[!d$cl %in% names(COL)] <- NA
d$shown[d$is_asv & d$reads < 50] <- ""      # label only abundant ASVs

p2 <- ggtree(tr, layout = "circular", size = 0.25, colour = INK2) %<+% d +
  geom_tippoint(aes(subset = is_asv, size = reads, fill = cl), shape = 21, colour = "white", stroke = 0.3) +
  geom_tippoint(aes(subset = is_ours, shape = who, colour = cl), size = 1.6) +
  scale_fill_manual(values = COL, na.value = "#9a9890", name = "Paper clade (Xolalpa-Aroche et al. 2026)",
                    breaks = names(COL)[1:3]) +
  scale_colour_manual(values = COL, na.value = INK, guide = "none") +
  scale_shape_manual(values = c("Our genome 16S" = 15, "16S recovered from our metagenomes" = 17), name = NULL) +
  scale_size_area(max_size = 5, breaks = c(10, 100, 1000, 10000), name = "M. favosa ASV reads") +
  ggnewscale::new_scale_colour() +
  geom_tiplab2(aes(label = shown, colour = who), size = 1.3, offset = 0.02, show.legend = FALSE) +
  scale_colour_manual(values = c("GTDB r226 reference" = "#9a9890", "Our genome 16S" = INK,
                                 "16S recovered from our metagenomes" = INK, "M. favosa ASV" = INK), guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 3)), size = guide_legend(override.aes = list(fill = INK))) +
  labs(title = "Placement of Melipona favosa honey Lactobacillaceae ASVs",
       subtitle = "16S tree constrained to the GTDB r226 genome topology; ASVs grafted at their EPA-ng best placement") +
  theme(legend.position = "right", legend.text = element_text(size = 7, colour = INK2),
        legend.title = element_text(size = 7.5, colour = INK),
        plot.title = element_text(face = "bold", size = 11, colour = INK),
        plot.subtitle = element_text(size = 8, colour = INK2))
p2 <- p2 + xlim(0, max(p2$data$x, na.rm = TRUE) * 1.3)
ggsave(file.path(out, "fig2_placement_tree.pdf"), p2, width = 11, height = 9.5)
ggsave(file.path(out, "fig2_placement_tree.png"), p2, width = 11, height = 9.5, dpi = 300, bg = "white")
message("figures in ", out)
