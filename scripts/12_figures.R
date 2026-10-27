#!/usr/bin/env Rscript
# Step 5 (Mac): figures for the chapter.
#  fig1_genus_by_method  whole 10MFH bacterial community by genus, one stacked bar per classification method,
#                        in the style and colours of chapter Fig. 6.4B (amplicon/scripts/microbiome_plots.R:
#                        RColorBrewer Paired in alphabetical order, Other grey65, Unclassified grey88)
#  fig2_placement_tree   circular 16S cladogram with genus sectors and the favosa ASVs grafted at their
#                        EPA-ng placement, in the style of amplicon/scripts/plot_lab_tree_circular.R
# Usage: Rscript scripts/12_figures.R   (needs results/11_compare and results/08_tree)
suppressPackageStartupMessages({ library(ape); library(ggplot2); library(ggtree); library(dplyr) })
out <- "results/12_figures"; dir.create(out, showWarnings = FALSE)

# colours shared by both figures for the new assignments (from the VIT tree palette)
CLADE_COL <- c("Clado 3 Acetilactobacillus-like" = "#D55E00", "Clado 4 Acetilactobacillus-like" = "#009E73",
               "Clado 2 Nicoliella-related" = "#CC79A7", "Acetilactobacillus" = "#E69F00")

# ---- Figure 1 ---------------------------------------------------------------------------------
METHODS <- c("RDP (AMPtk)" = "RDP (AMPtk)\nsubmitted, Fig. 6.4B", "SILVA 138.2" = "SILVA 138.2",
             "GTDB r226" = "GTDB r226", "GTDB + our 16S" = "GTDB r226\n+ our 16S",
             "Tree placement (EPA-ng)" = "Tree placement\n(EPA-ng)")
s <- read.delim("results/11_compare/summary_genus_all.tsv", stringsAsFactors = FALSE)
keep <- unique(s$genus[s$pct_reads >= 0.3])                  # shown in at least one bar, as in Fig. 6.4B
s$genus[!s$genus %in% keep] <- "Other"
s <- s |> group_by(method, genus) |> summarise(pct = sum(pct_reads) / 100, .groups = "drop")
s$method <- factor(METHODS[s$method], levels = METHODS)

special <- c("Other", "Unclassified")
genera  <- sort(setdiff(unique(s$genus), c(special, names(CLADE_COL))))
# Fig. 6.4B colours for the genera it shows; Paired order continues for the rest
fig64 <- c(Bacteroides = "#A6CEE3", Comamonas = "#1F78B4", Dysgonomonas = "#B2DF8A", Fusobacterium = "#33A02C",
           Haemophilus = "#FB9A99", Lactobacillus = "#E31A1C", Neisseria = "#FDBF6F", Rothia = "#FF7F00",
           Sphingobacterium = "#CAB2D6", Streptococcus = "#6A3D9A", Veillonella = "#FFFF99")
extra <- c("#B15928", "#8DD3C7", "#BC80BD", "#80B1D3", "#FCCDE5")
new_g <- setdiff(genera, names(fig64))
pal <- c(fig64[intersect(genera, names(fig64))], setNames(extra[seq_along(new_g)], new_g),
         CLADE_COL[intersect(names(CLADE_COL), s$genus)], Other = "grey65", Unclassified = "grey88")
lev <- c(names(CLADE_COL)[names(CLADE_COL) %in% s$genus], genera, special)
s$genus <- factor(s$genus, levels = rev(lev))
s <- s |> arrange(method, genus) |> group_by(method) |> mutate(mid = cumsum(pct) - pct / 2) |> ungroup()

p1 <- ggplot(s, aes(method, pct, fill = genus)) +
  geom_col(width = 0.75) +
  geom_text(data = subset(s, pct >= 0.05), aes(y = 1 - mid, label = sprintf("%.0f%%", 100 * pct)),
            size = 3, colour = "grey15") +
  scale_fill_manual(values = pal, breaks = lev, name = "Genus / paper clade") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_bw(base_size = 11) +
  theme(legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"),
        panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 8.5)) +
  labs(title = "16S Genus Composition: Melipona favosa honey (10MFH) by reference",
       subtitle = "Same ASVs; only the classification reference changes. Genera < 0.3% in every bar are grouped as Other.",
       x = NULL, y = "Relative Abundance")
ggsave(file.path(out, "fig1_genus_by_method.pdf"), p1, width = 9, height = 6)
ggsave(file.path(out, "fig1_genus_by_method.png"), p1, width = 9, height = 6, dpi = 300, bg = "white")

# ---- Figure 2 ---------------------------------------------------------------------------------
tree <- read.tree("results/08_tree/graft/epa_result.newick")
info <- read.delim("results/08_tree/tip_info.tsv", stringsAsFactors = FALSE)
cls  <- read.delim("results/11_compare/favosa_LAB_classification.tsv", stringsAsFactors = FALSE)

out_tips <- info$tip[info$family %in% c("Streptococcaceae", "Enterococcaceae")]
tree <- root(tree, outgroup = intersect(out_tips, tree$tip.label), resolve.root = TRUE)
# ASVs with an anomalous terminal branch (> 0.3 subst/site) are left out of the drawing, as in the VIT tree
tl   <- setNames(tree$edge.length[match(seq_len(Ntip(tree)), tree$edge[, 2])], tree$tip.label)
drop <- names(tl)[tl > 0.3 & grepl("^ASV", names(tl))]
if (length(drop) > 0) tree <- drop.tip(tree, drop)
tree <- ladderize(tree)

epa_last <- sub(".*;", "", cls$EPA_taxopath)
asv <- data.frame(tip = cls$ASV, type = "ASV", family = "Lactobacillaceae",
                  genus = sapply(strsplit(cls$EPA_taxopath, ";"), function(x) if (length(x) > 1) x[2] else ""),
                  clade = ifelse(grepl("^Clado", epa_last), epa_last, ""),
                  label = cls$ASV, pct = cls$pct_LAB_reads)
d <- bind_rows(mutate(info, pct = NA_real_), asv) |>
  mutate(clade = gsub("_", " ", clade),
         grupo = case_when(
           type == "GTDB"       ~ "GTDB r226 reference",
           type == "genome"     ~ "Our genome (published)",
           type == "metagenome" ~ "16S recovered from our metagenomes",
           clade != ""          ~ paste0("ASV: ", clade),
           genus == "Acetilactobacillus" ~ "ASV: Acetilactobacillus",
           TRUE                 ~ "ASV: other / family only"),
         tiplab = case_when(
           type == "ASV" & pct >= 0.5 ~ sprintf("%s (%.1f%%)", label, pct),
           type == "ASV"              ~ "",
           type == "GTDB"             ~ label,
           TRUE ~ paste0(sub("^(MG16S|OWN)_", "", tip),
                         ifelse(clade != "", paste0(" [", sub("^(Clado \\d).*", "\\1", clade), "]"), "")))) |>
  rename(name = label)

cols <- c(setNames(CLADE_COL, paste0("ASV: ", names(CLADE_COL))),
          "ASV: other / family only" = "#999999",
          "Our genome (published)" = "#0072B2", "16S recovered from our metagenomes" = "#56B4E9",
          "GTDB r226 reference" = "grey55")
cols <- cols[names(cols) %in% d$grupo]

# genus sectors: MRCA of the GTDB and genome references of each genus (constrained to the genome topology)
hl_genera <- c(Acetilactobacillus = "#F5DEB3", Nicoliella = "#F2D7E8", Apilactobacillus = "#D9EAD3",
               Philodulcilactobacillus = "#E8E8E8", Fructobacillus = "#DDEBF7", Bombilactobacillus = "#FFF2CC",
               Lactobacillus = "#EDE7F6", Fructilactobacillus = "#E2F0D9")
hl <- do.call(rbind, lapply(names(hl_genera), function(g) {
  # Acetilactobacillus: all our sequences of the genus too (MG16S and ASVs), as the VIT tree did with its own
  # genomes; the Clado 3-4 block falls outside the GTDB Acetilactobacillus group in the 16S tree
  types <- if (g == "Acetilactobacillus") c("GTDB", "genome", "metagenome", "ASV") else c("GTDB", "genome")
  tips <- intersect(d$tip[d$type %in% types & d$genus == g], tree$tip.label)
  if (length(tips) < 2) return(NULL)
  data.frame(node = getMRCA(tree, tips), genus_hl = g)
}))

p2 <- ggtree(tree, layout = "circular", branch.length = "none", size = 0.25) %<+% d +
  geom_hilight(data = hl, aes(node = node, fill = genus_hl), alpha = 0.6, extend = 6) +
  scale_fill_manual(values = hl_genera, name = "Genus (shading)") +
  geom_tippoint(aes(colour = grupo, size = ifelse(type == "ASV", pmax(pct, 0.05), 0.4)), alpha = 0.9) +
  geom_tiplab(aes(label = tiplab, colour = grupo), size = 1.7, offset = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = "Sequence") +
  scale_size_continuous(range = c(0.5, 6), name = "% of Lactobacillaceae reads (ASV)",
                        breaks = c(0.1, 1, 5, 20, 60)) +
  theme(legend.position = "right", legend.text = element_text(size = 8)) +
  labs(title = "16S Lactobacillaceae: Melipona favosa honey ASVs (10MFH) on GTDB r226 + our published genomes and recovered 16S",
       subtitle = paste("Reference tree constrained to the GTDB genome topology (bac120); ASVs placed with EPA-ng;",
                        "cladogram (no branch lengths). ASV colour = EPA-ng assignment; label if >= 0.5 %.",
                        if (length(drop) > 0) paste("\nLeft out for anomalous branch (> 0.3):", paste(drop, collapse = ", ")) else ""))
ggsave(file.path(out, "fig2_placement_tree.pdf"), p2, width = 18, height = 16, limitsize = FALSE)
ggsave(file.path(out, "fig2_placement_tree.png"), p2, width = 18, height = 16, dpi = 200, bg = "white")
message("figures in ", out)
