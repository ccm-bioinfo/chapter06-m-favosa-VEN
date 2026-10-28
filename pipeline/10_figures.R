# 10_figures.R - Figures: (A) ITS1 top-15 genera + bacterial 16S bars per classification;
# (B) placement cladogram (gappa graft) rooted on Fructilactobacillus; (C) V4 pairwise-difference heatmap.
suppressPackageStartupMessages({library(ape); library(ggtree); library(ggplot2); library(patchwork)
  library(RColorBrewer); library(dplyr)})
D <- Sys.getenv("CH06_WORK"); I <- Sys.getenv("CH06_INPUTS")
F <- file.path(D, "figures"); dir.create(F, showWarnings = FALSE)
pal_for <- function(labs) {           # Paired in alphabetical order; Other grey, Unclassified light grey
  g <- sort(setdiff(labs, c("Other", "Unclassified")))
  cols <- if (length(g) <= 12) brewer.pal(12, "Paired")[seq_along(g)] else colorRampPalette(brewer.pal(12, "Paired"))(length(g))
  c(setNames(cols, g), Other = "grey50", Unclassified = "grey85")
}
# ---- fungi ----
its <- read.delim(file.path(I, "its/ITS_genus_counts_Mfav.tsv"), check.names = FALSE)
its$lab <- ifelse(is.na(its$Genus) | its$Genus == "", "Unclassified", its$Genus)
fz <- its %>% group_by(lab) %>% summarise(n = sum(`B-VEN-D9`)) %>% mutate(p = 100 * n / sum(n))
top <- fz %>% filter(lab != "Unclassified") %>% arrange(desc(n)) %>% head(15) %>% pull(lab)
fz <- fz %>% mutate(lab = ifelse(lab %in% c(top, "Unclassified"), lab, "Other")) %>% group_by(lab) %>%
  summarise(p = sum(p)) %>% mutate(bar = "ITS1")
fz$lab <- factor(fz$lab, levels = c(sort(top), "Other", "Unclassified"))
pf <- ggplot(fz, aes(bar, p, fill = lab)) + geom_col(width = 0.6) + scale_fill_manual(values = pal_for(levels(fz$lab))) +
  labs(x = NULL, y = "% of reads", fill = "Fungal genus", title = "Fungal composition (ITS1)") + theme_bw()
# ---- bacteria ----
L <- read.delim(file.path(D, "results/09_asv_labels_full.tsv"), check.names = FALSE)
cls <- c("RDP", "SILVA", "GTDB", "GTDB+ours", "tree")
long <- do.call(rbind, lapply(cls, function(c) data.frame(bar = c, lab = L[[c]], reads = L$reads)))
long <- long %>% group_by(bar, lab) %>% summarise(n = sum(reads), .groups = "drop") %>% group_by(bar) %>%
  mutate(p = 100 * n / sum(n)) %>% ungroup()
keep <- long %>% group_by(lab) %>% summarise(m = max(p)) %>% filter(m >= 0.3 | lab == "Unclassified") %>% pull(lab)
long <- long %>% mutate(lab = ifelse(lab %in% keep, lab, "Other")) %>% group_by(bar, lab) %>% summarise(p = sum(p), .groups = "drop")
long$bar <- factor(long$bar, levels = cls, labels = c("RDP", "SILVA", "GTDB", "GTDB + our 16S", "Tree placement"))
bl <- c(sort(setdiff(unique(long$lab), c("Other", "Unclassified"))), "Other", "Unclassified")
long$lab <- factor(long$lab, levels = bl)
pb <- ggplot(long, aes(bar, p, fill = lab)) + geom_col(width = 0.7) + scale_fill_manual(values = pal_for(bl)) +
  labs(x = NULL, y = "% of reads", fill = "Bacterial label", title = "Bacterial composition (16S V4)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(F, "fig_bars_its_16S.pdf"), pf + pb + plot_layout(widths = c(1, 3.2)), width = 13, height = 7)
ggsave(file.path(F, "fig_bars_its_16S.png"), pf + pb + plot_layout(widths = c(1, 3.2)), width = 13, height = 7, dpi = 200)
write.table(long, file.path(D, "results/10_bar_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# ---- placement tree ----
tr <- read.tree(file.path(D, "place/graft/epa_result.newick"))
tips <- read.delim(file.path(D, "results/05_tree_tips.tsv"))
out <- tips$tip[tips$genus == "Fructilactobacillus"]
tr <- root(tr, outgroup = out, resolve.root = TRUE)
lac <- L[L$amptk_lacto == "yes", ]
meta <- data.frame(label = tr$tip.label) %>%
  left_join(tips %>% transmute(label = tip, genus, clade, shape), by = "label") %>%
  left_join(L %>% transmute(label = ASV, asvlab = tree, reads), by = "label") %>%
  mutate(is_asv = grepl("^ASV", label),
         group = ifelse(is_asv, asvlab, ifelse(clade != "" & !is.na(clade), clade, genus)),
         genus2 = ifelse(is_asv, ifelse(grepl("Clade 3-4", asvlab), "Acetilactobacillus",
                          ifelse(grepl("Clade", asvlab), "Nicoliella", asvlab)), genus),
         shape = ifelse(is_asv, "ASV", shape),
         size = ifelse(is_asv, 100 * reads / sum(lac$reads), NA))
p <- ggtree(tr, branch.length = "none") %<+% meta
gen <- c("Acetilactobacillus", "Nicoliella", "Philodulcilactobacillus", "Fructilactobacillus")
gcol <- setNames(brewer.pal(4, "Pastel1"), gen)
for (g in gen) {
  tl <- meta$label[meta$genus2 == g & !is.na(meta$genus2)]
  if (length(tl) > 1) { n <- getMRCA(tr, tl); p <- p + geom_hilight(node = n, fill = gcol[g], alpha = 0.5, extend = 2) }
}
grp <- sort(unique(meta$group)); gpal <- pal_for(grp)
p <- p + geom_tippoint(data = function(d) d[d$isTip & !d$is_asv, ], aes(shape = shape, colour = group), size = 2.5) +
  geom_tippoint(data = function(d) d[d$isTip & d$is_asv, ], aes(size = size, colour = group), shape = 16, alpha = 0.8) +
  geom_tiplab(data = function(d) d[d$isTip & (!d$is_asv | d$size >= 0.5), ], size = 2.3, offset = 0.3) +
  scale_shape_manual(values = c(triangle = 17, square = 15, circle = 16), name = "Reference") +
  scale_size_continuous(name = "% Lactobacillaceae reads", range = c(0.8, 8)) +
  scale_colour_manual(values = gpal, name = "Genus / clade") + xlim(0, 45) +
  ggtitle("16S placement tree (cladogram; shaded: genus; ASVs grafted at best placement)")
ggsave(file.path(F, "fig_placement_tree.pdf"), p, width = 11, height = 16)
ggsave(file.path(F, "fig_placement_tree.png"), p, width = 11, height = 16, dpi = 200)
# ---- V4 heatmap ----
M <- read.delim(file.path(D, "results/09_v4_differences.tsv"), row.names = 1, check.names = FALSE)
ord <- rev(get_taxa_name(p))
hm <- as.data.frame(as.table(as.matrix(M[ord, ord]))); names(hm) <- c("a", "b", "d")
hm$a <- factor(hm$a, levels = ord); hm$b <- factor(hm$b, levels = rev(ord))
ph <- ggplot(hm, aes(a, b, fill = d)) + geom_tile() + scale_fill_viridis_c(name = "V4 differences", option = "magma", direction = -1) +
  theme_minimal(base_size = 6) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x = NULL, y = NULL)
ggsave(file.path(F, "fig_v4_heatmap.pdf"), ph, width = 14, height = 13)
ggsave(file.path(F, "fig_v4_heatmap.png"), ph, width = 14, height = 13, dpi = 200)
cat("figures written to", F, "\n")
