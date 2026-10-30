#!/usr/bin/env Rscript
# Where the 10MFH 16S reads that RDP (AMPtk, submitted chapter) did not call Lactobacillaceae go in the final
# classification of 11_compare.py (organelles removed with SILVA, Lactobacillaceae with GTDB r226 + our 16S,
# the rest with SILVA 138.2): alluvial from the RDP genus (left) to the final genus (right), one flow per
# ASV pooled by the pair of labels. The RDP Lactobacillaceae without genus and Lactobacillus ASVs are left out.
# RDP labels from the AMPtk taxonomy: genus, else "<family> (no genus)", else "No hit" / "Unclassified".
# Labels with < MIN reads on a side are pooled as Other; flows coloured by the RDP label (Fig. 6.4B colours).
# Usage (Mac, after 11_compare.py): Rscript exploration/scripts/19_rdp_gtdb_alluvial.R
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(ggalluvial); library(ggrepel) })
out <- "results/19_rdp_gtdb_alluvial"; dir.create(out, showWarnings = FALSE)
MIN <- 40

a <- read.delim("results/11_compare/genus_by_asv.tsv", check.names = FALSE, stringsAsFactors = FALSE)
a <- a[!a[["RDP (AMPtk)"]] %in% c("Lactobacillaceae (no genus)", "Lactobacillus"), ]
rank_of <- function(tax, r) {
  m <- regmatches(tax, regexpr(paste0(r, ":[^,]+"), tax))
  if (length(m) == 0) NA else gsub('"', "", sub("^.:", "", m))
}
a$rdp <- sapply(a$RDP_taxonomy, function(t) {
  g <- rank_of(t, "g"); f <- rank_of(t, "f")
  if (!is.na(g)) g else if (!is.na(f)) paste0(f, " (no genus)") else if (t == "No hit") "No hit" else "Unclassified"
})
a$gtdb <- a[["Final (SILVA + GTDB + our 16S)"]]
fold <- function(x, w) { keep <- names(which(tapply(w, x, sum) >= MIN)); ifelse(x %in% keep, x, "Other") }
a$rdp <- fold(a$rdp, a$reads); a$gtdb <- fold(a$gtdb, a$reads)

f <- a |> group_by(rdp, gtdb) |> summarise(reads = sum(reads), n_ASVs = n(), .groups = "drop")
n_all <- sum(f$reads)
write.table(f |> arrange(desc(reads)), file.path(out, "rdp_to_final.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

special <- c("No hit", "Unclassified", "Organelles (mitochondria, chloroplast)", "Other")
ord <- function(x) { s <- tapply(f$reads, x, sum); c(setdiff(names(sort(-s)), special), intersect(special, names(s))) }
f$rdp <- factor(f$rdp, levels = ord(f$rdp))
f$gtdb <- factor(paste0(f$gtdb, "\u200b"), levels = paste0(ord(f$gtdb), "\u200b"))   # zero-width suffix: strata of
                                                                                     # the two axes stay apart

# Fig. 6.4B colours (12_figures.R) for the genera it has; Paired-like extras for the rest
fig64 <- c(Bacteroides = "#A6CEE3", Comamonas = "#1F78B4", Dysgonomonas = "#B2DF8A", Fusobacterium = "#33A02C",
           Haemophilus = "#FB9A99", Lactobacillus = "#E31A1C", Neisseria = "#FDBF6F", Rothia = "#FF7F00",
           Sphingobacterium = "#CAB2D6", Streptococcus = "#6A3D9A", Veillonella = "#FFFF99")
extra <- c("#B15928", "#8DD3C7", "#BC80BD", "#80B1D3", "#FCCDE5")
new_l <- setdiff(levels(f$rdp), c(names(fig64), special))
pal <- c(fig64, setNames(extra[seq_along(new_l)], new_l), "No hit" = "grey88", Unclassified = "grey80", Other = "grey65",
         "Organelles (mitochondria, chloroplast)" = "grey40")

nice <- function(x) sub(" \\((no genus|mitochondria)", "\n(\\1", sub("\u200b$", "", x))
# stratum label positions: ggalluvial stacks the first level on top
side <- function(col, x, h) {
  s <- tapply(f$reads, f[[col]], sum)[levels(f[[col]])]
  data.frame(x = x, y = n_all - (cumsum(s) - s / 2), hjust = h,
             label = sprintf("%s  %s", nice(names(s)), trimws(format(s, big.mark = ","))))
}
labs_df <- rbind(side("rdp", 1 - 0.1, 1), side("gtdb", 2 + 0.1, 0))
p <- ggplot(f, aes(axis1 = rdp, axis2 = gtdb, y = reads)) +
  geom_alluvium(aes(fill = rdp), width = 0.18, alpha = 0.85, knot.pos = 0.4, colour = "white", linewidth = 0.2) +
  geom_stratum(width = 0.18, fill = "white", colour = "grey30", linewidth = 0.3) +
  geom_text_repel(data = subset(labs_df, x < 1.5), aes(x, y, label = label), inherit.aes = FALSE, hjust = 1,
                  nudge_x = -0.12, direction = "y", size = 3.6, segment.colour = "grey50", min.segment.length = 0,
                  box.padding = 0.12, lineheight = 0.85) +
  geom_text_repel(data = subset(labs_df, x > 1.5), aes(x, y, label = label), inherit.aes = FALSE, hjust = 0,
                  nudge_x = 0.12, direction = "y", size = 3.6, segment.colour = "grey50", min.segment.length = 0,
                  box.padding = 0.12, lineheight = 0.85) +
  scale_x_discrete(limits = c("RDP (AMPtk)\nsubmitted", "SILVA + GTDB\n+ curated 16S\n(final)"), expand = c(0.45, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = pal, guide = "none") +
  theme_void(base_size = 13) +
  theme(axis.text.x = element_text(size = 13, face = "bold"), plot.title = element_text(face = "bold"),
        plot.margin = margin(8, 8, 8, 8)) +
  labs(title = "10MFH 16S reads outside Lactobacillaceae (RDP): RDP genus to final genus",
       subtitle = sprintf(paste0("%s reads (%.1f%% of the 23,491 16S reads); RDP Lactobacillaceae without genus and ",
                                 "Lactobacillus left out.\nFlow colour = RDP label; numbers = reads; labels with < %d reads ",
                                 "on a side pooled as Other."), format(n_all, big.mark = ","), 100 * n_all / 23491, MIN))
ggsave(file.path(out, "rdp_final_alluvial.pdf"), p, width = 11, height = 9)
ggsave(file.path(out, "rdp_final_alluvial.png"), p, width = 11, height = 9, dpi = 250, bg = "white")
message("written ", out)
