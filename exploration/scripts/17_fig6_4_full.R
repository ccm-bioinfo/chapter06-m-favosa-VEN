#!/usr/bin/env Rscript
# Full chapter Fig. 6.4 on one 15 x 11 in page (panels of the author's Illustrator draft fig4b_fullv2.ai):
#  (A) ITS1 genus bar and (B) 16S genus bars by classification reference, top left;
#  (C) circular placement tree of Fig. 4b (Acetilactobacillus at the top), full height on the right, because at
#      10 pt its 109 tip labels need a ~10 in circle; shading and size legends in its corners, the other legends
#      below A and B.
# All text is Arial 10 pt, genera and species in italics (GTDB species with the genus abbreviated, as the shading
# legend names the genera), little white space; editable SVG.
# Usage (Mac, after 12_figures.R and 14_fig6_4.R): Rscript exploration/scripts/17_fig6_4_full.R -> results/17_fig6_4_full/
suppressPackageStartupMessages({ library(ape); library(ggplot2); library(ggtree); library(patchwork) })
out <- "results/17_fig6_4_full"; dir.create(out, showWarnings = FALSE)
FS <- 10; FONT <- "Arial"
GENERA <- c("Acetilactobacillus", "Nicoliella", "Fructilactobacillus")

# legend label with the genus in italics, e.g. "ASV: Clade 3-4 Acetilactobacillus-like"
ital <- function(l) {
  g <- GENERA[sapply(GENERA, grepl, l)]
  if (length(g) == 0) return(l)
  parts <- strsplit(l, g[1], fixed = TRUE)[[1]]
  pre <- parts[1]; post <- if (length(parts) > 1) parts[2] else ""
  bquote(.(pre) * italic(.(g[1])) * .(post))
}
# every text in the panel at FS points (bar percentages included)
txt10 <- function(p) {
  for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomText")) p$layers[[i]]$aes_params$size <- FS / .pt
  p
}
th <- theme_bw(base_size = FS, base_family = FONT) +
  theme(text = element_text(size = FS), axis.text = element_text(size = FS, colour = "black"),
        legend.text = element_text(size = FS), legend.title = element_text(size = FS, face = "bold"),
        plot.title = element_text(size = FS, face = "bold"), plot.tag = element_text(size = FS, face = "bold"),
        legend.text.align = 0, legend.key.size = unit(0.35, "cm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"), plot.margin = margin(2, 2, 2, 2), panel.grid.major.x = element_blank())

# ---- A and B --------------------------------------------------------------------------------------
a <- txt10(readRDS("results/14_fig6_4/fig6_4A.rds")) + th +
  labs(tag = "A", title = "Fungi (ITS1)") +
  theme(legend.position = "right") + guides(fill = guide_legend(ncol = 1))
b <- txt10(readRDS("results/12_figures/fig1_genus_by_method.rds")) + th +
  labs(tag = "B", title = "Bacteria (16S rRNA) by reference", subtitle = NULL) +
  scale_x_discrete(labels = c("RDP\n(AMPtk)", "SILVA\n138.2", "GTDB\nr226", "GTDB\n+ our\n16S", "Tree\n(EPA-ng)")) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(legend.position = "right") + guides(fill = guide_legend(ncol = 2, byrow = FALSE))
sc <- b$scales$get_scales("fill")
sc$labels <- lapply(sc$breaks, function(l)
  if (l %in% c("Other", "Unclassified", "Lactobacillaceae (no genus)")) l
  else if (grepl("Clade", l)) ital(l) else bquote(italic(.(l))))

# ---- C ---------------------------------------------------------------------------------------------
k <- readRDS("results/12_figures/fig4b_parts.rds")
sub <- k$sub; acet4 <- k$acet4; c34 <- k$c34
TIP <- FS / .pt
p <- ggtree(sub, layout = "circular", branch.length = "none", size = 0.3) %<+% k$d
dep <- max(p$data$x)
p <- p +
  geom_hilight(data = k$hl4, aes(node = node, fill = genus_hl), alpha = 0.6, extend = 0.75 * dep) +
  scale_fill_manual(values = k$hl_genera, name = "Genus (shading)",
                    labels = function(l) lapply(l, function(x) bquote(italic(.(x))))) +
  geom_tippoint(aes(colour = grupo, shape = grupo, size = ifelse(type == "ASV", pmax(pct, 0.05), ifelse(type == "genome", 1.6, 0.8))),
                alpha = 0.9) +
  geom_text(data = function(x) {
    x <- subset(x, isTip)
    x$tl <- ifelse(x$type == "GTDB", sub("^(Acetilactobacillus|Fructilactobacillus|Nicoliella) ", "\\1.", x$name), x$tiplab)
    x$tl <- sub("^([A-Z])[a-z]+\\.", "\\1. ", x$tl)
    x$ff <- ifelse(x$type == "GTDB", "italic", "plain")          # GTDB species names in italics
    a <- x$angle %% 360
    flip <- x$label %in% c34 | (!x$label %in% acet4 & a > 90 & a < 270)
    transform(x, ang = ifelse(flip, angle + 180, angle), hj = ifelse(flip, 1, 0))
  }, aes(x = x + 0.4, label = tl, colour = grupo, angle = ang, hjust = hj, fontface = ff),
  size = TIP, family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = k$cols, name = "Sequence", labels = function(l) lapply(l, ital)) +
  scale_shape_manual(values = k$shapes, name = "Sequence", labels = function(l) lapply(l, ital)) +
  scale_size_continuous(range = c(0.6, 6), name = "% of Lactobacillaceae\nreads (ASV)", breaks = c(0.1, 1, 5, 20, 60)) +
  # coord_polar draws the radius limit at 80% of the panel; values past it (labels, shading) are kept and fill
  # the outer 20%, so the tree uses the whole panel
  scale_x_continuous(limits = c(0, 1.4 * dep), oob = scales::oob_keep, expand = c(0, 0)) +
  # shading and size legends in the empty corners of the circle; the sequence legend goes to the right column
  guides(fill = guide_legend(position = "inside", theme = theme(legend.position.inside = c(0, 1), legend.justification.inside = c(0, 1))),
         size = guide_legend(position = "inside", ncol = 1, theme = theme(legend.position.inside = c(1, 0), legend.justification.inside = c(1, 0))),
         colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(tag = "C") +
  theme(text = element_text(size = FS, family = FONT), legend.text = element_text(size = FS),
        legend.title = element_text(size = FS, face = "bold"), plot.tag = element_text(size = FS, face = "bold"),
        legend.position = "right", legend.key.size = unit(0.35, "cm"), legend.text.align = 0,
        legend.box.spacing = unit(0, "pt"), plot.margin = margin(0, 0, 0, 0))
acet_mid <- mean(range(subset(p$data, isTip & label %in% acet4)$angle))
p <- rotate_tree(p, 90 - acet_mid)

# ---- page ------------------------------------------------------------------------------------------
# left: A and B, then the legends of A, B and the tree sequences; right: the tree, full height
W <- 15; H <- 11
leg <- function(x) cowplot::get_plot_component(x, "guide-box-right")
legA <- leg(a + guides(fill = guide_legend(ncol = 2)))
legB <- leg(b + guides(fill = guide_legend(ncol = 2)))
legC <- leg(p + guides(fill = "none", size = "none"))
a <- a + theme(legend.position = "none")
b <- b + theme(legend.position = "none")
p <- p + guides(colour = "none", shape = "none")
left <- (a + b + plot_layout(widths = c(1, 3.6))) / wrap_elements(full = legA) / wrap_elements(full = legB) /
  wrap_elements(full = legC) + plot_layout(heights = c(4.4, 1.9, 2.1, 2.2))
fig <- wrap_plots(left, wrap_elements(full = p), widths = c(4.4, 10.6))
ggsave(file.path(out, "fig6_4_full.pdf"), fig, width = W, height = H, device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
ggsave(file.path(out, "fig6_4_full.png"), fig, width = W, height = H, dpi = 250, bg = "white")
ggsave(file.path(out, "fig6_4_full.svg"), fig, width = W, height = H, device = svglite::svglite, fix_text_size = FALSE)   # editable text
message("written ", out)
