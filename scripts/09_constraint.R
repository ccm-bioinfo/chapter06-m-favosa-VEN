#!/usr/bin/env Rscript
# Step 3b (Mac; alnitak has no ape): constraint for the 16S reference tree.
# GTDB r226 genome topology (GTDB_* tips of data/refs/constraint_gtdb.tree). Our genome 16S and the MG16S with a
# paper clade are added as a polytomy at the MRCA of their genus (free within the genus); MG16S without a
# clade link stay unconstrained.
# Usage: Rscript scripts/09_constraint.R   (needs results/08_tree/tip_info.tsv from alnitak)
suppressPackageStartupMessages(library(ape))
dir  <- "results/08_tree"
info <- read.delim(file.path(dir, "tip_info.tsv"), stringsAsFactors = FALSE)

ct <- read.tree("data/refs/constraint_gtdb.tree")
ct <- keep.tip(ct, intersect(ct$tip.label, info$tip[info$type == "GTDB"]))
genus_of <- setNames(info$genus, info$tip)

add <- info[(info$type == "genome") | (info$type == "metagenome" & info$clade != ""), ]
extra <- list()
for (i in seq_len(nrow(add))) {
  gtips <- ct$tip.label[genus_of[ct$tip.label] %in% add$genus[i]]
  if (length(gtips) == 0) { message("no GTDB reference for genus (left free): ", add$tip[i]); next }
  key <- as.character(if (length(gtips) == 1) which(ct$tip.label == gtips) else getMRCA(ct, gtips))
  extra[[key]] <- c(extra[[key]], add$tip[i])
}
nwk <- function(node) {
  kids <- ct$edge[ct$edge[, 1] == node, 2]
  base <- if (node <= Ntip(ct)) ct$tip.label[node] else paste0("(", paste(sapply(kids, nwk), collapse = ","), ")")
  more <- extra[[as.character(node)]]
  if (is.null(more)) base else if (node <= Ntip(ct)) paste0("(", paste(c(base, more), collapse = ","), ")")
                      else paste0("(", paste(c(sapply(kids, nwk), more), collapse = ","), ")")
}
txt <- paste0(nwk(Ntip(ct) + 1), ";")
writeLines(txt, file.path(dir, "constraint.tree"))
chk <- read.tree(text = txt)
message("constraint: ", Ntip(chk), " tips (", sum(chk$tip.label %in% add$tip), " of ours constrained; ",
        sum(!info$tip %in% chk$tip.label), " free)")
