#!/bin/bash
# 08_place.sh - Add the selected ASVs to the trimmed reference alignment (MAFFT --addfragments --keeplength),
# place them with EPA-ng (RAxML-NG model), assign genus/clade of the best placement with gappa (rooted on
# Fructilactobacillus) and graft them onto the tree.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$CH06_WORK"
python3 "$CH06_PIPELINE/07_labels.py"
mkdir -p place; cd place
mafft --addfragments ../work/place_query.fa --keeplength --thread 8 ../tree/refs.trim.fa > ref_plus_asv.fa 2> mafft_add.log
seqkit grep -r -p '^ASV' ref_plus_asv.fa > query.aln.fa
awk -F'\t' '$2 ~ /Fructilactobacillus/ {print $1}' ../tree/taxon_file.tsv > outgroup.txt
epa-ng --ref-msa ../tree/refs.trim.fa --tree ../tree/raxml.raxml.bestTree --query query.aln.fa \
  --model ../tree/raxml.raxml.bestModel --threads 8 --redo -w . > epa.log 2>&1
gappa examine assign --jplace-path epa_result.jplace --taxon-file ../tree/taxon_file.tsv \
  --root-outgroup outgroup.txt --best-hit --per-query-results --allow-file-overwriting --out-dir assign > assign.log 2>&1
gappa examine graft --jplace-path epa_result.jplace --allow-file-overwriting --out-dir graft > graft.log 2>&1
echo DONE
