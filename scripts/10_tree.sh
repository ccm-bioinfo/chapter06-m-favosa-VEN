#!/usr/bin/env bash
# Step 3c: 16S reference tree constrained to the GTDB genome topology + EPA-ng placement of the
# B-VEN-D9 Lactobacillaceae ASVs (same recipe as /data2/MIEL/vit/scripts/31_tree_mg16S.sh).
# Needs results/08_tree/{refs_full.fa,queries.fa,taxon.tsv,constraint.tree}.
# Usage (alnitak, ~16 threads): setsid nohup bash scripts/10_tree.sh > logs/10_tree.log 2>&1 &
set -euo pipefail
cd "$(dirname "$0")/../results/08_tree"
export PATH=/data/nselem/envs/vit/bin:$PATH
T=16

mafft --thread $T --maxiterate 1000 --localpair refs_full.fa > refs_full.aln 2> mafft_refs.log
mafft --thread $T --addfragments queries.fa --keeplength refs_full.aln > all.aln 2> mafft_add.log
trimal -in all.aln -out all.trim.aln -gt 0.05
seqkit grep -r -v -p '^ASV' all.trim.aln > refs.trim.aln
seqkit grep -r -p '^ASV' all.trim.aln > query.trim.aln

iqtree3 -s refs.trim.aln -m GTR+F+G4 -g constraint.tree -B 1000 -T $T --prefix ref_constr -redo > iqtree.stdout 2>&1
raxml-ng --evaluate --msa refs.trim.aln --tree ref_constr.contree --model GTR+G+F --prefix refeval --threads $T --redo > /dev/null
epa-ng --ref-msa refs.trim.aln --tree refeval.raxml.bestTree --query query.trim.aln \
  --model refeval.raxml.bestModel -T $T --outdir . --redo > epa.log
gappa examine assign --jplace-path epa_result.jplace --taxon-file taxon.tsv \
  --per-query-results --best-hit --out-dir assign --allow-file-overwriting > /dev/null
# accumulated LWR per taxon (all placements, not only the best edge): support for genus vs clade
gappa examine assign --jplace-path epa_result.jplace --taxon-file taxon.tsv \
  --per-query-results --out-dir assign_all --allow-file-overwriting > /dev/null
gappa examine graft --jplace-path epa_result.jplace --out-dir graft --allow-file-overwriting > /dev/null
grep -E "Log-likelihood of the tree" ref_constr.iqtree

# cleanup: keep alignments, trees, jplace and assignments
rm -f ref_constr.ckp.gz ref_constr.mldist ref_constr.uniqueseq.phy ref_constr.splits.nex ref_constr.bionj \
      refeval.raxml.rba refeval.raxml.reduced.phy refeval.raxml.startTree all.aln
echo DONE
