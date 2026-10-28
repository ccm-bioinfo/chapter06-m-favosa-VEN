#!/bin/bash
# 06_build_tree.sh - Reference 16S tree: MAFFT L-INS-i, trimAl -gt 0.05, IQ-TREE (GTR+F+G4, 1000 UFBoot,
# constrained by the pruned GTDB bac120 tree + genus polytomies), RAxML-NG re-estimation of model parameters.
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$CH06_WORK/tree"
python3 "$CH06_PIPELINE/05_tree_inputs.py"
mafft --localpair --maxiterate 1000 --thread $THREADS refs.fa > refs.mafft.fa 2> mafft.log
trimal -in refs.mafft.fa -out refs.trim.fa -gt 0.05
iqtree3 -s refs.trim.fa -m GTR+F+G4 -B 1000 -g constraint.nwk -T 8 --seed 1 --prefix iqtree -redo > iqtree.stdout 2>&1
raxml-ng --evaluate --msa refs.trim.fa --tree iqtree.treefile --model GTR+G4+FC --threads 4 --seed 1 \
  --prefix raxml --redo > raxml.stdout 2>&1
echo DONE
