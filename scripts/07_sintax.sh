#!/usr/bin/env bash
# Step 2b: SINTAX of every B-VEN-D9 ASV against each reference (cutoff 0.8, fixed seed, 1 thread
# per run for reproducibility; the 4 runs go in parallel). RDP/AMPtk is the submitted classification.
# Usage (alnitak): bash scripts/07_sintax.sh   ->  results/06_sintax/<ref>.tsv
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=/data/nselem/envs/vit/bin:$PATH
R=results/06_sintax; Q=results/04_queries/favosa_all.fa

for ref in silva gtdb gtdb_ours gtdb_ours_sens; do
  vsearch --sintax $Q --db $R/refs/$ref.fa --tabbedout $R/$ref.tsv \
    --sintax_cutoff 0.8 --randseed 1 --threads 1 --quiet &
done
wait
wc -l $R/*.tsv

# cleanup: the formatted reference FASTAs are rebuilt by 06_sintax_refs.py in minutes
rm -f $R/refs/silva.fa $R/refs/gtdb.fa $R/refs/gtdb_ours_sens.fa
