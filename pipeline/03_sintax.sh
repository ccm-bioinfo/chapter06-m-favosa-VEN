#!/bin/bash
# 03_sintax.sh - SINTAX (cutoff 0.8, seed 1, 1 thread for reproducibility) of the 273 ASVs against
# SILVA 138.2, GTDB r226 (>=1200 bp), GTDB+ours and GTDB+ours without MAG_41/MAG_45 (sensitivity);
# exhaustive best V4 hit of the 69 Lactobacillaceae ASVs against our 16S (against GTDB: 04_besthit_gtdb.sh).
set -euo pipefail
source "$(dirname "$0")/config.sh"; I=$CH06_INPUTS; cd "$CH06_WORK"
python3 "$CH06_PIPELINE/02_build_ours_refs.py"
cat refs/gtdb_r226_ge1200.sintax.fa refs/ours.sintax.fa      > refs/gtdb_ours.sintax.fa
cat refs/gtdb_r226_ge1200.sintax.fa refs/ours_sens.sintax.fa > refs/gtdb_ours_sens.sintax.fa
for db in silva138.2_nr99 gtdb_r226_ge1200 gtdb_ours gtdb_ours_sens; do
  vsearch --sintax $I/asv/B-VEN-D9_ASVs.fa --db refs/$db.sintax.fa --sintax_cutoff 0.8 --randseed 1 \
    --threads 1 --no_progress --tabbedout results/03_sintax_${db}.tsv 2> logs/03_sintax_${db}.log
done
vsearch --usearch_global work/lacto_asvs.fa --db refs/ours.sintax.fa --id 0.5 --maxaccepts 0 --maxrejects 0 \
  --top_hits_only --threads $THREADS --no_progress --userfields query+target+id+alnlen+mism+opens \
  --userout results/03_besthit_ours.tsv --output_no_hits 2> logs/03_besthit_ours.log
rm -f refs/gtdb_ours.sintax.fa refs/gtdb_ours_sens.sintax.fa
echo DONE
