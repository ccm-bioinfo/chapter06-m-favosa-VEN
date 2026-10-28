#!/bin/bash
# 04_besthit_gtdb.sh - Exhaustive (maxaccepts 0, maxrejects 0) best V4 hit of the 69 Lactobacillaceae ASVs against
# GTDB r226 (>= 1200 bp) without sequences with ambiguous bases (an N counts as a match in vsearch).
set -euo pipefail
source "$(dirname "$0")/config.sh"; cd "$CH06_WORK"
vsearch --usearch_global work/lacto_asvs.fa --db refs/gtdb_r226_ge1200_noN.fa --id 0.5 --maxaccepts 0 \
  --maxrejects 0 --top_hits_only --threads $THREADS --no_progress --userfields query+target+id+alnlen+mism+opens \
  --userout results/04_besthit_gtdb_noN.tsv --output_no_hits 2> logs/04_besthit_gtdb.log
echo DONE
