#!/usr/bin/env bash
# run_all.sh - Run the whole pipeline (steps 00-09) on a Linux server; step 10 (figures, R with ggtree) can run on
# any machine that has $CH06_WORK. Long: use screen/tmux or setsid nohup bash pipeline/run_all.sh > run.log 2>&1 &
set -euo pipefail
D=$(cd "$(dirname "$0")" && pwd)
source "$D/config.sh"
bash "$D/00_get_inputs.sh"
python3 "$D/01_build_refs.py"
bash "$D/03_sintax.sh"          # runs 02_build_ours_refs.py first
bash "$D/04_besthit_gtdb.sh"
bash "$D/06_build_tree.sh"      # runs 05_tree_inputs.py first
bash "$D/08_place.sh"           # runs 07_labels.py first
python3 "$D/09_summarise.py" | tee "$CH06_WORK/results/09_summary.txt"
command -v Rscript > /dev/null && Rscript "$D/10_figures.R" || echo "run 10_figures.R where R + ggtree are installed"
# large formatted references are rebuilt by 01_build_refs.py in minutes
rm -f "$CH06_WORK"/refs/gtdb_r226_ge1200.sintax.fa "$CH06_WORK"/refs/gtdb_r226_ge1200_noN.fa "$CH06_WORK"/refs/silva138.2_nr99.sintax.fa
echo PIPELINE DONE
