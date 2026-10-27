#!/usr/bin/env bash
# Step 6: why our 16S resolve the favosa LAB: best-hit V4 identity of each Lactobacillaceae ASV to
# our 16S (OWN_*) versus GTDB r226 alone. Top hits only (an all-hits output fills the disk).
# Usage (alnitak, after 06_sintax_refs.py): bash scripts/13_identity.sh  ->  results/13_identity/
set -euo pipefail
cd "$(dirname "$0")/../results"
export PATH=/data/nselem/envs/vit/bin:$PATH
O=13_identity; mkdir -p $O; export TMPDIR=$PWD/$O
seqkit grep -r -p "^OWN_" 06_sintax/refs/gtdb_ours.fa > $O/ours.fa
seqkit grep -r -v -p "^OWN_" 06_sintax/refs/gtdb_ours.fa > $O/gtdb_only.fa
for db in ours gtdb_only; do
  vsearch --usearch_global 04_queries/favosa_LAB.fa --db $O/$db.fa --id 0.7 --maxaccepts 32 --maxrejects 2000 \
    --top_hits_only --userout $O/best_$db.tsv --userfields query+target+id --quiet --threads 32
done
rm $O/ours.fa $O/gtdb_only.fa
python3 -I - <<PY
import csv
reads = {r[0]: int(r[1]) for r in csv.reader(open("04_queries/favosa_counts.tsv"), delimiter="\t") if r[0] != "ASV"}
def best(f):
    d = {}
    for q, t, i in csv.reader(open(f), delimiter="\t"):
        if float(i) > d.get(q, (0,))[0]:
            d[q] = (float(i), t.split(";")[0], (t.split(",g:") + [""])[1].split(",")[0])
    return d
o, g = best("$O/best_ours.tsv"), best("$O/best_gtdb_only.tsv")
lab = sorted(reads.keys() & (o.keys() | g.keys()), key=lambda a: -reads[a])
with open("$O/identity_ours_vs_gtdb.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV", "reads", "best_gtdb_id", "best_gtdb_genus", "best_ours_id", "best_ours_ref"])
    for a in lab:
        w.writerow([a, reads[a], g.get(a, (0, "", ""))[0], g.get(a, (0, "", ""))[2], o.get(a, (0, ""))[0], o.get(a, (0, ""))[1]])
tot = sum(reads[a] for a in lab)
for thr in (100, 99, 97):
    print(f">={thr}% identity: ours {100 * sum(reads[a] for a in lab if o.get(a, (0,))[0] >= thr) / tot:.1f}% | "
          f"GTDB {100 * sum(reads[a] for a in lab if g.get(a, (0,))[0] >= thr) / tot:.1f}% of LAB reads")
PY
