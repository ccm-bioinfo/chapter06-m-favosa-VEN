#!/usr/bin/env python3
"""Step 3a: sequences and tip table for the Lactobacillaceae 16S tree.
  GTDB   the same GTDB r226 accessions as the previous constrained tree (tips GTDB_* of
         data/refs/constraint_gtdb.tree; rule in /data2/MIEL/vit/scripts/08_lab_tree_inputs.py), longest 16S
         each, minus 3 with defective 16S (saturated terminal branch)
  OWN_   our genome 16S with use=yes, Lactobacillaceae, >= 1200 bp (one per genome)
  MG16S_ the 23 Lactobacillaceae 16S recovered from our metagenomes (>= 1000 bp)
  ASV    the 69 B-VEN-D9 Lactobacillaceae ASVs (queries, placed with EPA-ng)
Writes results/08_tree/: refs_full.fa, queries.fa, tip_info.tsv, taxon.tsv (gappa)
Usage (alnitak): python3 scripts/08_tree_inputs.py"""
import csv, gzip, os, re

BAD = {"GTDB_GCF_016861915.1", "GTDB_GCF_016861695.1", "GTDB_GCF_005405165.1"}
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
data = os.path.join(root, "data")
out = os.path.join(root, "results", "08_tree")
os.makedirs(out, exist_ok=True)


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    h, s = None, []
    for line in op(path, "rt"):
        if line.startswith(">"):
            if h:
                yield h, "".join(s)
            h, s = line[1:].rstrip("\n"), []
        else:
            s.append(line.strip())
    if h:
        yield h, "".join(s)


tips = []  # (tip, type, family, genus, clade, label, seq)

want = set(re.findall(r"GTDB_[A-Z]{3}_[0-9.]+", open(f"{data}/refs/constraint_gtdb.tree").read())) - BAD
best = {}
for h, s in fasta(f"{data}/refs/ssu_all_r226.fna.gz"):
    tip = "GTDB_" + h.split("~")[0][3:]
    if tip in want and len(s) >= 1400 and len(s) > len(best.get(tip, ("", ""))[1]):
        best[tip] = (h.split(" ", 2)[1], s)
for tip, (lin, s) in sorted(best.items()):
    f, g, sp = (re.search(rf"{r}__([^;]*)", lin).group(1) for r in "fgs")
    tips.append((tip, "GTDB", f, g, "", sp, s))
print(f"GTDB: {len(best)} of {len(want)} accessions")

gtax = dict(l.rstrip("\n").split("\t") for l in open(f"{data}/ours/custom_taxonomy.tsv"))
labels = {r["id"]: r for r in csv.DictReader(open(f"{root}/results/06_sintax/refs/ours_labels.tsv"), delimiter="\t")}
seqs = dict((h.split(";")[0], s) for h, s in fasta(f"{data}/ours/custom_16S.fa"))
for sid, r in labels.items():
    if r["source"] != "genome":
        continue
    g = sid.split("__16S_")[0]
    if "f__Lactobacillaceae" not in gtax[g] or len(seqs[sid]) < 1200:
        continue
    clade = r["species_or_clade"] if r["species_or_clade"].startswith("Clado") else ""
    tips.append((f"OWN_{g}", "genome", "Lactobacillaceae", r["genus"], clade, g, seqs[sid]))
for tip, s in fasta(f"{data}/ours/mg16S_LAB.fa"):
    r = labels[tip]
    clade = r["species_or_clade"] if r["species_or_clade"].startswith("Clado") else ""
    tips.append((tip, "metagenome", "Lactobacillaceae", r["genus"], clade, tip.replace("MG16S_", ""), s))

with open(f"{out}/refs_full.fa", "w") as fo, open(f"{out}/tip_info.tsv", "w") as ti, open(f"{out}/taxon.tsv", "w") as tx:
    w = csv.writer(ti, delimiter="\t", lineterminator="\n")
    w.writerow(["tip", "type", "family", "genus", "clade", "label"])
    for tip, typ, f, g, clade, label, s in tips:
        fo.write(f">{tip}\n{s}\n")
        w.writerow([tip, typ, f, g, clade, label])
        tx.write(f"{tip}\t{f};{g or 'unresolved'}" + (f";{clade}" if clade else "") + "\n")

n = 0
with open(f"{out}/queries.fa", "w") as fo:
    for h, s in fasta(f"{root}/results/04_queries/favosa_LAB.fa"):
        fo.write(f">{h}\n{s}\n")
        n += 1
print("genomes:", sum(t[1] == "genome" for t in tips), " MG16S:", sum(t[1] == "metagenome" for t in tips), " ASV queries:", n)
