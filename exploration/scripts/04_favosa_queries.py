#!/usr/bin/env python3
"""Step 1: ASVs of the only M. favosa sample (B-VEN-D9, 10MFH), as query sequences.
LAB = AMPtk taxonomy contains f:Lactobacillaceae (Patricia's classification, independent of our 16S).
Usage: python3 exploration/scripts/04_favosa_queries.py   ->  results/04_queries/
  favosa_all.fa       every ASV with reads in B-VEN-D9 (SINTAX input)
  favosa_LAB.fa       Lactobacillaceae ASVs only (tree placement input)
  favosa_counts.tsv   ASV, reads, pct, is_LAB, AMPtk taxonomy"""
import csv, gzip, os

SAMPLE = "B-VEN-D9"
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
data = os.path.join(root, "data", "vit")
out = os.path.join(root, "results", "04_queries")
os.makedirs(out, exist_ok=True)

with open(os.path.join(data, "otu_table.txt")) as fh:
    otu = list(csv.reader(fh, delimiter="\t"))
col = otu[0].index(SAMPLE)
counts = {r[0]: int(float(r[col])) for r in otu[1:] if float(r[col]) > 0}
total = sum(counts.values())

with open(os.path.join(data, "taxonomy_amptk.txt")) as fh:
    amptk = {r[0]: r[1].split(";", 1)[-1].rstrip(";") for r in csv.reader(fh, delimiter="\t")}

seqs, sid = {}, None
for line in gzip.open(os.path.join(data, "ASVs.fa.gz"), "rt"):
    if line.startswith(">"):
        sid = line[1:].split(";")[0].strip()
        seqs[sid] = ""
    else:
        seqs[sid] += line.strip()

order = sorted(counts, key=lambda a: -counts[a])
is_lab = {a: "f:Lactobacillaceae" in amptk[a] for a in order}

with open(os.path.join(out, "favosa_all.fa"), "w") as fa, open(os.path.join(out, "favosa_LAB.fa"), "w") as fl:
    for a in order:
        rec = f">{a}\n{seqs[a]}\n"
        fa.write(rec)
        if is_lab[a]:
            fl.write(rec)

with open(os.path.join(out, "favosa_counts.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV", "reads", "pct", "is_LAB", "amptk_taxonomy"])
    for a in order:
        w.writerow([a, counts[a], f"{100 * counts[a] / total:.2f}", is_lab[a], amptk[a]])

lab = [a for a in order if is_lab[a]]
print(f"{SAMPLE}: {len(order)} ASVs, {total} reads")
print(f"Lactobacillaceae (AMPtk): {len(lab)} ASVs, {sum(counts[a] for a in lab)} reads "
      f"({100 * sum(counts[a] for a in lab) / total:.1f}%)")
print(f"lengths: {min(len(seqs[a]) for a in order)}-{max(len(seqs[a]) for a in order)} bp")
