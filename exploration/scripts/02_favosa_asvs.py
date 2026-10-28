#!/usr/bin/env python3
"""ASVs present in the only M. favosa sample (B-VEN-D9, 10MFH) with their previous
classifications and the recovered 16S / MAG clade they link to.
Usage: python3 exploration/scripts/02_favosa_asvs.py > results/favosa_ASVs.tsv  (Mac: needs Mac-only links)"""
import csv, os, sys

SAMPLE = "B-VEN-D9"
data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")


def rows(name):
    with open(os.path.join(data, name)) as fh:
        return list(csv.reader(fh, delimiter="\t"))


def last_rank(tax):
    ranks = [r.split("(")[0] for r in tax.replace(";", ",").split(",") if ":" in r]
    return ranks[-1].replace('"', "") if ranks else "No hit"


otu = rows("vit/otu_table_gtdbcustom_taxonomy.txt")
col = otu[0].index(SAMPLE)
counts = {r[0]: int(float(r[col])) for r in otu[1:] if float(r[col]) > 0}
prev = {r[0]: r[-1].split(";")[-1] for r in otu[1:]}
amptk = {r[0]: r[1].rstrip(";").split(";")[-1] for r in rows("vit/taxonomy_amptk.txt")[1:]}
gtdb = {r[0]: r[3] if len(r) > 3 else "" for r in rows("vit/sintax_gtdb.tsv")}  # cutoff-filtered
link = {r[0]: r for r in rows("ours/ASV_16S_recuperado_bin.tsv")[1:]}
total = sum(counts.values())

w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["ASV", "reads", "pct", "is_LAB", "amptk_last_rank", "gtdb_r226_sintax0.8",
            "previous_gtdb_custom", "recovered_16S", "linked_bin", "bin_gtdb", "clade_label"])
for asv, n in sorted(counts.items(), key=lambda x: -x[1]):
    lk = link.get(asv, [""] * 11)
    w.writerow([asv, n, f"{100 * n / total:.2f}", "Lactobacillaceae" in prev[asv],
                last_rank(amptk.get(asv, "")), last_rank(gtdb.get(asv, "")), last_rank(prev[asv]),
                lk[1], lk[3], lk[4], lk[10]])
