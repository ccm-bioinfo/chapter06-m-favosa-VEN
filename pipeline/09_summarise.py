#!/usr/bin/env python3
# 09_summarise.py - Combine labels with gappa placement, compute % reads per label per classification,
# Lactobacillaceae shares, SINTAX-vs-tree agreement, sensitivity changes, best-hit identities and the V4
# pairwise differences (substitutions + internal indel columns) among all sequences of the placement tree.
import csv, sys, re, os
from collections import defaultdict
D = os.environ["CH06_WORK"] + "/"
O = D + "results/"
GT = D + "results/04_besthit_gtdb_noN.tsv"

lab = list(csv.DictReader(open(D + "results/07_asv_labels.tsv"), delimiter="\t"))
pq = defaultdict(list)
for r in csv.DictReader(open(D + "place/assign/per_query.tsv"), delimiter="\t"):
    pq[r["name"]].append(r)
place = {}
for a, rows in pq.items():
    best = max(rows, key=lambda r: (float(r["LWR"]), len(r["taxopath"].split(";"))))
    place[a] = (best["taxopath"].split(";")[-1].replace("Lactobacillaceae", "Lactobacillaceae (no genus)"),
                best["LWR"])
for r in lab:
    r["tree"] = place[r["ASV"]][0] if r["ASV"] in place else r["GTDB+ours"]
    r["placed"] = "yes" if r["ASV"] in place else "no"
    r["tree_LWR"] = place[r["ASV"]][1] if r["ASV"] in place else ""

def best(p):
    b = {}
    for l in open(p):
        x = l.rstrip("\n").split("\t")
        if x[0] not in b or float(x[2]) > float(b[x[0]][1]):
            b[x[0]] = (x[1].split(";")[0], x[2])
    return b
bo, bg = best(D + "results/03_besthit_ours.tsv"), best(GT)
for r in lab:
    a = r["ASV"]
    r["best_ours"], r["id_ours"] = bo.get(a, ("", ""))
    r["best_gtdb"], r["id_gtdb"] = bg.get(a, ("", ""))

cols = ["RDP", "SILVA", "GTDB", "GTDB+ours", "tree", "final"]
with open(O + "09_asv_labels_full.tsv", "w") as o:
    keys = ["ASV", "reads", "amptk_lacto"] + cols + ["GTDB+ours_sens", "placed", "tree_LWR",
            "best_ours", "id_ours", "best_gtdb", "id_gtdb"]
    o.write("\t".join(keys) + "\n")
    for r in lab: o.write("\t".join(r[k] for k in keys) + "\n")

tot = sum(int(r["reads"]) for r in lab)
lac = [r for r in lab if r["amptk_lacto"] == "yes"]; ltot = sum(int(r["reads"]) for r in lac)
print(f"total reads {tot}; Lactobacillaceae (AMPtk) ASVs {len(lac)} reads {ltot}; placed {len(place)}")
with open(O + "09_percent_reads_per_label.tsv", "w") as o:
    o.write("classification\tset\tlabel\treads\tpercent\n")
    for c in cols:
        for setname, rows, T in (("all", lab, tot), ("Lactobacillaceae", lac, ltot)):
            d = defaultdict(int)
            for r in rows: d[r[c]] += int(r["reads"])
            for k, v in sorted(d.items(), key=lambda x: -x[1]):
                o.write(f"{c}\t{setname}\t{k}\t{v}\t{100*v/T:.2f}\n")
agree = sum(int(r["reads"]) for r in lac if r["tree"] == r["GTDB+ours"])
print(f"SINTAX(GTDB+ours) vs tree agreement: {100*agree/ltot:.2f}% of Lactobacillaceae reads")
ch = [r["ASV"] + f"({r['GTDB+ours']}->{r['GTDB+ours_sens']},{r['reads']})" for r in lab if r["GTDB+ours"] != r["GTDB+ours_sens"]]
print("sensitivity changes:", len(ch), ch)
org = [r for r in lab if r["final"].startswith("Organelles")]
print(f"organelle 16S (SILVA): {len(org)} ASVs, {sum(int(r['reads']) for r in org)} reads; "
      f"final = tree for every Lactobacillaceae ASV placed: {all(r['final'] == r['GTDB+ours'] for r in lab if r['placed'] == 'yes')}")
for nm, k in (("ours", "id_ours"), ("GTDB", "id_gtdb")):
    v = sum(int(r["reads"]) for r in lac if r[k] and float(r[k]) >= 99)
    print(f"best hit >=99% to {nm}: {100*v/ltot:.2f}% of Lactobacillaceae reads")
print("placed but not AMPtk-Lactobacillaceae:", [a for a in place if a not in {r['ASV'] for r in lac}])
print("AMPtk-Lactobacillaceae not placed:", [(r["ASV"], r["GTDB+ours"]) for r in lac if r["placed"] == "no"])
print("\nTop 9 Lactobacillaceae ASVs:")
for r in sorted(lac, key=lambda r: -int(r["reads"]))[:9]:
    print("\t".join([r["ASV"], r["reads"]] + [r[c] for c in cols] + [r["best_ours"], r["id_ours"], r["best_gtdb"].split("~")[0], r["id_gtdb"]]))

# ---- V4 differences ----
seqs, h = {}, None
for l in open(D + "place/ref_plus_asv.fa"):
    l = l.strip()
    if l.startswith(">"): h = l[1:].split()[0]; seqs[h] = []
    else: seqs[h].append(l.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
def span(s):
    idx = [i for i, c in enumerate(s) if c not in "-."]
    return idx[0], idx[-1]
asvs = [k for k in seqs if k.startswith("ASV")]
v0 = min(span(seqs[a])[0] for a in asvs); v1 = max(span(seqs[a])[1] for a in asvs)
print(f"\nV4 window: alignment columns {v0+1}-{v1+1} of {len(next(iter(seqs.values())))}")
V = {k: s[v0:v1+1] for k, s in seqs.items()}
def diff(a, b):
    sa, sb = span(a) if a.strip("-.") else (0, -1), span(b) if b.strip("-.") else (0, -1)
    lo, hi = max(sa[0], sb[0]), min(sa[1], sb[1]); n = 0
    for i in range(lo, hi + 1):
        x, y = a[i], b[i]
        if x in "-." and y in "-.": continue
        if x == "N" or y == "N": continue
        if x != y: n += 1
    return n
names = list(V)
with open(O + "09_v4_differences.tsv", "w") as o:
    o.write("seq\t" + "\t".join(names) + "\n")
    for a in names: o.write(a + "\t" + "\t".join(str(diff(V[a], V[b])) for b in names) + "\n")
refs = [n for n in names if not n.startswith("ASV")]
for a in ("ASV1", "ASV2"):
    d = sorted((diff(V[a], V[r]), r) for r in refs)
    rec = [x for x in d if re.match(r"[PT]16S_", x[1])]
    print(a, "closest ref:", d[:3], "closest recovered:", rec[:3])
print("ASV1 vs ASV2:", diff(V["ASV1"], V["ASV2"]))
print("T16S_1 vs T16S_2:", diff(V["T16S_1"], V["T16S_2"]))
