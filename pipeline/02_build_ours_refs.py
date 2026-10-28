#!/usr/bin/env python3
# 02_build_ours_refs.py - Build "our 16S" SINTAX references (14 genome 16S of Table S1 + 23 recovered
# Lactobacillaceae 16S of Table S2), the GTDB+ours DB, the sensitivity DB (without MAG_41_Melli_001 and
# MAG_45_Melli_020) and the Lactobacillaceae ASV query file (AMPtk family Lactobacillaceae, n=69).
import csv, os, re
W = os.environ["CH06_WORK"] + "/"; I = os.environ["CH06_INPUTS"] + "/"
HI = "d:Bacteria,p:Bacillota,c:Bacilli,o:Lactobacillales,f:Lactobacillaceae"
CLADE_GENUS = {"Clade_2": "Nicoliella", "Clade_1": "Nicoliella", "Clade_3-4": "Acetilactobacillus"}

def fasta(p):
    d, h = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): h = l[1:].split()[0]; d[h] = []
        elif h: d[h].append(l)
    return {k: "".join(v).upper() for k, v in d.items()}

def clean(x): return re.sub(r"[,;:\s]+", "_", x.strip())

def clade_lab(c):
    if c in ("-", ""): return None
    if c.startswith("Clade 3") or c.startswith("Clade 4"): return "Clade_3-4"
    return "Clade_" + c.split()[1]

gen16 = fasta(I + "genome_16S/genome_16S_barrnap.fa")
rec16 = fasta(I + "recovered_16S/recovered_16S.fa")
ours = []  # (id, seq, tax, origin)
for r in csv.DictReader(open(I + "tables/Table_S1_genomes.tsv"), delimiter="\t"):
    if r["source"].startswith("GTDB") or r["16S_used"] in ("-", ""): continue
    cl = clade_lab(r["paper_clade_Xolalpa2026"])
    g = CLADE_GENUS[cl] if cl else r["genus"]
    sp = cl if cl else clean(r["species"])
    ours.append((r["genome"], gen16[r["16S_used"]], f"{HI},g:{g},s:{sp}", "genome"))
for r in csv.DictReader(open(I + "tables/Table_S2_recovered_16S.tsv"), delimiter="\t"):
    if r["Lactobacillaceae_reference"] != "yes": continue
    cl = clade_lab(r["clade_assigned"])
    tax = HI
    if cl: tax += f",g:{CLADE_GENUS[cl]},s:{cl}"
    elif r["genus"] not in ("-", ""): tax += f",g:{r['genus']}"
    ours.append((r["16S"], rec16[r["16S"]], tax, "recovered"))

with open(W + "refs/ours.sintax.fa", "w") as o, open(W + "results/02_ours_reference_taxonomy.tsv", "w") as t:
    t.write("id\torigin\tlength\ttax\n")
    for i, s, tax, og in ours:
        o.write(f">{i};tax={tax};\n{s}\n"); t.write(f"{i}\t{og}\t{len(s)}\t{tax}\n")
with open(W + "refs/ours_sens.sintax.fa", "w") as o:
    for i, s, tax, og in ours:
        if i in ("MAG_41_Melli_001", "MAG_45_Melli_020"): continue
        o.write(f">{i};tax={tax};\n{s}\n")
print("ours references:", len(ours), "genomes:", sum(x[3] == "genome" for x in ours),
      "recovered:", sum(x[3] == "recovered" for x in ours))

asv = fasta(I + "asv/B-VEN-D9_ASVs.fa")
lac = [r["ASV"] for r in csv.DictReader(open(I + "asv/B-VEN-D9_ASV_table.tsv"), delimiter="\t")
       if "f:Lactobacillaceae" in r["amptk_taxonomy"]]
with open(W + "work/lacto_asvs.fa", "w") as o:
    for a in lac: o.write(f">{a}\n{asv[a]}\n")
print("Lactobacillaceae ASVs (AMPtk):", len(lac))
