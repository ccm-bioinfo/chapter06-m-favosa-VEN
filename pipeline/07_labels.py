#!/usr/bin/env python3
# 07_labels.py - Per-ASV labels for RDP (AMPtk, from the package), SINTAX SILVA, GTDB, GTDB+ours and the
# sensitivity run (clade if any at species rank, else genus, else "Lactobacillaceae (no genus)", else
# "Unclassified"); writes the FASTA of ASVs to place (GTDB+ours = Nicoliella/Acetilactobacillus/clade or
# Lactobacillaceae without genus).
import csv, os
W = os.environ["CH06_WORK"] + "/"; I = os.environ["CH06_INPUTS"] + "/"
NOGEN = ("uncultured", "unknown", "unidentified", "metagenome", "incertae")

def lab_from(ranks):
    s, g, f = ranks.get("s", ""), ranks.get("g", ""), ranks.get("f", "")
    if s.startswith("Clade"): return s.replace("_", " ")
    if g and not g.lower().startswith(NOGEN): return g
    if f == "Lactobacillaceae": return "Lactobacillaceae (no genus)"
    return "Unclassified"

def sintax(p):
    out = {}
    for l in open(p):
        x = l.rstrip("\n").split("\t")
        r = dict(t.split(":", 1) for t in x[3].split(",") if ":" in t) if len(x) > 3 and x[3] else {}
        full = x[1] if len(x) > 1 else ""
        out[x[0]] = (lab_from(r), full)
    return out

tab = list(csv.DictReader(open(I + "asv/B-VEN-D9_ASV_table.tsv"), delimiter="\t"))
runs = {k: sintax(W + f"results/03_sintax_{v}.tsv") for k, v in
        [("SILVA", "silva138.2_nr99"), ("GTDB", "gtdb_r226_ge1200"), ("GTDB+ours", "gtdb_ours"),
         ("GTDB+ours_sens", "gtdb_ours_sens")]}
asv = {}; h = None
for l in open(I + "asv/B-VEN-D9_ASVs.fa"):
    l = l.strip()
    if l.startswith(">"): h = l[1:].split()[0].split(";")[0]; asv[h] = ""
    else: asv[h] += l
with open(W + "results/07_asv_labels.tsv", "w") as o, open(W + "work/place_query.fa", "w") as q:
    o.write("ASV\treads\tamptk_lacto\tRDP\tSILVA\tGTDB\tGTDB+ours\tGTDB+ours_sens\tsintax_gtdb_ours_full\n")
    nq = 0
    for r in tab:
        t = r["amptk_taxonomy"].replace('"', "")
        rr = dict(x.split(":", 1) for x in t.split(",") if ":" in x)
        rdp = lab_from({"g": rr.get("g", ""), "f": rr.get("f", "")})
        a = r["ASV"]
        L = [runs[k][a][0] for k in ("SILVA", "GTDB", "GTDB+ours", "GTDB+ours_sens")]
        o.write(f"{a}\t{r['reads']}\t{'yes' if rr.get('f') == 'Lactobacillaceae' else 'no'}\t{rdp}\t" +
                "\t".join(L) + f"\t{runs['GTDB+ours'][a][1]}\n")
        if L[2] in ("Nicoliella", "Acetilactobacillus", "Lactobacillaceae (no genus)") or L[2].startswith("Clade"):
            q.write(f">{a}\n{asv[a]}\n"); nq += 1
print("ASVs to place:", nq)
