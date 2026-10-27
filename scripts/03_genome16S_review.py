#!/usr/bin/env python3
"""One row per barrnap 16S from the custom genomes, for manual review before building references.
Decisions per genome (published, paper clade, use, note) come from docs/genome16S_decisions.tsv
and override the default (use = longest validated copy).
Usage (alnitak): python3 scripts/03_genome16S_review.py > results/03_genome16S_review.tsv"""
import csv, os, re, sys

data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "ours")


def table(name):
    with open(os.path.join(data, name)) as fh:
        return list(csv.reader(fh, delimiter="\t"))


def fasta(name):
    seqs, sid = {}, None
    for line in open(os.path.join(data, name)):
        if line.startswith(">"):
            sid = line[1:].strip().split(";")[0]
            seqs[sid] = ""
        else:
            seqs[sid] += line.strip()
    return seqs


def group(genome):
    if genome.startswith("UNU_"):
        return "UNU (excluded: no cooperation letter)"
    if genome.startswith("MAG_"):
        return "own MAG"
    if re.match(r"(N_)?barcode\d+", genome):
        return "own isolate"
    if genome.startswith("CR_"):
        return "Costa Rica MAG"
    return "public genome"


def rank(tax, r):
    m = re.search(rf"{r}[:_]{{1,2}}([^,;(]+)", tax)
    return m.group(1).replace("_", " ") if m else ""


seqs = fasta("custom_16S.fa")
clean = set(fasta("custom_16S.clean.noUNU.sintax.fa"))
gtax = {r[0]: r[1] for r in table("custom_taxonomy.tsv")}
check = {r[0]: r[5] for r in table("custom_16S_check.tsv")[1:]}
root = os.path.join(data, "..", "..")
with open(os.path.join(root, "docs", "genome16S_decisions.tsv")) as fh:
    decision = {r["genome"]: r for r in csv.DictReader(fh, delimiter="\t")}
sintax = {r[0].split(";")[0]: (r[3] if len(r) > 3 else "") for r in table("custom_16S_vs_gtdb.tsv")}

# longest validated copy per genome = the one the tree uses; sorted so that ties (copies of equal length)
# always go to the lowest copy number (iterating the set directly changed the pick between runs)
longest = {}
for sid in sorted(clean):
    g = sid.split("__16S_")[0]
    if g not in longest or len(seqs[sid]) > len(seqs[longest[g]]):
        longest[g] = sid

w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["seq", "genome", "group", "length", "longest_valid_copy", "genome_genus", "genome_species",
            "16S_sintax_genus", "check", "in_clean_noUNU", "flag", "published", "clade_paper", "use", "note"])
for sid in sorted(seqs, key=lambda x: (group(x.split("__16S_")[0]), x)):
    g = sid.split("__16S_")[0]
    t = gtax.get(g, "")
    is_ref = longest.get(g) == sid
    flag = "short (<1000 bp)" if is_ref and len(seqs[sid]) < 1000 else ""
    d = decision.get(g, {})
    pub = d.get("published") or ("public" if group(g) == "public genome" else "no" if g.startswith("UNU_") else "")
    use = "yes" if is_ref else "no"
    if is_ref and d.get("use"):
        use = d["use"]
    w.writerow([sid, g, group(g), len(seqs[sid]), is_ref, rank(t, "g"), rank(t, "s"),
                rank(sintax.get(sid, ""), "g"), check.get(sid, "not checked"), sid in clean, flag,
                pub, d.get("clade_paper", ""), use, d.get("note", "") if is_ref or g not in longest else ""])
