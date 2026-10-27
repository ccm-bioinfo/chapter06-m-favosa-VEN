#!/usr/bin/env python3
"""Tables for the Zenodo deposit (zenodo/tables/):
  genomes_tree.tsv   reference genomes of the Nicoliella-Acetilactobacillus subtree (Fig. 4/4b, rooted on
                     Fructilactobacillus): GTDB r226 representatives, Xolalpa 2026 MAGs and isolates and public genomes
                     (docs/public_genomes.tsv), with the 16S used and its length. ASVs and recovered 16S are not genomes
                     and are left out (the recovered 16S are in recovered_16S.tsv)
  recovered_16S.tsv  16S recovered from the Xolalpa 2026 metagenomes (phyloFlash; P16S general, T16S targeted to
                     ASV1/ASV2): strongest MarkerMAG link, the Xolalpa 2026 MAG of the same species by skani
                     (ANI >= 95 %) that names the clade, and the best link to a different paper clade when there is one.
                     Bin species labels from unpublished genomes are not included.
Usage (Mac, after 12_figures.R): python3 scripts/15_zenodo_tables.py"""
import csv, os, re

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
out = os.path.join(root, "zenodo", "tables")
os.makedirs(out, exist_ok=True)


def rows(path):
    with open(os.path.join(root, path)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fasta_len(path):
    n, h = {}, None
    for line in open(os.path.join(root, path)):
        if line.startswith(">"):
            h = line[1:].split()[0]
            n[h] = 0
        else:
            n[h] += len(line.strip())
    return n


CLADE = {"1": "Clade 1 Nicoliella-related", "2": "Clade 2 Nicoliella-related",
         "3": "Clade 3-4 Acetilactobacillus-like", "4": "Clade 3-4 Acetilactobacillus-like"}

# ---- genomes of the Fig. 4 subtree -----------------------------------------------------------------------
info = {r["tip"]: r for r in rows("results/08_tree/tip_info.tsv")}
fig4 = [t.strip() for t in open(os.path.join(root, "results/12_figures/fig4_tips.txt")) if t.strip()]
length = fasta_len("results/08_tree/refs_full.fa")
review = {r["genome"]: r for r in rows("results/03_genome16S_review.tsv") if r["use"] == "yes"}
public = {r["genome"]: r for r in rows("docs/public_genomes.tsv")}
outgroup = {t for t in fig4 if info.get(t, {}).get("genus") == "Fructilactobacillus"}
with open(f"{out}/genomes_tree.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["tip", "source", "genome_type", "record", "genus", "species", "clade_Xolalpa2026",
                "16S_used", "16S_length_bp", "role"])
    for t in fig4:
        r = info.get(t)
        if not r or r["type"] not in ("GTDB", "genome"):
            continue
        role = "outgroup" if t in outgroup else "reference"
        if r["type"] == "GTDB":
            w.writerow([t, "GTDB r226 representative", "-", t.replace("GTDB_", ""), r["genus"], r["label"], "-",
                        "longest GTDB r226 SSU >= 1400 bp", length[t], role])
            continue
        g = r["label"]
        v = review[g]
        if v["group"] == "public genome":
            p = public[g]
            src, gtype, rec = "public genome", p["type"], p["assembly"]
        else:
            src = "Xolalpa 2026"
            gtype = "MAG" if v["group"] == "own MAG" else "isolate"
            rec = v["published"]
        w.writerow([t, src, gtype, rec, r["genus"], v["genome_species"], CLADE.get(v["clade_paper"], "-"),
                    v["seq"], length[t], role])

# ---- recovered 16S ---------------------------------------------------------------------------------------
labels = {r["cluster"]: r for r in rows("docs/mg16S_labels.tsv")}
notes = {r["16S"]: r["note"] for r in rows("docs/mg16S_notes.tsv")}
links = {}
for r in rows("data/ours/ligas_16S_MAG.tsv"):
    links.setdefault(r["16S"], []).append(r)


def paper_clade(r):  # clade of the Xolalpa 2026 MAG that names the bin (1-4, before merging 3 and 4)
    m = re.match(r"Clado (\d)", r["clado"])
    return f"Clade {m.group(1)}" if m else ""


RELATION = {"es el MAG": "bin is the MAG", "misma especie, otro bin": "same species, other bin"}
with open(f"{out}/recovered_16S.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["16S", "assembly", "length_bp", "family_SINTAX", "genus_SINTAX", "samples_with_reads",
                "Lactobacillaceae_reference", "in_tree", "n_bins_linked", "strongest_bin", "strongest_pairs",
                "Xolalpa2026_MAG_same_species", "ANI_bin_vs_MAG", "relation", "paper_clade_strongest_link",
                "clade_assigned", "best_other_paper_clade", "best_other_clade_bin", "best_other_clade_pairs", "note"])
    for r in rows("data/ours/mg16S_info.tsv"):
        c = r["cluster"]
        lk = sorted(links.get(c, []), key=lambda x: -int(x["pares"]))
        top = lk[0] if lk else None
        lab = labels.get(c)
        clade = lab["clade"].replace("_", " ") if lab and lab["clade"] else "-"
        top_clade = paper_clade(top) if top else ""
        other = next((x for x in lk if paper_clade(x) and paper_clade(x) != top_clade), None) if top_clade else None
        in_tree = "yes" if r["tip"] in info else ("no (saturated branch)" if lab else "no")
        w.writerow([c, "targeted to ASV1/ASV2" if c.startswith("T16S") else "general", r["longitud"], r["familia"],
                    r["genero"], r["muestras_origen"], "yes" if lab else "no", in_tree,
                    len({x["bin"] for x in lk}), top["bin"] if top else "-", top["pares"] if top else "0",
                    (top["MAG_articulo"] or "-") if top else "-", (top["ANI"] or "-") if top else "-",
                    RELATION.get(top["relacion"], "no Xolalpa 2026 MAG >= 95 % ANI") if top else "no link",
                    top_clade or "-", clade, paper_clade(other) if other else "-",
                    other["bin"] if other else "-", other["pares"] if other else "-", notes.get(c, "")])
print("written", out)
