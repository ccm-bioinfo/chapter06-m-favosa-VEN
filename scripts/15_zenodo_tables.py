#!/usr/bin/env python3
"""Tables for the Zenodo deposit (results/15_zenodo/):
  genomes_used.tsv   genomes whose barrnap 16S entered the references (use=yes in results/03_genome16S_review.tsv;
                     every use=no genome, including UNU, is left out), with the copy used and whether it is in the tree
  recovered_16S.tsv  16S recovered from the Xolalpa 2026 metagenomes (phyloFlash; P16S general, T16S targeted to
                     ASV1/ASV2): strongest MarkerMAG link, the Xolalpa 2026 MAG of the same species by skani
                     (ANI >= 95 %) that names the clade, and the best link to a different clade when there is one.
                     Bin species labels from unpublished genomes are not included.
Usage (Mac, after 03, 05, 08): python3 scripts/15_zenodo_tables.py"""
import csv, os, re

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
out = os.path.join(root, "results", "15_zenodo")
os.makedirs(out, exist_ok=True)


def rows(path):
    with open(os.path.join(root, path)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


tips = {r["tip"] for r in rows("results/08_tree/tip_info.tsv")}
CLADE = {"1": "Clade 1 Nicoliella-related", "2": "Clade 2 Nicoliella-related",
         "3": "Clade 3-4 Acetilactobacillus-like", "4": "Clade 3-4 Acetilactobacillus-like"}

# ---- genomes ---------------------------------------------------------------------------------------------
SOURCE = {"own MAG": "Xolalpa 2026 MAG", "own isolate": "Xolalpa 2026 isolate", "public genome": "public isolate genome"}
with open(f"{out}/genomes_used.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome", "source", "record", "genus_GTDB-Tk", "species", "clade_Xolalpa2026", "16S_copy_used",
                "16S_length_bp", "in_SINTAX_reference", "in_tree", "not_in_tree_because"])
    for r in rows("results/03_genome16S_review.tsv"):
        if r["use"] != "yes":
            continue
        g = r["genome"]
        acc = re.search(r"(GC[AF])_(\d+)_(\d+)$", g)
        record = r["published"] if r["published"] != "public" else (f"{acc[1]}_{acc[2]}.{acc[3]}" if acc else "NCBI (public isolate)")
        in_tree = f"OWN_{g}" in tips
        why = "" if in_tree else "16S < 1200 bp" if int(r["length"]) < 1200 else "not Lactobacillaceae"  # tree rule of 08
        w.writerow([g, SOURCE[r["group"]], record, r["genome_genus"] or "-", r["genome_species"],
                    CLADE.get(r["clade_paper"], "-"), r["seq"], r["length"], "yes", "yes" if in_tree else "no", why])

# ---- recovered 16S ---------------------------------------------------------------------------------------
labels = {r["cluster"]: r for r in rows("docs/mg16S_labels.tsv")}
links = {}
for r in rows("data/ours/ligas_16S_MAG.tsv"):
    links.setdefault(r["16S"], []).append(r)
BAD_TREE = {"P16S_15", "P16S_62"}  # saturated terminal branch (08_tree_inputs.py)


def paper_clade(r):  # clade of the Xolalpa 2026 MAG that names the bin (1-4, before merging 3 and 4)
    m = re.match(r"Clado (\d)", r["clado"])
    return f"Clade {m.group(1)}" if m else ""


RELATION = {"es el MAG": "bin is the MAG", "misma especie, otro bin": "same species, other bin"}


with open(f"{out}/recovered_16S.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["16S", "assembly", "length_bp", "family_SINTAX", "genus_SINTAX", "samples_with_reads",
                "Lactobacillaceae_reference", "in_tree", "n_bins_linked", "strongest_bin", "strongest_pairs",
                "Xolalpa2026_MAG_same_species", "ANI_bin_vs_MAG", "relation", "paper_clade_strongest_link",
                "clade_assigned", "best_other_paper_clade", "best_other_clade_bin", "best_other_clade_pairs"])
    for r in rows("data/ours/mg16S_info.tsv"):
        c = r["cluster"]
        lk = sorted(links.get(c, []), key=lambda x: -int(x["pares"]))
        top = lk[0] if lk else None
        lab = labels.get(c)
        clade = lab["clade"].replace("_", " ") if lab and lab["clade"] else "-"
        top_clade = paper_clade(top) if top else ""
        other = next((x for x in lk if paper_clade(x) and paper_clade(x) != top_clade), None) if top_clade else None
        w.writerow([c, "targeted to ASV1/ASV2" if c.startswith("T16S") else "general", r["longitud"], r["familia"],
                    r["genero"], r["muestras_origen"], "yes" if lab else "no",
                    "no (saturated branch)" if c in BAD_TREE else ("yes" if lab else "no"),
                    len({x["bin"] for x in lk}), top["bin"] if top else "-", top["pares"] if top else "0",
                    (top["MAG_articulo"] or "-") if top else "-", (top["ANI"] or "-") if top else "-",
                    RELATION.get(top["relacion"], "no Xolalpa 2026 MAG >= 95 % ANI") if top else "no link",
                    top_clade or "-", clade, paper_clade(other) if other else "-", other["bin"] if other else "-",
                    other["pares"] if other else "-"])
print("written", out)
