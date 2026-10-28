#!/usr/bin/env python3
"""Tables for the Zenodo deposit (zenodo/tables/):
  genomes_tree.tsv   genomes of the Nicoliella-Acetilactobacillus subtree (Fig. 4/4b, rooted on Fructilactobacillus):
                     GTDB r226 representatives (with the public strain they correspond to) and public genomes
                     (docs/public_genomes.tsv), plus every Xolalpa 2026 clade genome (docs/xolalpa2026_genomes.tsv),
                     including those without a usable 16S, with barrnap copies, 16S used and MarkerMAG links received
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

# ---- genomes: the Fig. 4 subtree plus every Xolalpa 2026 clade genome ----------------------------------------
# Xolalpa 2026 genomes come from docs/xolalpa2026_genomes.tsv (in_table = yes): all were searched with barrnap and
# were MarkerMAG link targets, with or without a 16S of their own
info = {r["tip"]: r for r in rows("results/08_tree/tip_info.tsv")}
fig4 = [t.strip() for t in open(os.path.join(root, "results/12_figures/fig4_tips.txt")) if t.strip()]
length = fasta_len("results/08_tree/refs_full.fa")
all16S = rows("results/03_genome16S_review.tsv")
review = {r["genome"]: r for r in all16S if r["use"] == "yes"}
copies = {}
for r in all16S:
    copies[r["genome"]] = copies.get(r["genome"], 0) + 1
public = {r["genome"]: r for r in rows("docs/public_genomes.tsv")}
strain_of_gtdb = {r["same_as_GTDB_tip"]: f"{r['strain']} ({r['genome']})" for r in public.values() if r["same_as_GTDB_tip"] != "-"}
xol = {r["genome"]: r for r in rows("docs/xolalpa2026_genomes.tsv") if r["in_table"] == "yes"}
mm = {}
for r in rows("data/ours/ligas_16S_MAG.tsv"):
    if r["MAG_articulo"]:
        k = (r["MAG_articulo"], r["16S"])
        mm[k] = max(mm.get(k, 0), int(r["pares"]))
linked = {}
for (g, c), n in sorted(mm.items(), key=lambda x: -x[1]):
    linked.setdefault(g, []).append(f"{c} ({n})")
outgroup = {t for t in fig4 if info.get(t, {}).get("genus") == "Fructilactobacillus"}
with open(f"{out}/genomes_tree.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome", "source", "genome_type", "record", "genus", "species", "same_strain_as",
                "paper_clade_Xolalpa2026", "barrnap_16S_copies", "16S_used", "16S_length_bp", "MarkerMAG_16S_linked",
                "in_tree", "role", "note"])
    done = set()
    for t in fig4:
        r = info.get(t)
        if not r or r["type"] not in ("GTDB", "genome"):
            continue
        role = "outgroup" if t in outgroup else "reference"
        if r["type"] == "GTDB":
            acc = t.replace("GTDB_", "")
            w.writerow([acc, "GTDB r226 representative", "-", acc, r["genus"], r["label"], strain_of_gtdb.get(acc, "-"),
                        "-", "-", "longest GTDB r226 SSU >= 1400 bp", length[t], "-", "yes", role, ""])
            continue
        g = r["label"]
        v = review[g]
        if g in xol:
            continue  # written below with the Xolalpa 2026 genomes
        p = public[g]
        w.writerow([g, "public genome", p["type"], p["assembly"], r["genus"], v["genome_species"], "-", "-",
                    copies.get(g, 0), v["seq"], length[t], "-", "yes", role, p["note"]])
    for g, x in xol.items():
        v = review.get(g)
        t = f"OWN_{g}"
        rec = "Zenodo 10.5281/zenodo.18446951" if x["type"] == "MAG" else (v["published"] if v else "PRJNA1235305")
        w.writerow([g, "Xolalpa 2026", x["type"], rec, x["GTDB_genus_paper"], v["genome_species"] if v else "-", "-",
                    f"Clade {x['clade_paper']}", copies.get(g, 0), v["seq"] if v else "-",
                    length.get(t, "-") if t in info else (v["length"] if v else "-"), "; ".join(linked.get(g, [])) or "-",
                    "yes" if t in info else "no", "reference" if t in info else "searched, no usable 16S", x["decision"]])

# ---- recovered 16S ---------------------------------------------------------------------------------------
labels = {r["cluster"]: r for r in rows("docs/mg16S_labels.tsv")}
notes = {r["16S"]: r["note"] for r in rows("docs/mg16S_notes.tsv")}
links = {}
for r in rows("data/ours/ligas_16S_MAG.tsv"):
    links.setdefault(r["16S"], []).append(r)


def paper_clade(r):  # clade of the Xolalpa 2026 MAG that names the bin (1-4, before merging 3 and 4)
    m = re.match(r"Clado (\d)", r["clado"])
    return f"Clade {m.group(1)}" if m else ""


CLADE_GENUS = {v: ("Nicoliella" if "Nicoliella" in v else "Acetilactobacillus") for v in CLADE.values()}  # genus of the clade, else SINTAX genus
RELATION = {"es el MAG": "bin is the MAG", "misma especie, otro bin": "same species, other bin"}
with open(f"{out}/recovered_16S.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["16S", "assembly", "length_bp", "family_SINTAX", "genus_SINTAX", "samples_with_reads",
                "Lactobacillaceae_reference", "in_tree", "n_bins_linked", "strongest_bin", "strongest_pairs",
                "Xolalpa2026_MAG_same_species", "ANI_bin_vs_MAG", "relation", "paper_clade_strongest_link",
                "clade_assigned", "genus", "best_other_paper_clade", "best_other_clade_bin", "best_other_clade_pairs", "note"])
    for r in rows("data/ours/mg16S_info.tsv"):
        c = r["cluster"]
        lk = sorted(links.get(c, []), key=lambda x: -int(x["pares"]))
        top = lk[0] if lk else None
        lab = labels.get(c)
        clade = lab["clade"].replace("_", " ") if lab and lab["clade"] else "-"
        top_clade = paper_clade(top) if top else ""
        other = next((x for x in lk if paper_clade(x) and paper_clade(x) != top_clade), None) if top_clade else None
        genus = CLADE_GENUS.get(clade, r["genero"])
        # tree of the methods: recovered Lactobacillaceae 16S whose genus is Nicoliella or Acetilactobacillus
        in_tree = "yes" if lab and genus in ("Nicoliella", "Acetilactobacillus") else "no"
        w.writerow([c, "targeted to ASV1/ASV2" if c.startswith("T16S") else "general", r["longitud"], r["familia"],
                    r["genero"], r["muestras_origen"], "yes" if lab else "no", in_tree,
                    len({x["bin"] for x in lk}), top["bin"] if top else "-", top["pares"] if top else "0",
                    (top["MAG_articulo"] or "-") if top else "-", (top["ANI"] or "-") if top else "-",
                    RELATION.get(top["relacion"], "no Xolalpa 2026 MAG >= 95 % ANI") if top else "no link",
                    top_clade or "-", clade, CLADE_GENUS.get(clade, r["genero"]), paper_clade(other) if other else "-",
                    other["bin"] if other else "-", other["pares"] if other else "-", notes.get(c, "")])
print("written", out)
