#!/usr/bin/env python3
"""Step 4: compare the classifications of the B-VEN-D9 Lactobacillaceae ASVs.
Methods: RDP (AMPtk, submitted chapter), SILVA 138.2, GTDB r226, GTDB + our 16S (SINTAX 0.8), the
sensitivity run without the MAG_41/45 16S, and EPA-ng placement on the constrained tree (gappa best hit).
Writes results/11_compare/:
  favosa_LAB_classification.tsv  one row per ASV (label = deepest resolved rank; conf = SINTAX confidence)
  summary_by_method.tsv          % of LAB reads per method and label (input of the figures)
Usage: python3 scripts/11_compare.py   (Mac or alnitak, after steps 2 and 3)"""
import csv, os, re

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
res = os.path.join(root, "results")
out = os.path.join(res, "11_compare")
os.makedirs(out, exist_ok=True)
METHODS = ["RDP (AMPtk)", "SILVA 138.2", "GTDB r226", "GTDB + our 16S", "Tree placement (EPA-ng)"]


def tsv(path):
    with open(path) as fh:
        return list(csv.reader(fh, delimiter="\t"))


def label_from(ranks):
    """ranks: list of (rank, name) to the deepest resolved level -> display label."""
    names = dict(ranks)
    if names.get("s", "").startswith("Clado"):
        return names["s"].replace("_", " ", 2).replace("_", " ")
    if "g" in names:
        return names["g"]
    if names.get("f") == "Lactobacillaceae":
        return "Lactobacillaceae (family only)"
    return "Not resolved to Lactobacillaceae"


def sintax(path):
    d = {}
    for r in tsv(path):
        cut = r[3] if len(r) > 3 else ""
        ranks = [tuple(x.split(":", 1)) for x in cut.split(",") if ":" in x]
        probs = dict((m.group(1), float(m.group(2))) for m in re.finditer(r"(\w):[^(,]+\(([\d.]+)\)", r[1]))
        conf = probs.get(ranks[-1][0], "") if ranks else ""
        d[r[0]] = (label_from(ranks), ",".join(f"{k}:{v}" for k, v in ranks), conf)
    return d


counts = {r[0]: (int(r[1]), r[3] == "True", r[4]) for r in tsv(f"{res}/04_queries/favosa_counts.tsv")[1:]}
lab = sorted((a for a in counts if counts[a][1]), key=lambda a: -counts[a][0])
total = sum(counts[a][0] for a in lab)

rdp = {}
for a in lab:
    ranks = [tuple(x.split(":", 1)) for x in counts[a][2].replace('"', "").split(",") if ":" in x]
    lbl = label_from(ranks)
    rdp[a] = (lbl + " (RDP)" if lbl == "Lactobacillus" else lbl, counts[a][2], "")
runs = {k: sintax(f"{res}/06_sintax/{k}.tsv") for k in ("silva", "gtdb", "gtdb_ours", "gtdb_ours_sens")}

epa = {}
pq = f"{res}/08_tree/assign/per_query.tsv"
if os.path.exists(pq):
    for r in tsv(pq)[1:]:
        path = r[5].split(";")
        ranks = [("f", path[0])] + ([("g", path[1])] if len(path) > 1 and path[1] != "unresolved" else []) \
                + ([("s", path[2])] if len(path) > 2 else [])
        epa[r[0]] = (label_from(ranks), r[5], r[1])

alwr = {}  # accumulated LWR of the best-hit taxopath (sum over all placements)
pa = f"{res}/08_tree/assign_all/per_query.tsv"
if os.path.exists(pa):
    for r in tsv(pa)[1:]:
        alwr[(r[0], r[5])] = r[3]

per_method = {"RDP (AMPtk)": rdp, "SILVA 138.2": runs["silva"], "GTDB r226": runs["gtdb"],
              "GTDB + our 16S": runs["gtdb_ours"], "Tree placement (EPA-ng)": epa}

with open(f"{out}/favosa_LAB_classification.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["ASV", "reads", "pct_LAB_reads", "RDP_AMPtk", "SILVA", "SILVA_conf", "GTDB", "GTDB_conf",
                "GTDB_ours", "GTDB_ours_conf", "GTDB_ours_sens", "same_as_sens", "EPA_taxopath", "EPA_LWR",
                "EPA_aLWR", "SINTAX_vs_EPA"])
    for a in lab:
        s, g, o, ss = runs["silva"][a], runs["gtdb"][a], runs["gtdb_ours"][a], runs["gtdb_ours_sens"][a]
        e = epa.get(a, ("", "", ""))
        w.writerow([a, counts[a][0], f"{100 * counts[a][0] / total:.2f}", rdp[a][0], s[0], s[2], g[0], g[2],
                    o[0], o[2], ss[0], o[0] == ss[0], e[1], e[2], alwr.get((a, e[1]), ""),
                    "same" if o[0] == e[0] else "differ"])

with open(f"{out}/summary_by_method.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["method", "label", "n_ASVs", "reads", "pct_LAB_reads"])
    for m in METHODS:
        if not per_method[m]:
            continue
        agg = {}
        for a in lab:
            agg.setdefault(per_method[m][a][0], []).append(a)
        for lbl, asvs in sorted(agg.items(), key=lambda x: -sum(counts[a][0] for a in x[1])):
            n = sum(counts[a][0] for a in asvs)
            w.writerow([m, lbl, len(asvs), n, f"{100 * n / total:.2f}"])

# whole community at genus level, as in chapter Fig. 6.4B (paper clade replaces the genus when assigned;
# no genus = "Lactobacillaceae (no genus)" or "Unclassified"); the tree only places Lactobacillaceae, the rest keeps GTDB + our 16S
def genus_label(ranks):
    names = dict(ranks)
    if names.get("s", "").startswith("Clado"):
        return names["s"].replace("_", " ", 2).replace("_", " ")
    g = names.get("g", "")
    if not g or "Incertae_Sedis" in g or "uncultured" in g.lower():
        return "Lactobacillaceae (no genus)" if names.get("f") == "Lactobacillaceae" else "Unclassified"
    return g


def ranks_of(tax):
    return [tuple(x.split(":", 1)) for x in tax.replace('"', "").split(",") if ":" in x]


raw = {k: {r[0]: (r[3] if len(r) > 3 else "") for r in tsv(f"{res}/06_sintax/{k}.tsv")}
       for k in ("silva", "gtdb", "gtdb_ours")}
genus = {"RDP (AMPtk)": {a: genus_label(ranks_of(counts[a][2])) for a in counts},
         "SILVA 138.2": {a: genus_label(ranks_of(raw["silva"][a])) for a in counts},
         "GTDB r226": {a: genus_label(ranks_of(raw["gtdb"][a])) for a in counts},
         "GTDB + our 16S": {a: genus_label(ranks_of(raw["gtdb_ours"][a])) for a in counts}}
genus["Tree placement (EPA-ng)"] = {a: (genus_label([(r, n) for r, n in [("f", "")] + [
    (("g", "s")[i], p) for i, p in enumerate(epa[a][1].split(";")[1:3]) if p != "unresolved"]])
    if a in epa else genus["GTDB + our 16S"][a]) for a in counts}
all_reads = sum(c[0] for c in counts.values())
with open(f"{out}/summary_genus_all.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["method", "genus", "n_ASVs", "reads", "pct_reads"])
    for m in METHODS:
        agg = {}
        for a in counts:
            agg.setdefault(genus[m][a], []).append(a)
        for g, asvs in sorted(agg.items(), key=lambda x: -sum(counts[a][0] for a in x[1])):
            n = sum(counts[a][0] for a in asvs)
            w.writerow([m, g, len(asvs), n, f"{100 * n / all_reads:.2f}"])

sens_diff = [a for a in lab if runs["gtdb_ours"][a][0] != runs["gtdb_ours_sens"][a][0]]
print(f"{len(lab)} LAB ASVs, {total} reads; EPA placements: {len(epa)}; "
      f"ASVs changing without MAG_41/45 16S: {len(sens_diff)} {sens_diff}")
