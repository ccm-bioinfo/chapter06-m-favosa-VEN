#!/usr/bin/env python3
# 05_tree_inputs.py - Collect the 41 reference 16S of the placement tree (14 GTDB reps: longest SSU >= 1400 bp,
# 7 public genomes, 7 Xolalpa genomes, 13 recovered Nicoliella/Acetilactobacillus 16S), write tip metadata and the
# IQ-TREE constraint: GTDB bac120 r226 tree pruned to the 14 reps, our genome 16S + clade-linked recovered 16S
# added as a polytomy at the MRCA of their genus, clade-less recovered 16S left out (unconstrained).
import csv, gzip, os, re, sys
sys.setrecursionlimit(100000)
W = os.environ["CH06_WORK"] + "/"; I = os.environ["CH06_INPUTS"] + "/"
S1 = list(csv.DictReader(open(I + "tables/Table_S1_genomes.tsv"), delimiter="\t"))
S2 = list(csv.DictReader(open(I + "tables/Table_S2_recovered_16S.tsv"), delimiter="\t"))

def fasta(p):
    d, h = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): h = l[1:].split()[0]; d[h] = []
        elif h: d[h].append(l)
    return {k: "".join(v).upper() for k, v in d.items()}

def clade_lab(c):
    if c in ("-", ""): return ""
    if c.startswith("Clade 3") or c.startswith("Clade 4"): return "Clade 3-4"
    return "Clade " + c.split()[1]

tips = []  # name, seq, genus, clade, type, shape
reps = {r["genome"]: r for r in S1 if r["source"].startswith("GTDB")}
best = {}
with gzip.open(os.environ["GTDB_SSU"], "rt") as fh:
    h = None; seq = []
    def done(h, seq):
        if h is None: return
        acc = h.split("~")[0][3:]
        if acc in reps:
            s = "".join(seq).upper()
            if len(s) >= 1400 and (acc not in best or len(s) > len(best[acc])): best[acc] = s
    for l in fh:
        if l.startswith(">"): done(h, seq); h = l[1:].split()[0]; seq = []
        else: seq.append(l.strip())
    done(h, seq)
for acc, r in reps.items():
    tips.append((acc, best[acc], r["genus"], "", "GTDB representative", "circle"))
g16 = fasta(I + "genome_16S/genome_16S_barrnap.fa")
for r in S1:
    if r["source"].startswith("GTDB") or r["in_tree"] != "yes": continue
    cl = clade_lab(r["paper_clade_Xolalpa2026"])
    shape = "triangle" if r["genome_type"] == "MAG" else "square"
    typ = "Xolalpa genome" if r["source"].startswith("Xolalpa") else "public genome"
    tips.append((r["genome"], g16[r["16S_used"]], r["genus"], cl, typ, shape))
rec = fasta(I + "recovered_16S/recovered_16S.fa")
for r in S2:
    if r["genus"] in ("Nicoliella", "Acetilactobacillus"):
        tips.append((r["16S"], rec[r["16S"]], r["genus"], clade_lab(r["clade_assigned"]), "recovered 16S", "circle"))

with open(W + "tree/refs.fa", "w") as o, open(W + "results/05_tree_tips.tsv", "w") as t:
    t.write("tip\tgenus\tclade\ttype\tshape\tlength\n")
    for n, s, g, c, ty, sh in tips:
        o.write(f">{n}\n{s}\n"); t.write(f"{n}\t{g}\t{c}\t{ty}\t{sh}\t{len(s)}\n")
# gappa taxonomy file: Lactobacillaceae;genus[;clade]
with open(W + "tree/taxon_file.tsv", "w") as o:
    for n, s, g, c, ty, sh in tips:
        o.write(f"{n}\tLactobacillaceae;{g}" + (f";{c}" if c else "") + "\n")
from collections import Counter
print("tips:", len(tips), Counter(x[4] for x in tips))

# ---- prune GTDB tree ----
txt = open(os.environ["GTDB_TREE"]).read().strip()
def parse(txt):  # iterative Newick parser (deep tree), keeps topology and leaf names only
    root = {"c": [], "name": ""}; stack = [root]; cur = None; i = 0; n = len(txt)
    lab = re.compile(r"'[^']*'|[^,():;]*"); bl = re.compile(r":[^,();]*")
    while i < n:
        ch = txt[i]
        if ch == "(":
            nd = {"c": [], "name": ""}; stack[-1]["c"].append(nd); stack.append(nd); i += 1; cur = None; continue
        if ch == ",": i += 1; cur = None; continue
        if ch == ")": cur = stack.pop(); i += 1
        elif ch == ";": break
        else: cur = {"c": [], "name": ""}; stack[-1]["c"].append(cur)
        m = lab.match(txt, i); cur["name"] = m.group(0).strip("'"); i = m.end()
        m = bl.match(txt, i)
        if m: i = m.end()
    return root["c"][0]
root = parse(txt)
want = {("RS_" + a): a for a in reps} | {("GB_" + a): a for a in reps}
def prune(root):  # iterative post-order pruning, suppressing unary nodes
    order, st = [], [root]
    while st:
        x = st.pop(); order.append(x); st.extend(x["c"])
    res = {}
    for x in reversed(order):
        if not x["c"]:
            res[id(x)] = {"c": [], "name": want[x["name"]]} if x["name"] in want else None
        else:
            kids = [res[id(c)] for c in x["c"] if res[id(c)]]
            res[id(x)] = None if not kids else (kids[0] if len(kids) == 1 else {"c": kids, "name": ""})
    return res[id(root)]
pr = prune(root)
def leaves(n): return [n["name"]] if not n["c"] else [x for c in n["c"] for x in leaves(c)]
assert len(leaves(pr)) == len(reps), leaves(pr)
def nwk(n): return n["name"] if not n["c"] else "(" + ",".join(nwk(c) for c in n["c"]) + ")"
open(W + "tree/gtdb_pruned.nwk", "w").write(nwk(pr) + ";\n")
# add ours at genus MRCA
genus_of = {a: r["genus"] for a, r in reps.items()}
add = {"Nicoliella": [], "Acetilactobacillus": []}
for n, s, g, c, ty, sh in tips:
    if ty in ("public genome", "Xolalpa genome") or (ty == "recovered 16S" and c):
        add[g].append(n)
def mrca(n, g):
    L = leaves(n); gl = [x for x in L if genus_of.get(x) == g]
    for c in n["c"]:
        if all(x in leaves(c) for x in gl): return mrca(c, g)
    return n
for g, extra in add.items():
    m = mrca(pr, g)
    if not m["c"]:  # single representative: turn the leaf into a polytomy
        m["c"] = [{"c": [], "name": m["name"]}]; m["name"] = ""
    m["c"] += [{"c": [], "name": x} for x in extra]
open(W + "tree/constraint.nwk", "w").write(nwk(pr) + ";\n")
print("pruned:", nwk(pr))
