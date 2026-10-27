#!/usr/bin/env python3
"""Step 2a: SINTAX reference databases (d:,p:,c:,o:,f:,g:,s:), all 16S >= 1200 bp except our own.
  silva.fa            SILVA 138.2 NR99, Bacteria + Archaea, to genus (SILVA species = organism names, not curated)
  gtdb.fa             GTDB r226 SSU
  gtdb_ours.fa        GTDB + our 16S: genomes with use=yes (results/03_genome16S_review.tsv) and the 23 MG16S
                      (docs/mg16S_labels.tsv); paper clade as species, e.g. s:Clado_3-4_Acetilactobacillus-like (Clades 3 and 4 merged)
  gtdb_ours_sens.fa   same without the divergent Clade 2 16S of MAG_41_Melli_001 / MAG_45_Melli_020
  ours_labels.tsv     id, source, genus, species/clade of each of our references
Usage (alnitak): python3 scripts/06_sintax_refs.py   ->  results/06_sintax/refs/"""
import csv, gzip, os, re

MINLEN = 1200
SENS_DROP = {"MAG_41_Melli_001", "MAG_45_Melli_020"}
LAB = "d:Bacteria,p:Bacillota,c:Bacilli,o:Lactobacillales,f:Lactobacillaceae"
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
data = os.path.join(root, "data")
out = os.path.join(root, "results", "06_sintax", "refs")
os.makedirs(out, exist_ok=True)


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    h, s = None, []
    for line in op(path, "rt"):
        if line.startswith(">"):
            if h:
                yield h, "".join(s)
            h, s = line[1:].rstrip("\n"), []
        else:
            s.append(line.strip())
    if h:
        yield h, "".join(s)


def clean(name):
    return re.sub(r"[ ,;:]", "_", name.strip())


# SILVA: 'acc Bacteria;Phylum;Class;Order;Family;Genus;organism' -> first 6 ranks
n = 0
with open(f"{out}/silva.fa", "w") as fo:
    for h, s in fasta(f"{data}/refs/silva_138.2_NR99.fasta"):
        sid, lin = h.split(" ", 1)
        ranks = lin.split(";")
        if ranks[0] not in ("Bacteria", "Archaea") or len(s) < MINLEN:
            continue
        tax = ",".join(f"{r}:{clean(v)}" for r, v in zip("dpcofg", ranks[:6]) if v and "uncultured" not in v.lower())
        fo.write(f">{sid};tax={tax};\n{s.upper().replace('U', 'T')}\n")
        n += 1
print("SILVA:", n)

# GTDB: 'acc~contig d__..;..;s__Species [..]'
n = 0
with open(f"{out}/gtdb.fa", "w") as fo:
    for h, s in fasta(f"{data}/refs/ssu_all_r226.fna.gz"):
        if len(s) < MINLEN:
            continue
        sid, lin = h.split(" ", 1)
        lin = lin.split(" [")[0]  # lineage up to the [location=..] fields; species names contain spaces
        tax = ",".join(f"{p.split('__')[0]}:{clean(p.split('__')[1])}" for p in lin.split(";") if p.split("__")[1])
        fo.write(f">{sid};tax={tax};\n{s}\n")
        n += 1
print("GTDB:", n)

# our references
gtax = dict(l.rstrip("\n").split("\t") for l in open(f"{data}/ours/custom_taxonomy.tsv"))
with open(os.path.join(root, "docs", "genome16S_decisions.tsv")) as fh:
    clade_of = {r["genome"]: r["clade_paper"] for r in csv.DictReader(fh, delimiter="\t")}
CLADE = {"1": "Clado_1_Nicoliella-related", "2": "Clado_2_Nicoliella-related",
         "3": "Clado_3-4_Acetilactobacillus-like", "4": "Clado_3-4_Acetilactobacillus-like"}
with open(os.path.join(root, "results", "03_genome16S_review.tsv")) as fh:
    use = {r["seq"] for r in csv.DictReader(fh, delimiter="\t") if r["use"] == "yes"}
seqs = dict((h.split(";")[0], s) for h, s in fasta(f"{data}/ours/custom_16S.fa"))

ours = []  # (id, source, tax, genus, species, seq)
for sid in sorted(use):
    g = sid.split("__16S_")[0]
    lin = gtax[g]
    sp = CLADE.get(clade_of.get(g, ""), "")
    parts = [p for p in lin.split(";") if p.split("__")[1]]
    if sp:
        parts = [p for p in parts if not p.startswith("s__")] + [f"s__{sp}"]
    tax = ",".join(f"{p.split('__')[0]}:{clean(p.split('__')[1])}" for p in parts)
    genus = re.search(r"g__([^;]*)", lin).group(1)
    ours.append((sid, "genome", tax, genus, sp or re.search(r"s__(.*)$", lin).group(1), seqs[sid]))
with open(os.path.join(root, "docs", "mg16S_labels.tsv")) as fh:
    mg = {r["tip"]: r for r in csv.DictReader(fh, delimiter="\t")}
for tip, s in fasta(f"{data}/ours/mg16S_LAB.fa"):
    r = mg[tip]
    genus = "" if r["genus"] == "-" else r["genus"]
    sp = r["clade"] or (f"{genus}_sp._{tip}" if genus else "")
    tax = LAB + (f",g:{genus}" if genus else "") + (f",s:{sp}" if sp else "")
    ours.append((tip, "metagenome", tax, genus, sp, s))

with open(f"{out}/ours_labels.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "source", "genus", "species_or_clade", "length"])
    for sid, src, tax, genus, sp, s in ours:
        w.writerow([sid, src, genus, sp, len(s)])

for name, drop in (("gtdb_ours", set()), ("gtdb_ours_sens", SENS_DROP)):
    with open(f"{out}/{name}.fa", "w") as fo:
        for line in open(f"{out}/gtdb.fa"):
            fo.write(line)
        for sid, src, tax, genus, sp, s in ours:
            if sid.split("__16S_")[0] not in drop:
                fo.write(f">OWN_{sid};tax={tax};\n{s}\n")
print("ours:", len(ours), "(genomes", sum(o[1] == "genome" for o in ours), "+ MG16S", sum(o[1] == "metagenome" for o in ours), ")")
