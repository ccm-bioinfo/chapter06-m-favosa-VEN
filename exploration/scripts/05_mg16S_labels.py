#!/usr/bin/env python3
"""Genus and paper clade of each Lactobacillaceae 16S recovered from our metagenomes (MG16S, >=1000 bp).
Rule fixed before any classification:
  clade = paper clade only if the STRONGEST MarkerMAG link of that 16S goes to a paper MAG (etiqueta_liga);
  genus = genus of that clade (Clades 1-2 Nicoliella, 3-4 Acetilactobacillus); otherwise the 16S's own
          SINTAX genus (not the bin's: some bins are mixed); '-' = family level only.
Usage (Mac, needs Mac-only link tables): python3 exploration/scripts/05_mg16S_labels.py > docs/mg16S_labels.tsv"""
import csv, os, re, sys

data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "ours")
CLADE_GENUS = {"1": "Nicoliella", "2": "Nicoliella", "3": "Acetilactobacillus", "4": "Acetilactobacillus"}
CLADE_NAME = {"1": "Clade_1_Nicoliella-related", "2": "Clade_2_Nicoliella-related",
              "3": "Clade_3-4_Acetilactobacillus-like", "4": "Clade_3-4_Acetilactobacillus-like"}
# Clades 3 and 4 merged: their 16S are not separable (T16S_1 links to Clade 3 and Clade 4 bins, no clade 3/4
# genome 16S in the references); the genus is the firm call


def rows(name):
    with open(os.path.join(data, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


tips = [l[1:].strip() for l in open(os.path.join(data, "mg16S_LAB.fa")) if l.startswith(">")]
info = {r["tip"]: r for r in rows("mg16S_info.tsv")}
link = {r["16S"]: r["etiqueta_liga"] for r in rows("ligas_16S_resumen.tsv")}

w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["tip", "cluster", "length", "genus", "clade", "strongest_link"])
for tip in tips:
    r = info[tip]
    lk = link.get(r["cluster"], "")
    m = re.match(r"Clado (\d)", lk)
    clade = CLADE_NAME[m.group(1)] if m else ""
    genus = CLADE_GENUS[m.group(1)] if m else r["genero"]
    w.writerow([tip, r["cluster"], r["longitud"], genus, clade, lk or "no link"])
