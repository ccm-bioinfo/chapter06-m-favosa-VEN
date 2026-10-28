#!/usr/bin/env python3
"""Count records and size of every input linked in data/ (run on each host).
Usage: python3 exploration/scripts/01_inputs_overview.py <host> > overview_<host>.tsv"""
import csv, gzip, os, sys

host = sys.argv[1]
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
print("link\thost\tn_records\tunit\tsize_MB")
with open(os.path.join(root, "data", "manifest.tsv")) as fh:
    for row in csv.reader((l for l in fh if not l.startswith("#")), delimiter="\t"):
        path = os.path.join(root, "data", row[0])
        if not os.path.exists(path):
            continue
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            if ".fa" in path or ".fna" in path:
                n, unit = sum(1 for l in f if l.startswith(">")), "sequences"
            elif path.endswith(".tree"):
                n, unit = f.read().count(",") + 1, "tips"
            else:
                n, unit = sum(1 for l in f if l.strip() and not l.startswith("#")), "rows"
        size = os.path.getsize(os.path.realpath(path)) / 1e6
        print(f"{row[0]}\t{host}\t{n}\t{unit}\t{size:.2f}")
