#!/usr/bin/env python3
# 01_build_refs.py - Format GTDB r226 SSU (>= 1200 bp) and SILVA 138.2 NR99 as vsearch SINTAX references, plus the
# GTDB set without ambiguous bases (N) used for the best-hit search (vsearch counts an N as a match).
# Output: refs/gtdb_r226_ge1200.sintax.fa, refs/gtdb_r226_ge1200_noN.fa, refs/silva138.2_nr99.sintax.fa
import gzip, os, re
W = os.environ["CH06_WORK"] + "/"
GTDB, SILVA = os.environ["GTDB_SSU"], os.environ["SILVA_NR99"]


def clean(x):
    return re.sub(r"[,;:\s]+", "_", x.strip())


def fasta(fh):
    h, s = None, []
    for l in fh:
        l = l.rstrip("\n")
        if l.startswith(">"):
            if h is not None:
                yield h, "".join(s)
            h, s = l[1:], []
        else:
            s.append(l)
    if h is not None:
        yield h, "".join(s)


n = k = kn = 0
with gzip.open(GTDB, "rt") as fh, open(W + "refs/gtdb_r226_ge1200.sintax.fa", "w") as out, \
        open(W + "refs/gtdb_r226_ge1200_noN.fa", "w") as outn:
    for h, s in fasta(fh):
        n += 1
        if int(re.search(r"\[ssu_len=(\d+)\]", h).group(1)) < 1200:
            continue
        acc, tax = h.split(" ", 1)
        ranks = []
        for t in tax.split(" [")[0].split(";"):
            r, name = t.split("__", 1)
            if name:
                ranks.append(f"{r}:{clean(name)}")
        s = s.upper().replace("U", "T")
        rec = f">{acc};tax={','.join(ranks)};\n{s}\n"
        out.write(rec); k += 1
        if "N" not in s:
            outn.write(rec); kn += 1
print("GTDB total", n, "kept >= 1200 bp", k, "of them without N", kn)

n = 0
R = ["d", "p", "c", "o", "f", "g"]
with open(SILVA) as fh, open(W + "refs/silva138.2_nr99.sintax.fa", "w") as out:
    for h, s in fasta(fh):
        acc, tax = h.split(" ", 1)
        f = tax.split(";")
        lin = f[:-1][:6] if f[0] in ("Bacteria", "Archaea") else f[:1]  # drop organism; Eukaryota at domain only
        ranks = [f"{R[i]}:{clean(x)}" for i, x in enumerate(lin) if x]
        out.write(f">{acc};tax={','.join(ranks)};\n{s.upper().replace('U', 'T')}\n"); n += 1
print("SILVA written", n)
