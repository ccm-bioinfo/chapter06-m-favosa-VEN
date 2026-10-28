# Clean pipeline

Reproduces the analysis exactly as written in [`../docs/METHODS.md`](../docs/METHODS.md), starting from the input
package in Zenodo ([10.5281/zenodo.23690055](https://doi.org/10.5281/zenodo.23690055)) and public references. It was
written by an independent rerun that only had the methods text and the package (no access to `exploration/`).

| Step | What it does |
|---|---|
| `config.sh` | paths (`CH06_WORK`, `CH06_INPUTS`, `GTDB_SSU`, `GTDB_TREE`, `SILVA_NR99`) and `THREADS` |
| `00_get_inputs.sh` | downloads the Zenodo package (md5 check), GTDB r226 SSU and the bac120 tree; SILVA 138.2 NR99 must be downloaded by hand from arb-silva.de |
| `01_build_refs.py` | SINTAX-formatted GTDB (>= 1,200 bp), GTDB without N, SILVA |
| `02_build_ours_refs.py` | our 16S references from Tables S1 and S2 (+ sensitivity set without MAG_41/MAG_45), Lactobacillaceae ASVs |
| `03_sintax.sh` | SINTAX of the 273 ASVs against SILVA, GTDB, GTDB + ours and the sensitivity set; best hit against our 16S |
| `04_besthit_gtdb.sh` | exhaustive best hit against GTDB without N |
| `05_tree_inputs.py`, `06_build_tree.sh` | reference 16S tree (Nicoliella–Acetilactobacillus + Fructilactobacillus), constraint from the pruned GTDB bac120 tree, MAFFT, trimAl, IQ-TREE, RAxML-NG |
| `07_labels.py`, `08_place.sh` | labels per ASV and method; EPA-ng placement and gappa assignment of the 68 ASVs |
| `09_summarise.py` | % of reads per label and method, agreement, sensitivity, identities, V4 differences |
| `10_figures.R` | bar figure (ITS + 16S), placement tree, V4 heatmap |

```bash
CH06_WORK=/data2/<PROJECT>/run SILVA_NR99=/path/to/SILVA_138.2_NR99.fasta source pipeline/config.sh
setsid nohup bash pipeline/run_all.sh > $CH06_WORK/run.log 2>&1 &
```

Tools: vsearch 2.32, MAFFT 7.526, trimAl 1.5, IQ-TREE 3.1, RAxML-NG 2.0, EPA-ng 0.3.8, gappa 0.9, seqkit 2;
R 4.5 with ape, ggtree, ggplot2, patchwork, RColorBrewer, dplyr.
