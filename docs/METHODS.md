# Detailed methods: Lactobacillaceae of *Melipona favosa* honey reclassified with stingless bee 16S

This is the complete description of the analysis; the chapter text is a shorter version of it. The scripts in
[`pipeline/`](../pipeline) implement exactly what is written here, starting from the input package in Zenodo
([10.5281/zenodo.23690055](https://doi.org/10.5281/zenodo.23690055)). They were written by an independent rerun that
only had this text and the package. The scripts in [`exploration/`](../exploration) are the earlier analysis that led to
the decisions recorded in Tables S1 and S2.

## 1. Inputs

| Input | Content | Source |
|---|---|---|
| ASVs | 273 V4 ASVs (515F/806R, ~253 bp) of honey sample B-VEN-D9, reads per ASV and AMPtk (RDP) taxonomy | AMPtk/DADA2 run of the submitted chapter (Vit and Stajich 2026, 10.5281/zenodo.20481566); package `asv/` |
| ITS1 | genus counts of the same sample | same run; package `its/` |
| Table S1 | genomes of the 16S tree and every Xolalpa-Aroche et al. (2026) clade genome; barrnap 16S copies, copy used | package `tables/`, `genome_16S/` |
| Table S2 | 16S recovered from the 17 Xolalpa-Aroche et al. (2026) metagenomes, MarkerMAG link, clade, genus | package `tables/`, `recovered_16S/` |
| GTDB r226 SSU | `ssu_all_r226.fna.gz` | data.gtdb.ecogenomic.org, release 226 |
| GTDB r226 bac120 tree | `bac120_r226.tree` | same |
| SILVA 138.2 NR99 | SSU Ref NR99 (headers `accession Domain;...;Genus;organism`) | arb-silva.de (not redistributed) |

## 2. Steps done before the package (not repeated by the pipeline)

**16S recovered from metagenomes (Table S2).**
- phyloFlash (Gruber-Vodicka et al., 2020) recovered the 16S reads of each of the 17 shotgun metagenomes where the MAGs were assembled, assembled them with SPAdes (Bankevich et al., 2012) and reconstructed them with EMIRGE (Miller et al., 2011).
- Retained: sequences longer than 1,200 bp (P16S, general assembly), and sequences longer than 1,000 bp from an assembly targeted to the reads of ASV1 and ASV2 (T16S). All were clustered with vsearch (Rognes et al., 2016).
- MarkerMAG (Song et al., 2022) was run in each metagenome against all the bins of that sample. It counts read pairs with one read in the 16S and its mate in a contig of the bin.
- skani (Shaw and Yu, 2023) named each bin after the Xolalpa-Aroche et al. (2026) MAG of the same species (ANI ≥ 95 %). Each 16S received the clade of its strongest link only.
- Clades 3 and 4 were merged (Clade 3-4) because no genome of these clades has a usable 16S, and the 16S that represent them link to bins of both clades (T16S_1: 6,392 pairs to a Clade 3 bin and 4,623 to Clade 4 bins).
- Column `genus` of Table S2 is the genus of the clade (Nicoliella for Clades 1 and 2, Acetilactobacillus for Clade 3-4), otherwise the SINTAX genus against GTDB.
- Column `Lactobacillaceae_reference` marks the 23 Lactobacillaceae 16S used as references. P16S_24 keeps a Clade 3-4 label that is flagged as doubtful (weak link, 450 pairs).

**16S of genomes (Table S1).**
- barrnap (Seemann, 2018) found 67 16S copies in the 25 non-GTDB genomes: 15 Xolalpa-Aroche MAGs, 3 isolates and 7 public genomes of Nicoliella and Acetilactobacillus. Public genomes of the same strain as a GTDB representative were not included.
- Copies whose genus (SINTAX against GTDB) disagrees with the genome's GTDB-Tk genus (Chaumeil et al., 2022) were discarded. Of the rest, the longest was used; ties go to the lowest copy number.
- Column `16S_used` of Table S1 is the copy used and is the input of the pipeline. 14 genomes have one: 7 Xolalpa-Aroche (Clade 2) and 7 public.

## 3. SINTAX classification (pipeline 01–03)

**Reference formats** (vsearch SINTAX, `tax=d:,p:,c:,o:,f:,g:,s:`):
- *GTDB r226*: sequences ≥ 1,200 bp (`ssu_len` in the header), all ranks down to species. Characters `, ; : space` in names are replaced by `_`.
- *SILVA 138.2 NR99*: Bacteria and Archaea to genus (the organism field is dropped); Eukaryota at domain only. No length filter.
- *Our 16S* (37 sequences):
  - The 14 genome 16S of Table S1 and the 23 Lactobacillaceae 16S of Table S2.
  - Ranks `d:Bacteria,p:Bacillota,c:Bacilli,o:Lactobacillales,f:Lactobacillaceae`, then genus and species.
  - Sequences linked to a clade carry `g:<genus of the clade>,s:Clade_2` or `s:Clade_3-4`.
  - Other genome 16S carry their genus and species from Table S1.
  - Other recovered 16S carry their genus from Table S2, or only the family when the genus is `-`.
- *GTDB + ours* = GTDB r226 ≥ 1,200 bp plus our 16S.
- *Sensitivity* = the same without the 16S of MAG_41_Melli_001 and MAG_45_Melli_020 (35 of ours).

**Run.** All 273 ASVs were run against each reference with `vsearch --sintax --sintax_cutoff 0.8 --randseed 1 --threads 1`. One thread with a fixed seed makes the result reproducible.

**Label of each ASV per method:**
- the clade, if the species rank passing the cutoff is a clade;
- otherwise the genus passing the cutoff. Genus names starting with uncultured, unknown, unidentified, metagenome or incertae count as no genus;
- otherwise "Lactobacillaceae (no genus)" if the family is Lactobacillaceae, else "Unclassified".

RDP labels come from the AMPtk taxonomy of the package with the same rule. *Lactobacillaceae ASVs* are the 69 that AMPtk assigned to that family (20,998 of 23,491 reads).

**Organelle 16S (pipeline 07).** GTDB has no mitochondrial or chloroplast 16S, so SINTAX against GTDB puts organelle reads on the closest bacterium (473 reads of two ASVs became *Microbacterium*; by BLAST they are 93 % identical to the *Aspergillus niger* mitochondrion, see `docs/exploratory_figures.md`). Before any other label, ASVs whose SILVA classification passing the cutoff has `f:Mitochondria` or `o:Chloroplast` are labelled "Organelles (mitochondria, chloroplast)" in the SILVA, GTDB, GTDB + ours and sensitivity columns (21 ASVs, 589 reads). RDP keeps the submitted labels.

**Final classification (Table 6.5; pipeline 07, column `final`):**
1. organelle 16S as above;
2. ASVs whose GTDB + ours classification passing the cutoff is in Lactobacillaceae (the same 69 ASVs) take the GTDB + ours label;
3. every other ASV takes the SILVA label.

GTDB + ours is used only for Lactobacillaceae, where placement and identity back it up. For the other bacteria SILVA agrees with RDP and with BLAST, while GTDB sent *Haemophilus* (ASV59, 100 % to *H. parainfluenzae*) to *Escherichia* and *Fusobacterium* (ASV112) to *Rothia*. No Lactobacillaceae ASV is an organelle; the pipeline stops if one is.

## 4. Best-hit identity (pipeline 03–04)

Each of the 69 Lactobacillaceae ASVs was compared with our 16S and with GTDB r226 alone using an exhaustive search: `vsearch --usearch_global --maxaccepts 0 --maxrejects 0 --top_hits_only --id 0.5`.

The GTDB set excludes every sequence with an ambiguous base (N). vsearch counts an N as a match, and some GTDB SSU sequences are 19–32 % N, so they appear 100 % identical to any ASV.

## 5. Phylogenetic placement (pipeline 05–08)

**References (41 tips):**
- *GTDB*: the GTDB r226 representatives of Acetilactobacillus, Philodulcilactobacillus and Nicoliella (n = 5) and of Fructilactobacillus, the outgroup (n = 9). Each is represented by its longest SSU ≥ 1,400 bp.
- *Genomes*: the 14 genome 16S of Table S1.
- *Recovered*: the 13 recovered 16S whose `genus` in Table S2 is Nicoliella or Acetilactobacillus.

**Alignment and tree:**
1. The references were aligned with MAFFT L-INS-i (`--localpair --maxiterate 1000`; Katoh and Standley, 2013).
2. Columns with more than 95 % gaps were removed with trimAl (`-gt 0.05`; Capella-Gutiérrez et al., 2009).
3. The tree was inferred with IQ-TREE (`-m GTR+F+G4 -B 1000 --seed 1`; Minh et al., 2020) under a constraint tree:
   - the GTDB r226 bac120 tree pruned to the 14 GTDB representatives, keeping their topology;
   - the genome 16S and the recovered 16S with a clade added as an unresolved group (polytomy) at the most recent common ancestor of the GTDB representatives of their genus. A genus with a single representative becomes a polytomy with that leaf;
   - the recovered 16S without a clade left out of the constraint (unconstrained).
4. RAxML-NG `--evaluate` (Kozlov et al., 2019) re-estimated the model parameters (GTR+G4+FC) on the IQ-TREE tree.

**Placement:**
- *Which ASVs*: every ASV labelled by SINTAX with GTDB + ours as Nicoliella, Acetilactobacillus, one of their clades, or Lactobacillaceae (no genus). These are 68 ASVs: all the Lactobacillaceae ASVs except ASV404, which SINTAX calls Limosilactobacillus.
- *Alignment*: the ASVs were added to the trimmed alignment with `mafft --addfragments --keeplength`.
- *Placement*: EPA-ng (Barbera et al., 2019) placed them on the RAxML-NG tree with its model.
- *Assignment*: gappa (Czech et al., 2020) used `examine assign --best-hit`, rooted on Fructilactobacillus, with taxonomy `Lactobacillaceae;genus[;clade]`. Each ASV gets the deepest label of its best placement (likelihood weight).
- *Grafting*: `gappa examine graft` grafted the ASVs onto the tree.
- *ASVs not placed* keep their SINTAX (GTDB + ours) label.

## 6. V4 differences (pipeline 09)

- *Window*: the alignment columns from the first to the last column covered by any placed ASV.
- *Sequences compared*: the 41 references and the 68 ASVs, all against each other.
- *What counts*: a column counts 1 if both sequences have different bases, or one has a base and the other an internal gap.
- *What does not count*: terminal gaps (outside the overlap of the two sequences) and N.
- *Limitation*: insertions present only in ASVs are lost because of `--keeplength`.

## 7. Figures (pipeline 10)

**Fungi.** The 15 most abundant ITS1 genera are shown; the rest are grouped as Other, and reads without genus as Unclassified.

**Bacteria.**
- One stacked bar per classification (RDP, SILVA, GTDB, GTDB + ours, tree placement, and the final SILVA + GTDB + curated), with all 23,491 reads.
- Organelle 16S is its own dark grey segment (all bars except RDP).
- Labels below 0.3 % of the reads in every bar are grouped as Other; Unclassified stays separate.
- Colours are RColorBrewer Paired in alphabetical order, Other grey and Unclassified light grey.

**Placement tree.**
- A cladogram (no branch lengths) of the grafted tree, rooted on Fructilactobacillus, with the area of each genus shaded. Clade 3-4 is shaded with Acetilactobacillus and Clade 2 with Nicoliella.
- MAGs are triangles, isolate genomes squares, GTDB representatives, recovered 16S and ASVs circles.
- ASV point size is its share of the Lactobacillaceae reads.

**V4 heatmap.** The V4 differences are shown in the order of the tree.

## 8. Expected results (checks)

| Quantity | Value |
|---|---|
| Clade 3-4 (GTDB + ours), % of all 16S reads | 88.84 |
| Organelle 16S (SILVA), ASVs / reads / % of all 16S reads | 21 / 589 / 2.51 |
| Final classification, % of all 16S reads | Clade 3-4 88.84, *Dysgonomonas* 1.50, *Haemophilus* 0.99, *Comamonas* 0.60, *Acetilactobacillus* 0.52, Unclassified 0.60 (3 Enterobacteriaceae ASVs, 13 reads, change genus between SILVA runs; 0.55 in the exploration) |
| Clade 3-4, % of Lactobacillaceae reads (SINTAX and tree) | 99.39 |
| SINTAX vs tree, same label | 99.98–99.99 % of Lactobacillaceae reads |
| Sensitivity run | only ASV83 (2 reads) changes |
| Lactobacillaceae reads ≥ 99 % identical to our 16S / to GTDB without N | 99.4 % / 0.01 % (best GTDB hits 93–94 %) |
| V4 differences ASV1–ASV2, ASV1–T16S_2, ASV2–T16S_1, T16S_1–T16S_2 | 6, 0, 1, 5 |

## 9. Software

vsearch 2.32.0, MAFFT 7.526, trimAl 1.5, IQ-TREE 3.1.4, RAxML-NG 2.0.3, EPA-ng 0.3.8, gappa 0.9.0, seqkit 2.14; R 4.5.1 with
ape 5.8, ggtree 3.16, ggplot2 4.0, patchwork, RColorBrewer, dplyr. References for each tool are in the chapter.

## 10. Authorship and use of AI

Claude Code (Anthropic; model Claude Opus 5.5) was used as a programming and writing assistant.

- **N. Sélem-Mojica**: designed the analysis; provided the 16S recovered from the metagenomes with their MarkerMAG and
  skani links; decided which genomes, 16S and clades are used (Tables S1 and S2, `docs/*.tsv`); reviewed each step and
  figure; wrote the methods.
- **Claude**: wrote and ran the scripts in `exploration/` (classification, placement, comparisons, figures), built
  Tables S1 and S2 and the Zenodo package, and drafted parts of the methods, results and discussion, revised by the authors.
- **Reproducibility check**: a separate Claude session with access only to this text and the Zenodo package repeated the
  analysis and obtained the expected results of section 8. Its scripts, with parametrised paths, are `pipeline/`.

Percentages of 16S reads in the chapter (Fig. 6.4B, Table 6.5) are of all 23,491 reads, including the 614 reads without
an RDP hit.
