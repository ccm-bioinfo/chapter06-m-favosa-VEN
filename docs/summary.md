# Chapter 6 (*Melipona favosa*, Venezuela): Lactobacillaceae reclassified with stingless-bee 16S

Repository: github.com/ccm-bioinfo/chapter06-m-favosa-VEN · Alnitak: `/data2/MIEL/chapter06-m-favosa-VEN` ·
Data: Zenodo 10.5281/zenodo.23690055 (draft; nothing is committed to GitHub)

## Question
The submitted chapter classified the 16S ASVs of the only *M. favosa* honey sample (B-VEN-D9, 10MFH) with
AMPtk/RDP: most Lactobacillaceae reads had no genus and part was called *Lactobacillus*. Does adding the full-length
16S we recovered from stingless-bee MAGs, isolates and metagenomes (Xolalpa-Aroche et al. 2026, ISME Commun,
10.1093/ismeco/ycag063) to the references resolve them?

## Answer
Yes, to genus. With GTDB r226 + our 16S, **99.4 % of the Lactobacillaceae reads (88.8 % of all 16S reads) are
*Acetilactobacillus* of the stingless-bee lineage of Clades 3–4** of Xolalpa-Aroche et al. 2026. SINTAX and
phylogenetic placement agree on 99.98 % of the reads. The abundant ASVs are 99.2–100 % identical over V4 to our 16S,
but at most 94 % identical to any GTDB r226 16S. That is why RDP, SILVA and GTDB cannot assign them a genus.
Clades 3 and 4 are reported together ("Clade 3-4") because the available 16S do not separate them (see Caveats).

## Data
| Input | Content | Source |
|---|---|---|
| ASVs | 273 V4 ASVs (515F/806R, ~253 bp), 23,491 reads in B-VEN-D9; 69 Lactobacillaceae by AMPtk (20,998 reads, 89.4 %) | VIT/APIBA (stajichlab), AMPtk run of the submitted chapter |
| ITS1 | genus table of the same sample, 44,408 reads | VIT/APIBA tables (not recalculated) |
| Our genome 16S | barrnap 16S of 28 genomes used as references: 7 Xolalpa 2026 MAGs, 3 published *Nicoliella* isolates (barcode03/04/20, PRJNA1235305), 18 public bee-LAB isolate genomes | Zenodo 10.5281/zenodo.18446951 and NCBI |
| Our metagenome 16S | 23 Lactobacillaceae 16S (>= 1,000 bp) recovered from 17 Mexican stingless-bee metagenomes (phyloFlash), linked to MAGs with MarkerMAG | this work, from the Xolalpa 2026 metagenomes |
| GTDB r226 SSU | 747,353 sequences >= 1,200 bp; for the tree, 139 Lactobacillaceae genomes of the earlier constrained tree | gtdb.ecogenomic.org |
| SILVA 138.2 NR99 | 356,078 Bacteria/Archaea >= 1,200 bp | arb-silva.de |

The selection of genomes is decided in one place, `docs/genome16S_decisions.tsv`. Every later step, SINTAX and
the tree alike, uses the same set.
- **Excluded genomes:**
  - unpublished isolates (barcode07–11, 15, 17);
  - the Costa Rica MAG CR_A551_006;
  - MAG_42_Melli_004 (16S of only 569 bp);
  - MAG_89_Scapto_001 (its only 16S is a *Nicoliella* contaminant);
  - *F. fructivorans* DSM 20203 (its 16S is *Wolbachia*);
  - Es01, SGEP1A5 and *A. jinshanensis* HSLZ-75, which are the same assemblies as GTDB references GCF_037546185.1, GCF_023380205.1 and GCF_004359375.1;
  - all UNU isolates (no cooperation letter).
- **One 16S per genome:** the longest copy validated against the genome's GTDB-Tk genus. Ties go to the lowest copy number, so reruns are reproducible.
- **Excluded from the tree:** two metagenome 16S (P16S_15, P16S_62), because their terminal branches were saturated.

## Methods
1. **Clade labels of the metagenome 16S.**
   - A recovered 16S gets a paper clade only if its strongest MarkerMAG link (read pairs) goes to a bin of the same species (ANI >= 95 %) as a Xolalpa 2026 MAG.
   - This rule was fixed before any classification.
   - Clades 3 and 4 are merged as "Clade 3-4 *Acetilactobacillus*-like".
   - Labelled 16S: P16S_19 (Clade 2); T16S_1, T16S_2 and P16S_24 (Clade 3-4).
2. **SINTAX.**
   - Tool: vsearch 2.32.0, cutoff 0.8, `--randseed 1`, one thread per run.
   - All 273 ASVs were classified against four references:
     - SILVA 138.2;
     - GTDB r226;
     - GTDB r226 + our 16S (paper clade written as the species rank);
     - the same without the 16S of MAG_41/MAG_45, as a sensitivity test.
3. **Phylogenetic placement.**
   - Alignment and trimming:
     - MAFFT 7.526 L-INS-i aligned the 184 references (139 GTDB, 24 genomes, 21 metagenome 16S).
     - The 69 Lactobacillaceae ASVs were added with `--addfragments --keeplength`.
     - trimAl 1.5 (`-gt 0.05`) trimmed the alignment.
   - Reference tree:
     - IQ-TREE 3.1.4, GTR+F+G4, 1,000 ultrafast bootstraps.
     - It was constrained to the GTDB bac120 genome topology, with our 16S free inside the MRCA of their genus.
     - RAxML-NG 2.0.3 re-estimated the model parameters.
   - Placement: EPA-ng 0.3.8, then gappa 0.9.0 (taxonomic assignment by best hit, plus the support of each assignment summed over all placements, and the grafted tree).
4. **Identity.** Best V4 hit of each Lactobacillaceae ASV against our 16S and against GTDB r226 alone (vsearch `--usearch_global`).
5. **V4 differences.**
   - Pairwise differences between the 101 sequences of the *Nicoliella*–*Acetilactobacillus* subtree.
   - Only the alignment window covered by the ASVs is used (253 columns).
   - Substitutions and internal indels are counted; N and ambiguity codes are not.
6. **Figures.**
   - R 4.5.1, ggplot2 4.0.3, ggtree 3.16.3, ape 5.8.1.
   - Genus palette of chapter Fig. 6.4: RColorBrewer Paired, Other grey65, Unclassified grey88.

## Results
### Whole 16S community (% of the 23,491 reads)
| Reference | Main label | % |
|---|---|---|
| RDP (AMPtk, submitted) | Lactobacillaceae without genus 77.9 · *Lactobacillus* 11.5 | 89.4 |
| SILVA 138.2 | Lactobacillaceae without genus 88.9 · *Lactobacillus* 0.5 | 89.4 |
| GTDB r226 | Unclassified 89.9 | – |
| GTDB r226 + our 16S | **Clade 3-4 *Acetilactobacillus*-like 88.8** | 88.8 |
| Tree placement (EPA-ng) | **Clade 3-4 *Acetilactobacillus*-like 88.8** | 88.8 |

### The nine Lactobacillaceae ASVs with >= 100 reads
| ASV | % of Lactobacillaceae reads | RDP (submitted) | GTDB + our 16S (confidence) | Tree (support summed over placements) | Best hit, ours / GTDB |
|---|---|---|---|---|---|
| ASV2 | 66.5 | family only | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 99.6 % T16S_1 / 93.7 % |
| ASV1 | 10.9 | family only | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 100 % T16S_2 / 93.7 % |
| ASV66 | 9.3 | *Lactobacillus* | Clade 3-4 (1.00) | Clade 3-4 (0.97) | 99.6 % T16S_2 / 93.3 % |
| ASV37 | 3.7 | family only | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 99.2 % T16S_1 / 93.3 % |
| ASV168 | 2.2 | *Lactobacillus* | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 99.2 % T16S_2 / 93.3 % |
| ASV161 | 2.0 | family only | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 99.6 % T16S_2 / 93.3 % |
| ASV81 | 1.6 | family only | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 99.2 % T16S_1 / 93.3 % |
| ASV56 | 1.2 | *Lactobacillus* | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 100 % T16S_1 / 93.3 % |
| ASV118 | 0.5 | family only | Clade 3-4 (1.00) | Clade 3-4 (1.00) | 99.2 % T16S_1 / 94.1 % |

- **Identity:** 99.4 % of the Lactobacillaceae reads are >= 99 % identical to one of our 16S. None is >= 99 % identical to a GTDB r226 16S, and only 0.5 % are >= 97 % identical.
- **Clade 2 (*Nicoliella*-related):** only ASV10 and ASV83, 3 reads in total.
- **Sensitivity:** removing the 16S of MAG_41/MAG_45 changes only ASV83 (2 reads).

### V4 differences (Fig. 5)
- **Two V4 variants inside Clade 3-4.**
  - ASV1, ASV66 and relatives sit 0–1 differences from T16S_2.
  - ASV2, ASV37 and relatives sit 1 difference from T16S_1.
  - The two variants are 5–6 differences apart.
- **Clade 3-4 against GTDB *Acetilactobacillus*:** 21–27 differences from *A. jinshanensis*.
- **Clade 2 against public *Nicoliella*:** 0–7 differences, so V4 does not separate *Nicoliella* species.

### Fungi (ITS1, Table 6.4)
- **Denominator:** Table 6.4 percentages use the sum of the 15 listed genera (41,679 reads), not all 44,408 ITS reads.
- **Typo in the text:** *Debaryomyces* is 26.8 % of the 15 genera (Table 6.4: 27 %) and 25.1 % of all reads. The text's 37 % is a typo for 27 %.

## Figures (`results/12_figures/`, `results/14_fig6_4/`)
- **Fig. 6.4 (revised):**
  - A: ITS1 genera (VIT data);
  - B: the same 16S ASVs by classification reference.
- **fig1:** genus composition by reference; this is panel B.
- **fig2:** circular placement tree of all Lactobacillaceae references and ASVs.
- **fig3:** zoom of the *Nicoliella* sector, every tip with its origin.
- **fig4 / fig4b:** the *Nicoliella*–*Acetilactobacillus* subtree, rooted on *Fructilactobacillus*, rectangular and circular. Proposed for the article.
- **fig5:** heatmap of pairwise V4 differences in that subtree, with the full matrix as TSV.

The figures use these symbols:
- triangles: Xolalpa 2026 MAGs;
- squares: isolates (Xolalpa 2026 and public);
- circles: GTDB references, metagenome 16S and ASVs.

## Caveats
1. **Clades 3 and 4 are not separated.**
   - No Clade 3 or 4 genome contributes a usable 16S, so these clades are represented only by metagenome 16S labelled by MarkerMAG.
   - T16S_1, the reference of ASV2, links to a Clade 3 bin (6,392 pairs) and also to Clade 4 bins (4,623 pairs).
   - Separately, the tree's support for Clade 4 alone was only 0.56–0.62. Merged, the support is 0.97–1.0.
   - The two V4 variants exist, but which variant belongs to which paper clade is unresolved. The genus call is firm.
2. **P16S_24 carries a Clade 3-4 label from a weak link (450 pairs)** but is 19–24 V4 differences from T16S_1/T16S_2 and is placed with the common *Acetilactobacillus*. Its label is probably wrong. It does not change any abundant ASV.
3. **MAG_41/MAG_45 16S differ from the other Clade 2 sequences.** They are kept because they passed the a priori rule; the sensitivity test shows they do not drive the result.
4. **One sample.** The whole 16S result rests on a single honey sample (B-VEN-D9).

## Reproduce
- **Step order:** 00 → 03 → 04 → 05 (Mac) → 06 → 07 → 13 → 08 → 09 (Mac) → 10 → 11 → 12 (Mac) → 14 (Mac).
- **Starting from Zenodo:** use `scripts/00_download_zenodo.sh` instead of `00_link_data.sh`.
- **Environment:** heavy steps run on Alnitak in `/data/nselem/envs/vit`.
