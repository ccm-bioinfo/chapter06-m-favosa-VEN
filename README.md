# chapter06-m-favosa-VEN

Línea: MIEL (mieles)
Título original: Capítulo6_M_favosa_VEN
Objetivo: Revisar Lactobacillaceae y ampliar el análisis

Alnitak: /data2/MIEL/chapter06-m-favosa-VEN
Local: ~/Documents/MIEL/chapter06-m-favosa-VEN
Drive: MIEL/chapter06-m-favosa-VEN

Insumos compartidos: ../m-favosa-VEN-microbiome/ (local) y /data2/MIEL/vit/ (Alnitak); no copiar datos.

## Question
Do the Lactobacillaceae ASVs of the only *M. favosa* sample (B-VEN-D9, 10MFH) get a better
classification when the references include our full-length 16S recovered from MAGs/genomes
(Clades 1–4, Xolalpa-Aroche et al. 2026) instead of GTDB/AMPtk alone?

## Data policy
No data are committed to GitHub: `data/` holds only symlinks (see `data/manifest.tsv`).
Inputs and final results will be deposited in Zenodo record
[10.5281/zenodo.23690055](https://doi.org/10.5281/zenodo.23690055) (DOI reserved, draft; metadata in
`zenodo/metadata.json`) as `chapter06-inputs.tar.gz`, keeping the `data/` layout.
Outside the Mac/alnitak, start with `scripts/00_download_zenodo.sh` instead of `00_link_data.sh`;
later steps are identical. Public references (GTDB r226 SSU) are downloaded from their source.
Related records: MAGs 10.5281/zenodo.18446951 (Xolalpa-Aroche et al. 2026, 10.1093/ismeco/ycag063);
full APIBA 16S dataset 10.5281/zenodo.20481566.

## Layout
- `data/`     symlinks to shared inputs (`manifest.tsv`; recreate with `scripts/00_link_data.sh`)
- `docs/inputs_overview.tsv`  what each input is, where it lives (Mac/alnitak), size and record count
- `scripts/`  numbered steps
- `results/`  outputs (not tracked)

## Steps
| Script | Host | Output |
|---|---|---|
| `00_link_data.sh` | both | `data/` symlinks |
| `00_download_zenodo.sh` | anywhere | `data/` from Zenodo (once the record is published) |
| `01_inputs_overview.py <host>` | both | input counts (merged into `docs/inputs_overview.tsv`) |
| `02_favosa_asvs.py` | Mac | `results/favosa_ASVs.tsv`: first look at B-VEN-D9 ASVs and earlier classifications |
| `03_genome16S_review.py` | alnitak | `results/03_genome16S_review.tsv`: barrnap 16S per genome + decisions in `docs/genome16S_decisions.tsv` |
| `04_favosa_queries.py` | both | `results/04_queries/`: B-VEN-D9 ASVs (all 273; 69 Lactobacillaceae by AMPtk) |
| `05_mg16S_labels.py` | Mac | `docs/mg16S_labels.tsv`: genus and paper clade of each recovered 16S (strongest MarkerMAG link) |
| `06_sintax_refs.py` | alnitak | `results/06_sintax/refs/`: SILVA 138.2, GTDB r226, GTDB + our 16S, sensitivity (no MAG_41/45 16S) |
| `07_sintax.sh` | alnitak | `results/06_sintax/*.tsv`: SINTAX (cutoff 0.8, seed 1) of the 273 ASVs per reference |
| `08_tree_inputs.py` | alnitak | `results/08_tree/`: GTDB + our genome 16S + MG16S references, ASV queries, tip table |
| `09_constraint.R` | Mac | `results/08_tree/constraint.tree`: GTDB genome topology, our 16S at their genus MRCA |
| `10_tree.sh` | alnitak | constrained IQ-TREE, EPA-ng placement, gappa assignment and grafted tree |
| `11_compare.py` | both | `results/11_compare/`: per-ASV classification by method, LAB summary, whole-community genus summary |
| `13_identity.sh` | alnitak | `results/13_identity/`: best-hit identity of each LAB ASV to our 16S vs GTDB |
| `14_fig6_4.R` | Mac | `results/14_fig6_4/`: chapter Fig. 6.4 (A ITS1, B 16S) redrawn side by side from the VIT tables |
| `12_figures.R` | Mac | `results/12_figures/`: fig1 genus composition per method (style of chapter Fig. 6.4B), fig2 circular placement tree (style of VIT lab16S_circular) |

Order: 00 → 03 → 04 → 05 (Mac) → 06 → 07 → 08 → 09 (Mac) → 10 → 11 → 12 (Mac).
