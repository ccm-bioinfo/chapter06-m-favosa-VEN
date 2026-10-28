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

Methods: [`docs/METHODS.md`](docs/METHODS.md). Summary of the exploration: [`docs/summary.md`](docs/summary.md).

## Layout
| Folder | What it is |
|---|---|
| [`pipeline/`](pipeline) | **clean pipeline**: reproduces the analysis from the Zenodo input package, following [`docs/METHODS.md`](docs/METHODS.md) |
| [`docs/METHODS.md`](docs/METHODS.md) | complete methods (the chapter has a shorter version) |
| [`exploration/`](exploration) | the exploratory analysis that led to the decisions (genome selection, clade labels, duplicates, figures tried) |
| `docs/` | decision tables used by the exploration (`genome16S_decisions.tsv`, `public_genomes.tsv`, `xolalpa2026_genomes.tsv`, `mg16S_labels.tsv`, `mg16S_notes.tsv`) and the summary of the exploration (`summary.md`) |
| `zenodo/tables/` | Tables S1 and S2 as deposited in Zenodo |

## Data policy
No data are committed to GitHub. The inputs of the clean pipeline are in Zenodo
[10.5281/zenodo.23690055](https://doi.org/10.5281/zenodo.23690055) (`chapter06-inputs.tar.gz`); public references
(GTDB r226, SILVA 138.2) are downloaded from their source. Related records: MAGs 10.5281/zenodo.18446951
(Xolalpa-Aroche et al. 2026, 10.1093/ismeco/ycag063); isolate genomes NCBI BioProject PRJNA1235305; full APIBA 16S/ITS
dataset 10.5281/zenodo.20481566. The exploration reads the shared inputs through symlinks in `data/`
(`data/manifest.tsv`, `exploration/scripts/00_link_data.sh`).

## Run the clean pipeline
```bash
CH06_WORK=/data2/<PROJECT>/run SILVA_NR99=/path/to/SILVA_138.2_NR99.fasta source pipeline/config.sh
setsid nohup bash pipeline/run_all.sh > $CH06_WORK/run.log 2>&1 &
```
See [`pipeline/README.md`](pipeline/README.md) for the steps.

## Exploration steps
Scripts in `exploration/scripts/`, run from the repository root. Order:
00 → 03 → 04 → 05 (Mac) → 06 → 07 → 13 → 08 → 09 (Mac) → 10 → 11 → 12 (Mac) → 14 (Mac) → 15 (Mac) → 16.
