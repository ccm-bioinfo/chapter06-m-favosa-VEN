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
Inputs and final results will be deposited in a new Zenodo record (draft metadata:
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
| `00_download_zenodo.sh` | anywhere | `data/` from Zenodo (skeleton until the record exists) |
| `01_inputs_overview.py <host>` | both | input counts (merged into `docs/inputs_overview.tsv`) |
| `02_favosa_asvs.py` | Mac | `results/favosa_ASVs.tsv`: B-VEN-D9 ASVs with AMPtk / GTDB / previous custom taxonomy and linked 16S-MAG clade |
