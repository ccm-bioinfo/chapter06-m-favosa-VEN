# data/

Symlinks only — no data are stored in this repo. Recreate them with
`bash scripts/00_link_data.sh` (auto-detects Mac vs alnitak). Sources and descriptions: `manifest.tsv`.

- `vit/`  Patricia Vit (APIBA-ULA) 16S V4 metabarcoding. Tables cover all 39 samples;
  this chapter uses only **M. favosa = B-VEN-D9 (10MFH, honey, Paraguaná, VEN)**.
  5MMH is *M. mimetica* (ECU), not *M. favosa*.
- `ours/` Our full-length 16S: recovered from 17 MEX metagenomes (phyloFlash + MarkerMAG links to MAGs;
  Xolalpa-Aroche et al. 2026) and barrnap 16S from MAGs/isolates/public genomes (UNU excluded).
- `refs/` GTDB r226 SSU and the GTDB genome topology used to constrain the 16S tree.

Files available on only one side (marked `-` in the manifest):
- alnitak only: `ours/custom_16S.clean.noUNU.sintax.fa`, `refs/ssu_all_r226.fna.gz`
- Mac only: `ours/ASV_16S_recuperado_bin.tsv`, `ours/ligas_16S_*.tsv`
