#!/usr/bin/env bash
# Zenodo input package (chapter06-inputs.tar.gz): everything the methods start from, so the analysis can be
# repeated without the slow steps (AMPtk, phyloFlash/EMIRGE, MarkerMAG, skani). Public references (GTDB r226,
# SILVA 138.2) are not redistributed; README.txt gives their versions and URLs. No UNU rows or sequences.
# Usage (alnitak, after 15_zenodo_tables.py and with the ITS tables copied to data/vit/):
#   bash scripts/16_zenodo_package.sh   ->  zenodo_package/chapter06-inputs.tar.gz
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=/data/nselem/envs/vit/bin:$PATH
P=zenodo_package/chapter06-inputs
rm -rf "$P" && mkdir -p $P/{tables,asv,its,recovered_16S,genome_16S}

cp zenodo/tables/genomes_tree.tsv $P/tables/Table_S1_genomes.tsv
cp zenodo/tables/recovered_16S.tsv $P/tables/Table_S2_recovered_16S.tsv

# 16S ASVs of the M. favosa honey sample B-VEN-D9: sequences, reads and AMPtk taxonomy (submitted chapter)
cp results/04_queries/favosa_all.fa $P/asv/B-VEN-D9_ASVs.fa
cut -f1,2,5 results/04_queries/favosa_counts.tsv > $P/asv/B-VEN-D9_ASV_table.tsv

# ITS1 genus tables of the same sample (VIT)
cp data/vit/ITS_genus_counts_Mfav.tsv data/vit/ITS_genus_relabund_Mfav.tsv $P/its/

# 16S recovered from the Xolalpa 2026 metagenomes (phyloFlash): general P16S and targeted T16S
cat /data2/MIEL/vit/16S_link/pool/pooled_16S_99.fa /data2/MIEL/vit/16S_link/target/target_16S.fa > $P/recovered_16S/recovered_16S.fa

# every barrnap 16S copy of the non-GTDB genomes of Table S1 (the copy used is in the table)
awk -F'\t' 'NR>1 && $2!="GTDB r226 representative"{print $1}' $P/tables/Table_S1_genomes.tsv > $P/genome_ids.tmp
seqkit fx2tab data/ours/custom_16S.fa | awk -F'\t' 'NR==FNR{k[$1]=1; next} {id=$1; sub(/;.*/, "", id); g=id; sub(/__16S_.*/, "", g); if (g in k) print ">"id"\n"$2}' \
  $P/genome_ids.tmp - > $P/genome_16S/genome_16S_barrnap.fa
rm $P/genome_ids.tmp

cat > $P/README.txt <<'EOF'
Chapter 6 (Melipona favosa honey, Venezuela): inputs for the Lactobacillaceae reclassification
Repository: https://github.com/ccm-bioinfo/chapter06-m-favosa-VEN

tables/Table_S1_genomes.tsv            genomes of the Nicoliella-Acetilactobacillus 16S tree (GTDB r226 representatives,
                                       public genomes and every Xolalpa-Aroche et al. 2026 clade genome), barrnap copies,
                                       16S copy used and MarkerMAG links received
tables/Table_S2_recovered_16S.tsv      16S recovered from the 17 Xolalpa-Aroche et al. 2026 metagenomes: strongest
                                       MarkerMAG link, MAG of the same species (skani, ANI >= 95 %), clade and genus
asv/B-VEN-D9_ASVs.fa                   273 16S V4 ASVs (515F/806R) of sample B-VEN-D9 (AMPtk/DADA2)
asv/B-VEN-D9_ASV_table.tsv             reads per ASV in B-VEN-D9 and AMPtk taxonomy (RDP)
its/ITS_genus_*_Mfav.tsv               ITS1 genus counts and relative abundance of the same sample
recovered_16S/recovered_16S.fa         the 16S of Table S2 (P16S general assembly, T16S assembly targeted to ASV1/ASV2)
genome_16S/genome_16S_barrnap.fa       all barrnap 16S copies of the non-GTDB genomes of Table S1 (<genome>__16S_<n>)

Public references, not redistributed (download from the source):
  GTDB r226 SSU           https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_all/ssu_all_r226.fna.gz
  GTDB r226 bac120 tree   https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_r226.tree
  SILVA 138.2 NR99        https://www.arb-silva.de/ (SILVA_138.2_SSURef_NR99)

Related records: MAGs 10.5281/zenodo.18446951; isolate genomes NCBI BioProject PRJNA1235305;
full APIBA 16S/ITS dataset 10.5281/zenodo.20481566.
EOF
( cd $P && find . -type f ! -name md5sums.txt | sort | xargs md5sum > md5sums.txt )
tar -czf zenodo_package/chapter06-inputs.tar.gz -C zenodo_package chapter06-inputs
grep -c ">" $P/asv/B-VEN-D9_ASVs.fa $P/recovered_16S/recovered_16S.fa $P/genome_16S/genome_16S_barrnap.fa
ls -la zenodo_package/chapter06-inputs.tar.gz
