#!/usr/bin/env bash
# 00_get_inputs.sh - Download the input package (Zenodo 10.5281/zenodo.23690055), GTDB r226 SSU and the GTDB r226
# bac120 tree into $CH06_WORK/inputs. SILVA 138.2 NR99 cannot be redistributed: download it from arb-silva.de and
# set SILVA_NR99 to its path.
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "$CH06_WORK/inputs"
[ -d chapter06-inputs ] || { curl -L --fail -o chapter06-inputs.tar.gz \
  "https://zenodo.org/records/23690055/files/chapter06-inputs.tar.gz?download=1"
  tar -xzf chapter06-inputs.tar.gz && rm chapter06-inputs.tar.gz
  python3 -c 'import hashlib,sys; [sys.exit("md5 mismatch: "+f) for m,f in (l.split() for l in open("chapter06-inputs/md5sums.txt")) if hashlib.md5(open("chapter06-inputs/"+f,"rb").read()).hexdigest()!=m]; print("md5 ok")'; }
G=https://data.gtdb.ecogenomic.org/releases/release226/226.0
[ -s "$GTDB_SSU" ]  || curl -L --fail -o "$GTDB_SSU"  $G/genomic_files_all/ssu_all_r226.fna.gz
[ -s "$GTDB_TREE" ] || curl -L --fail -o "$GTDB_TREE" $G/bac120_r226.tree
[ -s "$SILVA_NR99" ] || echo "WARNING: set SILVA_NR99 to the SILVA 138.2 NR99 FASTA (not downloaded automatically)"
ls -la
