#!/usr/bin/env bash
# Create data/ symlinks from data/manifest.tsv. Works on the Mac and on alnitak.
# Usage: bash exploration/scripts/00_link_data.sh [local|alnitak]   (auto-detected if omitted)
set -euo pipefail
cd "$(dirname "$0")/../.."

host=${1:-$( [[ -d /data2/MIEL/vit ]] && echo alnitak || echo local )}
case $host in
  local)   base=$HOME/Documents/MIEL/m-favosa-VEN-microbiome/amplicon; col=2 ;;
  alnitak) base=/data2/MIEL/vit; col=3 ;;
  *) echo "unknown host: $host" >&2; exit 1 ;;
esac

grep -v '^#' data/manifest.tsv | while IFS=$'\t' read -r link src_local src_alnitak desc _zenodo; do
  src=$( [[ $col == 2 ]] && echo "$src_local" || echo "$src_alnitak" )
  if [[ $src == "-" ]]; then echo "skip ($host): $link"; continue; fi
  [[ $src == /* ]] || src=$base/$src   # absolute paths (shared databases) are used as-is
  if [[ ! -e $src ]]; then echo "MISSING: $src" >&2; continue; fi
  mkdir -p "data/$(dirname "$link")"
  ln -sfn "$src" "data/$link"
done
echo "links in data/ ($host):"; find data -type l | sort
