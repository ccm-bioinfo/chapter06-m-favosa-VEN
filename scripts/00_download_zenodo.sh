#!/usr/bin/env bash
# Fill data/ from the Zenodo record instead of symlinks (for anyone outside alnitak).
# Same layout as scripts/00_link_data.sh, so every later step runs unchanged.
# Usage: bash scripts/00_download_zenodo.sh
# SKELETON: set ZENODO_RECORD once the record is published (DOI reserved in the draft).
set -euo pipefail
cd "$(dirname "$0")/.."

ZENODO_RECORD=${ZENODO_RECORD:-TODO}   # e.g. 23700000 -> https://zenodo.org/records/23700000
ARCHIVE=chapter06-inputs.tar.gz

if [[ $ZENODO_RECORD == TODO ]]; then
  echo "No Zenodo record yet; use scripts/00_link_data.sh on the Mac or alnitak." >&2; exit 1
fi

# 1. Deposited inputs: one archive keeping the vit/ ours/ refs/ layout, checked against Zenodo's md5
url=https://zenodo.org/records/$ZENODO_RECORD/files/$ARCHIVE
md5=$(curl -s https://zenodo.org/api/records/$ZENODO_RECORD |
      python3 -I -c "import json,sys; print([f['checksum'] for f in json.load(sys.stdin)['files'] if f['key']=='$ARCHIVE'][0].split(':')[1])")
curl -L --fail -o data/$ARCHIVE "$url?download=1"
echo "$md5  data/$ARCHIVE" | md5sum -c -
tar -xzf data/$ARCHIVE -C data && rm data/$ARCHIVE

# 2. Public references: downloaded from their source (manifest column 5 is a URL), not redeposited
grep -v '^#' data/manifest.tsv | while IFS=$'\t' read -r link _ _ _ zenodo; do
  [[ $zenodo == http* ]] || continue
  mkdir -p "data/$(dirname "$link")"
  [[ -s data/$link ]] || curl -L --fail -o "data/$link" "$zenodo"
done

# 3. Check every manifest entry is present
grep -v '^#' data/manifest.tsv | cut -f1 | while read -r link; do
  [[ -s data/$link ]] || echo "MISSING: data/$link" >&2
done
