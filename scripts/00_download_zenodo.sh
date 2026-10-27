#!/usr/bin/env bash
# Fill data/ from the Zenodo record instead of symlinks (for anyone outside alnitak).
# Same layout as scripts/00_link_data.sh, so every later step runs unchanged.
# Usage: bash scripts/00_download_zenodo.sh
# Record 23690055 (DOI 10.5281/zenodo.23690055) is a draft until the archive is uploaded.
set -euo pipefail
cd "$(dirname "$0")/.."

ZENODO_RECORD=${ZENODO_RECORD:-23690055}
ARCHIVE=chapter06-inputs.tar.gz

if [[ $(curl -s -o /dev/null -w '%{http_code}' https://zenodo.org/api/records/$ZENODO_RECORD) != 200 ]]; then
  echo "Zenodo record $ZENODO_RECORD not published yet; use scripts/00_link_data.sh on the Mac or alnitak." >&2; exit 1
fi

# 1. Deposited inputs: one archive keeping the vit/ ours/ refs/ layout, checked against Zenodo's md5
url=https://zenodo.org/records/$ZENODO_RECORD/files/$ARCHIVE
md5=$(curl -s https://zenodo.org/api/records/$ZENODO_RECORD |
      python3 -I -c "import json,sys; print([f['checksum'] for f in json.load(sys.stdin)['files'] if f['key']=='$ARCHIVE'][0].split(':')[1])")
curl -L --fail -o data/$ARCHIVE "$url?download=1"
if command -v md5sum > /dev/null; then echo "$md5  data/$ARCHIVE" | md5sum -c -
else [ "$(md5 -q data/$ARCHIVE)" = "$md5" ] || { echo "md5 mismatch: data/$ARCHIVE" >&2; exit 1; }; fi   # macOS has md5, not md5sum
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
