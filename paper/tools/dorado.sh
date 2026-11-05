#!/usr/bin/env bash
# paper/tools/dorado.sh VERSION: one Dorado release (0.9.2, 1.4.0 or 2.1.2) from ONT's Linux x64 archive (dorado.tsv) into
# TOOLS_DIR/dorado-VERSION-linux-x64/, and the models the paper used with it (dorado.tsv's `models`) into DORADO_MODELS,
# downloaded by that release's `dorado download` (the models' files do not depend on the release that fetches them).
# Checks, exit non-zero on any difference:
#   - the archive's SHA-256 against dorado.tsv (when recorded there);
#   - bin/dorado's SHA-256 against the binary the paper's runs used (dorado.tsv `binary_sha256`);
#   - every model file's SHA-256 against models.sha256.
# The archive (2.6-3.5 GB) goes to DL_DIR (default TOOLS_DIR/.dl, resumable) and is deleted once unpacked.
# Settings (paper/config.sh, tools/env.sh): TOOLS_DIR, DORADO_MODELS, CURL, DL_DIR.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
source "$HERE/env.sh"
V=${1:?usage: dorado.sh VERSION (0.9.2, 1.4.0 or 2.1.2)}
row=$(awk -F'\t' -v v="$V" 'NR > 1 && $1 == v' "$HERE/dorado.tsv")
[[ -n $row ]] || { echo "dorado.sh: no Dorado $V in $HERE/dorado.tsv" >&2; exit 2; }
IFS=$'\t' read -r _ URL ARCHIVE_SHA BINARY_SHA MODELS <<< "$row"
D=$TOOLS_DIR/dorado-$V-linux-x64
mkdir -p "$TOOLS_DIR" "$DORADO_MODELS"

if [[ ! -x $D/bin/dorado ]]; then
  dl=${DL_DIR:-$TOOLS_DIR/.dl}; mkdir -p "$dl"
  f=$dl/$(basename "$URL")
  echo "dorado.sh: downloading $URL"
  $CURL -C - -o "$f.part" "$URL"
  got=$(sha256sum "$f.part" | cut -d' ' -f1)
  if [[ $ARCHIVE_SHA == . ]]; then
    echo "dorado.sh: archive SHA-256 $got (none recorded in dorado.tsv; bin/dorado is checked below)"
  elif [[ $got != "$ARCHIVE_SHA" ]]; then
    echo "dorado.sh: $f: SHA-256 $got, expected $ARCHIVE_SHA" >&2; exit 1
  fi
  echo "dorado.sh: unpacking into $TOOLS_DIR"
  tar -xzf "$f.part" -C "$TOOLS_DIR"
  rm -f -- "$f.part"
fi
got=$(sha256sum "$D/bin/dorado" | cut -d' ' -f1)
[[ $got == "$BINARY_SHA" ]] || { echo "dorado.sh: $D/bin/dorado: SHA-256 $got, expected $BINARY_SHA" >&2; exit 1; }
echo "dorado.sh: $D/bin/dorado is the paper's Dorado $V ($("$D/bin/dorado" --version 2>&1 | tail -1))"

# the models: dorado 0.9.2 names the target folder --directory, later releases --models-directory
opt=--models-directory
"$D/bin/dorado" download --help 2>&1 | grep -q -- '--models-directory' || opt=--directory
for m in ${MODELS//,/ }; do
  if [[ ! -s $DORADO_MODELS/$m/config.toml ]]; then
    echo "dorado.sh: downloading model $m"
    rm -rf -- "${DORADO_MODELS:?}/$m"
    "$D/bin/dorado" download --model "$m" "$opt" "$DORADO_MODELS"
  fi
  (cd "$DORADO_MODELS" && grep -v '^#' "$HERE/models.sha256" | grep -F "  $m/" | sha256sum -c --quiet) \
    || { echo "dorado.sh: model $m differs from the paper's (models.sha256)" >&2; exit 1; }
  echo "dorado.sh: model $m: every file equals the paper's"
done
