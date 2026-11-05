#!/usr/bin/env bash
# paper/check.sh KEY [OUT]: compare SQUASH's outputs of a paper set (OUT, default OUT_DIR/KEY) with the paper's,
# expected/KEY.sha256 (sha256sum lines; a .gz file by its decompressed text, as gzip bytes depend on the compressor).
# Exit 0 when every listed file is there and equal.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/config.sh"
KEY=${1:?usage: check.sh KEY [OUT]}
OUT=${2:-$OUT_DIR/$KEY}
EXP=$PAPER/expected/$KEY.sha256
[[ -s $EXP ]] || { echo "check.sh: no $EXP: the paper's outputs of $KEY are not recorded yet (PARITY.md)" >&2; exit 2; }
bad=0 n=0
while read -r sum name; do
  [[ -z $sum || $sum == \#* ]] && continue
  n=$((n + 1))
  f=$OUT/$name
  if [[ ! -f $f ]]; then echo "MISSING $name"; bad=$((bad + 1)); continue; fi
  if [[ $f == *.gz ]]; then got=$(gzip -dc -- "$f" | sha256sum | cut -d' ' -f1); else got=$(sha256sum -- "$f" | cut -d' ' -f1); fi
  if [[ $got == "$sum" ]]; then echo "same    $name"; else echo "DIFFERS $name"; bad=$((bad + 1)); fi
done < "$EXP"
if (( bad )); then echo "check.sh: $KEY: $bad of $n files differ or are missing ($OUT)"; exit 1; fi
echo "check.sh: $KEY: all $n files equal the paper's ($OUT)"
