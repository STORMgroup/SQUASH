#!/usr/bin/env bash
# paper/eval/strands.sh KEY [masked]: the strand of each unitig (minimap2 -cx asm10 --secondary=no of the spelled
# unitigs on the reference; the strand with more aligned bases, strands.py) and the spelled and placed FASTA split by
# strand -> OUT_DIR/eval/KEY/squash/strands/ (spelled.paf, squash_{spelled,placed}_{fwd,rev}.fa, summary.tsv), and
# the closed replicons (closed.py) -> squash/closed.tsv. `masked`: on the set's masked reference, in squash/masked/
# (the spelled unitigs copied there beside that placement).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: strands.sh KEY [masked]}"
D=$ES R=$REF
if [[ ${2:-} == masked ]]; then
  [[ $MASKED != - ]] || die "$KEY has no masked reference"
  D=$ES/masked R=$EVAL_REFS/$MASKED
  mkdir -p "$D/asm"; rm -f "$D/asm/squash_spelled.fa"; cp "$ES/asm/squash_spelled.fa" "$D/asm/squash_spelled.fa"
fi
need "$D/asm/squash_spelled.fa" "$D/asm/squash_placed.fa" "$R"
mkdir -p "$D/strands"
log "strands on $(basename "$R")"
"$MINIMAP2" -cx asm10 --secondary=no -t "$EVAL_THREADS" "$R" "$D/asm/squash_spelled.fa" > "$D/strands/spelled.paf" \
  2> "$D/strands/map.log"
"$REMORA_PYTHON" "$EVAL/strands.py" "$D/strands/spelled.paf" "$D/asm" "$D/strands"
[[ ${2:-} == masked ]] || "$REMORA_PYTHON" "$EVAL/closed.py" "$KEY" "$ES" > "$ES/closed.tsv"
column -t "$D/strands/summary.tsv" >&2
