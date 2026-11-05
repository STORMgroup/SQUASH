#!/usr/bin/env bash
# paper/eval/spell.sh KEY: SQUASH's unitigs spelled through the sup move tables (spell.py) ->
# OUT_DIR/eval/KEY/squash/asm/squash_spelled.fa and squash_spelled.pieces.tsv. Reads SQUASH's layout.tsv.gz and its
# kept work folder (KEEP_WORK=1: work/graph/input/, work/signal/spans.bin) and OUT_DIR/reads/KEY/sup.moves.bam.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: spell.sh KEY}"
squash_inputs
need "$MOVES"
mkdir -p "$ES/asm" "$ES/logs"
log "spelled sequences"
"$REMORA_PYTHON" "$EVAL/spell.py" --layout "$LAYOUT" --graph-input "$GRAPH_INPUT" --spans "$SPANS" --moves "$MOVES" \
  --output "$ES/asm/squash_spelled.fa" > "$ES/logs/spell.log" 2>&1 || { cat "$ES/logs/spell.log" >&2; die "spell.py failed"; }
cat "$ES/logs/spell.log" >&2
