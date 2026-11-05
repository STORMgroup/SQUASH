#!/usr/bin/env bash
# paper/eval/tier1.sh KEY: SQUASH's tier-1 contiguity from the raw signal alone (tier1.py squash: signal-estimated
# bases of its layout, G = twice the reference length) -> OUT_DIR/eval/KEY/squash/tier1.tsv.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: tier1.sh KEY}"
squash_inputs
need "$GENOME_REF"
mkdir -p "$ES/logs"
log "tier 1"
"$REMORA_PYTHON" "$EVAL/tier1.py" squash --key "$KEY" --gfa "$GFA" --layout "$LAYOUT" --graph-input "$GRAPH_INPUT" \
  --spans "$SPANS" --signal-lengths "$SIGNAL_LENGTHS" --genome "$GENOME_REF" --output "$ES/tier1.tsv" 2>&1 | tee "$ES/logs/tier1.log" >&2
