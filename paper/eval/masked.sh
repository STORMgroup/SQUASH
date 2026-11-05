#!/usr/bin/env bash
# paper/eval/masked.sh KEY: SQUASH's evaluation on the set's masked reference (eval/sets.tsv's masked column: D.
# melanogaster's repeat-rich regions, A. thaliana's tandem arrays): the placement redone on it (place.sh KEY masked),
# the strands split on it (strands.sh KEY masked), QUAST (quast.sh KEY masked) -> OUT_DIR/eval/KEY/squash/masked/ and
# squash/quast_masked/.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
bash "$HERE/place.sh" "${1:?usage: masked.sh KEY}" masked
bash "$HERE/strands.sh" "$1" masked
bash "$HERE/quast.sh" "$1" masked
