#!/usr/bin/env bash
# paper/tables/make.sh [NAME ...]: the paper's tables and data figures from the outputs under OUT_DIR (config.sh):
# OUT_DIR/tables/NAME.tex (bare table bodies, as the paper inputs them), OUT_DIR/tables/values.tsv (every number with
# its source file) and OUT_DIR/figures/runtime.pdf. NAME = a script of this folder without .py (e.g.
# datasets, fig_runtime); default: every table, then the figures. A value whose input is missing prints as ?? and
# its script lists the missing files; make.sh carries on and exits 1 if a script failed. PYTHON (config.sh) needs
# matplotlib for the figures (requirements.txt).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
ALL=(datasets contiguity_both correctness_both
     supp_commands supp_presets supp_versions supp_cont30 supp_corr30 supp_athal_masked supp_dmel_masked supp_hac
     supp_runtime_single supp_runtime supp_strands
     fig_runtime)
names=("$@")
(( ${#names[@]} )) || names=("${ALL[@]}")
failed=()
for n in "${names[@]}"; do
  [[ -f $HERE/$n.py ]] || { echo "make.sh: no script $HERE/$n.py" >&2; exit 2; }
  "$PYTHON" "$HERE/$n.py" || failed+=("$n")
done
echo "make.sh: tables in $OUT_DIR/tables, figures in $OUT_DIR/figures" >&2
if (( ${#failed[@]} )); then echo "make.sh: failed: ${failed[*]}" >&2; exit 1; fi
