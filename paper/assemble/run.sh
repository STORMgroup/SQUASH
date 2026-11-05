#!/usr/bin/env bash
# paper/assemble/run.sh KEY [TOOL ...]: the other assemblers on a paper set, one after another, each by its script:
# TOOL nanomdbg (nanomdbg.sh), flye (flye.sh), hifiasm (hifiasm.sh KEY sup), hifiasm_hac (hifiasm.sh KEY hac), rawsamble
# (rawsamble.sh). Default: the tools the paper ran on that set (assemble/options.tsv `paper_tools`: hifiasm on HAC reads
# for D1-D8 only, Rawsamble + miniasm on every set but D. melanogaster and A. thaliana). KEY may be main, supp or all.
# For the paper's timing, run it on a machine of its own with ASM_THREADS (64) CPUs, one tool at a time (as here).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
(( $# )) || { echo "usage: run.sh KEY|main|supp|all [TOOL ...]" >&2; exit 2; }
SEL=$1; shift
mapfile -t keys < <(awk -F'\t' -v s="$SEL" 'NR > 1 && (s == "all" || $1 == s || (s == "main" && $2 != "supp") ||
  (s == "supp" && $2 == "supp")) {print $1}' "$PAPER/sets.tsv")
(( ${#keys[@]} )) || { echo "run.sh: no set $SEL in $PAPER/sets.tsv" >&2; exit 2; }
for k in "${keys[@]}"; do
  if (( $# )); then tools=("$@"); else
    IFS=, read -r -a tools < <(awk -F'\t' -v k="$k" '$1 == k {print $4}' "$HERE/options.tsv"); fi
  for t in "${tools[@]}"; do
    case $t in
      nanomdbg|flye|rawsamble) bash "$HERE/$t.sh" "$k" ;;
      hifiasm) bash "$HERE/hifiasm.sh" "$k" sup ;;
      hifiasm_hac) bash "$HERE/hifiasm.sh" "$k" hac ;;
      *) echo "run.sh: unknown tool $t (nanomdbg, flye, hifiasm, hifiasm_hac, rawsamble)" >&2; exit 2 ;;
    esac
  done
done
