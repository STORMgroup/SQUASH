#!/usr/bin/env bash
# paper/basecall/run.sh KEY [KEY ...]: the basecalls of paper sets: sup.sh (move-table BAM + the assemblers' FASTQ) and
# hac.sh (hifiasm's HAC reads) for each KEY (a key of basecall/sets.tsv, or main = D1-D8, supp, all). TIMED=1 also runs
# the paper's three timed routes (timed.sh KEY sup gpu, KEY hac gpu, KEY hac cpu); TIMED_ONLY=1 runs only those. On one
# machine with a GPU (the CPU route wants BASECALL_CPU_THREADS CPUs). Stops at the first failure.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/../config.sh"
(( $# )) || { echo "usage: run.sh KEY|main|supp|all [KEY ...]" >&2; exit 2; }
keys=()
for a in "$@"; do
  case $a in
    main|supp|all) mapfile -t -O "${#keys[@]}" keys < <(awk -F'\t' -v s="$a" 'NR > 1 && (s == "all" ||
                     (s == "main" && $2 != "supp") || (s == "supp" && $2 == "supp")) {print $1}' "$PAPER/sets.tsv") ;;
    *) keys+=("$a") ;;
  esac
done
for k in "${keys[@]}"; do
  echo "=== basecall/run.sh: $k ($(date '+%F %T'))"
  if [[ ${TIMED_ONLY:-0} != 1 ]]; then
    bash "$HERE/sup.sh" "$k"
    bash "$HERE/hac.sh" "$k"
  fi
  if [[ ${TIMED:-0} == 1 || ${TIMED_ONLY:-0} == 1 ]]; then
    bash "$HERE/timed.sh" "$k" sup gpu
    bash "$HERE/timed.sh" "$k" hac gpu
    bash "$HERE/timed.sh" "$k" hac cpu
  fi
done
