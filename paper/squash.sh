#!/usr/bin/env bash
# paper/squash.sh KEY [SQUASH OPTION ...]: SQUASH on a paper set as the paper ran it: the set's preset (sets.tsv),
# --read-ids for the sets with an assembled list (data/sets/KEY/assembled.txt.gz), THREADS threads. Reads
# POD5_DIR/KEY.pod5 (data/extract.sh KEY makes it), writes RUN_OUT (default OUT_DIR/KEY, see config.sh) with
# RUN_THREADS threads (default THREADS); KEEP_WORK=1 keeps the work folder. Timed by /usr/bin/time -v when it exists
# (RUN_OUT.time.txt). Extra options go to squash as they are (e.g. --tmp-dir DIR).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/config.sh"
KEY=${1:?usage: squash.sh KEY [SQUASH OPTION ...] (keys: column 1 of sets.tsv)}
shift
preset=$(awk -F'\t' -v k="$KEY" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == k {print $(h["preset"])}' \
  "$PAPER/sets.tsv")
[[ -n $preset ]] || { echo "squash.sh: no set $KEY in $PAPER/sets.tsv" >&2; exit 2; }
pod5=$POD5_DIR/$KEY.pod5
[[ -s $pod5 ]] || { echo "squash.sh: no $pod5 (make -C paper data SET=$KEY)" >&2; exit 1; }
[[ -x $SQUASH ]] || { echo "squash.sh: no squash binary at $SQUASH (make at the repository's root)" >&2; exit 1; }
RUN_OUT=${RUN_OUT:-$OUT_DIR/$KEY}
args=(-x "$preset" -t "${RUN_THREADS:-$THREADS}" -o "$RUN_OUT")
ids=$PAPER/data/sets/$KEY/assembled.txt.gz
[[ -s $ids ]] && args+=(--read-ids "$ids")
[[ -n $TMP_DIR ]] && args+=(--tmp-dir "$TMP_DIR/$KEY")
[[ ${KEEP_WORK:-0} == 1 ]] && args+=(--keep-work)
mkdir -p "$(dirname "$RUN_OUT")"
timer=()
[[ -x /usr/bin/time ]] && timer=(/usr/bin/time -v -o "$RUN_OUT.time.txt")
echo "squash.sh: ${timer[*]} $SQUASH ${args[*]} $* $pod5" >&2
${timer[@]+"${timer[@]}"} "$SQUASH" "${args[@]}" "$@" "$pod5"
