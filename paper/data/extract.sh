#!/usr/bin/env bash
# paper/data/extract.sh KEY [KEY ...]: rebuild a paper set's POD5 from its public source and check it.
#   0. POD5_DIR/KEY.pod5 already there (a file or a link): only step 3;
#   1. the source's download script (download/SOURCE.sh, sets.tsv's `sources`), unless the source is complete or a
#      local copy is named (SRC_<SOURCE>, see download/lib.sh); streamed sources (human, D. melanogaster) keep only the
#      reads of their sets' lists;
#   2. pod5_cut.py: the reads of sets/KEY/reads.txt.gz into POD5_DIR/KEY.pod5, in the list's order (the order the
#      paper's run read), every run info first;
#   3. the check: the POD5's read ids in file order equal the list, and pod5_digest.py's signal digest equals
#      sets.tsv's `pod5_digest` (pod5_digest.py stops at the first read that differs from the list).
# Settings (paper/config.sh or the environment): SOURCES_DIR (default paper/data/sources), POD5_DIR (default
# paper/data/pod5), PYTHON (python with pod5 0.3.39), CUT_THREADS (cut workers, default 4), DROP_SOURCE=1 (delete a
# downloaded source after the cut; never a local copy). Exit non-zero on any failed step or check.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[[ -f $HERE/../config.sh ]] && source "$HERE/../config.sh"
SETS_TSV=${SETS_TSV:-$HERE/../sets.tsv}
SOURCES_DIR=${SOURCES_DIR:-$HERE/sources} POD5_DIR=${POD5_DIR:-$HERE/pod5}
PYTHON=${PYTHON:-python3} CUT_THREADS=${CUT_THREADS:-4}
export SOURCES_DIR PYTHON
(( $# >= 1 )) || { echo "usage: $0 KEY [KEY ...] (keys: column 1 of $SETS_TSV)" >&2; exit 2; }
"$PYTHON" -c 'import pod5, sys; sys.exit(0 if pod5.__version__ == "0.3.39" else 1)' 2>/dev/null \
  || echo "extract.sh: note: $PYTHON has no pod5 0.3.39 (the paper's POD5s were written with it)" >&2

col() {  # col KEY COLUMN: a value of sets.tsv
  awk -F'\t' -v k="$1" -v c="$2" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == k {print $(h[c]); f = 1}
    END {exit f ? 0 : 1}' "$SETS_TSV"
}
lists_of() {  # lists_of SOURCE: the read lists of every set that takes reads from SOURCE
  awk -F'\t' -v s="$1" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} index("," $(h["sources"]) ",", "," s ",") {print $1}' \
    "$SETS_TSV" | while read -r k; do echo "$HERE/sets/$k/reads.txt.gz"; done
}

for KEY in "$@"; do
  src=$(col "$KEY" sources) || { echo "extract.sh: no set $KEY in $SETS_TSV" >&2; exit 2; }
  digest=$(col "$KEY" pod5_digest) n=$(col "$KEY" reads)
  LIST=$HERE/sets/$KEY/reads.txt.gz OUT=$POD5_DIR/$KEY.pod5
  [[ -s $LIST ]] || { echo "extract.sh: no $LIST" >&2; exit 1; }
  echo "=== $KEY: $n reads from ${src//,/ + } ($(date '+%F %T'))"
  mkdir -p "$POD5_DIR"
  if [[ -s $OUT ]]; then
    echo "$OUT exists: checking it only (no download, no cut)"
  else
    dirs=()
    for s in ${src//,/ }; do
      case $s in
        hg002_giab2025) mapfile -t L < <(lists_of "$s"); bash "$HERE/download/$s.sh" "${L[@]}" ;;
        dmel_bkim2023) mapfile -t L < <(lists_of "$s"); KEEP="${L[*]}" bash "$HERE/download/$s.sh" ;;
        *) bash "$HERE/download/$s.sh" ;;
      esac
      dirs+=("$SOURCES_DIR/$s")
    done
    "$PYTHON" "$HERE/pod5_cut.py" cut "$LIST" "$OUT" "${dirs[@]}" --procs "$CUT_THREADS"
  fi
  if [[ $digest == . || -z $digest ]]; then
    echo "extract.sh: sets.tsv has no digest for $KEY: ids checked only" >&2
    "$PYTHON" "$HERE/pod5_digest.py" --ids "$LIST" "$OUT"
  else
    "$PYTHON" "$HERE/pod5_digest.py" --ids "$LIST" --expect "$digest" "$OUT"
  fi
  echo "=== $KEY: $OUT checked ($(date '+%F %T'))"
  if [[ ${DROP_SOURCE:-0} == 1 ]]; then
    for s in ${src//,/ }; do [[ -L $SOURCES_DIR/$s ]] || rm -rf -- "${SOURCES_DIR:?}/$s"; done
  fi
done
