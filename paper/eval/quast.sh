#!/usr/bin/env bash
# paper/eval/quast.sh KEY [masked]: QUAST 5.2.0 (defaults, --no-icarus) of SQUASH's spelled and placed unitigs and their
# per-strand splits (labels SQUASH (spelled), ...) -> OUT_DIR/eval/KEY/squash/quast/ against the set's reference,
# quast_hap/ against its hap reference (human: both haplotypes' q-arms); `masked`: the masked placement and split
# (squash/masked/) against the masked reference -> quast_masked/.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: quast.sh KEY [masked]}"
run() {  # run OUT REF DIR
  local files; read -r -a files <<< "$(sq_files "$3")"
  need "$2" "${files[@]}"
  log "QUAST $1 ($(basename "$2"))"
  rm -rf -- "$ES/.$1.tmp"
  "$QUAST" -r "$2" -o "$ES/.$1.tmp" -t "$EVAL_THREADS" --no-icarus -l "$SQ_LABELS" "${files[@]}" > "$ES/logs/$1.log" 2>&1 \
    || { tail -20 "$ES/logs/$1.log" >&2; die "QUAST failed ($ES/logs/$1.log)"; }
  rm -rf -- "${ES:?}/$1"; mv "$ES/.$1.tmp" "$ES/$1"
  grep -E "^(Assembly|# contigs	|NGA50|# misassemblies|Genome fraction)" "$ES/$1/report.tsv" | cut -c1-200 >&2
}
mkdir -p "$ES/logs"
if [[ ${2:-} == masked ]]; then
  [[ $MASKED != - ]] || die "$KEY has no masked reference"
  run quast_masked "$EVAL_REFS/$MASKED" "$ES/masked"
else
  run quast "$REF" "$ES"
  [[ $HAP == - ]] || run quast_hap "$EVAL_REFS/$HAP" "$ES"
fi
