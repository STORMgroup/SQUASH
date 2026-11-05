#!/usr/bin/env bash
# paper/eval/eval.sh KEY [STEP ...]: the paper's evaluation of one set, step by step (each skipped when its output is
# there; FORCE=1 runs it again). Steps, in order (default: all; `squash` = the first seven, `compare` = the last two):
#   refs           the evaluation references (refs.sh)                    -> EVAL_REFS/
#   spell          SQUASH's spelled unitigs (spell.sh)                    -> eval/KEY/squash/asm/squash_spelled.fa
#   place          SQUASH's placed unitigs (place.sh)                     -> eval/KEY/squash/asm/squash_placed.fa
#   strands        strand split and closed replicons (strands.sh)         -> eval/KEY/squash/strands/, closed.tsv
#   quast          QUAST of SQUASH's files (quast.sh)                     -> eval/KEY/squash/quast/ (+ quast_hap/)
#   tier1          SQUASH's tier 1 (tier1.sh)                             -> eval/KEY/squash/tier1.tsv
#   masked         the masked reference's placement and QUAST (masked.sh) -> eval/KEY/squash/quast_masked/
#   compare_quast  QUAST of the comparators (compare_quast.sh)            -> eval/KEY/compare/quast/ (+ hap, masked, hac)
#   compare_tier1  the comparators' tier 1 (compare_tier1.sh)             -> eval/KEY/compare/tier1.tsv
# Needs SQUASH's run of the set with its work folder (KEEP_WORK=1 bash squash.sh KEY), the sup move records
# (OUT_DIR/reads/KEY/sup.moves.bam), the set's POD5, and for the comparator steps their assemblies (OUT_DIR/assembly/KEY/).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/lib.sh"
set_paths "${1:?usage: eval.sh KEY [STEP ...]}"
shift
ALL=(refs spell place strands quast tier1 masked compare_quast compare_tier1)
steps=()
for s in "${@:-all}"; do
  case $s in
    all) steps+=("${ALL[@]}") ;;
    squash) steps+=("${ALL[@]:0:7}") ;;
    compare) steps+=("${ALL[@]:7}") ;;
    *) [[ " ${ALL[*]} " == *" $s "* ]] || die "no step $s (steps: ${ALL[*]}, squash, compare, all)"; steps+=("$s") ;;
  esac
done
done_file() {  # the output that marks a step done
  case $1 in
    refs) echo "-" ;; spell) echo "$ES/asm/squash_spelled.fa" ;; place) echo "$ES/asm/squash_placed.fa" ;;
    strands) echo "$ES/strands/summary.tsv" ;; quast) [[ $HAP == - ]] && echo "$ES/quast/report.tsv" || echo "$ES/quast_hap/report.tsv" ;;
    tier1) echo "$ES/tier1.tsv" ;; masked) echo "$ES/quast_masked/report.tsv" ;;
    compare_quast) [[ $(ecol "$KEY" hac) == 1 ]] && echo "$EC/quast_hac/report.tsv" || echo "$EC/quast/report.tsv" ;;
    compare_tier1) echo "$EC/tier1.tsv" ;;
  esac
}
for s in "${steps[@]}"; do
  [[ $s == masked && $MASKED == - ]] && continue
  f=$(done_file "$s")
  if [[ ${FORCE:-0} != 1 && $f != - && -s $f ]]; then echo "eval.sh: $KEY $s done ($f)" >&2; continue; fi
  case $s in
    refs) bash "$HERE/refs.sh" "$KEY" ;;
    *) bash "$HERE/$s.sh" "$KEY" ;;
  esac
done
echo "eval.sh: $KEY: ${steps[*]} done ($E)" >&2
