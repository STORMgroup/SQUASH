# paper/eval/lib.sh: sourced by every script of the evaluation: paper/config.sh, the evaluation's tools and settings,
# a set's row of eval/sets.tsv and the paths of one set's inputs and outputs (set_paths KEY).
EVAL=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$EVAL/../config.sh"
# The evaluation's tools come from tools/env.sh (sourced by config.sh): QUAST (QUAST 5.2.0), MINIMAP2 (minimap2
# 2.30-r1287), SAMTOOLS, REMORA_PYTHON (a python with Remora 3.3.0, pysam, numpy and pod5: every Python step here) and
# REMORA_LEVELS (RawHash2's R10.4.1 9-mer levels, Remora's k-mer table).
# Threads of the evaluation's steps; its large temporary files (default OUT_DIR/eval/tmp); the evaluation references
# (refs.sh writes them).
: "${EVAL_THREADS:=$THREADS}"
: "${EVAL_TMP:=$OUT_DIR/eval/tmp}"
: "${EVAL_REFS:=$SOURCES_DIR/references/eval}"
export EVAL_THREADS EVAL_TMP EVAL_REFS LC_ALL=C PYTHONDONTWRITEBYTECODE=1

die() { echo "${0##*/}: $*" >&2; exit 1; }
log() { echo "=== ${KEY:-} ${0##*/} $* ($(date '+%F %T'))" >&2; }

ecol() {  # ecol KEY COLUMN: a value of eval/sets.tsv ('-' = none)
  awk -F'\t' -v k="$1" -v c="$2" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == k {print $(h[c]); f = 1}
    END {exit f ? 0 : 1}' "$EVAL/sets.tsv"
}

# set_paths KEY: the set's inputs (SQUASH's outputs and kept work folder, the sup move records, the POD5, the
# comparators' assemblies) and the evaluation's folders.
set_paths() {
  KEY=$1
  ecol "$KEY" reference > /dev/null || die "no set $KEY in $EVAL/sets.tsv"
  SQ=$OUT_DIR/$KEY                                   # SQUASH's outputs (squash.sh)
  GRAPH_INPUT=$SQ/work/graph/input                   # its graph's input reads (KEEP_WORK=1)
  SPANS=$SQ/work/signal/spans.bin                    # its events' POD5 samples (KEEP_WORK=1)
  SIGNAL_LENGTHS=$SQ/work/bp_scale/signal_lengths.tsv
  LAYOUT=$SQ/layout.tsv.gz GFA=$SQ/assembly_graph.gfa
  MOVES=$OUT_DIR/reads/$KEY/sup.moves.bam            # dorado sup with move tables
  POD5=$POD5_DIR/$KEY.pod5
  ASM=$OUT_DIR/assembly/$KEY                         # the comparators: ASM/TOOL/asm/
  E=$OUT_DIR/eval/$KEY
  ES=$E/squash EC=$E/compare
  local r p g
  r=$(ecol "$KEY" reference) p=$(ecol "$KEY" placement) g=$(ecol "$KEY" genome)
  REF=$EVAL_REFS/$r                                  # QUAST's reference and the strand split's
  PLACE_REF=$EVAL_REFS/$([[ $p == - ]] && echo "$r" || echo "$p")    # the placement's
  GENOME_REF=$EVAL_REFS/$([[ $g == - ]] && echo "$r" || echo "$g")   # tier 1's genome length
  HAP=$(ecol "$KEY" hap) MASKED=$(ecol "$KEY" masked)
  export KEY SQ GRAPH_INPUT SPANS SIGNAL_LENGTHS LAYOUT GFA MOVES POD5 ASM E ES EC REF PLACE_REF GENOME_REF HAP MASKED
}

need() {  # need FILE ...: each file exists and is not empty
  local f
  for f in "$@"; do [[ -s $f ]] || die "missing $f"; done
}

squash_inputs() {  # what the steps on SQUASH's assembly read
  [[ -s $SPANS && -s $GRAPH_INPUT/token_spans.u64le ]] ||
    die "$SQ has no work folder: run SQUASH with it kept (KEEP_WORK=1 bash squash.sh $KEY)"
  need "$LAYOUT" "$GFA" "$GRAPH_INPUT/reads.tsv.gz" "$GRAPH_INPUT/tokens.u32le" "$SIGNAL_LENGTHS"
}

# QUAST labels of SQUASH's files, in the paper's evaluation's order (QUAST 5.2.0 defaults, --no-icarus)
SQ_LABELS="SQUASH (spelled),SQUASH (placed),SQUASH (spelled; forward-strand unitigs),SQUASH (spelled; reverse-strand unitigs),SQUASH (placed; forward-strand unitigs),SQUASH (placed; reverse-strand unitigs)"
sq_files() {  # sq_files DIR: SQUASH's six FASTA files of an evaluation folder (asm/, strands/), in SQ_LABELS' order
  echo "$1/asm/squash_spelled.fa" "$1/asm/squash_placed.fa" "$1/strands/squash_spelled_fwd.fa" \
    "$1/strands/squash_spelled_rev.fa" "$1/strands/squash_placed_fwd.fa" "$1/strands/squash_placed_rev.fa"
}

workdir() {  # workdir NAME: a new temporary folder under EVAL_TMP (removed at exit)
  mkdir -p "$EVAL_TMP"
  TMPW=$(mktemp -d "$EVAL_TMP/$1.XXXXXX")
  trap 'rm -rf -- "$TMPW"' EXIT
}
