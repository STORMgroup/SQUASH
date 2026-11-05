#!/usr/bin/env bash
# paper/eval/place.sh KEY [masked]: SQUASH's unitigs placed on the reference: the events of the layout's reads given
# reference positions by Remora's signal refinement of their sup alignments (remora_evref.py --collapse-check: the sup
# reads with their move tables aligned by minimap2 -y -Y -ax map-ont --secondary=no --MD, Remora 3.3.0 with RawHash2's
# R10.4.1 9-mer levels), then each layout piece written as the reference sequence its events lie on (project.py).
# Reference: eval/sets.tsv's placement column (default its reference); `masked`: the set's masked reference, into
# squash/masked/. -> OUT_DIR/eval/KEY/squash[/masked]/asm/squash_placed.fa (+ .blocks.tsv), reads.txt (the layout's
# reads) and evref/ (KEY.segs.tsv, KEY.summary.txt and the logs).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: place.sh KEY [masked]}"
D=$ES R=$PLACE_REF
if [[ ${2:-} == masked ]]; then
  [[ $MASKED != - ]] || die "$KEY has no masked reference"
  D=$ES/masked R=$EVAL_REFS/$MASKED
fi
squash_inputs
need "$MOVES" "$POD5" "$R" "$REMORA_LEVELS"
mkdir -p "$D/asm" "$D/evref"
workdir "place_$KEY"
log "placement on $(basename "$R")"
gzip -cdf "$LAYOUT" | tail -n +2 | cut -f5 | grep -v '^-$' | sort -u > "$D/reads.txt"
awk '/^>/ {print $1; next} {print}' "$R" > "$TMPW/reference.fa"
"$SAMTOOLS" faidx "$TMPW/reference.fa"
"$REMORA_PYTHON" -c 'import sys; sys.path.insert(0, sys.argv[1]); import common
want = {l.strip() for l in open(sys.argv[2]) if l.strip()}
n = common.cut_spans(sys.argv[3], want, sys.argv[4])
print(f"{n} of {len(want)} layout reads have spans", file=sys.stderr)' "$EVAL" "$D/reads.txt" "$SPANS" "$TMPW/spans.bin"
"$SAMTOOLS" view -b -F 0x900 -N "$D/reads.txt" -o "$TMPW/moves.bam" "$MOVES"
"$SAMTOOLS" fastq -T mv,ts,ns,sp,pi "$TMPW/moves.bam" 2> /dev/null \
  | "$MINIMAP2" -t "$EVAL_THREADS" -y -Y -ax map-ont --secondary=no --MD "$TMPW/reference.fa" - 2> "$D/evref/mm2.log" \
  | "$SAMTOOLS" sort -@ 4 -o "$TMPW/aln.bam" -
"$SAMTOOLS" index "$TMPW/aln.bam"
TMPDIR=$TMPW "$REMORA_PYTHON" "$EVAL/remora_evref.py" --pod5 "$POD5" --bam "$TMPW/aln.bam" --spans "$TMPW/spans.bin" \
  --reference "$TMPW/reference.fa" --levels "$REMORA_LEVELS" --mode remora --collapse-check --out "$TMPW/$KEY" \
  --threads "$EVAL_THREADS" --samtools "$SAMTOOLS" > "$D/evref/remora_evref.log" 2>&1 \
  || { tail -20 "$D/evref/remora_evref.log" >&2; die "remora_evref.py failed"; }
cp "$TMPW/$KEY.segs.tsv" "$TMPW/$KEY.summary.txt" "$D/evref/"
cat "$D/evref/$KEY.summary.txt" >&2
rm -f "$D/asm/squash_placed.fa"
"$REMORA_PYTHON" "$EVAL/project.py" --layout "$LAYOUT" --evref "$TMPW/$KEY.evref.bin" --segs "$TMPW/$KEY.segs.tsv" \
  --reference "$TMPW/reference.fa" --output "$D/asm/squash_placed.fa" >&2
