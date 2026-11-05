#!/usr/bin/env bash
# paper/eval/compare_tier1.sh KEY: the comparators' tier-1 contiguity (tier1.py compare) -> OUT_DIR/eval/KEY/compare/
# tier1.tsv, one row per tool: Rawsamble + miniasm (the sets with rawsamble = 1: its GFA's read layouts in
# signal-estimated bases against twice the reference length, through the reads' events without event merging, which
# SQUASH's signal step writes with --event-merge 0), nanoMDBG, Flye, hifiasm (basecalled bases against the reference
# length).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: compare_tier1.sh KEY}"
NANO=$ASM/nanomdbg/asm/contigs.fasta.gz FLYE=$ASM/flye/asm HIFI=$ASM/hifiasm/asm/asm.bp.p_ctg.fa
RSG=$ASM/rawsamble/asm/ava.gfa
need "$NANO" "$FLYE/assembly.fasta" "$HIFI" "$GENOME_REF" "$SIGNAL_LENGTHS"
mkdir -p "$EC/logs"
workdir "ctier1_$KEY"
rs=()
if [[ $(ecol "$KEY" rawsamble) == 1 ]]; then
  need "$RSG" "$POD5" "$SQUASH"
  log "events without event merging of Rawsamble's reads"
  awk -F'\t' '$1 == "a" {sub(/:[0-9]+-[0-9]+$/, "", $4); print $4}' "$RSG" | sort -u > "$TMPW/reads.txt"
  "$SQUASH" -x small --event-merge 0 --stop-after signal --read-ids "$TMPW/reads.txt" -t "$EVAL_THREADS" \
    -o "$TMPW/signal" "$POD5" > "$EC/logs/spans_raw.log" 2>&1 || { tail "$EC/logs/spans_raw.log" >&2; die "squash failed"; }
  rs=(--rawsamble "$RSG" --spans-raw "$TMPW/signal/work/signal/spans.bin")
fi
log "tier 1"
"$REMORA_PYTHON" "$EVAL/tier1.py" compare --key "$KEY" --genome "$GENOME_REF" --signal-lengths "$SIGNAL_LENGTHS" \
  --nanomdbg "$NANO" --flye "$FLYE" --hifiasm "$HIFI" ${rs[@]+"${rs[@]}"} --output "$EC/tier1.tsv" 2>&1 | tee "$EC/logs/tier1.log" >&2
