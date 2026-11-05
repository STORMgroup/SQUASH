#!/usr/bin/env bash
# paper/assemble/flye.sh KEY: Flye 2.9.6-b1802 on a paper set's sup reads as the paper ran it (assemble/lib.sh's
# protocol), with Flye's folder first on PATH:
#   flye --nano-hq LOCAL/reads.fastq -g SIZE -t ASM_THREADS -o LOCAL/asm
# SIZE per set (assemble/options.tsv): the genome length in kb the paper's runs passed (5130k for E. coli; elsewhere G of
# sets.tsv in whole kb). Flye 2.9.6 uses -g only with --asm-coverage, so it does not change the assembly.
# Input SUP_FASTQ (default OUT_DIR/reads/KEY/sup.fastq.gz). Output OUT_DIR/assembly/KEY/flye/: asm/assembly.fasta,
# asm/assembly_graph.gfa, asm/assembly_info.txt, asm/flye.log; logs/flye.time.txt, logs/flye.log, logs/staging.time.txt;
# stats.tsv (Flye's length summary, SHA-256 of assembly.fasta); done.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
as_pin "$0" "$@"
as_start "${1:?usage: flye.sh KEY}" flye
GS=$(as_col "$KEY" flye_genome_size)
stage_reads "${SUP_FASTQ:-$OUT_DIR/reads/$KEY/sup.fastq.gz}"
PATH="$(dirname "$(command -v "$FLYE")"):$PATH" timed flye "$FLYE" --nano-hq "$LOCAL/reads.fastq" -g "$GS" \
  -t "$ASM_THREADS" -o "$LOCAL/asm" > "${O:?}/logs/flye.log" 2>&1
[[ $DRY_RUN == 1 ]] && exit 0
for f in assembly.fasta assembly_graph.gfa assembly_info.txt flye.log; do cp "$LOCAL/asm/$f" "$O/asm/" 2> /dev/null || true; done
grep -E "Total length|Fragments|Fragments N50|Largest frg" "$O/asm/flye.log" | tail -4 > "$O/stats.tsv" || true
printf 'sha256_assembly\t%s\n' "$(sha256sum < "$O/asm/assembly.fasta" | cut -c1-64)" >> "$O/stats.tsv"
as_done
