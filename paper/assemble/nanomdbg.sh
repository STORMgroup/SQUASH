#!/usr/bin/env bash
# paper/assemble/nanomdbg.sh KEY: nanoMDBG (metaMDBG 1.2) on a paper set's sup reads as the paper ran it (assemble/lib.sh's
# protocol), with minimap2 2.30-r1287 first on its PATH (MINIMAP2's folder):
#   metaMDBG asm --out-dir LOCAL/asm --in-ont LOCAL/reads.fastq --threads ASM_THREADS
# Input SUP_FASTQ (default OUT_DIR/reads/KEY/sup.fastq.gz). Output OUT_DIR/assembly/KEY/nanomdbg/: asm/contigs.fasta.gz,
# asm/metaMDBG.log; logs/nanomdbg.time.txt, logs/nanomdbg.log, logs/staging.time.txt; stats.tsv (contigs, bases, SHA-256
# of the contigs' text); done. metaMDBG is multi-threaded and not reproducible run to run (contig counts and most contigs
# agree; a few differ by some bases).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
as_pin "$0" "$@"
as_start "${1:?usage: nanomdbg.sh KEY}" nanomdbg
stage_reads "${SUP_FASTQ:-$OUT_DIR/reads/$KEY/sup.fastq.gz}"
PATH="$(dirname "$(command -v "$MINIMAP2")"):$PATH" timed nanomdbg "$METAMDBG" asm --out-dir "$LOCAL/asm" \
  --in-ont "$LOCAL/reads.fastq" --threads "$ASM_THREADS" > "${O:?}/logs/nanomdbg.log" 2>&1
[[ $DRY_RUN == 1 ]] && exit 0
cp "$LOCAL/asm/contigs.fasta.gz" "$O/asm/"; cp "$LOCAL/asm/metaMDBG.log" "$O/asm/" 2> /dev/null || true
gzip -cd "$O/asm/contigs.fasta.gz" | awk '/^>/ {n++} !/^>/ {b += length($0)} END {printf "contigs\t%d\nbases\t%d\n", n, b}' > "$O/stats.tsv"
printf 'sha256_contigs_text\t%s\n' "$(gzip -cd "$O/asm/contigs.fasta.gz" | sha256sum | cut -c1-64)" >> "$O/stats.tsv"
as_done
