#!/usr/bin/env bash
# paper/eval/compare_quast.sh KEY: QUAST 5.2.0 (defaults, --no-icarus) of the basecalled assemblers' runs on the set's
# sup reads (OUT_DIR/assembly/KEY/: nanoMDBG's contigs.fasta.gz, Flye's assembly.fasta, hifiasm's primary contigs
# asm.bp.p_ctg.fa; labels nanoMDBG, Flye, hifiasm) -> OUT_DIR/eval/KEY/compare/quast/ against the set's reference,
# quast_hap/ and quast_masked/ against its hap and masked references; hifiasm on the HAC reads (hifiasm_hac, the sets
# with hac = 1; label hifiasm (HAC)) -> quast_hac/ against the reference.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
set_paths "${1:?usage: compare_quast.sh KEY}"
NANO=$ASM/nanomdbg/asm/contigs.fasta.gz FLYE=$ASM/flye/asm/assembly.fasta HIFI=$ASM/hifiasm/asm/asm.bp.p_ctg.fa
HAC=$ASM/hifiasm_hac/asm/asm.bp.p_ctg.fa
need "$NANO" "$FLYE" "$HIFI" "$REF"
[[ $(ecol "$KEY" hac) != 1 ]] || need "$HAC"
mkdir -p "$EC/logs"
workdir "compare_$KEY"
gzip -cd "$NANO" > "$TMPW/nanoMDBG.fa"
run() {  # run OUT REF LABELS FILE ...
  local out=$1 ref=$2 labels=$3; shift 3
  log "QUAST $out ($(basename "$ref"))"
  rm -rf -- "$EC/.$out.tmp"
  "$QUAST" -r "$ref" -o "$EC/.$out.tmp" -t "$EVAL_THREADS" --no-icarus -l "$labels" "$@" > "$EC/logs/$out.log" 2>&1 \
    || { tail -20 "$EC/logs/$out.log" >&2; die "QUAST failed ($EC/logs/$out.log)"; }
  rm -rf -- "${EC:?}/$out"; mv "$EC/.$out.tmp" "$EC/$out"
  grep -E "^(Assembly|# contigs	|NGA50|# misassemblies|Genome fraction)" "$EC/$out/report.tsv" >&2
}
run quast "$REF" "nanoMDBG,Flye,hifiasm" "$TMPW/nanoMDBG.fa" "$FLYE" "$HIFI"
[[ $HAP == - ]] || run quast_hap "$EVAL_REFS/$HAP" "nanoMDBG,Flye,hifiasm" "$TMPW/nanoMDBG.fa" "$FLYE" "$HIFI"
[[ $MASKED == - ]] || run quast_masked "$EVAL_REFS/$MASKED" "nanoMDBG,Flye,hifiasm" "$TMPW/nanoMDBG.fa" "$FLYE" "$HIFI"
[[ $(ecol "$KEY" hac) != 1 ]] || run quast_hac "$REF" "hifiasm (HAC)" "$HAC"
