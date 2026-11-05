#!/usr/bin/env bash
# paper/assemble/hifiasm.sh KEY [sup|hac]: hifiasm 0.25.0-r726 (ONT mode) on a paper set's sup reads (default) or HAC
# reads as the paper ran it (assemble/lib.sh's protocol), in LOCAL/asm:
#   hifiasm --ont -t ASM_THREADS [-l0] -o asm LOCAL/reads.fastq
# -l0 (no duplicate purging, the hifiasm README's setting for haploid and inbred genomes) on every set but the human one
# (D6: defaults; assemble/options.tsv). Input SUP_FASTQ or HAC_FASTQ (default OUT_DIR/reads/KEY/sup.fastq.gz,
# hac.fastq.gz). Output OUT_DIR/assembly/KEY/hifiasm/ (sup) or hifiasm_hac/ (HAC): asm/asm.bp.p_ctg.fa (the primary
# contigs: the S lines of asm.bp.p_ctg.gfa as FASTA, what the evaluation scores), asm/asm.bp.p_ctg.gfa, the haplotype
# contigs asm.bp.hap1.p_ctg.fa / hap2 when hifiasm writes them, asm/*.noseq.gfa (the graphs without sequence), scaffold
# files when written; logs/hifiasm.time.txt, logs/hifiasm.log, logs/version.txt, logs/outputs.txt,
# logs/staging.time.txt; stats.tsv (version, options, contigs, bases and SHA-256 per contig set); done.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
as_pin "$0" "$@"
K=${1:?usage: hifiasm.sh KEY [sup|hac]} READS=${2:-sup}
case $READS in
  sup) as_start "$K" hifiasm; SRC=${SUP_FASTQ:-$OUT_DIR/reads/$K/sup.fastq.gz} ;;
  hac) as_start "$K" hifiasm_hac; SRC=${HAC_FASTQ:-$OUT_DIR/reads/$K/hac.fastq.gz} ;;
  *) echo "usage: hifiasm.sh KEY [sup|hac]" >&2; exit 2 ;;
esac
OPTS=(); o=$(as_col "$KEY" hifiasm_options); [[ $o == - ]] || read -r -a OPTS <<< "$o"
[[ $DRY_RUN == 1 ]] || "$HIFIASM" --version > "$O/logs/version.txt"
stage_reads "$SRC"
mkdir -p "$LOCAL/asm"
(cd "$LOCAL/asm" && timed hifiasm "$HIFIASM" --ont -t "$ASM_THREADS" "${OPTS[@]}" -o asm \
  "$LOCAL/reads.fastq" > "${O:?}/logs/hifiasm.log" 2>&1)
[[ $DRY_RUN == 1 ]] && exit 0
printf 'hifiasm\t%s\noptions\t--ont -t %s %s\n' "$(cat "$O/logs/version.txt")" "$ASM_THREADS" "${OPTS[*]:-}" > "$O/stats.tsv"
for g in p_ctg hap1.p_ctg hap2.p_ctg; do
  f=$LOCAL/asm/asm.bp.$g.gfa
  [[ -s $f ]] || continue
  awk -F'\t' '$1 == "S" {print ">" $2; print $3}' "$f" > "$O/asm/asm.bp.$g.fa"
  awk -F'\t' 'BEGIN {OFS = "\t"} $1 == "S" {$3 = "*"} {print}' "$f" > "$O/asm/asm.bp.$g.noseq.gfa"
  awk -v g="$g" '/^>/ {n++} !/^>/ {b += length($0)} END {printf "%s_contigs\t%d\n%s_bases\t%d\n", g, n, g, b}' "$O/asm/asm.bp.$g.fa" >> "$O/stats.tsv"
  printf '%s_sha256\t%s\n' "$g" "$(sha256sum < "$O/asm/asm.bp.$g.fa" | cut -c1-64)" >> "$O/stats.tsv"
done
cp "$LOCAL/asm/asm.bp.p_ctg.gfa" "$O/asm/" 2> /dev/null || true
ls -l "$LOCAL/asm" > "$O/logs/outputs.txt"
for f in "$LOCAL"/asm/*scaf*; do
  [[ -s $f ]] || continue
  case $f in
    *.gfa) awk -F'\t' '$1 == "S" {print ">" $2; print $3}' "$f" > "$O/asm/$(basename "${f%.gfa}").fa" ;;
    *.fa|*.fasta) cp "$f" "$O/asm/" ;;
  esac
done
[[ -s $O/asm/asm.bp.p_ctg.fa ]] || { echo "hifiasm.sh: no primary contigs" >&2; exit 1; }
as_done
