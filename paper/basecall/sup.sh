#!/usr/bin/env bash
# paper/basecall/sup.sh KEY: the sup basecalls of a paper set as the paper made them (basecall/sets.tsv):
#   OUT_DIR/reads/KEY/sup.moves.bam   dorado basecaller --emit-moves on POD5_DIR/KEY.pod5: unaligned records with move
#                                     tables (the evaluation's input: spelled and placed sequences);
#   OUT_DIR/reads/KEY/sup.fastq.gz    the assemblers' reads: the BAM's primary records (samtools fastq -F 0x900) as
#                                     `@id`, sequence, `+`, qualities.
# Dorado and model per set: D1-D5, D8 and the bacterial / A. rabiei supplement sets Dorado 1.4.0 sup@v5.2.0 with read
# splitting off (--disable-read-splitting); D. melanogaster (D7, dmk60; 4 kHz, which Dorado 1.0 and later reject) Dorado
# 0.9.2 sup@v4.1.0 with its default read splitting (a split read's pieces carry their own ids and pi:Z:<read id>; the
# FASTQ keeps them, the BAM too); human (D6) the move tables from Dorado 1.4.0 sup@v5.0.0 (the model of ONT's basecalls)
# with read splitting off, and the reads from ONT's own sup v5.0.0 basecalls of the two flow cells (giab_2025.01
# calls.sorted.bam, data/download/manifests/hg002_giab2025_bams.tsv): their chr13 primary records (samtools view -F 0x904,
# region chr13, through the BAM index over HTTPS: ~4% of the 329 GB) of the set's reads, PAW70337 first, in the BAMs'
# order; the FASTQ's text is then checked against the paper's (sets.tsv `sup_fastq_sha256`). Dorado's records come in its
# output order. Dorado needs a GPU (its default device).
# STEP: moves (the BAM), fastq (the FASTQ; from the BAM, which must exist, except for human), default both. Steps already
# done (their output exists) are skipped. Settings: POD5_DIR, OUT_DIR, LOCAL_DIR (temporary files), BC_THREADS (samtools /
# pigz threads, 8), DORADO_* and DORADO_MODELS, SAMTOOLS (tools/env.sh), CURL. Logs in OUT_DIR/reads/KEY/logs/
# (sup.dorado.log, sup.time.txt); OUT_DIR/reads/KEY/sup.stats.tsv: reads, bases of the FASTQ.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
KEY=${1:?usage: sup.sh KEY [moves|fastq] (keys: column 1 of basecall/sets.tsv)} STEP=${2:-both}
[[ $STEP =~ ^(both|moves|fastq)$ ]] || { echo "sup.sh: STEP must be moves or fastq" >&2; exit 2; }
V=$(bc_col "$KEY" sup_dorado) || { echo "sup.sh: no set $KEY in $BC_HERE/sets.tsv" >&2; exit 2; }
MODEL=$(model_dir "$(bc_col "$KEY" sup_model)")
SPLIT=$(bc_col "$KEY" sup_split) FASTQ_FROM=$(bc_col "$KEY" sup_fastq)
DORADO=$(dorado_bin "$V")
O=$OUT_DIR/reads/$KEY
mkdir -p "$O/logs"
args=(basecaller --emit-moves)
[[ $SPLIT == off ]] && args+=(--disable-read-splitting)
if [[ $STEP != fastq && ! -s $O/sup.moves.bam ]]; then
  POD5=$(set_pod5 "$KEY")
  echo "sup.sh: $KEY: $DORADO ${args[*]} $MODEL $POD5"
  [[ $DRY_RUN != 1 ]] || exit 0
  /usr/bin/time -v -o "$O/logs/sup.time.txt" "$DORADO" "${args[@]}" "$MODEL" "$POD5" > "$O/sup.moves.bam.part" \
    2> "$O/logs/sup.dorado.log"
  mv -- "$O/sup.moves.bam.part" "$O/sup.moves.bam"
fi
[[ ! -s $O/sup.moves.bam ]] || echo "sup.sh: $KEY: $("$SAMTOOLS" view -c "$O/sup.moves.bam") records in $O/sup.moves.bam"
[[ $STEP != moves ]] || exit 0

if [[ ! -s $O/sup.fastq.gz ]]; then
  case $FASTQ_FROM in
    dorado)
      [[ -s $O/sup.moves.bam ]] || { echo "sup.sh: no $O/sup.moves.bam (sup.sh $KEY moves)" >&2; exit 1; }
      "$SAMTOOLS" fastq -@ "$BC_THREADS" -F 0x900 "$O/sup.moves.bam" 2> /dev/null | fastq_lines | gz_out "$O/sup.fastq.gz" ;;
    ont_giab2025)
      # ONT's basecalls: the chr13 primary records of the set's reads, flow cell by flow cell (the BAMs are sorted)
      mkdir -p "$LOCAL_DIR"; W=$(mktemp -d "$LOCAL_DIR/sup_$KEY.XXXXXX"); trap 'rm -rf -- "$W"' EXIT
      gzip -dc "$PAPER/data/sets/$KEY/reads.txt.gz" > "$W/ids.txt"
      for fc in PAW70337 PAW71238; do
        url=https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/$fc/calls.sorted.bam
        $CURL -o "$W/$fc.bam.bai" "$url.bai"
        echo "sup.sh: $KEY: chr13 primary records of $url"
        "$SAMTOOLS" view -@ "$BC_THREADS" -b -F 0x904 -X -o "$W/chr13.$fc.bam" "$url" "$W/$fc.bam.bai" chr13
      done
      "$SAMTOOLS" cat "$W/chr13.PAW70337.bam" "$W/chr13.PAW71238.bam" | "$SAMTOOLS" view -b -N "$W/ids.txt" - \
        | "$SAMTOOLS" fastq - 2> /dev/null | fastq_lines | gz_out "$O/sup.fastq.gz" ;;
    *) echo "sup.sh: unknown sup_fastq $FASTQ_FROM for $KEY" >&2; exit 2 ;;
  esac
fi
read -r n b < <(gzip -dc "$O/sup.fastq.gz" | awk 'NR % 4 == 2 {n++; b += length($0)} END {print n + 0, b + 0}')
printf 'reads\t%s\nbases\t%s\n' "$n" "$b" > "$O/sup.stats.tsv"
echo "sup.sh: $KEY: $n reads, $b bases in $O/sup.fastq.gz"
want=$(bc_col "$KEY" sup_fastq_sha256)
if [[ $want != . ]]; then   # human: ONT's basecalls, so the FASTQ is the paper's byte for byte
  got=$(gzip -dc "$O/sup.fastq.gz" | sha256sum | cut -d' ' -f1)
  [[ $got == "$want" ]] || { echo "sup.sh: $KEY: the FASTQ's text has SHA-256 $got, the paper's $want" >&2; exit 1; }
  echo "sup.sh: $KEY: the FASTQ equals the paper's (SHA-256 of its text $want)"
fi
