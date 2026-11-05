#!/usr/bin/env bash
# paper/basecall/hac.sh KEY: the HAC reads of a paper set for hifiasm's HAC run, as the paper made them: dorado basecaller
# --emit-fastq with Dorado's default read splitting on POD5_DIR/KEY.pod5, Dorado 2.1.2 hac@v5.2.0 (D. melanogaster, 4 kHz:
# Dorado 0.9.2 hac@v4.1.0; basecall/sets.tsv), the FASTQ as Dorado writes it -> OUT_DIR/reads/KEY/hac.fastq.gz (+
# hac.stats.tsv: reads, bases). The same command as the timed GPU HAC route (timed.sh KEY hac gpu), here untimed.
# Needs a GPU. Settings: POD5_DIR, OUT_DIR, LOCAL_DIR (the FASTQ before compression), BC_THREADS, DORADO_*,
# DORADO_MODELS. Logs in OUT_DIR/reads/KEY/logs/ (hac.dorado.log, hac.time.txt).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
KEY=${1:?usage: hac.sh KEY (keys: column 1 of basecall/sets.tsv)}
V=$(bc_col "$KEY" hac_dorado) || { echo "hac.sh: no set $KEY in $BC_HERE/sets.tsv" >&2; exit 2; }
MODEL=$(model_dir "$(bc_col "$KEY" hac_model)")
DORADO=$(dorado_bin "$V")
POD5=$(set_pod5 "$KEY")
O=$OUT_DIR/reads/$KEY
mkdir -p "$O/logs"
if [[ ! -s $O/hac.fastq.gz ]]; then
  echo "hac.sh: $KEY: $DORADO basecaller --emit-fastq $MODEL $POD5"
  [[ $DRY_RUN != 1 ]] || exit 0
  mkdir -p "$LOCAL_DIR"; W=$(mktemp -d "$LOCAL_DIR/hac_$KEY.XXXXXX"); trap 'rm -rf -- "$W"' EXIT
  /usr/bin/time -v -o "$O/logs/hac.time.txt" "$DORADO" basecaller --emit-fastq "$MODEL" "$POD5" > "$W/hac.fastq" \
    2> "$O/logs/hac.dorado.log"
  gz_out "$O/hac.fastq.gz" < "$W/hac.fastq"
fi
read -r n b < <(gzip -dc "$O/hac.fastq.gz" | awk 'NR % 4 == 2 {n++; b += length($0)} END {print n + 0, b + 0}')
printf 'reads\t%s\nbases\t%s\n' "$n" "$b" > "$O/hac.stats.tsv"
echo "hac.sh: $KEY: $n reads, $b bases in $O/hac.fastq.gz"
