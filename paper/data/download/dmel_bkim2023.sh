#!/usr/bin/env bash
# Drosophila melanogaster ISO1 (Bernard Kim, Petrov lab; ONT open data contribution melanogaster_bkim_2023.01):
# s3://ont-open-data/contrib/melanogaster_bkim_2023.01/flowcells/D.melanogaster.R1041.400bps/, six runs on flow cells
# FAV70669 and FAV81161 (R10.4.1, 400 bps, 4 kHz), 562 FAST5 files (245,841,648,088 bytes) + the six sequencing summaries
# (587,080,927 bytes; only the curation reads them); file list and sizes in manifests/dmel_bkim2023.tsv (S3 gives no
# MD5: sizes only). CC BY-NC 4.0 (ONT open data). Sets: dmk60_both60, dmk100_both100.
# Each FAST5 is downloaded, converted (`pod5 convert fast5 --threads 1`, pod5 0.3.39: the paper's conversion) and deleted:
#   default            -> SOURCES_DIR/dmel_bkim2023/RUN/NAME.pod5 (every read, ~244 GB)
#   KEEP=LIST [LIST..] -> SOURCES_DIR/dmel_bkim2023/RUN/NAME.pod5 holds only the listed reads (pod5_cut.py part; the
#                         full conversion deleted), so the disk peak is PAR files in flight + the kept reads
# A local copy instead: SRC_DMEL_BKIM2023=FOLDER (the converted runs). PAR=N files at once (4). SUMMARIES=1: also the
# sequencing summaries (-> SOURCES_DIR/dmel_bkim2023/summaries/).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
S=dmel_bkim2023
d=$(source_dir $S) || { log "$S: done before ($d)"; exit 0; }
keep=""
if [[ -n ${KEEP:-} ]]; then
  keep=$d/.keep_ids.txt
  for l in $KEEP; do if [[ $l == *.gz ]]; then gzip -dc "$l"; else cat "$l"; fi; done | sort -u > "$keep.tmp"
  mv "$keep.tmp" "$keep"; log "$S: keeping $(wc -l < "$keep") listed reads"
fi
export -f fetch log die; export CURL d keep PYTHON PAPER_DATA
one() {
  local key=$1 bytes=$2 run=$3 kind=$4 n out
  n=$(basename "$key"); n=${n%.fast5}
  [[ -e $d/.ok/$n ]] && return 0
  mkdir -p "$d/.ok" "$d/.dl" "$d/$run"
  if [[ $kind == summary ]]; then
    [[ ${SUMMARIES:-0} == 1 ]] || return 0
    mkdir -p "$d/summaries"; fetch "https://ont-open-data.s3.amazonaws.com/$key" "$d/summaries/$n" "$bytes"
    touch "$d/.ok/$n"; return 0
  fi
  fetch "https://ont-open-data.s3.amazonaws.com/$key" "$d/.dl/$n.fast5" "$bytes"
  rm -f -- "$d/.dl/$n.pod5"
  "$PYTHON" -m pod5.tools.main convert fast5 --threads 1 --output "$d/.dl/$n.pod5" "$d/.dl/$n.fast5" > "$d/.dl/$n.log" 2>&1 \
    || { cat "$d/.dl/$n.log" >&2; die "convert $n"; }
  rm -f -- "$d/.dl/$n.fast5"
  out=$d/$run/$n.pod5
  rm -f -- "$out"
  if [[ -n $keep ]]; then
    "$PYTHON" "$PAPER_DATA/pod5_cut.py" part "$keep" "$d/.dl/$n.pod5" "$out"
    rm -f -- "$d/.dl/$n.pod5"
  else
    mv "$d/.dl/$n.pod5" "$out"
  fi
  rm -f -- "$d/.dl/$n.log"; touch "$d/.ok/$n"
}
export -f one
tail -n +2 "$MANIFESTS/$S.tsv" | xargs -P "${PAR:-4}" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r k b r t <<< "{}"; one "$k" "$b" "$r" "$t"'
want=$(awk -F'\t' 'NR > 1 && $4 == "fast5"' "$MANIFESTS/$S.tsv" | wc -l)
got=$(awk -F'\t' 'NR > 1 && $4 == "fast5" {n = $1; sub(/.*\//, "", n); sub(/\.fast5$/, "", n); print n}' "$MANIFESTS/$S.tsv" \
  | while read -r n; do [[ -e $d/.ok/$n ]] && echo "$n"; done | wc -l)
(( got == want )) || die "$S: $got of $want FAST5 files converted"
rmdir "$d/.dl" 2>/dev/null || true
mark_complete "$d" $S
