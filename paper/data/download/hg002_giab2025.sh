#!/usr/bin/env bash
# Human HG002 (ONT open data giab_2025.01): flow cells PAW70337 and PAW71238 (R10.4.1, 5 kHz, POD5 with pass and fail
# reads), s3://ont-open-data/giab_2025.01/flowcells/HG002/{PAW70337,PAW71238}/pod5/, 408 POD5 files, 3,284,066,513,176
# bytes; file list and sizes in manifests/hg002_giab2025.tsv (S3 gives no MD5: sizes only). Terms: ONT open data
# (sources.tsv). Set: h2k60_both60 (its reads, a chr13 q-arm cohort, are in
# paper/data/sets/h2k60_both60/reads.txt.gz; making that cohort needs ONT's sup BAMs, manifests/hg002_giab2025_bams.tsv,
# 0.33 TB).
# Streamed: each POD5 is downloaded, the listed reads cut out (pod5_cut.py part), and the file deleted, so the disk
# peak is PAR files in flight (~13 GB each) + the kept reads (~62 GB for h2k60).
# Usage: hg002_giab2025.sh LIST [LIST ...]   (reads.txt.gz of the set(s) to keep; required)
# -> SOURCES_DIR/hg002_giab2025/parts/NAME.pod5. A local copy of the run's POD5s instead: SRC_HG002_GIAB2025=FOLDER (the
# extraction then cuts from it directly). PAR=N files at once (4).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
S=hg002_giab2025
(( $# >= 1 )) || die "usage: $0 LIST [LIST ...] (the read lists of the sets to keep)"
d=$(source_dir $S) || { log "$S: done before ($d)"; exit 0; }
keep=$d/.keep_ids.txt
for l in "$@"; do if [[ $l == *.gz ]]; then gzip -dc "$l"; else cat "$l"; fi; done | sort -u > "$keep.tmp"
if [[ -s $keep ]] && ! cmp -s "$keep" "$keep.tmp"; then
  rm -f -- "$keep.tmp"; die "$d was cut for another read list: remove it (or use another SOURCES_DIR)"
fi
mv "$keep.tmp" "$keep"; log "$S: keeping $(wc -l < "$keep") listed reads"
export -f fetch log die; export CURL d keep PYTHON PAPER_DATA
one() {
  local key=$1 bytes=$2 n
  n=$(basename "$key" .pod5)
  [[ -e $d/.ok/$n ]] && return 0
  mkdir -p "$d/.ok" "$d/.dl" "$d/parts"
  fetch "https://ont-open-data.s3.amazonaws.com/$key" "$d/.dl/$n.pod5" "$bytes"
  rm -f -- "$d/parts/$n.pod5"
  "$PYTHON" "$PAPER_DATA/pod5_cut.py" part "$keep" "$d/.dl/$n.pod5" "$d/parts/$n.pod5"
  rm -f -- "$d/.dl/$n.pod5"; touch "$d/.ok/$n"
}
export -f one
tail -n +2 "$MANIFESTS/$S.tsv" | xargs -P "${PAR:-4}" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r k b f <<< "{}"; one "$k" "$b"'
want=$(( $(wc -l < "$MANIFESTS/$S.tsv") - 1 )) got=$(ls "$d/.ok" | wc -l)
(( got == want )) || die "$S: $got of $want POD5 files cut"
rmdir "$d/.dl" 2>/dev/null || true
mark_complete "$d" $S
