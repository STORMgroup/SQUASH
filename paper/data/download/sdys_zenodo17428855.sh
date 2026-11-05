#!/usr/bin/env bash
# S. dysgalactiae subsp. equisimilis UT9728 and UT10237 (Schababerle et al. 2025): Zenodo 10.5281/zenodo.17428855, the
# MinION run FBA61885 (SQK-RBK114-24, 2025-01-11, every read: POD5s not split by pass / fail), 95 files (88 POD5, 7 zips
# of two POD5 each), 199,930,325,908 bytes; MD5s as Zenodo lists them (manifests/sdys_zenodo17428855.tsv). CC BY-SA 4.0.
# Sets: u9k60_both60, u9k100_both100 (UT9728; their reads, from the paper's demultiplexing, are in
# paper/data/sets/KEY/reads.txt.gz).
# -> SOURCES_DIR/sdys_zenodo17428855/*.pod5. A local copy instead: SRC_SDYS_ZENODO17428855=FOLDER. PAR=N files at once (4).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
S=sdys_zenodo17428855
d=$(source_dir $S) || { log "$S: done before ($d)"; exit 0; }
export -f fetch log die; export CURL d
tail -n +2 "$MANIFESTS/$S.tsv" | xargs -P "${PAR:-4}" -d '\n' -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r f m b u <<< "{}"
  [[ -e $d/.ok/$f ]] && exit 0
  mkdir -p "$d/.ok" "$d/.dl"
  if [[ $f == *.zip ]]; then
    fetch "$u" "$d/.dl/$f" "$b" "$m"; (cd "$d" && unzip -o -q "$d/.dl/$f"); rm -f -- "$d/.dl/$f"
  else
    fetch "$u" "$d/$f" "$b" "$m"
  fi
  touch "$d/.ok/$f"'
n=$(ls "$d/.ok" | wc -l); want=$(( $(wc -l < "$MANIFESTS/$S.tsv") - 1 ))
(( n == want )) || die "$S: $n of $want files"
rmdir "$d/.dl" 2>/dev/null || true
mark_complete "$d" $S
