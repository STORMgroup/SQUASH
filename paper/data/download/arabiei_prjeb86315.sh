#!/usr/bin/env bash
# Ascochyta rabiei ME14, seven strains (Curtin University CCDM, ENA PRJEB86315, "Boscalid - a new selectable marker for
# Ascochyta lentis and Ascochyta rabiei"): runs ERR14664019-ERR14664025, one POD5 tarball per strain (PromethION flow
# cell PAY64335, SQK-NBD114-24, pass and fail reads), 167,001,284,436 bytes; MD5s as ENA lists them
# (manifests/arabiei_prjeb86315.tsv, ENA filereport submitted_md5). Sets: ark60_both60, ark100_both100 (every strain pooled).
# One connection per file (concurrent ranged ENA transfers can give corrupt files; the MD5 check catches them).
# -> SOURCES_DIR/arabiei_prjeb86315/STRAIN/... (each tarball unpacked into its own folder). A local copy instead:
# SRC_ARABIEI_PRJEB86315=FOLDER (holding the seven unpacked tarballs).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
S=arabiei_prjeb86315
d=$(source_dir $S) || { log "$S: done before ($d)"; exit 0; }
while IFS=$'\t' read -r acc title ftp md5 bytes; do
  strain=$(basename "$ftp" .tar.gz)
  [[ -e $d/.ok.$acc ]] && continue
  rm -rf -- "${d:?}/$strain" "$d/.unpack"
  tar_source "$d" "https://$ftp" "$bytes" "$md5" -z
  mkdir -p "$d/$strain"; find "$d/.unpack" -mindepth 1 -maxdepth 1 -exec mv -t "$d/$strain/" {} +
  rm -rf -- "$d/.unpack"; touch "$d/.ok.$acc"
  log "$acc $strain: $(find "$d/$strain" -name '*.pod5' | wc -l) POD5 files"
done < <(tail -n +2 "$MANIFESTS/$S.tsv")
mark_complete "$d" $S
