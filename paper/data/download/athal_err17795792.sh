#!/usr/bin/env bash
# Arabidopsis thaliana Col-0 (Earlham Institute, ENA PRJEB79853 "Nanopore-only assembly and methylation profiling of
# Arabidopsis thaliana across chemistries, basecalling models and assembly workflows"): run ERR17795792, flow cell
# PBE26394 (FLO-PRO114M, SQK-LSK114, 5 kHz), ATRun3_pod5.tar.gz, 163,758,872,949 bytes, MD5
# f19752ffe88f54ec4b97ea6d5e1a4d46 (ENA submitted_md5); 75 POD5 files (its pod5_skip/ folder: every read, live
# basecalling skipped). INSDC terms. Set: awk60_both60 (its reads are in paper/data/sets/awk60_both60/reads.txt.gz).
# -> SOURCES_DIR/athal_err17795792/ (the tarball's tree). A local copy instead: SRC_ATHAL_ERR17795792=FOLDER.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
S=athal_err17795792
d=$(source_dir $S) || { log "$S: done before ($d)"; exit 0; }
tar_source "$d" "https://ftp.sra.ebi.ac.uk/vol1/run/ERR177/ERR17795792/ATRun3_pod5.tar.gz" 163758872949 \
  f19752ffe88f54ec4b97ea6d5e1a4d46 -z
find "$d/.unpack" -mindepth 1 -maxdepth 1 -exec mv -t "$d/" {} +
rm -rf -- "$d/.unpack"
mark_complete "$d" $S
