#!/usr/bin/env bash
# C. jejuni ATCC 33560 (ATCC_33560__202309), NanoVarBench 202309 (Hall et al. 2024): figshare 10.26188/25495054 (file 45334624, one tar of POD5s,
# 41652019200 bytes, MD5 c4c74bbc91c853303fe86291cf17f730 as figshare lists it), CC BY 4.0. Sets: cjk60_both60, cjk100_both100.
# -> SOURCES_DIR/nvb_cj/*.pod5. A local copy instead: SRC_NVB_CJ=FOLDER. STREAM=1: unpack while downloading.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
nvb_source nvb_cj 45334624 41652019200 c4c74bbc91c853303fe86291cf17f730
