#!/usr/bin/env bash
# V. parahaemolyticus ATCC 17802 (ATCC_17802__202309), NanoVarBench 202309 (Hall et al. 2024): figshare 10.26188/25495063 (file 45334792, one tar of POD5s,
# 84951562240 bytes, MD5 d4e37c56998981ba4303fec0f3557918 as figshare lists it), CC BY 4.0. Sets: vpk60_both60, vpk100_both100.
# -> SOURCES_DIR/nvb_vp/*.pod5. A local copy instead: SRC_NVB_VP=FOLDER. STREAM=1: unpack while downloading.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
nvb_source nvb_vp 45334792 84951562240 d4e37c56998981ba4303fec0f3557918
