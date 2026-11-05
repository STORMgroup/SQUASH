#!/usr/bin/env bash
# E. coli ATCC 25922 (ATCC_25922__202309), NanoVarBench 202309 (Hall et al. 2024): figshare 10.26188/25521892 (file 45408628, one tar of POD5s,
# 124970854400 bytes, MD5 a2c74a3928ba91dd55e0654b0720561c as figshare lists it), CC BY 4.0. Sets: eck60_both60, eck100_both100.
# -> SOURCES_DIR/nvb_ec/*.pod5. A local copy instead: SRC_NVB_EC=FOLDER. STREAM=1: unpack while downloading.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
nvb_source nvb_ec 45408628 124970854400 a2c74a3928ba91dd55e0654b0720561c
