# cjk100_both100 (D2): Campylobacter jejuni ATCC 33560, 100x of G (~50x per strand)

- Source: `nvb_cj` (`../../sources.tsv`, `../../download/nvb_cj.sh`).
- Pool: the sample's NanoVarBench deposit; md5-order prefix of its >= 10-kb reads to 100x of G (k60) / 160x (k100).
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1024 -> cut 0.1178); reads in md5 order of their id until 100 x G = 100 x 1,768,448 bp signal-estimated bases (the read crossing the target kept).
- Counts: 11,790 pool reads seen in md5 order, 1,505 dropped by the noise rule, 10,285 kept: 176,859,588 signal-estimated bases (100.01x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text 1bf2443903ddfa082d7c444a5acd323f11fc2b4401295367ad1a8532e0c89a79.
- POD5 signal digest (`../../pod5_digest.py`): a716202daf7818da4d46a5d79cf8c733e9a9fe42ea70a1edac02c716a34c307d (10,285 reads, 2,210,745,028 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 a88da5f2312d06071f3f77b18eb606b59fbe36fa4566b50a657f87d9870a5375 (10,285 reads).
- QUAST reference: CP149499.1 (NanoVarBench ATCC_33560__202309 assembly).
