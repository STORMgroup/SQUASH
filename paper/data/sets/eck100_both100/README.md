# eck100_both100 (D1): Escherichia coli ATCC 25922, 100x of G (~50x per strand)

- Source: `nvb_ec` (`../../sources.tsv`, `../../download/nvb_ec.sh`).
- Pool: the sample's NanoVarBench deposit; md5-order prefix of its >= 10-kb reads to 150x of G.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1090 -> cut 0.1254); reads in md5 order of their id until 100 x G = 100 x 5,132,220 bp signal-estimated bases (the read crossing the target kept).
- Counts: 28,914 pool reads seen in md5 order, 2,495 dropped by the noise rule, 26,419 kept: 513,238,897 signal-estimated bases (100.00x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text 9ad4a9e51e66cb51dcc752d6b31f1f3f5c8af2d34e521dc7028c388a49650686.
- POD5 signal digest (`../../pod5_digest.py`): 99705e41ebba0663114ce8dca02bf303198379b9191d1c0f6c030d575aa0bcfe (26,419 reads, 6,415,485,924 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 c08736afca74ee8a12df45685b9300302fc6e6981e5a641f2abc5defc78f9071 (26,419 reads).
- QUAST reference: CP149500.1 (NanoVarBench ATCC_25922__202309 assembly; chromosome).
