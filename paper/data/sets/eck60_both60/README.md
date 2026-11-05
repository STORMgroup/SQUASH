# eck60_both60 (supplement): Escherichia coli ATCC 25922, 60x of G (~30x per strand)

- Source: `nvb_ec` (`../../sources.tsv`, `../../download/nvb_ec.sh`).
- Pool: the sample's NanoVarBench deposit; md5-order prefix of its >= 10-kb reads to 150x of G.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1090 -> cut 0.1254); reads in md5 order of their id until 60 x G = 60 x 5,132,220 bp signal-estimated bases (the read crossing the target kept).
- Counts: 17,325 pool reads seen in md5 order, 1,485 dropped by the noise rule, 15,840 kept: 307,960,759 signal-estimated bases (60.01x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text 450b3665c544a69f7cfa57e315952931f48a06ba7b17f0e824bae819350832fb.
- POD5 signal digest (`../../pod5_digest.py`): 68b50d1ea0d22e24128fe8a461d224076769132b85a851de59872272e72016c4 (15,840 reads, 3,849,509,312 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 8857740ff2f9a1eaab9235023fb75b2f4d31c9e6491ee3397340c9db5c1c98e6 (15,840 reads).
- QUAST reference: CP149500.1 (NanoVarBench ATCC_25922__202309 assembly; chromosome).
