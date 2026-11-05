# cjk60_both60 (supplement): Campylobacter jejuni ATCC 33560, 60x of G (~30x per strand)

- Source: `nvb_cj` (`../../sources.tsv`, `../../download/nvb_cj.sh`).
- Pool: the sample's NanoVarBench deposit; md5-order prefix of its >= 10-kb reads to 100x of G (k60) / 160x (k100).
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1024 -> cut 0.1178); reads in md5 order of their id until 60 x G = 60 x 1,768,448 bp signal-estimated bases (the read crossing the target kept).
- Counts: 7,099 pool reads seen in md5 order, 922 dropped by the noise rule, 6,177 kept: 106,116,428 signal-estimated bases (60.01x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text e4eb791bf3f62fd21d5f553b932f788c38c991ceda9fda2b3c0346ddff2565d1.
- POD5 signal digest (`../../pod5_digest.py`): eec3123c16c72cfa9235db9fad67be454d23aedb53f72014c96207c58ec39381 (6,177 reads, 1,326,455,443 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 26098dae926627d8b14a433cd3316f6579642fc74ad0ca6497769918d8e2f8fa (6,177 reads).
- QUAST reference: CP149499.1 (NanoVarBench ATCC_33560__202309 assembly).
