# vpk60_both60 (supplement): Vibrio parahaemolyticus ATCC 17802, 60x of G (~30x per strand)

- Source: `nvb_vp` (`../../sources.tsv`, `../../download/nvb_vp.sh`).
- Pool: the sample's NanoVarBench deposit; md5-order prefix of its >= 10-kb reads to 100x of G (k60) / 160x (k100).
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1063 -> cut 0.1223); reads in md5 order of their id until 60 x G = 60 x 5,147,091 bp signal-estimated bases (the read crossing the target kept).
- Counts: 17,789 pool reads seen in md5 order, 1,596 dropped by the noise rule, 16,193 kept: 308,836,160 signal-estimated bases (60.00x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text af9ab18139093cd3f1ed66b2e232dccf8d033cf928e7e64217ed5a2f863865eb.
- POD5 signal digest (`../../pod5_digest.py`): 5ec32d9e25b6c357a5c6cc410ca71b75f54e0e4a2f111a8261290b93469bc505 (16,193 reads, 3,860,451,932 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 eb57f399fb18401b3c2168e118f99ec70210d49d439fb876800d92937d05f3b2 (16,193 reads).
- QUAST reference: CP149505.1 + CP149506.1 (NanoVarBench ATCC_17802__202309 assembly; both chromosomes).
