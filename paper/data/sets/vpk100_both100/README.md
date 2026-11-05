# vpk100_both100 (D3): Vibrio parahaemolyticus ATCC 17802, 100x of G (~50x per strand)

- Source: `nvb_vp` (`../../sources.tsv`, `../../download/nvb_vp.sh`).
- Pool: the sample's NanoVarBench deposit; md5-order prefix of its >= 10-kb reads to 100x of G (k60) / 160x (k100).
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1064 -> cut 0.1223); reads in md5 order of their id until 100 x G = 100 x 5,147,091 bp signal-estimated bases (the read crossing the target kept).
- Counts: 29,917 pool reads seen in md5 order, 2,743 dropped by the noise rule, 27,174 kept: 514,718,265 signal-estimated bases (100.00x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text 77a9e55f49f57ce15c1d2dd150563ecc39463f1f9ead4030728dda7e41952601.
- POD5 signal digest (`../../pod5_digest.py`): 0ab601933549d5ef6b397e0a87291f5846c9773eb2da0a8a0b81b39a0a0a26ea (27,174 reads, 6,433,978,523 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 b96e20b4a484d4c6f24e722accbca8d33525f7b0a1a15d344f78f8dc71ec148c (27,174 reads).
- QUAST reference: CP149505.1 + CP149506.1 (NanoVarBench ATCC_17802__202309 assembly; both chromosomes).
