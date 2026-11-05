# u9k100_both100 (D4): Streptococcus dysgalactiae subsp. equisimilis UT9728, 100x of G (~50x per strand)

- Source: `sdys_zenodo17428855` (`../../sources.tsv`, `../../download/sdys_zenodo17428855.sh`).
- Pool: the read-id list ../../pools/sdys/ (made with basecalls: ../../../curation/sdys_pool.sh): the paper's demultiplexing of run FBA61885 (dorado 1.4.0 hac v5.2.0 --kit-name SQK-RBK114-24 --disable-read-splitting), barcodes 01-06 + 09-14 = UT9728.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1042 -> cut 0.1198); reads in md5 order of their id until 100 x G = 100 x 2,093,215 bp signal-estimated bases (the read crossing the target kept).
- Counts: 10,985 pool reads seen in md5 order, 723 dropped by the noise rule, 10,262 kept: 209,347,891 signal-estimated bases (100.01x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text 43bbaf39972f45fa66afab119b1718657bfff4da2957e24e1f4113499ed2f142.
- POD5 signal digest (`../../pod5_digest.py`): ffb11b8d2a5045eba95aa9572f679d29606e5e77dd696872723f552c13dadd62 (10,262 reads, 2,616,848,187 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 85911cdbb7fde43455162240ee5aace7a6641445eff84e6b47fa76c34367ecf3 (10,262 reads).
- QUAST reference: UT9728_IPa_Hyb (Zenodo 15558489, the depositors' hybrid assembly).
