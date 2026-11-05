# u9k60_both60 (supplement): Streptococcus dysgalactiae subsp. equisimilis UT9728, 60x of G (~30x per strand)

- Source: `sdys_zenodo17428855` (`../../sources.tsv`, `../../download/sdys_zenodo17428855.sh`).
- Pool: the read-id list ../../pools/sdys/ (made with basecalls: ../../../curation/sdys_pool.sh): the paper's demultiplexing of run FBA61885 (dorado 1.4.0 hac v5.2.0 --kit-name SQK-RBK114-24 --disable-read-splitting), barcodes 01-06 + 09-14 = UT9728.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1042 -> cut 0.1198); reads in md5 order of their id until 60 x G = 60 x 2,093,215 bp signal-estimated bases (the read crossing the target kept).
- Counts: 6,597 pool reads seen in md5 order, 425 dropped by the noise rule, 6,172 kept: 125,627,262 signal-estimated bases (60.02x of G). Preset `small`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: one block: any order gives the same outputs); SHA-256 of the text 487d0f09af2c35915ec83af64e0c0ae12f249cfafcd5b5e2b580a6537e0c09c6.
- POD5 signal digest (`../../pod5_digest.py`): bd008f6811c048b44d4e2c9d8f74d4c0fc2b1df20a4bc89f03487e688edb49cd (6,172 reads, 1,570,340,770 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 e7befa8c9fb5d2e6229088a3fc6cf5941812b604959cc95550cc2b0e4ef461ca (6,172 reads).
- QUAST reference: UT9728_IPa_Hyb (Zenodo 15558489, the depositors' hybrid assembly).
