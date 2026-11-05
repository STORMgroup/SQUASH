# ark60_both60 (supplement): Ascochyta rabiei ME14 (seven strains pooled), 60x of G (~30x per strand)

- Source: `arabiei_prjeb86315` (`../../sources.tsv`, `../../download/arabiei_prjeb86315.sh`).
- Pool: the seven strain deposits ERR14664019-25 pooled; every read >= 10 kb.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1008 -> cut 0.1159); reads in md5 order of their id until 60 x G = 60 x 40,846,718 bp signal-estimated bases (the read crossing the target kept).
- Counts: 133,832 pool reads seen in md5 order, 17,740 dropped by the noise rule, 116,092 kept: 2,450,822,674 signal-estimated bases (60.00x of G). Preset `large`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: 3 blocks); SHA-256 of the text 1d7247e3a6499867cde440664a13604b9c3635d49d6cf8b1656e436630adfc84.
- POD5 signal digest (`../../pod5_digest.py`): f55c94356a7fc3c0b421cd579eb3a3a0ab510c41ce88ab0e1498e9fcb09b63a9 (116,092 reads, 30,635,283,890 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 f1d8018bf2012fadd2285a2c3dd4beb8879b5ac6d6946a9ff9a101c6eeff25f3 (116,092 reads).
- QUAST reference: Me14 GCA_004011695.2 (22 chromosomes + mitochondrion).
