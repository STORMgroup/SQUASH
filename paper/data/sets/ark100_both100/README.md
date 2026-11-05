# ark100_both100 (D5): Ascochyta rabiei ME14 (seven strains pooled), 100x of G (~50x per strand)

- Source: `arabiei_prjeb86315` (`../../sources.tsv`, `../../download/arabiei_prjeb86315.sh`).
- Pool: the seven strain deposits ERR14664019-25 pooled; every read >= 10 kb.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.1008 -> cut 0.1159); reads in md5 order of their id until 100 x G = 100 x 40,846,718 bp signal-estimated bases (the read crossing the target kept).
- Counts: 222,812 pool reads seen in md5 order, 29,477 dropped by the noise rule, 193,335 kept: 4,084,683,749 signal-estimated bases (100.00x of G). Preset `large`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: 5 blocks); SHA-256 of the text 78ee4ab2584740e7f98f5f553d7f3ae4ed2b49f519c0af6a3bb90d1ec8574626.
- POD5 signal digest (`../../pod5_digest.py`): 5c60f84cfc366bcb15e4c710554d54e92bffc9525e479a9272c444fb1c2f7569 (193,335 reads, 51,058,548,929 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 d14dfe444f06c9276ddf43ef9bcb26775abc98c60bddce56b91b9d57fbc13a75 (193,335 reads).
- QUAST reference: Me14 GCA_004011695.2 (22 chromosomes + mitochondrion).
