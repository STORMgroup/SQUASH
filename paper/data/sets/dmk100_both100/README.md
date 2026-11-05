# dmk100_both100 (D7): Drosophila melanogaster ISO1, 100x of G (~50x per strand)

- Source: `dmel_bkim2023` (`../../sources.tsv`, `../../download/dmel_bkim2023.sh`).
- Pool: every read of the six runs, pre-filtered on the runs' MinKNOW sequencing summaries (duration x 400 >= 9 kb), then >= 10 kb by POD5.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.0989 -> cut 0.1138); reads in md5 order of their id until 100 x G = 100 x 137,547,960 bp signal-estimated bases (the read crossing the target kept).
- Counts: 627,727 pool reads seen in md5 order, 70,473 dropped by the noise rule, 557,254 kept: 13,754,808,873 signal-estimated bases (100.00x of G). Preset `large`.
- `assembled.txt.gz`: the 550,874 reads the paper's runs assembled, those with a basecalled record under their own id (SQUASH's `--read-ids`; SHA-256 of the text 06d06931dba54b8d8298ad620b55e684ad633a83bf33dcf932979203b5c0f2be). The POD5 keeps all 557,254 reads: SQUASH fixes the read order over every read of its input before it drops the others.
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: 14 blocks); SHA-256 of the text bdf5c53ac1ff30b93a3debdc307c66a45ccf19ef77004cecc9d41d72cf676ddf.
- POD5 signal digest (`../../pod5_digest.py`): 9204053c4031194f442cda20f1776a15b3c7444620253775e99c9b50f9a340a2 (557,254 reads, 137,548,089,056 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 f9ecfaf5693ef38516d79dca57a72c94fcae4f915355388351bcabcad9cfcb52 (550,874 reads).
- QUAST reference: release 6 GCF_000001215.4: X, 2L, 2R, 3L, 3R, 4, Y.
