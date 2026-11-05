# dmk60_both60 (supplement): Drosophila melanogaster ISO1, 60x of G (~30x per strand)

- Source: `dmel_bkim2023` (`../../sources.tsv`, `../../download/dmel_bkim2023.sh`).
- Pool: every read of the six runs, pre-filtered on the runs' MinKNOW sequencing summaries (duration x 400 >= 9 kb), then >= 10 kb by POD5.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.0989 -> cut 0.1138); reads in md5 order of their id until 60 x G = 60 x 137,547,960 bp signal-estimated bases (the read crossing the target kept).
- Counts: 376,269 pool reads seen in md5 order, 42,193 dropped by the noise rule, 334,076 kept: 8,252,886,347 signal-estimated bases (60.00x of G). Preset `large`.
- `assembled.txt.gz`: the 330,203 reads the paper's runs assembled, those with a basecalled record under their own id (SQUASH's `--read-ids`; SHA-256 of the text 32c59dcc9850c4d5f47130745199834405154813aac9bad91df57aa7d0ba87ad). The POD5 keeps all 334,076 reads: SQUASH fixes the read order over every read of its input before it drops the others.
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: 9 blocks); SHA-256 of the text 8418ac07bcc199498e2604bc834cced88386e67946e99b797d301dbb84ee9643.
- POD5 signal digest (`../../pod5_digest.py`): 271f19a0c395457aeeaa6744aac9195dc1e4e63d6d2b7bd47a2a5a2bd1944eef (334,076 reads, 82,528,860,885 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 7a1d028af286ec4f6265c3afc82289083fb7ff85b87fe8ffc5741a88e46999ed (330,203 reads).
- QUAST reference: release 6 GCF_000001215.4: X, 2L, 2R, 3L, 3R, 4, Y.
