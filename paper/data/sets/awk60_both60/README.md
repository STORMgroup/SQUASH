# awk60_both60 (D8): Arabidopsis thaliana Col-0, whole genome, 60x of G (~30x per strand)

- Source: `athal_err17795792` (`../../sources.tsv`, `../../download/athal_err17795792.sh`).
- Pool: the whole run ERR17795792; md5-order prefix of its >= 10-kb reads to 100x of G.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.0885 -> cut 0.1018); reads in md5 order of their id until 60 x G = 60 x 142,479,245 bp signal-estimated bases (the read crossing the target kept).
- Counts: 279,508 pool reads seen in md5 order, 28,596 dropped by the noise rule, 250,912 kept: 8,548,808,017 signal-estimated bases (60.00x of G). Preset `large`.
- `assembled.txt.gz`: the 250,544 reads the paper's runs assembled, those with a basecalled record under their own id (SQUASH's `--read-ids`; SHA-256 of the text a58203cff198adb75520a1202331a911a3aed9e47b128b53700c0a9f2b038f33). The POD5 keeps all 250,912 reads: SQUASH fixes the read order over every read of its input before it drops the others.
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: 7 blocks); SHA-256 of the text 844062395cd8b66878e530b74a6acf72d22bc9301d4cbb8a0cedd31edd634efa.
- POD5 signal digest (`../../pod5_digest.py`): d9587b45f048d04c8880fccd31ceb3408a3a1bf8951b091fa7f7212410d76971 (250,912 reads, 106,860,100,815 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 50ba0eb40602fb7ca2840e3a3fb24619c04603a4a3c9cb7a8b4ee6d96cf89e57 (250,544 reads).
- QUAST reference: Col-CC GCA_028009825.2 (five chromosomes) + chloroplast + mitochondrion.
