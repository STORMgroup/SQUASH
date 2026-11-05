# h2k60_both60 (D6): Homo sapiens HG002, chromosome 13 q-arm, 60x of G (~30x per strand)

- Source: `hg002_giab2025` (`../../sources.tsv`, `../../download/hg002_giab2025.sh`).
- Pool: the read-id list ../../pools/hs13q/ (made with ONT's basecalls and a reference: ../../../curation/hs13q_pool.sh): unsplit chr13 primary records of ONT's sup v5.0.0 GRCh38 BAMs (duration x 400 >= 9 kb), remapped to HG002 Q100: MAPQ >= 10 on chr13_MATERNAL from 17.24 Mb.
- Selection, from the raw signal alone (`../../../curation/`): signal-estimated length (samples x 400 / sample rate) >= 10 kb; noise (median absolute first difference / sd of the calibrated signal) <= 1.15 x the pool's median (0.0921 -> cut 0.1059); reads in md5 order of their id until 60 x G = 60 x 96,096,965 bp signal-estimated bases (the read crossing the target kept).
- Counts: 225,001 pool reads seen in md5 order, 18,192 dropped by the noise rule, 206,809 kept: 5,765,832,274 signal-estimated bases (60.00x of G). Preset `large`.
- The paper's runs assembled every read of the set (no `assembled.txt.gz`).
- `reads.txt.gz`: the reads in the order of the POD5 the paper's run read (SQUASH reads a POD5 in blocks of 40,000 reads: 6 blocks); SHA-256 of the text 955ad1922e0a3911253ae988137b86284f324175aa4502cd656e83ad36b6d1c1.
- POD5 signal digest (`../../pod5_digest.py`): ee70d87bf6182677cd2cddce3596fc1661e092835e17c664c337e0f4745c58ec (206,809 reads, 72,072,905,980 samples).
- SQUASH's q-symbols (`work/signal/qsymbols.tsv`): SHA-256 10ff2544f73ee7573f716f60934c48997ce67f02f464c6c55fca590dc2f8691e (206,809 reads).
- QUAST reference: HG002 Q100 v1.2 chr13_MATERNAL from 17,240,000 (chr13q_MAT; also both haplotypes).
