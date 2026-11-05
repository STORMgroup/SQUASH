# paper/data/sets: one folder per read set

- `reads.txt.gz`: the set's read ids, one per line (`\n`), in the file order of the POD5 the paper's run read; gzip -9
  without name or time stamp (`../../sets.tsv` gives the SHA-256 of the text and of the file). `../pod5_cut.py` writes
  the POD5 in this order and `../pod5_digest.py --ids` checks it.
- `assembled.txt.gz` (awk60_both60, dmk60_both60, dmk100_both100): the subset the paper's runs assembled
  (SQUASH's `--read-ids`): the reads with a basecalled record under their own id (the basecaller split or
  dropped the rest). The other 11 sets ran on every read of their POD5.
- `README.md`: source, pool, selection rule and counts, digests.
