# paper/data: the paper's read sets

The reads are not re-hosted. Each of the 14 sets is a list of read ids (`sets/KEY/reads.txt.gz`) taken from a public
run (`sources.tsv`); `extract.sh KEY` downloads the run if needed, cuts the listed reads into one POD5 in the list's
order and checks it against the POD5 the paper's run read. From the repository's root, with the Python of
[../README.md](../README.md#before-you-start) active, either line (the second does the same without make):

```sh
make -C paper data SET=cjk60_both60        # one set (a key of ../sets.tsv); SET=main, supp or all for several
bash paper/data/extract.sh cjk60_both60    # the same, without make
```

The POD5 lands in `POD5_DIR/KEY.pod5` (default `paper/data/pod5/`), the download in `SOURCES_DIR` (default
`paper/data/sources/`); both are set in `../config.sh`.

| Path | What |
|---|---|
| `../sets.tsv` | one row per set: key, paper label (D1-D8 main, `supp`), genome, source(s), pool, G, coverage, preset, bases per second, reads, assembled reads, the list's SHA-256, the POD5 signal digest, the q-symbols table's SHA-256, the QUAST reference |
| `sets/KEY/reads.txt.gz` | the set's read ids, one per line, in the file order of the POD5 the paper's run read |
| `sets/KEY/assembled.txt.gz` | awk60, dmk60, dmk100 only: the reads the runs assembled (SQUASH's `--read-ids`) |
| `sets/KEY/README.md` | source, pool, selection rule and counts, digests |
| `pools/sdys/`, `pools/hs13q/` | the pools of the three sets whose pool took basecalls or a reference to make, as read-id lists: the UT9728 reads of the paper's demultiplexing (u9k60, u9k100), the chr13 q-arm cohort (h2k60); each with a README. How they and every set's list are remade: [../curation/](../curation/README.md) |
| `sources.tsv` | every source run and reference: location, files, bytes, checksum and where it comes from, terms, citation, script |
| `download/SOURCE.sh` | one script per source (`lib.sh`: resumable, size / MD5 checked, local copies); `references.sh` the evaluation references |
| `pod5_cut.py` | the listed reads into one POD5 in the list's order, every run info first (`cut`; `part` + `join` for streamed sources) |
| `pod5_digest.py` | the signal digest of a POD5 (read id, calibration offset and scale, sample rate, sample count, int16 samples, in file order: one SHA-256); `--ids LIST` first checks the read ids in file order and stops at the first differing read |
| `extract.sh KEY [KEY ...]` | download if missing, cut, check ids and digest (exit non-zero on any difference) |

## Why the order matters

SQUASH reads a POD5 in blocks of 40,000 reads of file order, visits the blocks in the string order of `"{c}_0.fast5"`
(block 10 before block 2) and sorts the reads by id inside a block; only then are the reads outside `--read-ids`
dropped. So a set's POD5 must hold the same reads in the same order, and a POD5 cut down to the assembled reads would
move the block boundaries (awk60, dmk60 and dmk100 keep two lists for this). The bacterial sets have at most 40,000
reads, one block: any order gives the same outputs. `pod5 filter` and `pod5 subset` keep no order, hence `pod5_cut.py`.

## How the lists were made

Each list holds the reads of the set's selection (`../curation/README.md`: from the raw signal alone, applied to a pool
of the genome's reads), in the order of the POD5 the paper's run read: the source files in path order and each file's
reads in file order, as `pod5_cut.py` writes them. `extract.sh` checks a rebuilt POD5 against the list (read ids in
file order) and against `sets.tsv`'s signal digest (every sample of every read).

## Disk and time

Downloads for all 14 sets: ~4.31 TB (human 3.28 TB, streamed: only the listed reads are kept, a peak of ~65 GB plus the
files in flight; D. melanogaster's FAST5s converted and cut as they stream). Set POD5s: ~450 GB (the main sets ~336 GB).
Both `pod5_cut.py` and `pod5_digest.py` read every sample of the set, so their time is the storage's read speed
(dmk100's POD5 is 122 GB). pod5 converts older POD5 files (the A. thaliana run's) in a temporary folder in the working
directory when it opens them: run the scripts from a writable folder.

## Existing copies

A copy of a source that is already on disk can be named instead of downloading it (`SRC_<SOURCE>=FOLDER` in
`paper/config.sh` or `paper/config.local.sh`, see `download/lib.sh`). A folder with POD5s of other samples works too
(only the listed reads are taken), but a read present in two files is taken from the first in path order. A set POD5
that already exists in `POD5_DIR` (as `KEY.pod5`, a file or a link) is only checked, not cut again.
