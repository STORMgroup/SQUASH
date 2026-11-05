# paper/curation: how the paper's read lists were selected

Each of the 14 sets (`../data/sets/KEY/reads.txt.gz`) was selected from the raw signal alone: POD5 signal and metadata,
no basecalls and no reference, applied to a pool of the genome's reads. For 11 sets the pool is a public run or deposit
as it is; for u9k60 and u9k100 (which barcodes are strain UT9728's) and h2k60 (which reads come from the chr13 q-arm)
the pool took basecalls and a reference to make, so it is committed as a read-id list (`../data/pools/`) with the
scripts that made it. Rebuilding a set's POD5 from its committed list (`../data/extract.sh`) needs nothing here. Each
line below is a separate example (its comment says what it does): run the ones you need.

```sh
make -C paper curate SET=cjk60_both60      # remake the set's list from its source and compare it (from the repository's root)
bash paper/curation/sdys_pool.sh           # the UT9728 pool list (needs a GPU and the tools of make -C paper tools)
bash paper/curation/hs13q_pool.sh          # the chr13 q-arm pool list (minimap2 and samtools of make -C paper tools)
```

| File | What |
|---|---|
| `select_reads.py` | the selection, every value an option; reads POD5s (pod5 0.3.39) and writes `reads.txt.gz` (the kept reads in md5 order), `selection.tsv` (every pool read: length, noise, open-pore run, kept or why not), `params.tsv` (every option and count) |
| `remake.sh KEY ...` | a set's list remade from its source with its pool and options (below) and checked against `../data/sets/KEY/` (`make -C paper curate SET=KEY`) |
| `sdys_pool.sh` | the UT9728 pool list (`../data/pools/sdys/`): dorado demultiplexing on a GPU, barcode -> strain by mapping to the two strains' assemblies |
| `hs13q_pool.sh` | the chr13 q-arm pool list (`../data/pools/hs13q/`): ONT's sup BAMs, chr13 records remapped to HG002 Q100 v1.2 |
| `pools.py` | the steps of those two scripts (sampling, barcode -> strain, strain reads, cohort) |
| `test_select.py` | tests of `select_reads.py` on synthetic POD5s (`make -C paper test`) |

## The rules

`select_reads.py`'s steps, every value an option (the paper's in brackets):
1. the pool: every read of the source, or the reads of a pool list (`--pool-ids`), or the reads a run's sequencing
   summaries give a duration x bases per second of at least `--summary-min-length` (`--summary`, [9,000]);
2. signal-estimated length = samples x `--bases-per-second` [400] / sample rate, to the nearest base; reads of at least
   `--min-length` [10,000] in md5 order of the read id (the MD5 hex digest of the id's text); `--pool-coverage P`: the
   prefix of that order that reaches P x G signal-estimated bases (the read reaching it included);
3. noise = median absolute first difference / standard deviation of the calibrated signal (pA, in float64), to 6
   significant digits; reads whose noise exceeds `--noise-max` [1.15] x the pool's median are dropped (`--open-pore N`
   would also drop reads with a run of more than N samples at >= 0.85 x the read's `median_before`, 2,000 samples off
   each end; no paper set uses it);
4. the remaining reads in md5 order until `--coverage` [60 or 100] x G signal-estimated bases (`-G`; the read that
   reaches the target kept).

`select_reads.py --table OTHER/selection.tsv` selects again from a table written before (another coverage, or a shorter
prefix of the same pool) without reading the signal. The source folder must hold the run's or deposit's files and
nothing else: other samples' reads would enter the pool. Run it from a writable folder: pod5 migrates older POD5
tables (the A. thaliana run's) into a temporary folder in the working directory.

## The sets

| Sets | Source (`../data/download/`) | Pool | `select_reads.py` options besides `-G`, `--coverage` (`sets.tsv`) | Needs |
|---|---|---|---|---|
| eck60, eck100 | `nvb_ec` (125 GB) | the deposit; prefix to 150 x G | `--pool-coverage 150` | |
| cjk60, cjk100 | `nvb_cj` (41.7 GB) | the deposit; prefix to 100 x G (k60), 160 x G (k100) | `--pool-coverage 100` / `160` | |
| vpk60, vpk100 | `nvb_vp` (85 GB) | as C. jejuni | `--pool-coverage 100` / `160` | |
| u9k60, u9k100 | `sdys_zenodo17428855` (200 GB) | `../data/pools/sdys/ut9728.txt.gz`; prefix to 100 x G (k60), 160 x G (k100) | `--pool-ids ../data/pools/sdys/ut9728.txt.gz --pool-coverage 100` / `160` | the pool list (`sdys_pool.sh`: a GPU, Zenodo 15558489) |
| ark60, ark100 | `arabiei_prjeb86315` (167 GB) | the seven strain deposits; every >= 10-kb read | (none) | |
| awk60 | `athal_err17795792` (164 GB) | the run; prefix to 100 x G | `--pool-coverage 100` | |
| dmk60, dmk100 | `dmel_bkim2023` (246 GB; `SUMMARIES=1`, every read) | the six runs' reads with a summary duration x 400 >= 9 kb; every >= 10-kb read | `--summary SOURCES_DIR/dmel_bkim2023/summaries/*sequencing_summary* --summary-min-length 9000 --duration-digits 6` | the runs' MinKNOW summaries (public) |
| h2k60 | `hg002_giab2025` (3.28 TB streamed; `hg002_giab2025.sh ../data/pools/hs13q/hg002_chr13q.txt.gz` keeps ~100 GB) | `../data/pools/hs13q/hg002_chr13q.txt.gz`; every >= 10-kb read | `--pool-ids ../data/pools/hs13q/hg002_chr13q.txt.gz --duration-digits 6` | the pool list (`hs13q_pool.sh`: ONT's BAMs, 329 GB; HG002 Q100 v1.2) |

For example (`remake.sh` runs these, and checks the result):

    python3 select_reads.py --pod5 SOURCES_DIR/nvb_ec -G 5132220 --pool-coverage 150 --coverage 100 --out eck100
    python3 select_reads.py --table eck100/selection.tsv -G 5132220 --pool-coverage 150 --coverage 60 --out eck60

Two details the lists depend on:
- **Lengths of the D. melanogaster and human sets** are computed from the read's duration in seconds written to 6
  significant digits (`--duration-digits 6`), which moves a length by a base at most for reads under 1,000 s. Exact
  lengths select the same D. melanogaster reads, with slightly different base counts (98 fewer for dmk60, 75 more for
  dmk100, 17 more for h2k60 than their READMEs).
- **Noise** is computed in float64 from pod5's float32 calibrated signal and compared at 6 significant digits.

## Resources

`select_reads.py` reads every source file's read table (sample counts; fast), then the signal of the pool's reads only.
Its memory is the read index (~300 bytes per source read) plus one read's signal per worker: 1-2 GB. With 4-16 workers
it takes minutes on a bacterial deposit and about half an hour on the A. rabiei, A. thaliana and D. melanogaster runs,
mostly reading signal. The pool scripts: the S. dysgalactiae demultiplexing takes about 1.7 GPU-hours (one RTX A5000);
the human cohort reads 329 GB of BAMs and maps ~350,000 reads with minimap2 (about 1.5 h on 8 threads, 10 GB).

## Reproducibility

`remake.sh KEY` remakes a set's list from its source and compares it with the committed list (as a set: the list is in
POD5 order) and with the README's counts (pool reads seen, dropped by the noise rule, kept, signal-estimated bases, the
noise median and cut). For every set but h2k60, `select_reads.py` on the set's pool gives the committed list and the
README's counts exactly. h2k60's pool comes from the 3.28-TB human run: its committed list is the paper's selection, and
on the set's own reads the lengths, the stopping read and every read's noise agree with the README and the rule.
`pools.py` gives the paper's barcode -> strain assignment and UT9728 reads from the demultiplexing table; the GPU
demultiplexing and the human cohort's mapping are documented in `sdys_pool.sh` and `hs13q_pool.sh` with the tool
versions the lists were made with.
