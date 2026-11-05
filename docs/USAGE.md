# Using SQUASH

Build SQUASH first ([README.md](../README.md#install)); `squash` below is `target/release/squash`, or the command itself
once `make install` put it on your `PATH`.

```
squash -x small|large -o OUT_DIR [-t THREADS] [options] READS.pod5 [MORE.pod5 | POD5_DIR ...]
```

## Quick usage

Each line is a separate example (its comment says what it does): run the ones you need.

```sh
squash -x small -o ecoli ecoli.pod5                                  # a bacterium, every core
squash -x large -t 64 --tmp-dir /scratch/$USER -o dmel dmel_pod5/    # a eukaryote, a folder of POD5 files
squash -x large --rounds 20 --tip-max 30k -o dmel_r20 dmel_pod5/     # the preset with two values changed
squash -x large --dry-run -o dmel dmel_pod5/                         # print every value and stop
squash -x large --resume -o dmel dmel_pod5/                          # continue a stopped run
```

## Input

POD5 files, or folders of them (a folder's `*.pod5` files in name order). Reads are taken in file order, in blocks of
40,000 reads, each block sorted by read id: the same files in the same order give the same result. R10.4.1 at 4 or 5
kHz. `--read-ids FILE` (one id per line, plain or gzipped) assembles only the listed reads, after the read order has
been fixed over every read of the files.

## The preset and the main options

| Option | What it does |
|---|---|
| `-x`, `--preset small\|large` | Required. `small` for genomes under 10 Mb, `large` from 10 Mb; sets every default below |
| `-o`, `--out-dir DIR` | Required. The output folder (created if missing) |
| `-t`, `--threads N` | Threads (default: every core this process may use); outputs do not depend on it |
| `-g`, `--genome-size SIZE` | Expected genome size (e.g. `4.6m`): only for the coverage estimate and its warnings |
| `--read-ids FILE` | Assemble only these reads |
| `--tmp-dir DIR` | Large temporary files (default `OUT/work/tmp`); a local disk helps |
| `--keep-work` | Keep `OUT/work`, the intermediate files ([OUTPUTS.md](OUTPUTS.md)) |
| `--resume` | Continue a stopped run: stages already run with the same values are kept |
| `--stop-after STAGE` | Stop after `signal`, `correct`, `assemble`, `thread` or `layout` (keeps `OUT/work`) |
| `--dry-run` | Print every resolved value and stop |
| `--no-threading`, `--no-tangles`, `--no-joins` | Leave out read threading, its tangle rule, or the contig joins |
| `-h`, `--help`; `-V`, `--version` | The main options (`-h`) or every option by step (`--help`); the version |
| | A run without `--resume` starts over: an earlier `OUT/work` and earlier outputs in `OUT` are removed first |

## Every tunable value

The defaults are the paper's values under each preset. `-h` shows the main options only, `--help` all of them by step.

### Signal: events and q-symbols

| Option | small | large | What it sets |
|---|---|---|---|
| `--evdet-win` | 3,6 | 3,6 | Event detector window lengths, short,long (RawHash2 window_length1,2) |
| `--evdet-thr` | 5.5,3.5 | 5.5,3.5 | Event detector t-statistic thresholds, short,long (RawHash2 threshold1,2) |
| `--evdet-peak` | 0.2 | 0.2 | Event detector peak height (RawHash2 peak_height) |
| `--event-merge` | 0.30 | 0.30 | Merge consecutive events whose means stay within this many z-units (0 = off) |
| `-q`, `--quant-bits` | 3 | 3 | Bits per q-symbol (2^q quantization levels of event means) |
| `-e`, `--events` | 9 | 9 | q-symbols per token (each token hashes e consecutive q-symbols) |
| `-w`, `--window` | 16 | 16 | Minimizer window: one token per window of w hashes |
| `--bases-per-second` | 400 | 400 | Translocation speed of the chemistry; converts signal length to bp for every bp limit (R10.4.1 400 bps) |

### Correction in token space

| Option | small | large | What it sets |
|---|---|---|---|
| `--low-complexity` | 0.2 | 0.2 | Leave out of correction (as targets and donors) a read whose most frequent token is at least this fraction of its tokens |
| `--min-anchors` | 4 | 4 | Shared tokens a partner's chain needs |
| `--min-overlap` | 500 | 500 | Overlap a partner needs, in q-symbols |
| `--min-cover` | 0.5 | 0.5 | Fraction of the implied overlap a partner's chain must cover |
| `--coverage` | 20 | 20 | Coverage cap: chains are taken best first while some target token in their span has fewer than N |
| `--seed-repeat-fraction` | 0.0001 | 0.0001 | Fraction of the most frequent token values not used as seeds (nanoMDBG's repeat rule) |
| `--rounds` | 12 | 16 | Correction rounds at most |
| `--stop-below` | 0.0005 | 0.0005 | Stop early when fewer than this fraction of tokens change in a round |
| `--band` | 256 | 256 | Band of the token-level POA alignment |
| `--hops-from` | off | 2 | From this round, look for partners two hops away (partners of partners) |
| `--changed-from` | off | 2 | From this round, search new partners only for reads that changed |
| `--skip-from` | off | 8 | From this round, correct only the segments that changed |
| `--skip-pad` | 150 | 150 | Q-symbols of padding around a changed segment (with --skip-from) |
| `--qsym-check` | off | 0.667 | Q-symbol alignment check of partners: largest weighted edit cost per alignment column |

### Assembly (multi-k minimizer-space de Bruijn graph)

| Option | small | large | What it sets |
|---|---|---|---|
| `--density` | 8 | 8 | Assemble with 1 in N corrected tokens, chosen by value hash |
| `--repeat-ppm` | 100 | 100 | Leave out the most frequent token values, in parts per million (metaMDBG's repetitive-minimizer cut) |
| `--k-min` | 4 | 4 | Smallest k of the multi-k graph |
| `--k-max` | 32 | 32 | Largest k of the multi-k graph (the final graph's k) |
| `--abundance-cutoffs` | 2,3 | 2,3 | Progressive abundance cutoffs of the graph cleaning |
| `--tip-max` | 50000 | 50000 | Longest tip removed, in bp |
| `--superbubble-max` | 50000 | 50000 | Longest superbubble popped, in bp |
| `--max-path-rows` | 100000000 | 1000000000 | Resource cap of the graph construction (no effect unless reached) |

### Read threading and contig joins

| Option | small | large | What it sets |
|---|---|---|---|
| `--min-support` | 3 | 3 | Reads a threaded pairing of unitigs needs |
| `--max-conflict` | 0.2 | 0.2 | Largest fraction of conflicting reads of a threaded pairing |
| `--anchor-nodes` | 200 | 200 | Atomic nodes a unitig needs to anchor a tangle (200 nodes ~ 20 kb) |
| `--anchor-edge` | 40 | 40 | Nodes at an anchor's end where reads are counted |
| `--tangle-dominance` | 3 | 3 | A tangle's pairing needs this many times the reads of the next best |
| `--join-paths` | 2 | 2 | Independent read walks a contig join needs |
| `--join-threshold` | 0.62 | 0.62 | Largest weighted edit cost per column of a join hop's q-symbol alignment |
| `--join-min-nodes` | 200 | 200 | Nodes a unitig needs for its ends to be joined, or to be landed on |
| `--join-head` | 200 | 200 | A join lands only in the first N tokens of the other unitig |

### Signal layout

| Option | small | large | What it sets |
|---|---|---|---|
| `--layout-max-gap` | 32 | 32 | Unmatched windows allowed inside a read's piece of a unitig |
| `--layout-own-gap` | 2000 | 2000 | Largest gap, in q-symbols, between a piece's tokens read from its own signal |

## Fixed parts of the method

Not options (the paper's method; changing them changes what SQUASH is): the read order (blocks of 40,000 reads);
the q-symbol filter of the tokens (runs of equal q-symbols collapsed); the partner chains' gap (120 tokens), drift
(6 + 0.25 x gap) and look-back (60); the correction rounds' tie-break and anchored alignment (pad 8); guarded
re-recruitment; the q-symbol alignment check's window, band, path, gap cost (2), whole-overlap rule, anchored seed and
first-search chain cover (0.25), and its lazy form (checked only for the chains the coverage cap keeps); metaMDBG's tip
rule and the carry of each order's contigs into the next; the conversion of every bp value to tokens (the median
signal-estimated bp per token of the assembled tokens, rounded half to even); the joins' word size and hop rules.

## Diagnostics

Environment variables that print profiles to the steps' logs (`OUT/work/logs/`, kept with `--keep-work`); none of
them changes an output. The partner search: `SQUASH_PARTNERS_PROFILE=1` (phase times and counts),
`SQUASH_PARTNERS_COVER_PROBE=1` (the donors a bound on their anchors could skip), `SQUASH_PARTNERS_MEMTRACK=1` (memory,
once a second), `SQUASH_QSYM_PROFILE=1` (the q-symbol check). The correction rounds: `SQUASH_ROUNDS_TIMING=1`,
`SQUASH_ROUNDS_MEM=1`. The graph: `SQUASH_GRAPH_TIMING=1`, `SQUASH_GRAPH_MEM=1`, and `SQUASH_GRAPH_SERIAL=1` (each order
computed after the previous one is written, so that the phases' times do not overlap).

## Output

See [OUTPUTS.md](OUTPUTS.md).
