# Parity: SQUASH's outputs and the paper's runs

SQUASH 1.0.0 gives the same results as the runs the paper reports. On every paper set, a SQUASH run was compared with
the paper's run of that set:

1. **Intermediate data** (a run with `--keep-work`): every data file the paper's run kept and SQUASH also makes (the
   signal lengths, the uncorrected tokens, the final round's corrected tokens, the graph's input and its k = 32 files,
   the threaded graph, its unitig map, joins and stats), file by file (gzip files by their text, a first line that names
   the format left out).
2. **Graph and layout**: `assembly_graph.gfa`, `layout.tsv.gz` and `read_paths.tsv.gz` against the paper's evaluation
   of that run (its graph and the signal layout's tables), and `joins.tsv` and `stats.tsv` against the paper's run.
3. **Outputs**: the six output files were then recorded in `expected/KEY.sha256` (a `.gz` file by its decompressed
   text, since gzip bytes depend on the compressor), which `check.sh` and `make check` compare a run with.

| Set | Preset | Intermediate data | Graph and layout | Threads |
|---|---|---|---|---|
| `eck100_both100` (D1) | small | 17 / 17 identical | identical | 12, 64 |
| `cjk100_both100` (D2) | small | 17 / 17 identical | identical | 12, 64 |
| `vpk100_both100` (D3) | small | 17 / 17 identical | identical | 12, 64 |
| `u9k100_both100` (D4) | small | 17 / 17 identical | identical | 12, 64 |
| `ark100_both100` (D5) | large | | identical | 64 |
| `h2k60_both60` (D6) | large | | identical | 64 |
| `dmk100_both100` (D7) | large | | identical | 64 |
| `awk60_both60` (D8) | large | | identical | 64 |
| `eck60_both60` | small | 17 / 17 identical | identical | 12, 64 |
| `cjk60_both60` | small | 17 / 17 identical | identical | 12, 64 |
| `vpk60_both60` | small | 17 / 17 identical | identical | 12, 64 |
| `u9k60_both60` | small | 17 / 17 identical | identical | 12, 64 |
| `ark60_both60` | large | 14 / 14 identical | identical | 16, 64 |
| `dmk60_both60` | large | | identical | 64 |

The intermediate data of the large main sets were not compared (their runs were timed, without `--keep-work`); A.
rabiei's 30x run kept the graph's input but not its k = 32 files, and A. thaliana's run not its joins and stats. The outputs do not depend on the thread count: sets run
at several thread counts (the paper's runs used 64) gave the same outputs at each.

**From a fresh clone.** A clone built as the README says (rustup installing the pinned toolchain, `make`, `make test`),
with the Python of `environment.yml`, the default `config.sh` and the real download of C. jejuni's deposit: `make -C
paper check SET=cjk60_both60` (download, cut, digest check, SQUASH, outputs), `make -C paper threads SET=cjk60_both60`
(4 threads) and `make -C paper curate SET=cjk60_both60` (the read list remade from the deposit) all equal the paper's.

**The whole chain on D2.** From a clone, through `make tools`, `data`, `check` (`KEEP_WORK=1`), `basecall`, `assemble`,
`eval` and `tables`: SQUASH's outputs, every sup read with its move table and the HAC reads equal the paper's; hifiasm
and Rawsamble + miniasm give the paper's assemblies byte for byte and Flye and hifiasm on HAC reads the paper's circular
contig from another start; the evaluation (spelled and placed sequences, strands, QUAST, contiguity) and D2's rows of
the contiguity, correctness, HAC and strand tables equal the paper's. nanoMDBG differs from run to run (28 or 29 contigs
on D2, the paper's own binary too), and the times depend on the machine.
