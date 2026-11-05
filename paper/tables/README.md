# paper/tables/: the paper's tables and figures

One script per table of the paper and its supplement and per data figure, each reading the runs and evaluations under
`OUT_DIR` (`config.sh`); `make.sh` runs them all.

From the repository's root, with the Python of [../README.md](../README.md#before-you-start) active. Each line is a
separate example (its comment says what it does): run the ones you need.

```sh
make -C paper tables                       # every table, then the figure (= bash paper/tables/make.sh)
bash paper/tables/make.sh supp_runtime     # one of them (a script's name without .py)
python3 paper/tables/test_tables.py        # the tests of lib.py (no data needed)
```

Outputs:

- `OUT_DIR/tables/NAME.tex`: the bare table body (`\begin{tabular}` ... `\end{tabular}` and the commands it defines),
  as the paper `\input`s it. The paper's preamble defines `\sysname`, `\sysrelease` and `\best` (booktabs, xcolor).
- `OUT_DIR/tables/values.tsv`: every number of every table and figure as printed, with its source file. Columns
  `table`, `key` (the set), `row`, `column`, `value`, `source` (a path under `OUT_DIR`, or `paper/...`); one file per
  table or figure in `OUT_DIR/tables/values/`, values.tsv remade from them by each script.
- `OUT_DIR/figures/runtime.pdf`.

A value whose input is missing prints as `??` and its script lists the missing files on stderr; a tool that has no run
folder for a set prints as `not run` (in the paper: Rawsamble + miniasm on D7 and D8). `make.sh` exits 1 if a script
failed.

## The scripts

| Script | In the paper | What | Reads (under `OUT_DIR`; KEY = a set of `sets.tsv`) |
|---|---|---|---|
| `datasets.py` | datasets table | D1-D8: genome length, reads, average read length, coverage | `sets.tsv` |
| `contiguity_both.py` | contiguity table | unitigs, longest unitig, NG50, largest component of D1-D8 | `eval/KEY/squash/tier1.tsv`, `eval/KEY/compare/tier1.tsv`, `assembly/KEY/rawsamble/` |
| `correctness_both.py` | correctness table | QUAST of D1-D8: NGA50, misassemblies, genome fraction | `eval/KEY/{squash,compare}/quast/report.tsv` |
| `supp_runtime.py` | Supp. runtime and memory | elapsed time, CPU time, peak memory of the five pipelines on D1-D8 | `KEY.time.txt`, `assembly/KEY/*/logs/`, `timing/KEY/` |
| `supp_commands.py` | Supp. parameters | each tool's parameters (SQUASH's preset per set) | `sets.tsv` |
| `supp_presets.py` | Supp. presets | the values each preset sets | `docs/USAGE.md` |
| `supp_versions.py` | Supp. versions | each tool's version and source | `versions.tsv` |
| `supp_cont30.py`, `supp_corr30.py` | Supp. 30x per strand | the six supplement sets, as the contiguity and correctness tables | as those |
| `supp_athal_masked.py` | Supp. A. thaliana masked | D8 against Col-CC as is and with its tandem arrays masked | `eval/KEY/{squash,compare}/quast/`, `quast_masked/` |
| `supp_dmel_masked.py` | Supp. D. melanogaster masked | D7 against release 6 as is and with its repeat-rich regions masked | as above |
| `supp_hac.py` | Supp. HAC reads | hifiasm on the sup and the HAC reads | `eval/KEY/compare/quast/`, `quast_hac/` |
| `supp_runtime_single.py` | Supp. each tool alone | elapsed time and peak memory of each tool | as `supp_runtime.py` |
| `supp_strands.py` | Supp. strands | the strand of each SQUASH unitig | `eval/KEY/squash/asm/`, `strands/` |
| `fig_runtime.py` | runtime figure | the runtime table's pipelines as bars | as `supp_runtime.py` |
| `lib.py` | | what the scripts share: readers, formats, ranking, the two-panel layout, values.tsv | |

The main tables put the bacteria (`small` preset, D1-D4) and the eukaryotes (`large`, D5-D8) side by side; the
supplement's six sets (D1-D5, D7 at ~30x per strand) take the label of the main set of the same genome. The paper's
overview figure and its table of changes to metaMDBG's assembly rules are drawn and written by hand: not here.

## Inputs

| Path under `OUT_DIR` | What |
|---|---|
| `KEY.time.txt` | SQUASH's run (`squash.sh`): `/usr/bin/time -v` of the whole run |
| `assembly/KEY/TOOL/logs/STEP.time.txt` | the other tools' timed runs (`/usr/bin/time -v` per step, staging not counted): TOOL `nanomdbg`, `flye`, `hifiasm` (sup reads), `hifiasm_hac` (HAC reads), `rawsamble`; STEP `nanomdbg`, `flye`, `hifiasm`, or `rawsamble_index` + `rawsamble_overlap` + `miniasm` (summed) |
| `timing/KEY/basecall_sup_gpu.time.txt`, `basecall_hac_gpu.time.txt` | dorado on the GPU |
| `timing/KEY/basecall_hac_cpu.time.txt`, `basecall_hac_cpu.tsv` | dorado HAC on the CPU: the timed run and its rate (below) |
| `eval/KEY/squash/` | SQUASH's evaluation: `quast/report.tsv` (spelled and placed, forward- and reverse-strand columns), `quast_masked/` (D7, D8), `tier1.tsv` (contiguity from the signal alone), `asm/squash_spelled.fa`, `strands/squash_spelled_{fwd,rev}.fa` |
| `eval/KEY/compare/` | the other assemblers' evaluation: `quast/` (nanoMDBG, Flye, hifiasm on sup reads), `quast_masked/` (D7, D8), `quast_hac/` (hifiasm on HAC reads, D1-D8), `tier1.tsv` (one row per tool: Rawsamble + miniasm, nanoMDBG, Flye, hifiasm) |

Static facts: `sets.tsv` (the sets, their genome length G and signal-estimated bases), `versions.tsv`, `docs/USAGE.md`.
SQUASH's columns and rows in the evaluation's files are named `SQUASH (spelled)`, `SQUASH (spelled; forward-strand
unitigs)` ... (eval/README.md).

**Times.** Every tool ran with 64 threads on a node of its own with two AMD EPYC 7502 CPUs; dorado on one NVIDIA RTX
A5000, or on 64 of those cores. SQUASH's time is its whole run (the POD5 files to the signal layout and the output
files). Elapsed time as h:mm:ss (rounded to the second), CPU time = user + system, memory = the peak resident set of
`/usr/bin/time -v` in GB (10^6 kB). A pipeline that basecalls adds dorado's elapsed time to hifiasm's and takes the
larger peak; GPU basecalling has no CPU time (N/A).

### basecall_hac_cpu.tsv

HAC basecalling on the CPU is timed on the whole set for the bacteria and on a subset of about one hour for the large
sets; the set's time is its signal samples over the samples per second dorado reached. One header line, one row:

| Column | |
|---|---|
| `key` | the set |
| `samples_timed` | the signal samples of the timed reads |
| `dorado_s` | dorado's own time for them (its log's `Finished in (ms)`, in s) |
| `set_samples` | the set's signal samples |
| `set_est_s` | the set's time: `set_samples` x `dorado_s` / `samples_timed`, to the second |

The tables take `set_est_s` as the elapsed time, `set_est_s` x the timed run's CPU share (`Percent of CPU this job got`
of `basecall_hac_cpu.time.txt`) as the CPU time, and the timed run's peak memory.

## Requirements

Python 3.10 or later; the tables need its standard library only, the figures matplotlib 3.10.9 (`requirements.txt`).
The figures use the Liberation Sans font where it is installed (else matplotlib's default font).
