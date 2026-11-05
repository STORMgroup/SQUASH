# paper/: reproducing the paper

The paper's reads are not re-hosted: each of its 14 read sets is a list of read ids from a public nanopore run. The
scripts here download a set's run, cut the listed reads into one POD5 in the order the paper's run read them, check it
against the paper's POD5 (read ids and a digest of every signal sample), run SQUASH on it as the paper did and compare
the outputs with the paper's. They also run the other tools of the paper (basecalling, nanoMDBG, Flye, hifiasm,
Rawsamble + miniasm), the evaluation and the tables.

Every command below is run from the repository's root folder (`SQUASH/`).

## Before you start

1. Install SQUASH: the [README's Install section](../README.md#install), steps 1-3.
2. Only if you do not have conda yet: install Miniforge (into your home folder):

   ```sh
   curl -L -O https://github.com/conda-forge/miniforge/releases/download/26.9.4-1/Miniforge3-26.9.4-1-Linux-x86_64.sh
   bash Miniforge3-26.9.4-1-Linux-x86_64.sh -b -p $HOME/miniforge3
   source $HOME/miniforge3/bin/activate
   ```

3. The Python packages of these scripts, at the paper's versions (once):

   ```sh
   conda env create -f paper/environment.yml
   conda activate squash-paper
   ```

   In every new terminal, `conda activate squash-paper` brings it back (after `source $HOME/miniforge3/bin/activate`
   if you installed Miniforge in step 2). With pip instead of conda: `pip install -r paper/requirements.txt` (Python
   3.10 or later).
4. Check the free disk space where the data goes (`df -h .`): by default the downloads, the set POD5s and the outputs
   go inside `paper/` (`paper/data/sources`, `paper/data/pod5`, `paper/out`). To put them elsewhere, edit these lines
   of `paper/config.sh`, for example:

   ```sh
   : "${SOURCES_DIR:=/big/disk/squash/sources}"
   : "${POD5_DIR:=/big/disk/squash/pod5}"
   : "${OUT_DIR:=/big/disk/squash/out}"
   ```

## Reproduce one dataset

```sh
make -C paper check SET=cjk60_both60
```

This downloads the C. jejuni deposit (42 GB, about an hour; an interrupted download continues where it stopped when you
run the command again), cuts the set's 6,177 reads into one POD5, checks it, runs SQUASH (minutes) and compares its six
output files with the paper's. It ends with `check.sh: cjk60_both60: all 6 files equal the paper's`. The download needs
~90 GB of free disk while the archive is unpacked; 42 GB stay.

Three more checks, each optional and independent of the others:

- SQUASH again with 4 threads, which gives the same outputs:

  ```sh
  make -C paper threads SET=cjk60_both60
  ```

- the set's read list remade from the raw signal of the deposit:

  ```sh
  make -C paper curate SET=cjk60_both60
  ```

- the tests of these scripts (no download):

  ```sh
  make -C paper test
  ```

## Reproduce everything

`SET` is one set's key (the table below), `main` (the paper's datasets D1-D8), `supp` (the supplement's six) or `all`.
Run the steps in this order. Finished downloads, set POD5s, tool installs, basecalls and assembler runs are not
repeated when a step runs again (SQUASH's runs are).

| Step | Command | Needs | Makes |
|---|---|---|---|
| 1. Data | `make -C paper data SET=main` | the downloads of the table below | `POD5_DIR/KEY.pod5`, each checked against the paper's |
| 2. SQUASH | `KEEP_WORK=1 make -C paper check SET=main` | up to 30 GB of memory per set (the README's performance table); with `KEEP_WORK=1` (the evaluation needs the work folder) about the set POD5's size again on disk | `OUT_DIR/KEY/`, compared with the paper's outputs |
| 3. Read lists | `make -C paper curate SET=main` | the sources' downloads (D. melanogaster and human: more of their reads than step 1 keeps, [curation/](curation/README.md)) | each set's list remade from the raw signal and compared |
| 4. Other tools | `make -C paper tools` | conda (Before you start), ~21 GB of disk | Dorado, nanoMDBG, Flye, hifiasm, Rawsamble + miniasm, minimap2, samtools, QUAST and Remora at the paper's versions in `TOOLS_DIR` ([tools/](tools/README.md)) |
| 5. Basecalls | `make -C paper basecall SET=main` | an NVIDIA GPU | the sup reads with move tables and the HAC reads in `OUT_DIR/reads/KEY/` ([basecall/](basecall/README.md)) |
| 6. Basecalling times | `make -C paper timing SET=main` | a GPU; 64 CPU cores and ~190 GB of memory for the CPU route | `OUT_DIR/timing/KEY/` |
| 7. Other assemblers | `make -C paper assemble SET=main` | 64 CPU cores, up to ~80 GB of memory | nanoMDBG, Flye, hifiasm and Rawsamble + miniasm in `OUT_DIR/assembly/KEY/`, timed ([assemble/](assemble/README.md)) |
| 8. Evaluation | `make -C paper eval SET=main` | steps 2, 5 and 7; QUAST and Remora of step 4; up to ~1 h per large set | SQUASH's unitigs spelled and placed, QUAST, contiguity from the signal alone, and the other assemblers' QUAST, in `OUT_DIR/eval/KEY/` ([eval/](eval/README.md)) |
| 9. Tables | `make -C paper tables` | steps 2-8 | the paper's tables and runtime figure in `OUT_DIR/tables/`, `OUT_DIR/figures/` ([tables/](tables/README.md)) |

The same steps as commands, in order (each line is one step of the table; replace `main` by one set's key to do one
set):

```sh
make -C paper data SET=main
KEEP_WORK=1 make -C paper check SET=main
make -C paper curate SET=main
make -C paper tools
make -C paper basecall SET=main        # on a machine with an NVIDIA GPU
make -C paper timing SET=main          # on a machine with an NVIDIA GPU (and 64 cores for the CPU route)
make -C paper assemble SET=main
make -C paper eval SET=main
make -C paper tables
```

Steps 5 and 6 need a GPU; the others run on CPUs. When they run on different machines, put the repository and
`OUT_DIR` on storage both can read, and run steps 4 and 7 on the same kind of CPU: RawHash2 (Rawsamble) is compiled
for the CPU it is built on.

The paper timed every tool alone on a machine of its own with 64 threads (two AMD EPYC 7502 CPUs; Dorado on one NVIDIA
RTX A5000). Other machines give other times, but the same SQUASH outputs.

## The sets

`sets.tsv` has one row per set: the source, the pool and rule it was selected by (from the raw signal alone: reads of
at least 10 kb by signal-estimated length, a noise cut, then reads in md5 order of their id until 60x or 100x of the
genome, so 30x or 50x per strand), the genome size G, the preset, the read counts and signal-estimated bases, the
SHA-256s of the lists, the POD5 signal digest, the SHA-256 of the q-symbols table SQUASH writes
(`work/signal/qsymbols.tsv`) and the QUAST reference. `data/sets/KEY/README.md` describes each set.

| Key | Paper | Genome | G (Mb) | Coverage | Preset | Reads | Assembled | Source | Download (GB) |
|---|---|---|---|---|---|---|---|---|---|
| `eck100_both100` | D1 | Escherichia coli ATCC 25922 | 5.13 | 100x | `small` | 26,419 | | nvb_ec | 125 |
| `cjk100_both100` | D2 | Campylobacter jejuni ATCC 33560 | 1.77 | 100x | `small` | 10,285 | | nvb_cj | 42 |
| `vpk100_both100` | D3 | Vibrio parahaemolyticus ATCC 17802 | 5.15 | 100x | `small` | 27,174 | | nvb_vp | 85 |
| `u9k100_both100` | D4 | Streptococcus dysgalactiae subsp. equisimilis UT9728 | 2.09 | 100x | `small` | 10,262 | | sdys_zenodo17428855 | 200 |
| `ark100_both100` | D5 | Ascochyta rabiei ME14 (seven strains pooled) | 40.85 | 100x | `large` | 193,335 | | arabiei_prjeb86315 | 167 |
| `h2k60_both60` | D6 | Homo sapiens HG002, chromosome 13 q-arm | 96.10 | 60x | `large` | 206,809 | | hg002_giab2025 | 3,284 (streamed) |
| `dmk100_both100` | D7 | Drosophila melanogaster ISO1 | 137.55 | 100x | `large` | 557,254 | 550,874 | dmel_bkim2023 | 246 (streamed) |
| `awk60_both60` | D8 | Arabidopsis thaliana Col-0, whole genome | 142.48 | 60x | `large` | 250,912 | 250,544 | athal_err17795792 | 164 |
| `eck60_both60` | supp | Escherichia coli ATCC 25922 | 5.13 | 60x | `small` | 15,840 | | nvb_ec | 125 |
| `cjk60_both60` | supp | Campylobacter jejuni ATCC 33560 | 1.77 | 60x | `small` | 6,177 | | nvb_cj | 42 |
| `vpk60_both60` | supp | Vibrio parahaemolyticus ATCC 17802 | 5.15 | 60x | `small` | 16,193 | | nvb_vp | 85 |
| `u9k60_both60` | supp | Streptococcus dysgalactiae subsp. equisimilis UT9728 | 2.09 | 60x | `small` | 6,172 | | sdys_zenodo17428855 | 200 |
| `dmk60_both60` | supp | Drosophila melanogaster ISO1 | 137.55 | 60x | `large` | 334,076 | 330,203 | dmel_bkim2023 | 246 (streamed) |
| `ark60_both60` | supp | Ascochyta rabiei ME14 (seven strains pooled) | 40.85 | 60x | `large` | 116,092 | | arabiei_prjeb86315 | 167 |

Sets sharing a source download it once. All 14 sets: ~4.31 TB of downloads (the human run 3.28 TB, streamed: each
file is cut to the listed reads and deleted), ~450 GB of set POD5s. **Assembled**: the paper's runs of these three sets
assembled only the reads with a basecalled record under their own id; the POD5 keeps every listed read (the signal step
reads POD5s in blocks of 40,000 reads of file order, so dropping reads from the POD5 would move the blocks), and
`squash.sh` passes the assembled list as `--read-ids`.

Every list can be remade from its public source (step 3): `curation/select_reads.py` applies the rule to the set's
pool, read from the POD5s. The pools of u9k60, u9k100 (which barcodes are UT9728's) and h2k60 (which reads are the
chr13 q-arm's) took basecalls or a reference to make; they are committed as read-id lists in `data/pools/` with the
scripts that made them.

## Settings

Every script reads `config.sh` (the only file to edit); a value set in the environment wins, for example
`THREADS=16 make -C paper check SET=cjk60_both60`. The main settings:

| Setting | Default | What |
|---|---|---|
| `SOURCES_DIR`, `POD5_DIR`, `OUT_DIR` | `paper/data/sources`, `paper/data/pod5`, `paper/out` | the downloads, the set POD5s, every output |
| `TMP_DIR` | empty (inside `OUT_DIR`) | SQUASH's large temporary files: a fast local disk helps on the large sets |
| `THREADS` | every core | SQUASH's threads (its outputs do not depend on them) |
| `PYTHON` | `python3` | a Python with `requirements.txt`'s packages (the conda environment above) |
| `SRC_<SOURCE>` | unset | a copy of a source you already have, instead of downloading it (e.g. `SRC_NVB_CJ=/data/ATCC_33560__202309`; [data/download/](data/download/README.md)) |
| `TOOLS_DIR` and one variable per tool | `paper/tools/opt` | where step 4 installs the other tools ([tools/](tools/README.md)) |

## Expected outputs and parity

`expected/KEY.sha256` holds SQUASH's output files of each set (a `.gz` file by its decompressed text, since gzip bytes
depend on the compressor), recorded from runs checked against the paper's runs (`PARITY.md`). The outputs do not
depend on the thread count.

## Files

| Path | What |
|---|---|
| `config.sh` | the settings (the only file to edit; `config.local.sh`, git-ignored, for one machine's) |
| `requirements.txt`, `environment.yml` | the pinned Python packages of the scripts (pip or conda) |
| `versions.tsv` | every tool, model and reference the paper used, with its version and source |
| `Makefile` | `data`, `squash`, `check`, `threads`, `curate`, `tools`, `basecall`, `timing`, `assemble`, `eval`, `eval-refs`, `tables`, `test` |
| `sets.tsv` | one row per set (above) |
| `data/` | the read lists, sources, download scripts, `extract.sh`, `pod5_cut.py`, `pod5_digest.py` and their tests ([data/README.md](data/README.md)) |
| `curation/` | how the read lists were selected: `select_reads.py` (the raw signal alone), `remake.sh` (a set's list remade from its source and checked), the scripts that made the two pool lists of `data/pools/`, tests ([curation/README.md](curation/README.md)) |
| `squash.sh` | SQUASH on one set as the paper ran it |
| `check.sh` | a set's outputs against `expected/` |
| `expected/` | the outputs of the sets checked against the paper's runs (`PARITY.md`) |
| `PARITY.md` | how SQUASH was matched to the paper's runs |
| `tools/` | the other tools' installs, pinned to the paper's versions, and `env.sh` (each tool's variable) ([tools/README.md](tools/README.md)) |
| `basecall/` | Dorado's sup and HAC reads of a set and the timed basecalling routes ([basecall/README.md](basecall/README.md)) |
| `assemble/` | nanoMDBG, Flye, hifiasm and Rawsamble + miniasm on a set, timed ([assemble/README.md](assemble/README.md)) |
| `eval/` | the evaluation: spelled and placed sequences, QUAST 5.2.0, contiguity from the signal alone, the references it needs ([eval/README.md](eval/README.md)) |
| `tables/` | one script per table and data figure of the paper, from the outputs under `OUT_DIR` ([tables/README.md](tables/README.md)) |
