# paper/tools: the other tools, pinned to the paper's versions

`install.sh` installs every tool that `basecall/`, `assemble/` and the evaluation run, at the versions of
`../versions.tsv`, into `TOOLS_DIR` (default `paper/tools/opt`), then `check.sh` checks that each one prints the
paper's version and that the downloaded models and tables equal the paper's files. `env.sh` (sourced by every script
after `config.sh`) names each tool by a variable: a value set in the environment, `config.sh` or `config.local.sh`
wins, else the copy in `TOOLS_DIR`, else the tool's name on `PATH`. So a tool you already have at the paper's version can
be used by setting its variable.

Before: conda (or mamba) on your `PATH`, e.g. the Miniforge of [../README.md](../README.md#before-you-start). Then run
them from the repository's root. Each line is a separate example (its comment says what it does): run the ones you need.

```
make -C paper tools                    # the same as bash paper/tools/install.sh
bash paper/tools/install.sh            # everything (about 21 GB; Dorado's archives are fetched, unpacked and deleted)
bash paper/tools/install.sh flye quast # only these
bash paper/tools/check.sh              # versions against the paper's; exit 1 when any differs
```

Needs bash, curl, git, tar, sha256sum and conda, mamba or micromamba (`CONDA_CMD`). Dorado runs on Linux x64 with an
NVIDIA GPU (and on the CPU for the CPU route); the conda environments are linux-64.

| Tool | Version (`check.sh` expects) | How | Variable | Used by |
|---|---|---|---|---|
| Dorado | 0.9.2, 1.4.0, 2.1.2 | ONT's Linux x64 archives (`dorado.sh`, `dorado.tsv`): archive SHA-256 (where recorded), `bin/dorado` SHA-256 equal to the paper's | `DORADO_0_9_2`, `DORADO_1_4_0`, `DORADO_2_1_2` | `basecall/` |
| Dorado models | sup@v5.2.0, hac@v5.2.0, sup@v5.0.0 (1.4.0); sup@v4.1.0, hac@v4.1.0 (0.9.2) | `dorado download` of the release that lists them; every file's SHA-256 equal to the paper's (`models.sha256`) | `DORADO_MODELS` | `basecall/` |
| nanoMDBG (metaMDBG) | 1.2 | source, tag `metaMDBG-1.2` (`metamdbg.sh`; bioconda's 1.2 package requires minimap2 2.28, the paper's runs had 2.30) | `METAMDBG` | `assemble/nanomdbg.sh` |
| Flye | 2.9.6-b1802 | bioconda (`envs/flye.yml`) | `FLYE` | `assemble/flye.sh` |
| hifiasm | 0.25.0-r726 | bioconda (`envs/hifiasm.yml`) | `HIFIASM` | `assemble/hifiasm.sh` |
| Rawsamble (RawHash2) | 2.5, commit 92366af | source, STORMgroup/RawHash2 at that commit, CMake Release, POD5 only (`rawhash2.sh`) | `RAWHASH2` | `assemble/rawsamble.sh` |
| R10.4.1 9-mer pore model | `extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt` of that checkout (SHA-256 e5e15451...) | with RawHash2 | `RAWHASH2_PORE` | Rawsamble's index (`-p`) |
| R10.4.1 9-mer levels | `extern/local_kmer_models/uncalled_r1041_model_only_means.txt` of that checkout (SHA-256 0b28b79a...) | with RawHash2 | `REMORA_LEVELS` | the evaluation (Remora's refinement) |
| miniasm | 0.3-r179 | bioconda (`envs/miniasm.yml`) | `MINIASM` | `assemble/rawsamble.sh` |
| minimap2 | 2.30-r1287 | bioconda (`envs/minimap2.yml`) | `MINIMAP2` | nanoMDBG's PATH; the evaluation; curation |
| samtools | 1.23.1 (htslib 1.23.1) | bioconda (`envs/samtools.yml`) | `SAMTOOLS` | `basecall/`; curation |
| QUAST | 5.2.0 | bioconda, with its minimap2 2.28 (`envs/quast.yml`) | `QUAST` | the evaluation |
| Remora | 3.3.0 | PyPI `ont-remora` in Python 3.10.14, the paper's package versions (`envs/remora.yml`; torch's CPU build) | `REMORA_PYTHON` | the evaluation |
| GCC 15.1.0, CMake 3.31.8 | | conda-forge (`envs/build.yml`): the compilers of the two source builds | | `rawhash2.sh`, `metamdbg.sh` |

Checked on a fresh install (`check.sh`): every tool prints the paper's version, the three `bin/dorado` binaries are the
paper's byte for byte, the 412 files of the five models and the two 9-mer tables equal the paper's. The conda files pin
each tool's version and build string as the paper's environments had them. RawHash2 builds with
`-march=native` (its CMake default): build it on the kind of machine that runs it.

## Files

| File | What |
|---|---|
| `install.sh [TOOL ...]` | the parent: every tool (or the named ones) into `TOOLS_DIR`, then `check.sh` |
| `dorado.sh VERSION` | one Dorado release and its models, checked against `dorado.tsv` and `models.sha256` |
| `conda.sh NAME` | the environment `envs/NAME.yml` into `TOOLS_DIR/envs/NAME` |
| `rawhash2.sh`, `metamdbg.sh` | the two source builds (after `conda.sh build`) |
| `check.sh [TOOL ...]` | each tool's printed version against the paper's, the models' and tables' SHA-256s; writes `TOOLS_DIR/versions.tsv` |
| `env.sh` | the tool variables (sourced after `config.sh`) |
| `dorado.tsv` | per Dorado release: archive URL and SHA-256, the SHA-256 of the paper's `bin/dorado`, its models |
| `models.sha256` | the SHA-256 of every file of the five models, as the paper's runs had them |
| `envs/*.yml` | the conda environments |
