# paper/data/download: one script per source

`SOURCE.sh` puts source SOURCE (a row of `../sources.tsv`) in `SOURCES_DIR/SOURCE/` (default `paper/data/sources`).
`../extract.sh` calls them; they also run alone, from the repository's root. Each line is a separate example (its
comment says what it does): run the ones you need.

```sh
bash paper/data/download/nvb_cj.sh             # one source (41.7 GB); run it again to resume a stopped download
STREAM=1 bash paper/data/download/nvb_cj.sh    # unpack while downloading: half the disk, no resume
bash paper/data/download/references.sh         # the references of the evaluation (QUAST)
```

Shared code in `lib.sh`:
- resumable: each file is fetched with `curl -C -` into `FILE.part`, checked (size; MD5 where the archive gives one) and
  renamed; a finished source has `.complete` and is skipped;
- an existing copy instead of a download: `SRC_<SOURCE>=FOLDER` (upper case, `-` -> `_`, e.g. `SRC_NVB_EC`), a folder of
  the source's POD5 files; the script links `SOURCES_DIR/SOURCE` to it and checks nothing (it is your copy);
- `STREAM=1` for the tar sources: unpack while downloading (MD5 and size taken on the stream; half the disk, no resume);
- `PAR=N` parallel files (Zenodo, D. melanogaster, human; default 4); `PYTHON` = a python with pod5 0.3.39.

| Script | Source | Download | Check |
|---|---|---|---|
| `nvb_ec.sh`, `nvb_cj.sh`, `nvb_vp.sh` | NanoVarBench 202309 samples (figshare), one tar each, POD5s flattened into the folder | 125.0, 41.7, 85.0 GB | figshare MD5 |
| `sdys_zenodo17428855.sh` | S. dysgalactiae run FBA61885 (Zenodo), 88 POD5 + 7 zips | 199.9 GB | Zenodo MD5 per file (`manifests/`) |
| `arabiei_prjeb86315.sh` | A. rabiei ME14, seven strain tarballs (ENA), one folder per strain | 167.0 GB | ENA MD5 per file (`manifests/`) |
| `athal_err17795792.sh` | A. thaliana run ERR17795792 (ENA), one tar.gz | 163.8 GB | ENA MD5 |
| `dmel_bkim2023.sh` | D. melanogaster, ONT open data: 562 FAST5, each converted (`pod5 convert fast5 --threads 1`, pod5 0.3.39) and deleted; `KEEP="LIST ..."` keeps only the listed reads per file (`../extract.sh` passes the dm sets' lists) | 245.8 GB | sizes (`manifests/`; S3 lists no MD5) |
| `hg002_giab2025.sh LIST ...` | human HG002, ONT open data: 408 POD5, streamed (each file cut to the listed reads, then deleted) | 3.28 TB | sizes (`manifests/`) |
| `references.sh [REF ...]` | the evaluation references (NCBI efetch / Datasets, Zenodo, HG002 Q100 v1.2), reduced to the scored records | small (HG002 3 GB) | every record's length and sequence SHA-256 (`manifests/references.tsv`, from the paper's files) |
