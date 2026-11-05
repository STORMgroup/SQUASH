# paper/basecall: the basecalls of the comparison

The other assemblers read basecalled reads, and the evaluation reads Dorado's move tables. The scripts here make them
from a set's POD5 (`POD5_DIR/KEY.pod5`, `make -C paper data SET=KEY`) with the Dorado release and model of the
paper's runs (`sets.tsv`; the tools: `../tools/`), and time the basecalling as the runtime tables did.

Before: the set's POD5 (`make -C paper data SET=KEY`), the tools (`make -C paper tools`), and an NVIDIA GPU with its
driver (`nvidia-smi` lists it). Run them from the repository's root. Each line is a separate example (its comment says
what it does): run the ones you need.

```
bash paper/basecall/sup.sh cjk100_both100            # OUT_DIR/reads/KEY/sup.moves.bam, sup.fastq.gz
bash paper/basecall/hac.sh cjk100_both100            # OUT_DIR/reads/KEY/hac.fastq.gz
bash paper/basecall/timed.sh cjk100_both100 sup gpu  # OUT_DIR/timing/KEY/basecall_sup_gpu.*
bash paper/basecall/run.sh main                      # sup.sh + hac.sh of D1-D8 (TIMED=1: and the three timed routes)
DRY_RUN=1 bash paper/basecall/sup.sh KEY             # print the Dorado command only
```

## Scripts

| Script | Does | Writes |
|---|---|---|
| `sup.sh KEY [moves\|fastq]` | `dorado basecaller --emit-moves` (read splitting off where the paper turned it off): unaligned records with move tables; their primary records (`samtools fastq -F 0x900`) as the assemblers' FASTQ (`@id`, sequence, `+`, qualities). Human (D6): the FASTQ is ONT's own sup v5.0.0 basecalls (below), only the move tables come from Dorado | `OUT_DIR/reads/KEY/sup.moves.bam` (the evaluation's input), `sup.fastq.gz` (the assemblers' input), `sup.stats.tsv`, `logs/sup.dorado.log`, `logs/sup.time.txt` |
| `hac.sh KEY` | `dorado basecaller --emit-fastq`, HAC model, default read splitting: hifiasm's HAC run | `OUT_DIR/reads/KEY/hac.fastq.gz`, `hac.stats.tsv`, `logs/hac.*` |
| `timed.sh KEY sup\|hac gpu\|cpu` | the runtime tables' basecalling: FASTQ only, Dorado / model / POD5 staged to `LOCAL_DIR` before the clock, `/usr/bin/time -v` of Dorado alone; `cpu`: `--device cpu` on the first `BASECALL_CPU_THREADS` (64) CPUs by `taskset`, under `timeout CPU_CAP` (3h), the large sets on a one-hour subset (below) | `OUT_DIR/timing/KEY/basecall_KIND_DEV.time.txt`, `.dorado.log`, `.tsv`; cpu also `KEY.cpu_subset.ids.txt`, `.tsv` |
| `cpu_subset.py KEY TARGET SET.pod5 OUT_DIR` | the CPU route's subset: the set's reads in MD5 order of their id until TARGET samples (the read crossing it kept); `all` = every read | the subset's ids and counts |
| `test_cpu_subset.py` | its tests, on a synthetic POD5 (`$PYTHON basecall/test_cpu_subset.py`) | |
| `run.sh KEY\|main\|supp\|all ...` | the parent: `sup.sh` and `hac.sh` per set; `TIMED=1` also `timed.sh` sup gpu, hac gpu, hac cpu (`TIMED_ONLY=1`: only those) | |
| `lib.sh` | shared settings and helpers | |
| `sets.tsv` | per set: sample rate; Dorado release, model, read splitting and FASTQ source of the sup reads; the HAC reads' release and model; the timed runs' release and models; the CPU subset's samples; the SHA-256 of the sup FASTQ's text where it is the paper's byte for byte (D6) | |

The paper's runtime tables use `basecall_sup_gpu` and `basecall_hac_gpu` (wall time of `.time.txt`, peak memory) and
`basecall_hac_cpu` (`set_est_s` of `.tsv`: the set's samples / the subset's samples per second, where the rate is the
subset's samples over Dorado's "Finished in" time; the bacteria are basecalled whole, so measured).

## Releases and models (`sets.tsv`)

| Sets | sup reads + move tables | HAC reads, timed GPU / CPU basecalling | CPU subset |
|---|---|---|---|
| D1-D4, eck60, cjk60, vpk60, u9k60 (5 kHz) | Dorado 1.4.0 `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, `--disable-read-splitting` | Dorado 2.1.2 `hac@v5.2.0` / `sup@v5.2.0` | whole set |
| D5, ark60, D8 (5 kHz) | the same | the same | 3.67e9 samples (~1 h) |
| D6 human (5 kHz) | reads: ONT's sup v5.0.0 basecalls (giab_2025.01 `calls.sorted.bam` of PAW70337 and PAW71238): the set's chr13 primary records (`samtools view -F 0x904`, region chr13, read through the BAM index over HTTPS, ~4% of the 329 GB), PAW70337 first, in the BAMs' order: the paper's FASTQ byte for byte (`sup.sh` checks it); move tables: Dorado 1.4.0 `sup@v5.0.0` (ONT's model), `--disable-read-splitting` | the same as D5 | 3.67e9 samples |
| D7, dmk60 D. melanogaster (4 kHz, which Dorado 1.0 and later reject) | Dorado 0.9.2 `sup@v4.1.0`, its default read splitting (a split read's pieces carry their own ids and `pi:Z:<read id>`; the FASTQ keeps them, as the paper's did) | Dorado 0.9.2 `hac@v4.1.0` / `sup@v4.1.0` | 2.33e9 samples (~1 h) |

## Compared with the paper's runs

The Dorado commands are the paper's, with these differences:
- sup: the paper basecalled each source's reads in larger batches than one set (a pool of candidate reads; for the
  bacteria and A. rabiei aligned to the genome during basecalling with `--reference`, then reset to unaligned records
  with `samtools reset`, which restores each read's sequence and move table) and took the set's records; here Dorado
  reads the set's POD5 alone, without `--reference`. The records keep Dorado's `TS` tag (which `samtools reset` removed)
  and lack the alignment's leftover tags (`ms`, `nn`, `de`, `tp`, `cm`, `s1`, `s2`, `rl`). Records come in Dorado's output
  order, which varies run to run (the paper's FASTQs too).
- The timed GPU basecalling of D6, D7, D8 and dmk60 read the set's reads from the source POD5s with `--read-ids`; here the
  set's POD5 is staged like the other sets'.

Checked on 200 reads of cjk100_both100 (one RTX A5000): every sup record equals the paper's in sequence, qualities,
move table (`mv`, `ts`, `ns`) and read tags (`qs`, `du`, `sm`, `sd`, `mx`, `ch`, `st`, `rn`), both for the 100 reads the paper basecalled on an RTX A5000 and for
the 100 it basecalled on an RTX A6000, and so do the 200 FASTQ records; the HAC reads (253 records, read splitting
included) equal the paper's HAC reads. Other GPU models were not checked.

The CPU route used 150-190 GB of memory with 64 threads in the paper's runs.

Settings (`config.sh`, environment): `POD5_DIR`, `OUT_DIR`, `LOCAL_DIR` (default `TMPDIR`), `BC_THREADS` (samtools /
pigz threads, 8), `BASECALL_CPU_THREADS` (64), `CPU_CAP` (3h), `CUT_THREADS`, `PYTHON` (pod5 0.3.39), `CURL`, and the
tool variables of `../tools/env.sh`.
