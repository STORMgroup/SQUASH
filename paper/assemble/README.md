# paper/assemble: the other assemblers

One script per assembler, each on one paper set, as the paper's benchmark ran it: the tool alone on a machine of its
own, pinned to the first `ASM_THREADS` (64) CPUs, its input staged to a local disk before the clock, its working and
output folders on that disk, `/usr/bin/time -v` of each step of the tool, the outputs copied afterwards. The reads come
from `../basecall/` (sup reads for every tool, HAC reads for hifiasm's HAC run), Rawsamble's signal from the set's POD5.

Before: the tools (`make -C paper tools`) and the set's reads (`make -C paper basecall SET=KEY`). Run them from the
repository's root. Each line is a separate example (its comment says what it does): run the ones you need.

```
bash paper/assemble/run.sh cjk100_both100              # every tool the paper ran on the set
bash paper/assemble/run.sh main flye                   # Flye on D1-D8
bash paper/assemble/hifiasm.sh h2k60_both60 hac        # hifiasm on D6's HAC reads
DRY_RUN=1 bash paper/assemble/rawsamble.sh KEY         # print the timed commands only
```

| Script | Tool (version) | Timed steps: command (LOCAL = the working folder, STAGE = the staged input's) | Output in `OUT_DIR/assembly/KEY/TOOL/` |
|---|---|---|---|
| `nanomdbg.sh KEY` | nanoMDBG = metaMDBG 1.2, minimap2 2.30-r1287 first on its PATH | `nanomdbg`: `metaMDBG asm --out-dir LOCAL/asm --in-ont LOCAL/reads.fastq --threads 64` | `nanomdbg/asm/contigs.fasta.gz`, `asm/metaMDBG.log` |
| `flye.sh KEY` | Flye 2.9.6-b1802 | `flye`: `flye --nano-hq LOCAL/reads.fastq -g SIZE -t 64 -o LOCAL/asm` | `flye/asm/assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt`, `flye.log` |
| `hifiasm.sh KEY [sup\|hac]` | hifiasm 0.25.0-r726 | `hifiasm` (in LOCAL/asm): `hifiasm --ont -t 64 -l0 -o asm LOCAL/reads.fastq`; D6 without `-l0` | `hifiasm/` or `hifiasm_hac/`: `asm/asm.bp.p_ctg.fa` (primary contigs, the S lines of `asm.bp.p_ctg.gfa`), `asm.bp.p_ctg.gfa`, `asm.bp.*.noseq.gfa`, haplotype / scaffold FASTA when written |
| `rawsamble.sh KEY` | Rawsamble = RawHash2 2.5 (commit 92366af) + miniasm 0.3-r179 | `rawsamble_index`: `rawhash2 -x r10ava -t 64 -p PORE -d LOCAL/ava.ind STAGE/in.pod5`; `rawsamble_overlap`: `rawhash2 -x r10ava -t 64 LOCAL/ava.ind STAGE/in.pod5 > STAGE/ava.paf`; `miniasm`: `miniasm STAGE/ava.paf > LOCAL/ava.gfa` | `rawsamble/asm/ava.gfa` (the graph), `asm/ava.paf.gz` (the overlaps) |
| `run.sh KEY\|main\|supp\|all [TOOL ...]` | the parent: the named tools (`nanomdbg`, `flye`, `hifiasm`, `hifiasm_hac`, `rawsamble`), default the ones the paper ran on the set | | |

Every run also writes `logs/STEP.time.txt` (one per timed step; `staging.time.txt` for the untimed staging), the tool's
log (`logs/nanomdbg.log`, `flye.log`, `hifiasm.log`, `rawsamble_*.err`, `miniasm.err`), `stats.tsv` (contigs or
segments, bases, SHA-256 of the main output) and `done`. A run with `done` is skipped; remove its folder to run it again.

`options.tsv`, per set: Flye's `-g` (the genome length in kb the paper's runs passed: 5130k for E. coli, elsewhere G of
`../sets.tsv` in whole kb; Flye 2.9.6 reads `-g` only with `--asm-coverage`, so it does not change the assembly),
hifiasm's options (`-l0`, no duplicate purging, as the hifiasm README advises for haploid and inbred genomes; none for
D6) and the tools the paper ran (hifiasm on HAC reads for D1-D8 only; Rawsamble + miniasm on every set but D7, dmk60 and
D8, where its all-vs-all overlaps would take days).

Settings (`config.sh`, environment): `ASM_THREADS` (64), `LOCAL_DIR` (the working folders; default `TMPDIR`),
`STAGE_DIR` (the staged input; default `LOCAL_DIR`; the paper staged the large sets' inputs in RAM, `/dev/shm`),
`SUP_FASTQ`, `HAC_FASTQ`, `SET_POD5` (the inputs; defaults `OUT_DIR/reads/KEY/sup.fastq.gz`, `hac.fastq.gz`,
`POD5_DIR/KEY.pod5`), `DRY_RUN=1`, and the tool variables of `../tools/env.sh`. The paper's runs peaked at 77 GB of
memory (Rawsamble's overlaps on D6; Flye 59-67 GB on D7 and D8, the others below 35 GB) and took up to 21 hours (the same
step); the local disk holds the staged input and the tool's working folder.

## Reproducibility

Rerun through these scripts from the paper's own reads of cjk100_both100 (D2; 64 threads on 16 cores):
- hifiasm (sup and HAC reads): `asm.bp.p_ctg.fa` and the GFAs equal the paper's byte for byte;
- Rawsamble + miniasm: `ava.gfa` equals the paper's byte for byte; the overlaps have the same 500,463 lines (PAF
  columns 1-12) in another order, with another mapping-time tag (`mt:f`);
- Flye: one circular contig of 1,768,446 bp as the paper's, `assembly_info.txt` equal; its sequence was the paper's
  started elsewhere in one run, and in another run reverse-complemented with 5 substitutions (Flye is not reproducible
  run to run at the base level);
- nanoMDBG: not reproducible run to run (three runs, one of them with the paper's own metaMDBG binary: 28 contigs each,
  the paper's run 29; 2.52-2.53 Mb in all; the chromosome contig 1,768,451-1,768,465 bp; 21-23 of the paper's 29 contigs
  identical).

Rawsamble reads the set's POD5 as `make data` rebuilds it (the paper's reads in the paper's order, checked by its signal
digest). The other tools read basecalled reads: basecalled again (`../basecall/`) these are the paper's records
(checked on 200 reads on the same GPU model; D6's FASTQ is ONT's, byte for byte), but in Dorado's output order, so an
assembly can differ from the paper's where the tool depends on the read order, besides the run-to-run variation above.
