# paper/eval: the paper's evaluation

The scripts here score SQUASH's assembly of a paper set and the comparators' assemblies of the same reads as the
paper does. SQUASH outputs no bases: its final graph is described by the signal that tiles each unitig (its signal
layout, `layout.tsv.gz`), and the evaluation turns that signal into sequence only to score it.

- **Tier 2, QUAST 5.2.0** (defaults, `--no-icarus`) of SQUASH's unitigs written two ways: **spelled** (the sup bases of
  each layout piece's signal, read through dorado's move tables) and **placed** (the reference sequence that each
  piece's events lie on, by Remora's reference-anchored refinement of the reads' sup alignments), both strands together
  and per strand (SQUASH assembles each strand apart: each unitig on the strand with more aligned bases). The
  comparators (nanoMDBG, Flye, hifiasm on the sup reads; hifiasm on the HAC reads) by the same QUAST. Variants: human
  (D6) also against both haplotypes' q-arms; D. melanogaster (D7) and A. thaliana (D8) also with their repeat-rich
  regions masked (the placement redone on the masked reference); A. thaliana's reference holds Col-CC's chromosomes and
  its chloroplast and mitochondrion.
- **Tier 1**, contiguity from the raw signal alone: SQUASH's unitigs and Rawsamble + miniasm's segments in
  signal-estimated bases (samples x 400 bases per second / sample rate) against twice the genome (both strands), the
  basecalled assemblers in bases against the genome.

`eval.sh KEY` runs every step of one set; each step is also a script of its own.

## Inputs

| What | Where | Made by |
|---|---|---|
| SQUASH's outputs and its work folder | `OUT_DIR/KEY/` (`layout.tsv.gz`, `assembly_graph.gfa`; `work/graph/input/`, `work/signal/spans.bin`, `work/bp_scale/signal_lengths.tsv`) | `KEEP_WORK=1 bash squash.sh KEY` (the outputs are the same; the rest of `work/` may be deleted) |
| the set's POD5 | `POD5_DIR/KEY.pod5` | `data/extract.sh KEY` |
| the sup reads with move tables | `OUT_DIR/reads/KEY/sup.moves.bam` (dorado `--emit-moves`) | the basecalling |
| the comparators' assemblies | `OUT_DIR/assembly/KEY/TOOL/asm/`: `nanomdbg/asm/contigs.fasta.gz`, `flye/asm/assembly.fasta` and `assembly_graph.gfa`, `hifiasm/asm/asm.bp.p_ctg.fa`, `hifiasm_hac/asm/asm.bp.p_ctg.fa` (D1-D8), `rawsamble/asm/ava.gfa` (Rawsamble + miniasm: the bacteria, A. rabiei, human) | the comparators' runs |
| the references | `SOURCES_DIR/references/` | `data/download/references.sh` (run by `refs.sh` when missing) |

Why the work folder: a spelled piece (and a piece of tier 1) ends where the next piece takes over, at the first sample
of the take-over token in this read's signal. `layout.tsv.gz` gives each piece's own samples; the take-over token's
first q-symbol is in the graph's input tokens (`work/graph/input/`) and its sample in the events' spans
(`work/signal/spans.bin`), which the placement also reads for every event of the layout's reads.

## Tools

The tools are `tools/env.sh`'s (sourced by `config.sh`; `tools/install.sh` installs the paper's versions):

| Variable | What |
|---|---|
| `QUAST` | QUAST 5.2.0's `quast.py` |
| `MINIMAP2`, `SAMTOOLS` | minimap2 2.30-r1287 (the sup reads' alignments for Remora, the strand split) and samtools |
| `REMORA_PYTHON` | a python with Remora 3.3.0, pysam, numpy and pod5 (`tools/envs/remora.yml`): every Python step here |
| `REMORA_LEVELS` | RawHash2's R10.4.1 9-mer level table `uncalled_r1041_model_only_means.txt` (SHA-256 0b28b79a...ecf2): Remora's k-mer levels |
| `SQUASH` | the squash binary: its signal step gives Rawsamble's reads their events without event merging (tier 1) |

and the evaluation's settings (`config.sh` or the environment; `lib.sh` gives the defaults):

| Variable | What |
|---|---|
| `EVAL_THREADS` | threads of QUAST, minimap2 and Remora (default `THREADS`) |
| `EVAL_TMP` | large temporary files (default `OUT_DIR/eval/tmp`) |
| `EVAL_REFS` | the evaluation references (default `SOURCES_DIR/references/eval`) |

## Steps

| Step | Script | Output (`OUT_DIR/eval/KEY/`) |
|---|---|---|
| `refs` | `refs.sh KEY` (`mask_arrays.py`, `mask_dmel.py`) | the set's references in `EVAL_REFS`, each checked against the paper's (`references.tsv`) |
| `spell` | `spell.sh KEY` (`spell.py`) | `squash/asm/squash_spelled.fa`, `.pieces.tsv` |
| `place` | `place.sh KEY` (`remora_evref.py`, `project.py`) | `squash/asm/squash_placed.fa`, `.blocks.tsv`; `squash/reads.txt`, `squash/evref/` |
| `strands` | `strands.sh KEY` (`strands.py`, `closed.py`) | `squash/strands/` (`spelled.paf`, `squash_{spelled,placed}_{fwd,rev}.fa`, `summary.tsv`), `squash/closed.tsv` |
| `quast` | `quast.sh KEY` | `squash/quast/` (+ `quast_hap/` on D6) |
| `tier1` | `tier1.sh KEY` (`tier1.py squash`) | `squash/tier1.tsv` |
| `masked` | `masked.sh KEY` (D7, D8) | `squash/masked/` (placement and strands on the masked reference), `squash/quast_masked/` |
| `compare_quast` | `compare_quast.sh KEY` | `compare/quast/` (+ `quast_hap/`, `quast_masked/`, `quast_hac/`) |
| `compare_tier1` | `compare_tier1.sh KEY` (`tier1.py compare`) | `compare/tier1.tsv` |

`eval.sh KEY [STEP ...]` runs them in this order (`squash` = the first seven, `compare` = the last two), skipping a step
whose output is there (`FORCE=1` runs it again). Logs: `squash/logs/`, `compare/logs/`.

## Outputs

- `squash/quast/report.tsv`: six columns, `SQUASH (spelled)`, `SQUASH (placed)` (both strands) and
  `SQUASH (spelled; forward-strand unitigs)`, `... reverse-strand unitigs)`, `SQUASH (placed; forward-strand
  unitigs)`, `... reverse-strand unitigs)`; `quast_hap/`, `quast_masked/` the same against the variant references.
  Placed sequence is reference sequence: its identity, mismatches and indels mean nothing, only its structure
  (misassemblies, NGA50, genome fraction).
- `compare/quast*/report.tsv`: `nanoMDBG`, `Flye`, `hifiasm`; `quast_hac/`: `hifiasm (HAC)`. QUAST scores each assembly
  apart, so these equal the columns of one QUAST of every assembly.
- `squash/tier1.tsv`: one row (`key`, then the columns of the table below, `reference_length`, `gaps` = layout gaps,
  `pieces`, `sample_rate`); `compare/tier1.tsv`: one row per tool (Rawsamble + miniasm where it was run, nanoMDBG, Flye,
  hifiasm; then `reference_length`, `bases_per_second`, `sample_rate`). Columns: `assembler`, `unit`, `genome_size` (G),
  `segments`, `total_bp`, `largest_bp`, `N50`, `NG50`, `auN`, `auNG`, `graph`, `graph_segments`, `links`, `circular`,
  `components`, `graph_total_bp`, `largest_component_bp`, `largest_component_share`.
- `squash/strands/summary.tsv`: SQUASH's unitigs and spelled bases on each strand and on neither (no alignment: left out
  of the per-strand QUAST). `squash/closed.tsv`: per reference replicon, whether a circular unitig closes it on the
  forward and on the reverse strand (its alignments cover >= 95% of the replicon, its length 95-105% of it).

## The sets

`sets.tsv` gives each set's references (file names in `EVAL_REFS`): `reference` (QUAST and the strand split),
`placement` (the placement's, when not the reference: human both haplotypes, so reads of either haplotype are placed),
`genome` (tier 1's G, when not the reference: A. thaliana's five chromosomes), `hap` and `masked` (the variants), and
whether Rawsamble + miniasm (`rawsamble`) and hifiasm on HAC reads (`hac`) were run on it. `refs.sh` makes them from the
public references: the E. coli chromosome (CP149500.1, the first record of `ecoli_CP149500.fa`), the S. dysgalactiae
assembly with its record named `u9_UT9728`, the A. thaliana organelles (NCBI NC_000932.1, NC_037304.1) after Col-CC's
chromosomes, Col-CC's tandem arrays (each chromosome's 10-kb windows whose 21-mers are >= 30% high-copy, >= 20 copies,
merged across < 100 kb, clusters >= 100 kb: 24.58 Mb), and D. melanogaster's rDNA locus (UCSC dm6 rmsk rRNA elements
+- 20 kb), chrY and pericentromeric heterochromatin (the unbanded proximal arm ends of UCSC dm6 cytoBandIdeo): 18.0 Mb.

## Disk, memory, time

SQUASH's output folder with its work folder kept: 1.8-2.2 GB for a bacterium at 100x, 26 GB for A. rabiei at 60x (most
of it the correction rounds' tables, which the evaluation does not read). The evaluation streams the spans file and
reads only the layout's reads (a few hundred to ~25,000): a bacterium takes 3-5 min at 8-12 threads and ~1 GB of
memory, A. rabiei at 60x 15 min at 16 threads and 6.2 GB, most of it Remora's placement. QUAST of the larger
eukaryotes' assemblies takes 10-30 min per report at 8 threads.

## Parity with the paper's evaluation

These scripts were checked against the files of the paper's evaluation (QUAST labels and file names aside):
- every evaluation reference made from fresh public downloads (`references.sh`, `refs.sh`) holds the paper's records
  (`references.tsv`); the derived files (`ecoli_chr.fa`, `u9.fa`, `colcc_org.fa`, both masked files and their BEDs)
  are byte-identical to the paper's;
- D1-D4 (`eck100_both100`, `cjk100_both100`, `vpk100_both100`, `u9k100_both100`), from SQUASH's own runs
  (`KEEP_WORK=1`, outputs equal to `expected/`), every step: the spelled and placed FASTA, the strand split's PAF and
  four FASTA byte-identical; `squash/quast/report.tsv` byte-identical but its labels; `squash/tier1.tsv` and
  `closed.tsv` equal in every value; the comparators' QUAST (`quast/`, `quast_hac/`) equal value by value to the
  paper's one QUAST of every assembly, their tier 1 (Rawsamble + miniasm, nanoMDBG, Flye, hifiasm) equal in every
  column;
- A. rabiei at 60x (`ark60_both60`), from SQUASH's own run, every step: the same files byte-identical, QUAST, tier 1
  and the comparators' QUAST and tier 1 (Rawsamble + miniasm included) equal;
- D6, D7, D8 (the QUAST steps on the paper's spelled, placed and strand files; the comparators' runs):
  `squash/quast/`, `quast_hap/` (D6) and `quast_masked/` (D7, D8) equal value by value, `compare/quast*/` likewise,
  the comparators' tier 1 equal.

Remora's refinement is not byte-reproducible on D. melanogaster's reads (4 kHz): a second placement of D7 gives 6,536
of its 6,966 placed unitigs identical (+0.003% bases) and the same QUAST values but the placed genome fraction, which
moves in the third decimal. The bacteria's and A. rabiei's placements are byte-identical.

## Tests

`$REMORA_PYTHON test_eval.py`: the steps' functions on small made-up inputs (no Remora, no data).
