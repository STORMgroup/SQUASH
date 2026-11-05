# The method in one page, and where it is in the code

SQUASH assembles nanopore reads from their raw current, without basecalling. It follows nanoMDBG's design (correction
in minimizer space, then metaMDBG's multi-k assembly) on signal-derived tokens instead of base minimizers.

1. **Signal to tokens** (`src/signal.rs`, `src/events.rs`). Each read's POD5 signal, in pA, is cut into events by
   RawHash2's R10.4.1 event detector (windows 3 and 6, thresholds 5.5 and 3.5, peak 0.2); consecutive events whose
   means stay within 0.30 z-units are merged. Event means are quantized to 3-bit q-symbols (0-7). A token hashes 9
   consecutive q-symbols (runs of equal symbols collapsed); one token per window of 16 hashes is kept (a minimizer), so
   each read becomes an ordered list of tokens. Reads are processed in blocks of 40,000 of file order, each block by
   read id.
2. **Correction in token space** (`src/partners/`, `src/rounds/`). Reads whose most frequent
   token is at least 20% of their tokens are left out of correction (low complexity). For every read, partners are
   found by chaining shared tokens (at least 4 anchors, an implied overlap of at least 500 q-symbols, the chain covering
   half of it; the most frequent 0.01% of token values are not seeds, nanoMDBG's repeat rule), taken best first under
   a coverage cap of 20; the read's tokens are then replaced by a token-level partial-order consensus of its partners
   (band 256). Partners are searched again and the reads corrected again for up to 12 rounds (`-x small`) or 16
   (`-x large`), stopping when fewer than 0.05% of tokens change; a guard keeps a partner only when the uncorrected
   tokens still agree. For genomes of 10 Mb and more (`-x large`), each partner's q-symbols are also aligned to the
   read's over the whole overlap (weighted edit distance, gap 2) and the partner is kept at a cost of at most 0.667 per
   column, checked only for the chains the coverage cap would keep; from round 2 partners of partners are tried and
   only changed reads are searched again, and from round 8 only the changed segments are corrected.
3. **Assembly** (`src/sparse.rs`, `src/bp_scale.rs`, `src/graph/`). One in 8 corrected
   tokens (by value hash) is kept, without the most frequent 100 ppm of values (metaMDBG's repetitive-minimizer cut).
   The multi-k directed minimizer-space de Bruijn graph is built for k = 4 to 32 from the reads and the previous order's
   contigs, as metaMDBG does, and cleaned at every order with its rules: tips up to 50 kb, bubbles, superbubbles up to
   50 kb and progressive abundance filters (cutoffs 2, 3). The cleaned graph at k = 32 is the result. Lengths in bp
   come from the signal: a token's length is the median over reads of samples x 400 bases/s / sample rate / tokens.
4. **Read threading and contig joins** (`src/threading.rs`, `src/joins.rs`). Repeats the k = 32 graph leaves
   collapsed are resolved by the reads that cross them (pairings of unique anchors with at least 3 reads, at most 20%
   conflicting, 3x the next best in tangles), as Flye resolves repeats. A dead end of a unitig of at least 200 nodes is
   then joined to the start of another by two rules: first, when two walks that share no read, each hopping from read to
   read through first-round partners whose uncorrected q-symbols align (at most 0.62 per column), land in the first 200
   tokens of the same unitig; then, on the ends left, when at least 3 reads touch the end of one and the start of the
   other and each is the other's best partner by 3x the next.
5. **Signal layout** (`src/layout.rs`). Each unitig of the final graph is tiled with the stretches of read signal that
   make it up (read id and POD5 sample span), preferring pieces whose corrected tokens equal the read's own.

Raw signal has no reverse complement, so each strand assembles on its own: the graph holds each region once per
strand, and each strand needs its own coverage (about 20x per strand or more).

The driver (`src/pipeline.rs`) runs these steps as the paper's runs did, each in a child process (`squash __STEP`,
`src/main.rs`, `src/exec.rs`) with the paper's arguments, stage by stage: `src/prepare.rs` (the uncorrected tokens and
the low-complexity list), `src/correct.rs` (the partner search and the correction rounds), `src/assemble.rs` (the bp
scale, the graph, threading and joins, the layout and the output files). `src/presets.rs` holds every tunable value
and its default under each preset; `src/cli.rs` the command line; `src/io.rs` the file helpers the steps share. The
drivers build each step's arguments; the steps themselves are the modules named above.

In the code: `partners/` is the partner search (`index.rs` the token index, `chain.rs` the chaining, `search.rs` and
`targets.rs` the search per read, `signal.rs` with `signal/` the q-symbol alignment check and its SIMD kernels,
`heavy.rs` the targets with the most index hits, searched on threads of their own); `rounds/` the correction rounds (`align.rs` the token alignment,
`replay.rs` the consensus, `skip.rs` the segments left as they are); `graph/` the multi-k graph (`tokens.rs` its input,
`build.rs` each order, `clean.rs` tips, bubbles and filters, `contigs.rs` the unitigs, `orders.rs` the order loop and
its files).
