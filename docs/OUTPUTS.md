# Output files

Everything SQUASH writes to `OUT/`. Lengths are in tokens: a token is one assembly minimizer of the signal (1 in
`--density` corrected tokens); `approx_bp` converts tokens to bp with the run's median signal-estimated bp per token
(signal samples x `--bases-per-second` / sample rate, over the reads' tokens).

## assembly_graph.gfa

The final graph at k = `--k-max` (32), after read threading and contig joins, in GFA 1.0.

- `S <unitig> * LN:i:<tokens> IA:i:<abundance>`: one segment per unitig. There is no base sequence (`*`); `LN` is the
  unitig's length in tokens, `IA` its abundance (the read support of its k-min-mers that the graph cleaning uses,
  rounded).
- `L <from> + <to> + 31M`: links. The graph is directed (raw signal has no reverse complement), so every link joins the
  end of one unitig to the start of another in the forward orientation; linked unitigs overlap by k - 1 = 31 tokens.
- `L <u> + <u> + 0M`: a circular unitig's self-link (written after the other lines, as the paper's evaluation did).

Each region of the genome is in the graph once per strand; the two copies run in opposite directions and are not
paired.

## assembly_info.tsv

One row per unitig of the final graph:

| Column | |
|---|---|
| `unitig` | the unitig's id (the GFA segment name) |
| `tokens` | length in tokens |
| `approx_bp` | `tokens` x the run's bp per token, rounded: an estimate of its length in bases |
| `abundance` | the GFA's `IA` |
| `circular` | 1 for a circular unitig |
| `original_unitigs` | the unitigs of the graph before threading it was made from, comma-separated in path order (a repeat unitig copied by threading appears in each path through it) |

## layout.tsv.gz

The signal layout: for each unitig, the pieces of read signal that tile it, in order along the unitig (as miniasm lists
the reads of a unitig). Consecutive pieces share at least one k-min-mer, so each junction is placed in both reads'
signals.

| Column | |
|---|---|
| `unitig`, `circular`, `nodes` | the unitig, whether it is circular, its nodes (k-min-mers) |
| `rank` | the piece's place along the unitig (0 first) |
| `read_id`, `read_ordinal` | the read (POD5 read id) and its index in the run |
| `node_from`, `node_to` | the unitig nodes the piece contributes |
| `token_from`, `token_to` | the same in the unitig's tokens |
| `read_token_from`, `read_token_to` | where they are in the read's tokens |
| `q_from`, `q_to` | where they are in the read's q-symbols (= events) |
| `sample_from`, `sample_to` | where they are in the read's POD5 signal (samples, end exclusive): the stretch to cut out to basecall or polish this part of the unitig |
| `piece_node_from`, `piece_node_to` | the whole piece's nodes (beyond its contribution) |
| `verified` | 1 when the piece's corrected tokens equal the read's own uncorrected ones (read from its own signal); 0 for a fallback piece used where no verified piece reaches |
| `own_tokens` | the piece's own (uncorrected-equal) tokens |

Where no piece covers some nodes of a unitig, the layout has a gap there.

## read_paths.tsv.gz

Every read's path through the final graph: each run of the read's k-min-mers that follows one unitig collinearly, one
row per piece (`unitig`, `read_id`, `node_from`, `node_to`, `offset` = the piece's diagonal: the read's k-min-mer
window w lies on unitig node w + offset; `verified` = 1 for pieces read from the read's own signal, 0 for the others).

## joins.tsv, stats.tsv

`joins.tsv`: one row per contig join made (the rule, the two unitig ends, the reads of the walks). `stats.tsv`: counts
of the threading and join steps (reads, windows, threaded repeats and tangles, joins).

## params.tsv, squash.log

`params.tsv`: every value the run used (the preset's or yours), the version and the inputs; `squash --dry-run` prints
the same. `squash.log`: the run's log; every step's own log is in `OUT/work/logs/` while the work folder is kept.

## work/ (with --keep-work)

The intermediate files:

| Folder | Files |
|---|---|
| `signal/` | `qsymbols.tsv` (every read's q-symbols), `spans.bin` (every event's POD5 samples) |
| `tokens/` | the uncorrected tokens: `reads.tsv.gz` (per read: id, ordinal, q-symbols, first token, tokens), `tokens.u32le` (token values), `token_spans.u64le` (each token's first and end q-symbol) |
| `correction/` | `low_complexity.tsv` (the reads left out of correction), `partners.tsv.gz` (the first partner search), `qsym_verdicts.tsv` (`-x large`: the q-symbol check's verdicts), `roundN/tokens.tsv.gz` (the final round's corrected tokens), `roundN/sparse/tokens.tsv.gz` (1 in 8 of them, the graph's input) |
| `bp_scale/` | `signal_lengths.tsv`, `token_counts.tsv.gz`, `bp_per_token.tsv` (the bp scale of the sparse tokens) |
| `graph/` | `input/` (the sparse tokens as token arrays, the same three files as `tokens/`), `k32.unitigs.tsv.gz`, `k32.unitigs.gfa`, `k32.node_tokens.u32le` (the k = 32 graph; the smaller orders are removed) |
| `threaded/` | the same three `k32.*` files after read threading and joins, `unitig_map.tsv`, `joins.tsv`, `stats.tsv` |
| `layout/` | `layout.tsv`, `pile.tsv.gz`, `pile_unverified.tsv.gz` |
| `logs/`, `stages/` | each step's command line and log; one marker per finished stage (`--resume`) |
