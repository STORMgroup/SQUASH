# Frequently asked questions

**Which preset?** By genome size: `-x small` under 10 Mb (bacteria, archaea), `-x large` from 10 Mb (eukaryotes). The
paper tested `small` on bacteria of 1.8-5.2 Mb and `large` on genomes of 41-142 Mb. `large` adds a check that keeps
false correction partners out of larger genomes, and speed options that would cost quality on bacteria.

**How much coverage?** At least 20x per strand, i.e. 40x of the genome; the paper used 30x and 50x per strand. Raw
signal has no reverse complement, so the reads of each strand assemble on their own. Below 20x per strand the graph
fragments. `-g SIZE` (e.g. `-g 4.6m`) makes SQUASH estimate the coverage and warn when it is too low.

**Why does each region appear twice?** Once per strand (see above). The two copies run in opposite directions and are
not paired, so the graph holds about twice the genome length.

**Why no FASTA?** SQUASH never basecalls. Unitigs are sequences of signal tokens. `layout.tsv.gz` says which stretch of
which read's signal makes up each part of a unitig (read id and POD5 sample span), so a unitig can be turned into bases
by basecalling those stretches (for example with dorado on the listed reads) and joining them in layout order. The
paper does this only to evaluate the assembly.

**Which data?** ONT R10.4.1 POD5, 4 kHz or 5 kHz, as MinKNOW or a public deposit writes them; other chemistries are
untested. FAST5: convert it to POD5 first (`pip install pod5==0.3.39`, then `pod5 convert fast5 --output reads.pod5
fast5_folder/*.fast5`). If SQUASH cannot read a POD5 that was merged from several runs, give it the runs' own files.

**Are runs reproducible?** Yes: the same POD5 files in the same order with the same options give the same outputs, at
any thread count. Reordering the reads of a large set can change the result slightly.

**Memory, time, disk?** See the README's performance section. While it runs, SQUASH keeps its intermediate files in
`OUT/work` (about as large as the POD5s; deleted at the end unless `--keep-work`): for a large genome put the output
folder, or `--tmp-dir`, on a fast local disk. With `pigz` installed these files are compressed and smaller.

**A run stopped.** Every step's log is in `OUT/work/logs/` and the run's in `OUT/squash.log`. `--resume` continues a
stopped run from its first unfinished stage (stages already run with the same values are kept).
