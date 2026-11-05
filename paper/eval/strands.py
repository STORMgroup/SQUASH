#!/usr/bin/env python3
"""Strand of each unitig: SQUASH assembles each strand of the genome apart (a both-strand set gives every region twice,
once per strand), so the per-strand QUAST scores the unitigs of one strand against the reference. A unitig's strand is
the one with more aligned bases in minimap2's alignments of its spelled sequence to the reference (PAF column 11; equal:
forward); a unitig without an alignment is in neither file (and is counted in summary.tsv).
Writes OUT/squash_spelled_fwd.fa, _rev.fa and squash_placed_fwd.fa, _rev.fa (the records copied as they are) and
OUT/summary.tsv (unitigs and spelled bases of each class: all, forward, reverse, neither).
Usage: strands.py SPELLED.paf ASM_DIR OUT_DIR   (ASM_DIR holds squash_spelled.fa and squash_placed.fa)
"""
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

paf, asm, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
strand = defaultdict(lambda: [0, 0])
for line in open(paf):
    f = line.split("\t")
    strand[f[0]][f[4] == "-"] += int(f[10])
for kind in ("spelled", "placed"):
    fh = {s: open(out / f".squash_{kind}_{s}.fa.tmp", "w") for s in ("fwd", "rev")}
    cur = None
    for line in open(asm / f"squash_{kind}.fa"):
        if line.startswith(">"):
            u = line[1:].split()[0]
            cur = fh["rev" if strand[u][1] > strand[u][0] else "fwd"] if u in strand else None
        if cur:
            cur.write(line)
    for s, x in fh.items():
        x.close()
        os.replace(out / f".squash_{kind}_{s}.fa.tmp", out / f"squash_{kind}_{s}.fa")
allu = common.fasta_lengths(asm / "squash_spelled.fa")
fw = common.fasta_lengths(out / "squash_spelled_fwd.fa")
rv = common.fasta_lengths(out / "squash_spelled_rev.fa")
ne = {u: x for u, x in allu.items() if u not in fw and u not in rv}
with open(out / ".summary.tsv.tmp", "w") as fh:
    fh.write("class\tunitigs\tspelled_bp\n")
    for name, d in (("all", allu), ("forward", fw), ("reverse", rv), ("neither", ne)):
        fh.write(f"{name}\t{len(d)}\t{sum(d.values())}\n")
os.replace(out / ".summary.tsv.tmp", out / "summary.tsv")
