#!/usr/bin/env python3
"""Closed replicons: for each reference replicon and strand, the circular unitigs of SQUASH's assembly (spelled FASTA
header circular=1) whose alignments (the strand split's strands/spelled.paf, minimap2 -cx asm10) lie on that replicon in
that orientation for most of their aligned bases, cover >= 95% of the replicon, and whose spelled length is 95-105% of
the replicon's. Orientation + = forward strand, - = reverse (as the strand split). One row per replicon that a unitig
aligns to: key, replicon, its length, forward / reverse closed (yes / no), the closing unitigs (name:length).
Usage: closed.py KEY EVAL_DIR > closed.tsv   (EVAL_DIR = the set's eval/KEY/squash/)
"""
import sys
from collections import defaultdict
from pathlib import Path

key, ev = sys.argv[1], Path(sys.argv[2])
print("key\treplicon\treplicon_bp\tforward\treverse\tunitigs")
circ, ulen = set(), {}
for line in open(ev / "asm/squash_spelled.fa"):
    if line.startswith(">"):
        f = line[1:].split()
        ulen[f[0]] = int(next(x for x in f if x.startswith("LN:i:"))[5:])
        if "circular=1" in f:
            circ.add(f[0])
aln = defaultdict(lambda: defaultdict(int))    # unitig -> (target, strand) -> aligned query bases
iv = defaultdict(list)                         # (unitig, target) -> target intervals
tlen = {}
for line in open(ev / "strands/spelled.paf"):
    x = line.split("\t")
    q, qs, qe, s, t, tl, ts, te = x[0], int(x[2]), int(x[3]), x[4], x[5], int(x[6]), int(x[7]), int(x[8])
    tlen[t] = tl
    aln[q][(t, s)] += qe - qs
    iv[(q, t)].append((ts, te))
closed = defaultdict(lambda: {"+": [], "-": []})
for u in sorted(circ):
    if not aln[u]:
        continue
    (t, s), _ = max(aln[u].items(), key=lambda kv: kv[1])
    cov, end = 0, -1
    for a, b in sorted(iv[(u, t)]):
        if b > end:
            cov += b - max(a, end)
            end = b
    if cov >= 0.95 * tlen[t] and 0.95 * tlen[t] <= ulen[u] <= 1.05 * tlen[t]:
        closed[t][s].append(f"{u}:{ulen[u]}")
for t in sorted(tlen):
    c = closed[t]
    print(f"{key}\t{t}\t{tlen[t]}\t{'yes' if c['+'] else 'no'}\t{'yes' if c['-'] else 'no'}\t{','.join(c['+'] + c['-'])}")
