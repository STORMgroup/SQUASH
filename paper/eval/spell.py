#!/usr/bin/env python3
"""Spelled sequences: SQUASH's unitigs written in bases through the sup basecalls of the signal its layout assigns to
them, for QUAST. Evaluation only: SQUASH never sees bases; the basecalls are a read-out of the layout's signal, as
miniasm spells a unitig from its layout reads without consensus.

For each unitig, the rows of layout.tsv.gz in rank order:
- a piece gives the POD5 samples of its read from the start of its first unitig token to the start of the token where
  the next piece takes over (the first sample of that token's first event; both reads hold that token, so the cut is
  at the same place of the genome in both); the last piece of a linear unitig runs to the end of its last token, the
  last piece of a circular unitig to the start of the unitig's first token;
- those samples become the read's sup bases through dorado's move table (base k belongs to the stride block of its
  move: block = (sample - sp - ts) // stride), over every record of the read (dorado's split children carry pi:Z:parent
  and sp:i:offset);
- samples that no dorado record covers become Ns at the read's own bases per sample; a layout gap becomes Ns: its tokens
  times the median bases per token of the spelled pieces.

Inputs: SQUASH's layout.tsv.gz, its graph's input reads (work/graph/input/: each read's tokens' q-symbol spans), its
event spans (work/signal/spans.bin: each event's POD5 samples) and dorado's sup BAM with move tables (--emit-moves).
Writes FASTA (one record per unitig, `>unitig_ID LN:i:... circular=0|1 pieces=N gaps=N`, 80 bases a line) and
OUT.pieces.tsv (one row per piece).
Usage: spell.py --layout layout.tsv.gz --graph-input DIR --spans spans.bin --moves sup.moves.bam --output OUT.fa
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


class Moves:
    """One dorado record: sequence, move table, trimmed samples and its sample offset in the parent read."""

    def __init__(self, seq: str, stride: int, flags: np.ndarray, ts: int, sp: int):
        self.seq, self.stride, self.ts, self.sp = seq, stride, ts, sp
        # flags: one per block of `stride` samples, 1 = a new base starts in that block; cum[i] = bases started in
        # blocks 0..i
        self.cum = np.cumsum(flags.astype(np.int64))
        self.lo = sp + ts                                  # first sample with a move block
        self.hi = sp + ts + len(flags) * stride            # end of the last block

    def bases_before(self, sample: int) -> int:
        """Bases emitted in the blocks before the one holding `sample` (clamped to the record)."""
        blk = (sample - self.sp - self.ts) // self.stride
        if blk <= 0:
            return 0
        if blk >= len(self.cum):
            return int(self.cum[-1]) if len(self.cum) else 0
        return int(self.cum[blk - 1])

    def bases(self, s0: int, s1: int) -> str:
        """The bases whose blocks lie in samples [s0, s1)."""
        return self.seq[self.bases_before(s0):self.bases_before(s1)]


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def record_bases(r) -> str:
    """The record's bases in signal order: an aligned BAM stores a reverse-strand record's bases reverse complemented,
    while its move table stays in signal order (dorado's unaligned BAMs have no reverse records)."""
    seq = r.query_sequence or ""
    return seq.translate(COMPLEMENT)[::-1] if r.is_reverse else seq


def load_moves(path: Path, want: set[str]) -> dict[str, list[Moves]]:
    """parent read id -> its records sorted by sample offset (one record unless dorado split the read)."""
    import pysam
    out: dict[str, list[Moves]] = {}
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or not r.has_tag("mv"):   # one primary record per read (part)
                continue
            parent = r.get_tag("pi") if r.has_tag("pi") else r.query_name   # a split child counts as its parent read
            if parent not in want:
                continue
            mv = np.frombuffer(r.get_tag("mv"), dtype=np.int8)              # mv[0] = stride, mv[1:] = block flags
            out.setdefault(parent, []).append(Moves(record_bases(r), int(mv[0]), mv[1:], int(r.get_tag("ts")),
                                                    int(r.get_tag("sp")) if r.has_tag("sp") else 0))
    for v in out.values():
        v.sort(key=lambda m: m.sp)   # a split read's parts in signal order
    return out


def read_bases(records: list[Moves], s0: int, s1: int) -> str:
    """The read's bases of POD5 samples [s0, s1), through every record that overlaps them; samples no record covers
    become Ns, at the read's own bases per sample."""
    covered = sum(m.hi - m.lo for m in records)
    rate = sum(len(m.seq) for m in records) / covered if covered > 0 else 0.0
    out, pos = [], s0                     # pos = first sample not yet written
    for m in records:                     # records in signal order
        a, b = max(s0, m.lo), min(s1, m.hi)
        if b <= a:
            continue
        if a > pos:
            out.append("N" * round((a - pos) * rate))
        out.append(m.bases(a, b))
        pos = b
    if s1 > pos:
        out.append("N" * round((s1 - pos) * rate))
    return "".join(out)


def piece_samples(rows: list[dict], j: int, circular: bool, qs_of, st_of) -> tuple[int, int]:
    """Sample range [s0, s1) of layout row j of one unitig (rows in rank order, row j not a gap)."""
    r = rows[j]
    s0 = int(r["sample_from"])   # first sample of the piece's first token
    nxt = rows[j + 1] if j + 1 < len(rows) else (rows[0] if circular and len(rows) > 1 else None)
    if nxt is not None and nxt["read_id"] != "-":
        # end where the next piece takes over, in this read's signal: read token r1 is the take-over token (the first
        # token after this piece's contribution); qs[r1, 0] = its first q-symbol (event); st[q] = its first sample
        qs = qs_of(r["read_id"])
        r1 = int(r["read_token_to"])
        if r1 < len(qs):
            q = int(qs[r1, 0])
            st = st_of(r["read_id"])
            if 0 <= q < len(st):
                return s0, max(s0, int(st[q]))
    return s0, int(r["sample_to"])   # last piece of a linear unitig, or before a gap: to the end of its last token


def spell(layout_rows: list[dict], qs_of, st_of, moves: dict[str, list[Moves]]):
    """unitig -> (sequence, circular, pieces, gaps); the pieces table's rows; the median bases per token."""
    by_u: dict[int, list[dict]] = {}
    for r in layout_rows:
        by_u.setdefault(int(r["unitig"]), []).append(r)
    parts, table, per_token = {}, [], []
    for u, rows in by_u.items():                           # pass 1: every piece
        rows.sort(key=lambda x: int(x["rank"]))
        circ = rows[0]["circular"] == "1"
        seqs = []
        for j, r in enumerate(rows):
            if r["read_id"] == "-":                        # a layout gap: sized in pass 2
                seqs.append(None)
                continue
            s0, s1 = piece_samples(rows, j, circ, qs_of, st_of)
            b = read_bases(moves.get(r["read_id"], []), s0, s1)
            seqs.append(b)
            tok = int(r["token_to"]) - int(r["token_from"])
            if tok > 0 and b:
                per_token.append(len(b) / tok)
            table.append((u, int(r["rank"]), r["read_id"], s0, s1, len(b), b.count("N"), r.get("verified", "1")))
        parts[u] = (rows, circ, seqs)
    bpt = float(np.median(per_token)) if per_token else 0.0
    out = {}
    for u, (rows, circ, seqs) in parts.items():            # pass 2: the pieces joined in rank order, gaps as Ns
        s, gaps = [], 0
        for r, b in zip(rows, seqs):
            if b is None:
                gaps += 1
                s.append("N" * max(1, round((int(r["token_to"]) - int(r["token_from"])) * bpt)))
            else:
                s.append(b)
        out[u] = ("".join(s), circ, sum(1 for b in seqs if b is not None), gaps)
    return out, table, bpt


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--layout", type=Path, required=True, help="SQUASH's layout.tsv.gz")
    p.add_argument("--graph-input", type=Path, required=True, help="SQUASH's work/graph/input/")
    p.add_argument("--spans", type=Path, required=True, help="SQUASH's work/signal/spans.bin")
    p.add_argument("--moves", type=Path, required=True, help="dorado BAM with move tables (--emit-moves)")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    rows = common.layout_rows(a.layout)
    want = common.layout_reads(rows)
    qs = common.graph_reads(a.graph_input, want)
    spans = common.event_spans(a.spans, want)
    moves = load_moves(a.moves, want)
    missing = sorted(want - set(moves))
    out, table, bpt = spell(rows, qs.__getitem__, lambda r: spans[r][0], moves)
    tmp = a.output.with_name("." + a.output.name + ".tmp")
    kept = 0
    with tmp.open("w") as fh:
        for u in sorted(out):
            seq, circ, n, gaps = out[u]
            if not seq:
                continue
            kept += 1
            common.write_fasta_record(fh, f"unitig_{u} LN:i:{len(seq)} circular={int(circ)} pieces={n} gaps={gaps}", seq)
    with a.output.with_suffix(".pieces.tsv").open("w") as fh:
        fh.write("unitig\trank\tread_id\tsample_from\tsample_to\tbases\tn_bases\tverified\n")
        fh.writelines("\t".join(map(str, t)) + "\n" for t in table)
    tmp.rename(a.output)
    total = sum(len(v[0]) for v in out.values())
    print(f"{len(out)} unitigs ({kept} written), {total:,} bp, {len(table)} pieces, {sum(v[3] for v in out.values())} gaps, "
          f"{bpt:.2f} bases per token, reads without moves: {len(missing)}, "
          f"Ns where dorado gave no bases: {sum(t[6] for t in table):,}")


if __name__ == "__main__":
    main()
