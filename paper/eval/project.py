#!/usr/bin/env python3
"""Placed sequences: SQUASH's unitigs written as the reference sequence their signal lies on, for QUAST. Evaluation
only: none of the reads' bases are written.

A spelled unitig (spell.py) carries one read's basecalls per layout piece, so that read's basecalling errors count
against the assembly. Here each layout piece (a stretch of one read's signal) is placed instead:
- every event of the piece has a reference position (remora_evref.py's evref.bin: Remora's refinement of the read's sup
  alignment; statuses CONF and LOWQ placed, never UNPLACED);
- consecutive placed events that continue on one contig and strand (each within --merge bp of the previous, either
  direction) form a block; a jump or a change of contig or strand starts a new one. Blocks also continue across the
  pieces' hand-overs (both reads hold the hand-over token);
- the unitig is written as its blocks' reference sequence in order (reverse complement on the reverse strand), with
  Ns for what cannot be placed: a layout gap (its tokens x the median bp per token of the placed pieces), a fallback
  piece (layout column verified = 0: its read's own signal does not support the unitig there), runs of unplaced events
  between blocks that do not continue each other (their events x the piece's bp per event).

So QUAST scores where the unitig lies: a dropped tandem copy, a join of two loci or an inversion appear as they would in
a sequence assembly; a read's indels, stalls or base errors do not. Identity, mismatches and indels in QUAST's report
mean nothing here.

Usage: project.py --layout layout.tsv.gz --evref KEY.evref.bin --segs KEY.segs.tsv --reference evref's reference.fa
       --output OUT.fa [--merge 500] [--statuses 1,2,3] [--gap 10000000]
Writes OUT.fa (`>unitig_ID LN:i:... blocks=N placed=F`) and OUT.blocks.tsv (one row per block or N run).
"""
from __future__ import annotations

import argparse
import mmap
import statistics
import struct
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def offsets(names_lengths: list[tuple[str, int]], gap: int) -> list[tuple[int, str]]:
    """Global coordinates: contigs in header order, each after the previous one's length plus `gap`."""
    out, cur = [], 0
    for n, ln in names_lengths:
        out.append((cur, n))
        cur += ln + gap
    return out


def load_evref(path: Path, want: set[str]) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """read -> (global reference position (i64, -1 unplaced), segment (i32), status (u8)) per event."""
    out = {}
    with open(path, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        off = 0
        while off < len(mm):
            nl = struct.unpack_from("<I", mm, off)[0]
            name = mm[off + 4:off + 4 + nl].decode()
            off += 4 + nl
            n = struct.unpack_from("<I", mm, off)[0]
            off += 4
            if name in want:
                ref = np.frombuffer(mm, "<i8", n, off).copy()
                seg = np.frombuffer(mm, "<i4", n, off + 12 * n).copy()
                st = np.frombuffer(mm, "<u1", n, off + 16 * n).copy()
                out[name] = (ref, seg, st)
            off += 17 * n
    return out


def load_segs(path: Path) -> dict[int, str]:
    strand = {}
    with open(path) as fh:
        head = next(fh).rstrip("\n").split("\t")
        i, j = head.index("seg"), head.index("strand")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            strand[int(f[i])] = f[j]
    return strand


class Projector:
    """Blocks of one unitig, built piece by piece."""

    def __init__(self, merge: int):
        self.merge = merge
        self.items: list[list] = []   # ["B", contig, strand, first, last, lo, hi, events] or ["N", length, why]

    def _last_block(self):
        return self.items[-1] if self.items and self.items[-1][0] == "B" else None

    def add_event(self, contig: str, strand: str, pos: int, skipped_bp: float = 0.0) -> None:
        b = self._last_block()
        tol = self.merge + skipped_bp
        if b is not None and b[1] == contig and b[2] == strand and abs(pos - b[4]) <= tol:
            b[4] = pos
            b[5] = min(b[5], pos)
            b[6] = max(b[6], pos)
            b[7] += 1
            return
        self.items.append(["B", contig, strand, pos, pos, pos, pos, 1])

    def add_n(self, length: int, why: str) -> None:
        if length <= 0:
            return
        if self.items and self.items[-1][0] == "N" and self.items[-1][2] == why:
            self.items[-1][1] += length
        else:
            self.items.append(["N", length, why])

    def drop_bridged_ns(self) -> None:
        """A run of Ns between two blocks that continue each other (within --merge plus 1.5 x the run's length) is
        dropped and the blocks merged: what could not be placed lies where the flanks say."""
        out: list[list] = []
        pending: list[list] = []
        for it in self.items:
            if it[0] == "N":
                pending.append(it)
                continue
            a = out[-1] if out and out[-1][0] == "B" else None
            nlen = sum(n[1] for n in pending)
            if a is not None and a[1] == it[1] and a[2] == it[2] and abs(it[3] - a[4]) <= self.merge + 1.5 * nlen:
                a[4] = it[4]
                a[5] = min(a[5], it[5])
                a[6] = max(a[6], it[6])
                a[7] += it[7]
            else:
                out.extend(pending)
                out.append(it)
            pending = []
        out.extend(pending)
        self.items = out

    def sequence(self, ref: dict[str, str]) -> str:
        parts = []
        for it in self.items:
            if it[0] == "N":
                parts.append("N" * it[1])
            else:
                s = ref[it[1]][it[5]:it[6] + 1]
                parts.append(s if it[2] == "+" else revcomp(s))
        return "".join(parts)


def project(rows, ev, segs, places, merge: int, statuses: set[int], bp_token: float) -> Projector:
    """rows: one unitig's layout rows (dicts) in rank order -> Projector. places(global) -> (contig, pos)."""
    pr = Projector(merge)
    for r in rows:
        tokens = int(r["token_to"]) - int(r["token_from"])
        if r["read_id"] == "-":
            pr.add_n(int(round(tokens * bp_token)), "gap")
            continue
        q0, q1 = int(r["q_from"]), int(r["q_to"])
        if r.get("verified", "1") == "0":
            pr.add_n(int(round(tokens * bp_token)), "fallback")
            continue
        if r["read_id"] not in ev or q1 <= q0:
            pr.add_n(int(round(tokens * bp_token)), "no_placement")
            continue
        ref, seg, st = ev[r["read_id"]]
        q1 = min(q1, len(ref))
        ok = np.zeros(q1 - q0, dtype=bool)
        for s in statuses:
            ok |= st[q0:q1] == s
        ok &= ref[q0:q1] >= 0
        idx = np.flatnonzero(ok) + q0
        if len(idx) < 2:
            pr.add_n(int(round(tokens * bp_token)), "unplaced")
            continue
        # the piece's own bases per event, from its placed span
        span = abs(int(ref[idx[-1]]) - int(ref[idx[0]]))
        bpe = span / max(1, idx[-1] - idx[0]) if span > 0 else 1.35
        bpe = min(max(bpe, 0.5), 3.0)
        lead = idx[0] - q0
        if lead:
            pr.add_n(int(round(lead * bpe)), "unplaced")
        prev = None
        for q in idx.tolist():
            contig, pos = places(int(ref[q]))
            skipped = (q - prev - 1) * bpe if prev is not None else 0.0
            if skipped > 0:
                b = pr._last_block()
                if b is None or b[1] != contig or b[2] != segs.get(int(seg[q]), "+") or abs(pos - b[4]) > merge + skipped * 1.5:
                    pr.add_n(int(round(skipped)), "unplaced")
            pr.add_event(contig, segs.get(int(seg[q]), "+"), pos, skipped * 1.5)
            prev = q
        tail = q1 - 1 - idx[-1]
        if tail:
            pr.add_n(int(round(tail * bpe)), "unplaced")
    pr.drop_bridged_ns()
    return pr


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layout", type=Path, required=True, help="SQUASH's layout.tsv.gz")
    ap.add_argument("--evref", type=Path, required=True)
    ap.add_argument("--segs", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True, help="the evref's reference.fa (contig order = its axis)")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--merge", type=int, default=500,
                    help="largest jump (bp) that still continues a block (below QUAST's 1-kb misassembly threshold)")
    ap.add_argument("--statuses", default="1,2,3", help="evref statuses placed (1 CONF, 2 SINGLE, 3 LOWQ)")
    ap.add_argument("--gap", type=int, default=10_000_000, help="remora_evref.py --gap of the evref")
    a = ap.parse_args()
    if a.output.exists():
        raise SystemExit(f"refusing to overwrite {a.output}")
    ref = read_fasta(a.reference)
    offs = offsets([(n, len(s)) for n, s in ref.items()], a.gap)
    starts = np.array([o for o, _ in offs], dtype=np.int64)
    names = [n for _, n in offs]

    def places(g: int) -> tuple[str, int]:
        i = int(np.searchsorted(starts, g, side="right")) - 1
        return names[i], g - int(starts[i])

    rows = defaultdict(list)
    for r in common.layout_rows(a.layout):
        rows[int(r["unitig"])].append(r)
    want = {r["read_id"] for rs in rows.values() for r in rs if r["read_id"] != "-"}
    ev = load_evref(a.evref, want)
    segs = load_segs(a.segs)
    statuses = {int(x) for x in a.statuses.split(",")}
    # bp per token: placed pieces' reference span over their tokens (median), for gaps and fallback pieces
    rates = []
    for rs in rows.values():
        for r in rs:
            if r["read_id"] in ev and r.get("verified", "1") != "0":
                refv, _, st = ev[r["read_id"]]
                q0, q1 = int(r["q_from"]), min(int(r["q_to"]), len(refv))
                sel = refv[q0:q1][(st[q0:q1] == 1) & (refv[q0:q1] >= 0)]
                t = int(r["token_to"]) - int(r["token_from"])
                if len(sel) > 10 and t > 0:
                    rates.append((sel.max() - sel.min()) / t)
    bp_token = statistics.median(rates) if rates else 100.0
    a.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = a.output.with_name("." + a.output.name + ".tmp")
    blocks = a.output.with_suffix(".blocks.tsv")
    btmp = blocks.with_name("." + blocks.name + ".tmp")
    n_why = defaultdict(int)
    total = placed = 0
    with open(tmp, "w") as fa, open(btmp, "w") as bt:
        bt.write("unitig\titem\tkind\tcontig\tstrand\tstart\tend\tlength\tevents_or_reason\n")
        for u in sorted(rows):
            pr = project(rows[u], ev, segs, places, a.merge, statuses, bp_token)
            seq = pr.sequence(ref)
            nb = sum(1 for it in pr.items if it[0] == "B")
            pl = sum(it[6] - it[5] + 1 for it in pr.items if it[0] == "B")
            total += len(seq)
            placed += pl
            common.write_fasta_record(fa, f"unitig_{u} LN:i:{len(seq)} blocks={nb} placed={pl / max(1, len(seq)):.4f}", seq)
            for k, it in enumerate(pr.items):
                if it[0] == "B":
                    bt.write(f"{u}\t{k}\tblock\t{it[1]}\t{it[2]}\t{it[5]}\t{it[6] + 1}\t{it[6] - it[5] + 1}\t{it[7]}\n")
                else:
                    bt.write(f"{u}\t{k}\tN\t-\t-\t-\t-\t{it[1]}\t{it[2]}\n")
                    n_why[it[2]] += it[1]
    btmp.rename(blocks)
    tmp.rename(a.output)
    print(f"{len(rows)} unitigs, {total:,} bp written, {placed:,} placed ({100 * placed / max(1, total):.2f}%); Ns by reason: "
          + ", ".join(f"{k} {v:,}" for k, v in sorted(n_why.items())) + f"; {bp_token:.1f} bp per token")


if __name__ == "__main__":
    main()
