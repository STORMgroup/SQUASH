#!/usr/bin/env python3
"""Tier 1: contiguity from the raw signal alone (no basecalls for the raw-signal assemblers), as the paper reports it.

Length unit of a raw-signal assembly: bases = POD5 samples x bases per second (400, the R10.4.1 kit's speed) / sample
rate. Each segment's samples:
  SQUASH               its signal layout (layout.tsv.gz): each piece from its own first sample to where the next piece
                       takes over in this read's signal (spell.py's piece_samples: the spelled pieces' sample ranges), so
                       consecutive pieces join without a gap; a layout gap (no read) = its tokens x the median samples per
                       token of the pieces.
  Rawsamble + miniasm  each `a` line of miniasm's GFA: the read's Rawsamble units (its events without event merging) from
                       the first sample of its first event to the first sample of the event after its last (the read's
                       last sample at its end); a read without spans = its units x the median samples per unit.
  basecalled assemblers their contigs' lengths in bases (nanoMDBG, hifiasm: the contig list, a circular contig as a
                       self-link; Flye: its contigs, with its assembly_graph.gfa for the graph statistics).
Metrics: segments, total, largest, N50, NG50, auN, auNG; graph, graph segments, links (a link and its reverse complement
once), circular segments (a self-link or miniasm's `c` name), connected components and the largest one (its share of the
graph's total). NG50 / auNG against G: the reference length for a basecalled assembly, twice it for a raw-signal assembly
of a both-strand set (each strand is assembled apart). auN = sum L^2 / sum L; auNG = the area under the NGx curve.

Usage:
  tier1.py squash --key KEY --gfa assembly_graph.gfa --layout layout.tsv.gz --graph-input DIR --spans spans.bin
                  --signal-lengths signal_lengths.tsv --genome REF.fa --output tier1.tsv
  tier1.py compare --key KEY --genome REF.fa --signal-lengths signal_lengths.tsv --nanomdbg contigs.fasta.gz
                  --flye ASM_DIR --hifiasm asm.bp.p_ctg.fa [--rawsamble ava.gfa --spans-raw spans.bin] --output tier1.tsv
"""
from __future__ import annotations

import argparse
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from spell import piece_samples  # noqa: E402

COLUMNS = ["assembler", "unit", "genome_size", "segments", "total_bp", "largest_bp", "N50", "NG50", "auN", "auNG",
           "graph", "graph_segments", "links", "circular", "components", "graph_total_bp", "largest_component_bp",
           "largest_component_share"]
FLIP = {"+": "-", "-": "+"}
BASES_PER_SECOND = 400.0


# ---------------------------------------------------------------- length and graph statistics

def nx(lengths: list[float], target: float) -> float:
    """The length of the segment at which the decreasing cumulative length first reaches `target` (0 if never)."""
    cum = 0.0
    for x in sorted(lengths, reverse=True):
        cum += x
        if cum >= target:
            return x
    return 0.0


def aun(lengths: list[float]) -> float:
    total = sum(lengths)
    return sum(x * x for x in lengths) / total if total > 0 else 0.0


def aung(lengths: list[float], genome: float) -> float:
    """Area under the NGx curve: each segment's length times the part of [0, G] its cumulative interval covers, / G."""
    area, cum = 0.0, 0.0
    for x in sorted(lengths, reverse=True):
        covered = min(cum + x, genome) - min(cum, genome)
        area += x * covered
        cum += x
        if cum >= genome:
            break
    return area / genome if genome > 0 else 0.0


def read_gfa_graph(path: Path) -> tuple[list[str], set, set]:
    """(segment names in file order, canonical links, circular segments) of a GFA (S and L lines)."""
    names, links, circular = [], set(), set()
    with open(path) as fh:
        lines = fh.read().splitlines()
    for line in lines:
        f = line.split("\t")
        if f[0] == "S":
            names.append(f[1])
        elif f[0] == "L":
            a, oa, b, ob = f[1], f[2], f[3], f[4]
            rev = (b, FLIP[ob], a, FLIP[oa])
            links.add(min((a, oa, b, ob), rev))
            if a == b:
                circular.add(a)
    return names, links, circular


def components(names: list[str], links) -> list[list[str]]:
    """Connected components (links undirected)."""
    parent = {n: n for n in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, _, b, _ in links:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    groups = defaultdict(list)
    for n in names:
        groups[find(n)].append(n)
    return list(groups.values())


def metrics(assembler: str, unit: str, lengths: dict, genome: float, links=None, circular=None,
            graph: str = "the segments", graph_lengths: dict | None = None) -> dict:
    """Contig statistics of `lengths`; graph statistics of `graph_lengths` (default the same segments) with `links`."""
    vals = list(lengths.values())
    total = sum(vals)
    row = {"assembler": assembler, "unit": unit, "genome_size": int(genome), "segments": len(vals),
           "total_bp": round(total), "largest_bp": round(max(vals)) if vals else 0, "N50": round(nx(vals, total / 2)),
           "NG50": round(nx(vals, genome / 2)), "auN": round(aun(vals)), "auNG": round(aung(vals, genome)),
           "graph": "-", "graph_segments": "-", "links": "-", "circular": "-", "components": "-", "graph_total_bp": "-",
           "largest_component_bp": "-", "largest_component_share": "-"}
    if links is not None:
        glen = lengths if graph_lengths is None else graph_lengths
        gtotal = sum(glen.values())
        comps = components(list(glen), links)
        big = max((sum(glen[n] for n in c) for c in comps), default=0.0)
        row.update(graph=graph, graph_segments=len(glen), links=len(links), circular=len(circular or ()),
                   components=len(comps), graph_total_bp=round(gtotal), largest_component_bp=round(big),
                   largest_component_share=f"{big / gtotal:.4f}" if gtotal else "0")
    return row


def gfa_lengths(path: Path) -> dict[str, float]:
    """Segment lengths of a GFA: the LN:i: tag, else the sequence's length."""
    out = {}
    for line in open(path):
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            ln = next((int(x[5:]) for x in f[3:] if x.startswith("LN:i:")), None)
            out[f[1]] = float(ln if ln is not None else (len(f[2]) if f[2] != "*" else 0))
    return out


def contig_list(fa: Path) -> tuple[dict, set]:
    """Contig lengths and the circular ones (header `circular=yes` or `circular=1`) of a (gzipped) FASTA."""
    lens = common.fasta_lengths(fa)
    with common.open_text(fa) as fh:
        circ = {l[1:].split()[0] for l in fh if l.startswith(">") and ("circular=yes" in l or "circular=1" in l)}
    return lens, circ


# ---------------------------------------------------------------- raw-signal segments in samples

def squash_samples(layout: Path, graph_input: Path, spans: Path) -> tuple[dict[str, float], dict]:
    """unitig -> samples, from the signal layout as spell.py cuts its pieces; and counts."""
    rows = common.layout_rows(layout)
    want = common.layout_reads(rows)
    qs = common.graph_reads(graph_input, want)
    sp = common.event_spans(spans, want)
    by_u: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_u[r["unitig"]].append(r)
    pieces, per_token = {}, []
    for u, rs in by_u.items():
        rs.sort(key=lambda x: int(x["rank"]))
        circ = rs[0]["circular"] == "1"
        out = []
        for j, r in enumerate(rs):
            tokens = int(r["token_to"]) - int(r["token_from"])
            if r["read_id"] == "-" or r["read_id"] not in sp or r["read_id"] not in qs or int(r["sample_from"]) < 0:
                out.append((None, tokens))                 # a gap: sized below
                continue
            s0, s1 = piece_samples(rs, j, circ, qs.__getitem__, lambda rid: sp[rid][0])
            out.append((s1 - s0, tokens))
            if tokens > 0 and s1 > s0:
                per_token.append((s1 - s0) / tokens)
        pieces[u] = out
    spt = statistics.median(per_token) if per_token else 0.0
    samples = {u: sum(s if s is not None else t * spt for s, t in out) for u, out in pieces.items()}
    gaps = sum(1 for out in pieces.values() for s, _ in out if s is None)
    return samples, {"pieces": sum(len(v) for v in pieces.values()), "gaps": gaps, "samples_per_token": spt}


def read_layouts(gfa: Path):
    """miniasm's GFA: segment lengths (LN) and layouts {segment: [(offset, read, start, end, strand, length)]}."""
    lengths, layouts = {}, {}
    with open(gfa) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                lengths[f[1]] = int(next(x for x in f[3:] if x.startswith("LN:i:"))[5:])
            elif f[0] == "a":
                read, span = f[3].rsplit(":", 1)
                start, end = (int(x) for x in span.split("-"))
                layouts.setdefault(f[1], []).append((int(f[2]), read, start, end, f[4], int(f[5])))
    for seg in layouts:
        layouts[seg].sort()
    return lengths, layouts


def piece_units(start: int, end: int, strand: str, length: int) -> tuple[int, int]:
    """The read interval [lo, hi) in signal units an `a` line uses."""
    if strand == "+":
        return start, min(end, start + max(0, length))
    return max(start, end - max(0, length)), end


def rawsamble_samples(gfa: Path, spans_raw: Path) -> tuple[dict[str, float], dict]:
    """segment -> samples of its layout's read pieces (unmerged events -> samples); and counts."""
    lengths, layouts = read_layouts(gfa)
    want = {r for lay in layouts.values() for _, r, *_ in lay}
    sp = {k: (a.astype(np.int64), b.astype(np.int64)) for k, (a, b) in common.event_spans(spans_raw, want).items()}
    samples, missing, per_unit, pending = {}, 0, [], []
    for seg in lengths:
        total, unsized = 0.0, 0
        for _, read, start, end, strand, length in layouts.get(seg, []):
            lo, hi = piece_units(start, end, strand, length)
            if read not in sp or hi <= lo:
                missing += read not in sp
                unsized += max(0, hi - lo)
                continue
            st, en = sp[read]
            n = len(st)
            lo, hi = max(0, min(lo, n)), max(0, min(hi, n))
            if hi <= lo:
                continue
            s0 = int(st[lo])
            s1 = int(st[hi]) if hi < n else int(en[n - 1])
            total += s1 - s0
            per_unit.append((s1 - s0) / (hi - lo))
        if seg not in layouts:
            unsized += lengths[seg]
        samples[seg] = total
        pending.append((seg, unsized))
    spu = statistics.median(per_unit) if per_unit else 0.0
    for seg, unsized in pending:                        # units without spans: the median samples per unit
        samples[seg] += unsized * spu
    return samples, {"reads_without_spans": missing, "samples_per_unit": spu}


# ---------------------------------------------------------------- the two tables

def sample_rate(signal_lengths: Path) -> int:
    """The sample rate of the set's reads (SQUASH's work/bp_scale/signal_lengths.tsv); one rate expected."""
    rates = set()
    with open(signal_lengths) as fh:
        next(fh)
        for line in fh:
            rates.add(line.rstrip("\n").split("\t")[2])
    if len(rates) != 1:
        raise SystemExit(f"{signal_lengths}: sample rates {sorted(rates)}")
    return int(rates.pop())


def genome_length(ref: Path) -> int:
    return sum(common.fasta_lengths(ref).values())


def write(path: Path, cols: list[str], rows: list[dict]) -> None:
    tmp = path.with_name("." + path.name + ".tmp")
    tmp.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r.get(k, "")) for k in cols) + "\n" for r in rows))
    tmp.rename(path)


def squash_table(a) -> None:
    genome = genome_length(a.genome)
    v, f = BASES_PER_SECOND, sample_rate(a.signal_lengths)
    g_raw = 2 * genome                                  # both strands
    names, links, circ = read_gfa_graph(a.gfa)
    samples, info = squash_samples(a.layout, a.graph_input, a.spans)
    if set(samples) != set(names):
        raise SystemExit(f"layout unitigs and GFA segments differ ({len(samples)} vs {len(names)})")
    row = metrics("SQUASH", "signal-estimated bp", {n: samples[n] * v / f for n in names}, g_raw, links, circ)
    row.update(key=a.key, reference_length=genome, gaps=info["gaps"], pieces=info["pieces"], sample_rate=f)
    write(a.output, ["key"] + COLUMNS + ["reference_length", "gaps", "pieces", "sample_rate"], [row])
    print(f"{a.key}: SQUASH {row['segments']} unitigs, NG50 {row['NG50']}, auNG {row['auNG']}, largest component "
          f"{row['largest_component_bp']} (signal-estimated bp, G = {g_raw})")


def compare_table(a) -> None:
    genome = genome_length(a.genome)
    v, f = BASES_PER_SECOND, sample_rate(a.signal_lengths)
    scale, g_raw = v / f, 2 * genome
    rows = []
    if a.rawsamble:
        rnames, rlinks, rcirc = read_gfa_graph(a.rawsamble)
        rcirc |= {n for n in rnames if n.endswith("c")}   # miniasm's circular unitigs (utgNNNNNNc)
        rsamples, rinfo = rawsamble_samples(a.rawsamble, a.spans_raw)
        rows.append(metrics("Rawsamble + miniasm", "signal-estimated bp", {n: rsamples[n] * scale for n in rnames},
                            g_raw, rlinks, rcirc))
        print(f"{a.key}: Rawsamble reads without spans {rinfo['reads_without_spans']}, samples per unit "
              f"{rinfo['samples_per_unit']:.2f}")
    lens, ncirc = contig_list(a.nanomdbg)
    rows.append(metrics("nanoMDBG", "basecalled bp", lens, genome, {(x, "+", x, "+") for x in ncirc}, ncirc,
                        "contig list (no links)"))
    gfa = a.flye / "assembly_graph.gfa"
    flens = common.fasta_lengths(a.flye / "assembly.fasta")
    if gfa.is_file():
        _, glinks, gcirc = read_gfa_graph(gfa)
        rows.append(metrics("Flye", "basecalled bp", flens, genome, glinks, gcirc, "assembly_graph.gfa", gfa_lengths(gfa)))
    else:
        rows.append(metrics("Flye", "basecalled bp", flens, genome, set(), set(), "contigs (graph not kept)"))
    lens, hcirc = contig_list(a.hifiasm)
    rows.append(metrics("hifiasm", "basecalled bp", lens, genome, {(x, "+", x, "+") for x in hcirc}, hcirc,
                        "primary contig list (no links)"))
    for r in rows:
        r.update(key=a.key, reference_length=genome, bases_per_second=v, sample_rate=f)
    write(a.output, ["key"] + COLUMNS + ["reference_length", "bases_per_second", "sample_rate"], rows)
    print(f"{a.key}: " + "; ".join(f"{r['assembler']} NG50 {r['NG50']}" for r in rows))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="what", required=True)
    s = sub.add_parser("squash")
    s.add_argument("--gfa", type=Path, required=True)
    s.add_argument("--layout", type=Path, required=True)
    s.add_argument("--graph-input", type=Path, required=True)
    s.add_argument("--spans", type=Path, required=True)
    c = sub.add_parser("compare")
    c.add_argument("--nanomdbg", type=Path, required=True)
    c.add_argument("--flye", type=Path, required=True, help="Flye's output folder (assembly.fasta, assembly_graph.gfa)")
    c.add_argument("--hifiasm", type=Path, required=True)
    c.add_argument("--rawsamble", type=Path, help="miniasm's GFA of Rawsamble's overlaps (with `a` lines)")
    c.add_argument("--spans-raw", type=Path, help="with --rawsamble: the reads' event spans without event merging")
    for p in (s, c):
        p.add_argument("--key", required=True)
        p.add_argument("--genome", type=Path, required=True, help="the reference whose length is G")
        p.add_argument("--signal-lengths", type=Path, required=True, help="SQUASH's work/bp_scale/signal_lengths.tsv")
        p.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    if a.what == "compare" and bool(a.rawsamble) != bool(a.spans_raw):
        raise SystemExit("--rawsamble needs --spans-raw")
    (squash_table if a.what == "squash" else compare_table)(a)


if __name__ == "__main__":
    main()
