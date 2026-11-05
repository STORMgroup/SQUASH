#!/usr/bin/env python3
"""Reference positions of the layout reads' events, by Remora's reference-anchored signal refinement (Remora 3.3.0;
evaluation only). The input of project.py.

For each read (its sup basecalls aligned by minimap2 with the move table carried, `-y --MD`): Remora's refinement
(refine_signal_map.SigMapRefiner with the R10.4.1 9-mer level table) gives every reference base of the alignment its
first signal sample at single-sample resolution (Read.ref_to_signal after set_refine_signal_mapping(ref_mapping=True)).
Each event's first POD5 sample (the spans file) is placed at the reference base whose segment holds it: i with
ref_to_signal[i] <= sample - ts < ref_to_signal[i + 1]; position = alignment start + i, or alignment end - 1 - i on the
reverse strand. Events outside the alignment (clips, adapters) and reads Remora cannot refine are unplaced (status 0).
Every alignment record of the read is used (primary first, then supplementary by aligned length; minimap2 -Y keeps their
full query): an event takes the first record that places it. Status 1 (CONF) when that record has MAPQ >= --min-mapq,
else 3 (LOWQ). Dorado-split reads (pi tag) are not placed. Global coordinates: contigs in the reference's header order,
each after the previous one's length plus --gap.

--collapse-check: Remora's refined map sometimes puts most of a read's events on one reference base. Each record's
Remora placement is compared with the move table's for the same record: when its longest run of consecutive events on
one base is >= --collapse-min (500) and more than --collapse-factor (4) x the move table's (a pore stall puts many events
on one base in both), Remora is run again up to --collapse-retries (2) times, and if it still collapses the record takes
the move table's placement (segs.tsv aligner `remora>moves`).

Usage: remora_evref.py --pod5 P.pod5 --bam ALN.bam --spans spans.bin --reference REF.fa --levels LEVELS.txt
                       --mode remora|moves [--collapse-check] --out PREFIX [--threads N] [--samtools SAMTOOLS]
Writes PREFIX.evref.bin (per read: u32 name_len, name, u32 n, n x i64 ref, n x i32 base (-1), n x i32 seg, n x u8
status), PREFIX.segs.tsv, PREFIX.summary.txt."""
from __future__ import annotations

import argparse
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np

UNPLACED, CONF, LOWQ = 0, 1, 3


def reference_offsets(fai_or_fa: Path, gap: int) -> dict[str, int]:
    """Contig -> global offset (header order, previous lengths plus `gap`), as project.py."""
    lens, name = [], None
    fai = Path(str(fai_or_fa) + ".fai")
    if fai.is_file():
        lens = [(l.split("\t")[0], int(l.split("\t")[1])) for l in fai.read_text().splitlines()]
    else:
        n = 0
        for line in open(fai_or_fa):
            if line.startswith(">"):
                if name is not None:
                    lens.append((name, n))
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
        if name is not None:
            lens.append((name, n))
    out, cur = {}, 0
    for n, ln in lens:
        out[n] = cur
        cur += ln + gap
    return out


def place_events(starts: np.ndarray, ts: int, r2s: np.ndarray, aln_start: int, aln_end: int, reverse: bool) -> np.ndarray:
    """Reference position (alignment coordinates) of each event's first sample, -1 outside the alignment. r2s: the
    alignment's reference-to-signal map (len = reference bases + 1, along the read, samples after the trimmed ts)."""
    s = starts.astype(np.int64) - ts
    i = np.searchsorted(r2s, s, side="right") - 1
    ok = (i >= 0) & (i < len(r2s) - 1)
    pos = np.where(reverse, aln_end - 1 - i, aln_start + i)
    return np.where(ok, pos, -1).astype(np.int64)


def longest_run(pos: np.ndarray) -> int:
    """The longest run of consecutive placed events (pos >= 0, in event order) on one reference position."""
    x = pos[pos >= 0]
    if len(x) == 0:
        return 0
    idx = np.flatnonzero(np.concatenate(([True], x[1:] != x[:-1], [True])))
    return int(np.diff(idx).max())


def collapsed(pos: np.ndarray, control: np.ndarray, min_run: int, factor: float) -> bool:
    """Remora's placement `pos` collapsed: a run of >= min_run events on one base, > factor x the move table's."""
    run = longest_run(pos)
    return run >= min_run and run > factor * max(1, longest_run(control))


def write_record(fh, name: str, ref: np.ndarray, seg: np.ndarray, status: np.ndarray) -> None:
    n = len(ref)
    fh.write(struct.pack("<I", len(name)))
    fh.write(name.encode())
    fh.write(struct.pack("<I", n))
    fh.write(ref.astype("<i8").tobytes())
    fh.write(np.full(n, -1, dtype="<i4").tobytes())
    fh.write(seg.astype("<i4").tobytes())
    fh.write(status.astype("u1").tobytes())


def split_spans(path: Path, want: set[str] | None, chunks: int, outdir: Path) -> dict[str, int]:
    """One pass over the spans file: each wanted read's record appended to outdir/c{i}.spans (i = its index mod `chunks`
    in file order); returns read -> number of events."""
    counts = {}
    fhs = [open(outdir / f"c{i}.spans", "wb") for i in range(chunks)]
    k = 0
    with open(path, "rb", buffering=1 << 24) as src:
        while True:
            h = src.read(4)
            if not h:
                break
            nl = struct.unpack("<I", h)[0]
            name = src.read(nl)
            hd = src.read(8)
            n = struct.unpack("<II", hd)[1]
            body = src.read(8 * n)
            if want is None or name.decode() in want:
                fhs[k % chunks].write(h + name + hd + body)
                counts[name.decode()] = n
                k += 1
    for fh in fhs:
        fh.close()
    return counts


def read_starts(path: Path) -> dict[str, np.ndarray]:
    buf = path.read_bytes()
    off, out = 0, {}
    while off < len(buf):
        nl = struct.unpack_from("<I", buf, off)[0]
        name = buf[off + 4:off + 4 + nl].decode()
        off += 4 + nl
        _, n = struct.unpack_from("<II", buf, off)
        off += 8
        out[name] = np.frombuffer(buf, "<u4", n, off).astype(np.int64)
        off += 8 * n
    return out


def work(args):
    """One chunk: per read, its alignment records (primary first, then supplementary by aligned length) as
    (contig, strand, mapq, positions (i64, -1 outside the record) or None, error[, note])."""
    bam_path, spans_path, pod5_path, levels, mode = args[:5]
    min_run, factor, retries = args[5] if len(args) > 5 else (0, 0.0, 0)   # --collapse-check (min_run 0 = off)
    import pod5
    import pysam
    from remora import io, refine_signal_map
    starts = read_starts(Path(spans_path))
    refiner = refine_signal_map.SigMapRefiner(kmer_model_filename=levels, do_rough_rescale=True, scale_iters=0,
                                              do_fix_guage=True) if mode == "remora" else None
    recs: dict[str, list] = {}
    with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.has_tag("pi") or r.query_name not in starts:
                continue
            recs.setdefault(r.query_name, []).append(r)
    for v in recs.values():
        v.sort(key=lambda r: (r.is_supplementary, -r.query_alignment_length))
    out = []
    with pod5.Reader(pod5_path) as reader:
        for p in reader.reads(selection=list(recs), missing_ok=True):
            rid = str(p.read_id)
            rows = []
            for a in recs[rid]:
                strand = "-" if a.is_reverse else "+"
                try:
                    read = io.Read.from_pod5_and_alignment(p, a)
                    ts = int(a.get_tag("ts"))
                    control = None
                    if refiner is not None and min_run > 0:
                        control = place_events(starts[rid], ts, np.asarray(read.ref_to_signal, dtype=np.int64).copy(),
                                               a.reference_start, a.reference_end, a.is_reverse)
                    if refiner is not None:
                        read.set_refine_signal_mapping(refiner, ref_mapping=True)
                    r2s = np.asarray(read.ref_to_signal, dtype=np.int64)
                    pos = place_events(starts[rid], ts, r2s, a.reference_start, a.reference_end, a.is_reverse)
                    note = ""
                    if control is not None and collapsed(pos, control, min_run, factor):
                        note = "moves"
                        for _ in range(retries):
                            again = io.Read.from_pod5_and_alignment(p, a)
                            again.set_refine_signal_mapping(refiner, ref_mapping=True)
                            pos = place_events(starts[rid], ts, np.asarray(again.ref_to_signal, dtype=np.int64),
                                               a.reference_start, a.reference_end, a.is_reverse)
                            if not collapsed(pos, control, min_run, factor):
                                note = "retried"
                                break
                        if note == "moves":
                            pos = control
                    rows.append((a.reference_name, strand, a.mapping_quality, pos, "", note) if note
                                else (a.reference_name, strand, a.mapping_quality, pos, ""))
                except Exception as e:  # noqa: BLE001
                    rows.append((a.reference_name, strand, a.mapping_quality, None, str(e)[:120]))
            out.append((rid, rows))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pod5", type=Path, required=True)
    ap.add_argument("--bam", type=Path, required=True, help="minimap2 -y --MD alignment of the sup basecalls (move tags)")
    ap.add_argument("--spans", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True, help="the alignment's reference (contig order = the axis)")
    ap.add_argument("--levels", type=Path, required=True, help="Remora k-mer level table (R10.4.1 9-mer)")
    ap.add_argument("--mode", choices=("remora", "moves"), default="remora")
    ap.add_argument("--reads", type=Path, help="read ids to place (default every read of the spans file)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--gap", type=int, default=10_000_000)
    ap.add_argument("--collapse-check", action="store_true", help="check each Remora placement against the move table's")
    ap.add_argument("--collapse-min", type=int, default=500)
    ap.add_argument("--collapse-factor", type=float, default=4.0)
    ap.add_argument("--collapse-retries", type=int, default=2)
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()
    want = {l.strip() for l in open(a.reads) if l.strip()} if a.reads else None
    offs = reference_offsets(a.reference, a.gap)
    n = max(1, a.threads)
    nchunk = n * 4
    tmp = Path(tempfile.mkdtemp(prefix="remora_evref.", dir=os.environ.get("TMPDIR", "/tmp")))
    try:
        counts = split_spans(a.spans, want, nchunk, tmp)
        tasks = []
        for i in range(nchunk):
            names = read_starts(tmp / f"c{i}.spans").keys()
            if not names:
                continue
            ids = tmp / f"c{i}.txt"
            ids.write_text("".join(x + "\n" for x in names))
            sub = tmp / f"c{i}.bam"
            subprocess.run([a.samtools, "view", "-b", "-N", str(ids), "-o", str(sub), str(a.bam)], check=True)
            task = (str(sub), str(tmp / f"c{i}.spans"), str(a.pod5), str(a.levels), a.mode)
            if a.collapse_check and a.mode == "remora":
                task += ((a.collapse_min, a.collapse_factor, a.collapse_retries),)
            tasks.append(task)
        done, seg_id = set(), 0
        placed_reads = errors = ev_placed = 0
        err_examples = []
        notes = {"retried": 0, "moves": 0}
        with open(f"{a.out}.evref.bin", "wb") as fb, open(f"{a.out}.segs.tsv", "w") as fs, Pool(n) as pool:
            fs.write("seg\tread\taligner\tprimary\tcontig\tstrand\tmapq\n")
            for part in pool.imap_unordered(work, tasks):
                for rid, rows in part:
                    nev = counts[rid]
                    ref = np.full(nev, -1, np.int64)
                    seg = np.full(nev, -1, np.int64)
                    status = np.zeros(nev, np.uint8)
                    any_placed = False
                    for j, row in enumerate(rows):   # primary first: fills the events still free
                        contig, strand, mapq, pos, err = row[:5]
                        note = row[5] if len(row) > 5 else ""
                        if note:
                            notes[note] += 1
                        if pos is not None and contig in offs:
                            free = (pos >= 0) & (ref < 0)
                            ref[free] = pos[free] + offs[contig]
                            seg[free] = seg_id
                            status[free] = CONF if mapq >= a.min_mapq else LOWQ
                            any_placed = any_placed or bool(free.any())
                        elif err:
                            errors += 1
                            if len(err_examples) < 5:
                                err_examples.append(f"{rid}: {err}")
                        aligner = f"{a.mode}>moves" if note == "moves" else a.mode
                        fs.write(f"{seg_id}\t{rid}\t{aligner}\t{int(j == 0)}\t{contig}\t{strand}\t{mapq}\n")
                        seg_id += 1
                    placed_reads += int(any_placed)
                    ev_placed += int((ref >= 0).sum())
                    write_record(fb, rid, ref, seg, status)
                    done.add(rid)
            for rid, nev in counts.items():   # reads without a primary alignment: unplaced
                if rid not in done:
                    write_record(fb, rid, np.full(nev, -1, np.int64), np.full(nev, -1, np.int64), np.zeros(nev, np.uint8))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ev_total = sum(counts.values())
    lines = [f"mode {a.mode}: {len(counts)} reads with spans, {len(done)} with a primary alignment, {placed_reads} placed, "
             f"{errors} errors; events placed {ev_placed:,} of {ev_total:,} ({100 * ev_placed / max(1, ev_total):.2f}%)"]
    if a.collapse_check and a.mode == "remora":
        lines.append(f"collapse check (>= {a.collapse_min} events on one base, > {a.collapse_factor:g} x the move table's): "
                     f"{notes['retried']} records placed again by Remora, {notes['moves']} took the move table's placement")
    lines += [f"  error: {e}" for e in err_examples]
    Path(f"{a.out}.summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), file=sys.stderr)


if __name__ == "__main__":
    main()
