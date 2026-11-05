#!/usr/bin/env python3
"""Cut a set's reads out of its source POD5s into one POD5, in the order of the set's list (pod5 0.3.39).

SQUASH's signal step reads a POD5 in blocks of 40,000 reads of file order, so a set's POD5 must hold its reads in the
order of its list (paper/data/sets/KEY/reads.txt.gz); `pod5 filter` and `pod5 subset` keep no order. Steps (as the
paper's set POD5s were cut):
  1. index: the read ids of every source file (read tables only); every listed read must be found (a read in two
     source files is taken from the first, in path order);
  2. parts: the source files that hold listed reads, in path order, grouped into about 4 x PROCS parts of similar
     read counts; each part (PARTS/part_NNNNN.pod5, PROCS workers) has its files' listed reads, file after file, each
     in file order;
  3. join: one writer adds every run info, pore type and end reason first, so the read table needs no Arrow delta
     dictionary (SQUASH's POD5 reader, pod5lib 0.2.0, refuses those), then the reads in the list's order;
  4. check: the output's read ids in file order equal the list (then run pod5_digest.py against sets.tsv).
Usage:
  pod5_cut.py cut  LIST OUT.pod5 SOURCE [SOURCE ...] [--procs N] [--parts DIR]
      SOURCE: POD5 files or folders (a folder: every *.pod5 below it, sorted by path)
  pod5_cut.py part LIST SOURCE.pod5 PART.pod5
      the listed reads of one source file, in its file order (for downloads streamed file by file: cut a file, then
      delete it; no PART is written when the file holds no listed read); parts are POD5s, so `cut` takes them as sources
  pod5_cut.py join LIST OUT.pod5 PART [PART ...]
      parts made by `part` (in any order) into OUT, as step 3 + 4
LIST: one read id per line, plain or .gz. OUT is written as a hidden temporary file and renamed when complete.
"""
import argparse
import gzip
import shutil
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pod5
from pod5.api_utils import format_read_ids


def die(msg):
    print(f"pod5_cut.py: {msg}", file=sys.stderr)
    sys.exit(1)


def read_list(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        ids = [line.strip() for line in f if line.strip()]
    if len(set(ids)) != len(ids):
        die(f"{path}: duplicate read ids")
    return ids


def pod5_files(paths):
    files = []
    for p in map(Path, paths):
        if p.is_dir():
            files += sorted(p.rglob("*.pod5"))
        elif p.is_file():
            files.append(p)
        else:
            die(f"no such POD5 file or folder: {p}")
    return files


def ids_of(arr):
    buf = np.frombuffer(arr.buffers()[1], dtype=np.uint8)
    start = arr.offset * 16
    return format_read_ids(buf[start:start + 16 * len(arr)].reshape(-1, 16))


def file_ids(path):
    with pod5.Reader(path) as r:
        return r.read_ids


def tmp_of(dst):
    return dst.with_name("." + dst.name + ".tmp")


def write_part(sources, out):
    """the listed reads of each (source file, ids) in `sources`, file after file, each in its file order, into `out`;
    returns the count and the run infos, pore types and end reasons seen"""
    runs, pores, ends, n = {}, [], [], 0
    tmp = tmp_of(out)
    tmp.unlink(missing_ok=True)
    with pod5.Writer(tmp) as w:
        for src, ids in sources:
            k = 0
            with pod5.Reader(src) as r:
                for rec in r.reads(selection=ids, missing_ok=True):
                    ri = rec.run_info
                    runs.setdefault(ri.acquisition_id, ri)
                    p = str(rec.pore.pore_type)
                    if p not in pores:
                        pores.append(p)
                    if rec.end_reason not in ends:
                        ends.append(rec.end_reason)
                    w.add_read(rec.to_read())
                    k += 1
            if k != len(ids):
                raise RuntimeError(f"{src}: wrote {k} of the {len(ids)} listed reads it holds")
            n += k
    tmp.rename(out)
    return n, list(runs.values()), pores, ends


def _part_job(job):
    i, sources, out = job
    n, runs, pores, ends = write_part(sources, out)
    return i, n, runs, pores, ends


def groups(per_file, procs):
    """consecutive source files (path order) grouped into parts of about equal read counts: few part files even for
    sources of thousands of files (each POD5 has a few MB of overhead)"""
    total = sum(len(h) for _, _, h in per_file)
    target = max(1, -(-total // max(4 * procs, 1)))
    out, cur, n = [], [], 0
    for _, f, hits in per_file:
        cur.append((f, hits))
        n += len(hits)
        if n >= target:
            out.append(cur)
            cur, n = [], 0
    if cur:
        out.append(cur)
    return out


def part_metadata(part):
    """run infos, pore types and end reasons of a part, in order of first appearance (as pod5_rewrite.py)"""
    runs, pores, ends = {}, [], []
    with pod5.Reader(part) as r:
        for rec in r.reads():
            ri = rec.run_info
            runs.setdefault(ri.acquisition_id, ri)
            p = str(rec.pore.pore_type)
            if p not in pores:
                pores.append(p)
            if rec.end_reason not in ends:
                ends.append(rec.end_reason)
    return list(runs.values()), pores, ends


def join(want, dst, parts, meta=None):
    """parts -> dst: every run info, pore type and end reason first, then the reads in the order of `want`"""
    if meta is None:
        meta = [part_metadata(p) for p in parts]
    runs, pores, ends = {}, [], []
    for rs, ps, es in meta:
        for ri in rs:
            runs.setdefault(ri.acquisition_id, ri)
        pores += [p for p in ps if p not in pores]
        ends += [e for e in es if e not in ends]
    readers = [pod5.Reader(p) for p in parts]
    where = {}
    for pi, r in enumerate(readers):
        for bi in range(r.batch_count):
            for row, rid in enumerate(ids_of(r.get_batch(bi).read_id_column)):
                if rid in where:
                    die(f"read {rid} is in two parts")
                where[rid] = (pi, bi, row)
    missing = [x for x in want if x not in where]
    if missing:
        die(f"{len(missing):,} listed reads are in no part (first: {missing[0]})")
    if len(where) != len(want):
        die(f"the parts hold {len(where) - len(want):,} reads that are not listed")
    tmp = tmp_of(dst)
    tmp.unlink(missing_ok=True)
    cache = {}
    buf = []
    with pod5.Writer(tmp) as w:
        for ri in runs.values():
            w.add(ri)
        for p in pores:
            w.add(p)
        for e in ends:
            w.add(e)
        for rid in want:
            pi, bi, row = where[rid]
            got = cache.get(pi)
            if got is None or got[0] != bi:
                got = cache[pi] = (bi, readers[pi].get_batch(bi))
            rec = got[1].get_read(row)
            if str(rec.read_id) != rid:
                die(f"part {parts[pi]} batch {bi} row {row}: {rec.read_id}, expected {rid}")
            buf.append(rec.to_read())
            if len(buf) >= 100:
                w.add_reads(buf)
                buf = []
        if buf:
            w.add_reads(buf)
    for r in readers:
        r.close()
    got = file_ids(tmp)
    if got != want:
        i = next((i for i, (a, b) in enumerate(zip(got, want)) if a != b), min(len(got), len(want)))
        die(f"{tmp}: read ids differ from the list at read {i + 1:,}")
    tmp.rename(dst)
    print(f"wrote {len(want):,} reads to {dst} ({len(runs)} run infos, {len(pores)} pore types, {len(ends)} end "
          f"reasons; from {len(parts)} parts)", file=sys.stderr)


def cmd_cut(a):
    want = read_list(a.list)
    dst = Path(a.out)
    if dst.exists():
        die(f"{dst} exists")
    files = pod5_files(a.sources)
    wanted = set(want)
    found = {}
    per_file = []
    dup = 0
    for fi, f in enumerate(files):
        hits = []
        for x in file_ids(f):
            if x in wanted:
                if x in found:  # a read in two source files: the first file (in path order) holds it
                    dup += 1
                    if dup <= 3:
                        print(f"note: read {x} is in {files[found[x]]} and {f}: taken from the first", file=sys.stderr)
                    continue
                found[x] = fi
                hits.append(x)
        if hits:
            per_file.append((fi, f, hits))
    if dup:
        print(f"note: {dup:,} listed reads are in more than one source file (the digest check tells whether that matters)",
              file=sys.stderr)
    if len(found) != len(want):
        miss = next(x for x in want if x not in found)
        die(f"{len(want) - len(found):,} listed reads are in no source file (first: {miss})")
    parts_dir = Path(a.parts) if a.parts else dst.with_name(dst.name + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(gi, g, parts_dir / f"part_{gi:05d}.pod5") for gi, g in enumerate(groups(per_file, a.procs))]
    print(f"{len(want):,} listed reads in {len(per_file)} of {len(files)} source files -> {len(jobs)} parts in "
          f"{parts_dir}", file=sys.stderr)
    with Pool(max(1, a.procs)) as pool:
        res = sorted(pool.imap_unordered(_part_job, jobs), key=lambda x: x[0])
    if sum(n for _, n, *_ in res) != len(want):
        die("the parts hold another number of reads than the list")
    parts = [parts_dir / f"part_{i:05d}.pod5" for i, *_ in res]
    join(want, dst, parts, meta=[(rs, ps, es) for _, _, rs, ps, es in res])
    if not a.parts:
        shutil.rmtree(parts_dir)


def cmd_part(a):
    want = read_list(a.list)
    src, out = Path(a.source), Path(a.part)
    if out.exists():
        die(f"{out} exists")
    wanted = set(want)
    hits = [x for x in file_ids(src) if x in wanted]
    if not hits:
        print(f"{src.name}: no listed read (no part written)", file=sys.stderr)
        return
    n, *_ = write_part([(src, hits)], out)
    print(f"{src.name}: {n:,} listed reads -> {out}", file=sys.stderr)


def cmd_join(a):
    want = read_list(a.list)
    dst = Path(a.out)
    if dst.exists():
        die(f"{dst} exists")
    parts = [Path(p) for p in a.parts]
    join(want, dst, [p for p in parts if file_ids(p)])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cut", help="index + parts + join")
    c.add_argument("list")
    c.add_argument("out")
    c.add_argument("sources", nargs="+")
    c.add_argument("--procs", type=int, default=4, help="parallel part writers (default 4)")
    c.add_argument("--parts", help="keep the parts in this folder (default: OUT.parts, removed at the end)")
    c.set_defaults(func=cmd_cut)
    p = sub.add_parser("part", help="the listed reads of one source file")
    p.add_argument("list")
    p.add_argument("source")
    p.add_argument("part")
    p.set_defaults(func=cmd_part)
    j = sub.add_parser("join", help="parts -> one POD5 in the list's order")
    j.add_argument("list")
    j.add_argument("out")
    j.add_argument("parts", nargs="+")
    j.set_defaults(func=cmd_join)
    a = ap.parse_args()
    if pod5.__version__ != "0.3.39":
        print(f"pod5_cut.py: note: pod5 {pod5.__version__}; the paper's POD5s were written with 0.3.39",
              file=sys.stderr)
    a.func(a)


if __name__ == "__main__":
    main()
