#!/usr/bin/env python3
"""The signal digest of a POD5: the check that a set's POD5 holds what the paper's run read (pod5 0.3.39).

File SHA-256s of POD5s are not comparable across writers (the same reads written twice give other bytes), so a set's
POD5 is checked by its content, in file order. digest = SHA-256 of the line "squash-pod5-digest-v1\\n" followed by,
for every read in file order (files in the order given):

    read id              36 ASCII bytes (the canonical lower-case UUID)
    calibration offset   float32, little endian (the read table's value)
    calibration scale    float32, little endian
    sample rate          uint16, little endian (from the read's run info)
    sample count         uint64, little endian
    samples              sample count x int16, little endian (the raw ADC signal, decompressed)

These are the values SQUASH reads from a POD5: the signal step reads the id, calibration and samples, the bp-scale step
the sample count and sample rate. Two POD5s with equal digests give SQUASH the same input in the same order.

Usage:
  pod5_digest.py [--ids LIST] [--expect HEX] POD5 [POD5 ...]
      LIST (one read id per line, plain or .gz): first check that the POD5's read ids in file order equal LIST (read
      table only); stop at the first differing read. Then the digest (every sample is read).
  pod5_digest.py --select LIST [--expect HEX] SOURCE [SOURCE ...]
      the digest of LIST's reads read straight from SOURCE POD5s (files or folders, folders as sorted(rglob('*.pod5'))),
      without writing a POD5: the sources must yield LIST's reads in LIST's order (each file in file order, the files in
      the order given), e.g. to check a set's reads in their sources before cutting them.
Prints one TSV line: reads, samples, digest. Exit 1: --expect differs; 2: read ids or order differ from LIST;
3: a malformed read.
"""
import argparse
import gzip
import hashlib
import queue
import struct
import sys
import threading
from pathlib import Path

import numpy as np
import pod5
from pod5.api_utils import format_read_ids

MAGIC = b"squash-pod5-digest-v1\n"
HEADER = struct.Struct("<ffHQ")


def die(msg, code):
    print(msg, file=sys.stderr)
    sys.exit(code)


def read_list(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def pod5_files(paths):
    files = []
    for p in map(Path, paths):
        files += sorted(p.rglob("*.pod5")) if p.is_dir() else [p]
    return files


def ids_of(arr):
    """read ids of a FixedSizeBinary(16) Arrow array as strings"""
    buf = np.frombuffer(arr.buffers()[1], dtype=np.uint8)
    start = arr.offset * 16
    return format_read_ids(buf[start:start + 16 * len(arr)].reshape(-1, 16))


def sample_rates(reader):
    rates = {}
    table = reader.run_info_table
    for i in range(table.num_record_batches):
        b = table.get_batch(i)
        for acq, rate in zip(b.column("acquisition_id").to_pylist(), b.column("sample_rate").to_pylist()):
            rates[acq] = int(rate)
    return rates


def file_read_ids(files):
    out = []
    for f in files:
        with pod5.Reader(f) as r:
            out += r.read_ids
    return out


def check_ids(files, want):
    got = file_read_ids(files)
    for i, (a, b) in enumerate(zip(got, want)):
        if a != b:
            die(f"read ids differ from the list at read {i + 1:,}: POD5 {a}, list {b}", 2)
    if len(got) != len(want):
        die(f"read count differs: POD5 {len(got):,}, list {len(want):,}", 2)
    print(f"# read ids in file order equal the list ({len(got):,} reads)", file=sys.stderr)


def batches(files, select=None):
    """(ids, offsets, scales, rates, counts, samples) per read batch, in file order"""
    for f in files:
        with pod5.Reader(f) as r:
            rates = sample_rates(r)
            if select is None:
                it = r.read_batches(preload={"samples"})
            else:
                it = r.read_batches(selection=select, missing_ok=True, preload={"samples"})
            for b in it:
                cols = b.columns
                rows = getattr(b, "_selected_batch_rows", None)  # pod5 0.3.39: the rows a selection keeps
                take = (lambda c: c.take(rows)) if rows is not None else (lambda c: c)
                ids = ids_of(b.read_id_column)
                if not ids:
                    continue
                off = take(cols.calibration_offset).to_numpy(zero_copy_only=False).astype("<f4")
                scale = take(cols.calibration_scale).to_numpy(zero_copy_only=False).astype("<f4")
                acq = take(cols.run_info).to_pylist()  # acquisition ids (a dictionary column)
                n_tab = take(cols.num_samples).to_numpy(zero_copy_only=False)
                samples = b.cached_samples_column
                if len(samples) != len(ids):
                    die(f"{f}: {len(samples)} signals for {len(ids)} reads in a batch", 3)
                yield f, ids, off, scale, [rates[a] for a in acq], n_tab, samples


def digest(files, select_list=None):
    h = hashlib.sha256(MAGIC)
    q = queue.Queue(maxsize=4)
    err = []

    def hasher():
        while True:
            item = q.get()
            if item is None:
                return
            for head, s in item:
                h.update(head)
                h.update(s)

    t = threading.Thread(target=hasher, daemon=True)
    t.start()
    n_reads = n_samples = 0
    pos = 0
    want = select_list
    try:
        for f, ids, off, scale, rates, n_tab, samples in batches(files, want):
            item = []
            for i, rid in enumerate(ids):
                if want is not None:
                    if pos >= len(want) or want[pos] != rid:
                        exp = want[pos] if pos < len(want) else "(end of list)"
                        die(f"--select: read {pos + 1:,} of the sources is {rid}, the list has {exp} ({f}): the list is "
                            "not in the sources' order", 2)
                    pos += 1
                s = samples[i]
                if s.dtype != np.int16:
                    die(f"{f}: read {rid}: signal is {s.dtype}, not int16", 3)
                n = len(s)
                if n != int(n_tab[i]):
                    die(f"{f}: read {rid}: {n} samples, the read table says {int(n_tab[i])}", 3)
                if sys.byteorder != "little":
                    s = s.astype("<i2")
                item.append((rid.encode("ascii") + HEADER.pack(float(off[i]), float(scale[i]), rates[i], n),
                             np.ascontiguousarray(s)))
                n_reads += 1
                n_samples += n
            q.put(item)
    except BaseException as e:  # stop the hashing thread before leaving
        err.append(e)
    q.put(None)
    t.join()
    if err:
        raise err[0]
    if want is not None and pos != len(want):
        die(f"--select: {len(want) - pos:,} listed reads not found in the sources", 2)
    return n_reads, n_samples, h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pod5", nargs="+", help="POD5 file(s) (or, with --select, source files / folders)")
    ap.add_argument("--ids", help="the set's read list (reads.txt.gz): the POD5's ids in file order must equal it")
    ap.add_argument("--select", help="digest of this list's reads read from the sources, in the list's order")
    ap.add_argument("--expect", help="the digest the paper's POD5 had (sets.tsv): exit 1 when different")
    a = ap.parse_args()
    if a.ids and a.select:
        ap.error("--ids and --select exclude each other")
    if pod5.__version__ != "0.3.39":
        print(f"# note: pod5 {pod5.__version__} (the digests were made with 0.3.39; the digest does not depend on it)",
              file=sys.stderr)
    files = pod5_files(a.pod5)
    if not files:
        sys.exit("no POD5 file found")
    if a.ids:
        check_ids(files, read_list(a.ids))
    sel = read_list(a.select) if a.select else None
    n_reads, n_samples, hexd = digest(files, sel)
    print(f"{n_reads}\t{n_samples}\t{hexd}")
    if a.expect and a.expect.strip().lower() != hexd:
        print(f"digest differs: got {hexd}, expected {a.expect}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
