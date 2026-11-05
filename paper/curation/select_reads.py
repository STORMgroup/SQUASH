#!/usr/bin/env python3
"""Select a paper read set from the raw signal alone: POD5 signal and metadata, no basecalls or reference (pod5 0.3.39).

The rules, every value an option:
  1. the pool: every read of the sources; with --pool-ids only the listed reads (a pool list, ../data/pools/); with
     --summary only the reads a run's sequencing summary gives a length (duration x bases per second) of at least
     --summary-min-length;
  2. signal-estimated length = samples x bases per second / sample rate, to the nearest base (ties to even); with
     --duration-digits D the duration in seconds (samples / sample rate) is first rounded to D significant digits;
     the reads of at least --min-length, in md5 order of the read id (the MD5 hex digest of the id's text);
     --pool-coverage P keeps the prefix of that order whose signal-estimated bases reach P x G (the read reaching it
     included);
  3. noise = median absolute first difference / standard deviation of the calibrated signal (pA, float64), rounded to
     6 significant digits; a read whose noise exceeds --noise-max x the median noise of the pool is dropped;
     --open-pore N also drops a read whose longest run of samples >= 0.85 x its median_before, 2,000 samples off each
     end, is longer than N samples;
  4. the remaining reads in md5 order until their signal-estimated bases reach --coverage x G (the read reaching the
     target kept).
Outputs in OUT: reads.txt.gz (the kept reads in md5 order, gzip -9 -n), selection.tsv (every pool read in md5 order:
read_id, length, noise, open_pore_run, status = kept, noise, open_pore or unused (after the target)) and params.tsv
(every option and count).

Usage:
  select_reads.py --pod5 SOURCE [SOURCE ...] -G BP --coverage C --out OUT [options]
      SOURCE: POD5 files or folders (a folder: every *.pod5 below it, sorted by path); a read in two files is taken
      from the first
  select_reads.py --table TABLE -G BP --coverage C --out OUT [options]
      lengths, noise and open-pore runs from TABLE (an earlier run's selection.tsv, or any table with the columns
      read_id, length, noise, open_pore_run; '.' = not measured) instead of the POD5s: another coverage, pool coverage
      or noise rule over reads measured before. TABLE must hold every read the new pool takes.
  select_reads.py --pod5 SOURCE [SOURCE ...] --pool-out LIST [rules 1-2]
      write the read ids of rules 1 and 2 (before --pool-coverage), sorted, and stop (no signal is read)
"""
import os
import sys

import argparse
import gzip
import hashlib
import statistics
from multiprocessing import Pool
from pathlib import Path

import numpy as np

END, OPEN_FRACTION, TASK_READS = 2000, 0.85, 500
WANT = None  # read ids the index keeps (None: every read); set before the workers start


def die(msg):
    print(f"select_reads.py: {msg}", file=sys.stderr)
    sys.exit(1)


def log(msg):
    print(f"# {msg}", file=sys.stderr, flush=True)


def read_list(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line.split("\t")[0].strip() for line in f if line.strip()]


def write_list(path, ids):
    data = ("".join(f"{r}\n" for r in ids)).encode()
    tmp = Path(path).with_name("." + Path(path).name + ".tmp")
    if str(path).endswith(".gz"):
        with open(tmp, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as f:
            f.write(data)
    else:
        tmp.write_bytes(data)
    tmp.replace(path)


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


def md5(rid):
    return hashlib.md5(rid.encode()).hexdigest()


def signal_length(samples, rate, bps, digits):
    """signal-estimated length in bases: samples x bps / rate, or the duration rounded to `digits` significant digits
    x bps; to the nearest base, ties to even"""
    if digits:
        return round(float(f"{samples / rate:.{digits}g}") * bps)
    return round(samples * bps / rate)


def ids_of(arr):
    """read ids of a FixedSizeBinary(16) Arrow array as strings"""
    from pod5.api_utils import format_read_ids
    buf = np.frombuffer(arr.buffers()[1], dtype=np.uint8)
    start = arr.offset * 16
    return format_read_ids(buf[start:start + 16 * len(arr)].reshape(-1, 16))


def index_file(path):
    """(read id, sample count, sample rate) of every read of one POD5 (kept by WANT), read table only"""
    import pod5
    out = []
    with pod5.Reader(path) as r:
        rates = {}
        table = r.run_info_table
        for i in range(table.num_record_batches):
            b = table.get_batch(i)
            for acq, rate in zip(b.column("acquisition_id").to_pylist(), b.column("sample_rate").to_pylist()):
                rates[acq] = int(rate)
        for b in r.read_batches():
            ids = ids_of(b.read_id_column)
            n = b.columns.num_samples.to_numpy(zero_copy_only=False)
            acq = b.columns.run_info.to_pylist()
            for i, rid in enumerate(ids):
                if WANT is None or rid in WANT:
                    out.append((rid, int(n[i]), rates[acq[i]]))
    return out


def longest_run(mask):
    if not mask.any():
        return 0
    d = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
    return int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max())


def features(rec):
    """(noise to 6 significant digits, longest open-pore run) of one read"""
    x = rec.signal_pa.astype(np.float64)
    n = len(x)
    noise = float(np.median(np.abs(np.diff(x)))) if n > 1 else 0.0
    sd = float(x.std()) or 1e-6
    mb = rec.median_before
    core = x[END:-END] if n > 2 * END + 10 else x
    opn = core >= OPEN_FRACTION * mb if mb == mb and mb > 0 else np.zeros(len(core), bool)
    return float(f"{noise / sd:.6g}"), longest_run(opn)


def feature_task(task):
    import pod5
    path, ids = task
    out = []
    with pod5.Reader(path) as r:
        for rec in r.reads(selection=ids, missing_ok=True):
            out.append((str(rec.read_id), *features(rec)))
    return out


def summary_reads(paths, bps, min_len):
    """read ids with a sequencing-summary row of duration x bps >= min_len"""
    keep = set()
    for p in paths:
        opener = gzip.open if str(p).endswith(".gz") else open
        with opener(p, "rt") as f:
            head = f.readline().rstrip("\n").split("\t")
            if "read_id" not in head or "duration" not in head:
                die(f"{p}: no read_id and duration columns")
            ir, idur = head.index("read_id"), head.index("duration")
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) > max(ir, idur) and float(c[idur]) * bps >= min_len:
                    keep.add(c[ir])
    return keep


def read_table(path):
    """read id -> (length, noise or None, open-pore run or None) from a table with those columns"""
    out = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        try:
            ir, il, iq, io = (head.index(c) for c in ("read_id", "length", "noise", "open_pore_run"))
        except ValueError:
            die(f"{path}: needs the columns read_id, length, noise, open_pore_run")
        for line in f:
            c = line.rstrip("\n").split("\t")
            out[c[ir]] = (int(c[il]), None if c[iq] == "." else float(c[iq]), None if c[io] == "." else int(c[io]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], epilog="\n".join(__doc__.split("\n")[1:]),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pod5", nargs="+", metavar="SOURCE", help="POD5 files or folders")
    src.add_argument("--table", help="lengths and noise measured before (see above) instead of --pod5")
    ap.add_argument("--out", help="output folder (reads.txt.gz, selection.tsv, params.tsv)")
    ap.add_argument("--pool-out", help="write the reads of rules 1-2 (before --pool-coverage) to this list and stop")
    ap.add_argument("-G", "--genome-size", type=int, help="G in bases")
    ap.add_argument("--coverage", type=float, help="the target: COVERAGE x G signal-estimated bases (60, 100)")
    ap.add_argument("--pool-ids", help="pool list: only these reads (one id per line, plain or .gz)")
    ap.add_argument("--missing-ok", action="store_true", help="listed pool reads missing from the sources are no error")
    ap.add_argument("--summary", nargs="+", default=[],
                    help="sequencing summaries: only reads with a duration x bases per second >= --summary-min-length")
    ap.add_argument("--summary-min-length", type=float, default=9000)
    ap.add_argument("--min-length", type=float, default=10000, help="minimum signal-estimated length (10000)")
    ap.add_argument("--bases-per-second", type=float, default=400, help="translocation speed for the length (400)")
    ap.add_argument("--duration-digits", type=int, default=0,
                    help="length from the duration in seconds rounded to this many significant digits (0: exact)")
    ap.add_argument("--pool-coverage", type=float, help="pool = md5-order prefix of the long reads to this x G")
    ap.add_argument("--noise-max", type=float, default=1.15, help="noise cut as a multiple of the pool's median (1.15)")
    ap.add_argument("--open-pore", type=int, default=0, help="longest open-pore run kept, in samples (0: off)")
    ap.add_argument("--threads", type=int, default=4, help="worker processes (4)")
    a = ap.parse_args()
    if a.pool_out is None and (a.out is None or a.genome_size is None or a.coverage is None):
        ap.error("--out, -G and --coverage are required (or --pool-out)")
    if a.pool_out and a.table:
        ap.error("--pool-out needs --pod5")
    if a.pool_coverage is not None and a.genome_size is None:
        ap.error("--pool-coverage needs -G")
    bps, digits, G = a.bases_per_second, a.duration_digits, a.genome_size
    srcs = a.pod5 or []
    p = {"pod5": " ".join(srcs) if len(srcs) <= 3 else f"{srcs[0]} ... {srcs[-1]} ({len(srcs)} paths)" if srcs else ".",
         "table": a.table or ".", "pool_ids": a.pool_ids or ".",
         "summary": " ".join(a.summary) or ".", "summary_min_length": f"{a.summary_min_length:g}" if a.summary else ".",
         "min_length": f"{a.min_length:g}", "bases_per_second": f"{bps:g}", "duration_digits": digits or ".",
         "pool_coverage": "." if a.pool_coverage is None else f"{a.pool_coverage:g}", "noise_max": f"{a.noise_max:g}",
         "open_pore": a.open_pore or ".", "genome_size": G if G is not None else ".",
         "coverage": "." if a.coverage is None else f"{a.coverage:g}"}

    # rule 1: the pool's reads, with their lengths (rule 2)
    global WANT
    listed = None
    if a.pool_ids:
        listed = set(read_list(a.pool_ids))
        WANT = set(listed)
        p["pool_ids_listed"] = len(listed)
    if a.summary:
        long_in_summary = summary_reads(a.summary, bps, a.summary_min_length)
        p["summary_reads"] = len(long_in_summary)
        WANT = long_in_summary if WANT is None else WANT & long_in_summary
    table, files, where = None, [], {}
    length = {}
    if a.table:
        table = read_table(a.table)
        p["table_reads"] = len(table)
        for rid, (n, _, _) in table.items():
            if WANT is None or rid in WANT:
                length[rid] = n
    else:
        import pod5
        p["pod5_version"], p["numpy_version"] = pod5.__version__, np.__version__
        files = pod5_files(a.pod5)
        if not files:
            die("no POD5 file found")
        log(f"indexing {len(files):,} POD5 files")
        dups = 0
        with Pool(a.threads) as pool:
            for fi, rows in enumerate(pool.imap(index_file, files, chunksize=4)):
                for rid, n, rate in rows:
                    if rid in where:
                        dups += 1
                        continue
                    where[rid] = fi
                    length[rid] = signal_length(n, rate, bps, digits)
        p["pod5_files"], p["duplicate_reads"] = len(files), dups
    p["candidates"] = len(length)
    if listed is not None:
        missing = len(listed - set(length))
        p["pool_ids_missing"] = missing
        if missing and table is None and not a.missing_ok:
            die(f"{missing:,} of the {len(listed):,} listed pool reads are not in the sources (--missing-ok to go on)")
    long_reads = sorted((r for r, n in length.items() if n >= a.min_length), key=md5)
    p["long_reads"] = len(long_reads)
    p["long_bases"] = sum(length[r] for r in long_reads)
    log(f"{len(length):,} reads, {len(long_reads):,} of >= {a.min_length:g} bases")
    if a.pool_out:
        write_list(a.pool_out, sorted(long_reads))
        log(f"{len(long_reads):,} read ids -> {a.pool_out}")
        return
    pool_ids, acc = [], 0
    for r in long_reads:
        if a.pool_coverage is not None and acc >= a.pool_coverage * G:
            break
        pool_ids.append(r)
        acc += length[r]
    p["pool_reads"], p["pool_bases"], p["pool_depth"] = len(pool_ids), acc, f"{acc / G:.2f}"
    if not pool_ids:
        die("the pool is empty")

    # rule 3: the pool's noise and open-pore runs
    feat = {}
    if table is not None:
        for r in pool_ids:
            _, q, o = table[r]
            if q is None or (a.open_pore and o is None):
                die(f"{a.table}: no noise / open-pore run for pool read {r}")
            feat[r] = (q, o)
    else:
        by_file = {}
        for r in pool_ids:
            by_file.setdefault(where[r], []).append(r)
        tasks = [(files[fi], ids[i:i + TASK_READS]) for fi, ids in sorted(by_file.items())
                 for i in range(0, len(ids), TASK_READS)]
        log(f"signal features of {len(pool_ids):,} pool reads from {len(by_file):,} files")
        with Pool(a.threads) as pool:
            for rows in pool.imap_unordered(feature_task, tasks):
                for rid, q, o in rows:
                    feat[rid] = (q, o)
        if len(feat) != len(pool_ids):
            die(f"{len(pool_ids) - len(feat):,} pool reads not read back from the sources")
    med = statistics.median(feat[r][0] for r in pool_ids)
    cut = a.noise_max * med
    p["noise_median"], p["noise_cut"] = repr(med), repr(cut)

    # rule 4: md5 order to the target
    target = a.coverage * G
    status, kept, acc = {}, [], 0
    seen = noisy = open_pore = 0
    for r in pool_ids:
        if acc >= target:
            status[r] = "unused"
            continue
        seen += 1
        q, o = feat[r]
        if q > cut:
            status[r] = "noise"
            noisy += 1
        elif a.open_pore and o > a.open_pore:
            status[r] = "open_pore"
            open_pore += 1
        else:
            status[r] = "kept"
            kept.append(r)
            acc += length[r]
    p.update({"target_bases": f"{target:.0f}", "seen": seen, "dropped_noise": noisy, "dropped_open_pore": open_pore,
              "kept": len(kept), "kept_bases": acc, "depth": f"{acc / G:.2f}"})
    if acc < target:
        log(f"warning: {acc / G:.2f}x after the rules, {a.coverage:g}x asked (the pool is exhausted)")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_list(out / "reads.txt.gz", kept)
    with open(out / ".selection.tsv", "w") as f:
        f.write("read_id\tlength\tnoise\topen_pore_run\tstatus\n")
        for r in pool_ids:
            q, o = feat[r]
            f.write(f"{r}\t{length[r]}\t{q:.6g}\t{'.' if o is None else o}\t{status[r]}\n")
    (out / ".selection.tsv").replace(out / "selection.tsv")
    with open(out / ".params.tsv", "w") as f:
        f.write("name\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in p.items()))
    (out / ".params.tsv").replace(out / "params.tsv")
    print(f"{len(kept)}\t{acc}\t{acc / G:.2f}\t{seen}\t{noisy}\t{open_pore}\t{med:.4f}\t{cut:.4f}")


if __name__ == "__main__":
    main()
