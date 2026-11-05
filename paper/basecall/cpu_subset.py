#!/usr/bin/env python3
"""The reads of the CPU basecalling route's timed subset (pod5 0.3.39).

The paper timed CPU basecalling of the large sets on a subset of about one hour of signal and scaled it to the set: the
set's reads in MD5 order of their read id (hex digest of the id's text), taken until their samples reach TARGET (the
read that crosses it is kept): 3.67e9 samples at 5 kHz, 2.33e9 at 4 kHz; TARGET `all` takes every read (the bacteria,
basecalled whole). The set's time is then its samples / (the subset's samples / Dorado's "Finished in" time).

Usage: cpu_subset.py KEY TARGET SET.pod5 OUT_DIR
Writes OUT_DIR/KEY.cpu_subset.ids.txt (the subset's read ids, MD5 order) and OUT_DIR/KEY.cpu_subset.tsv (key, set_reads,
set_samples, sample_rate, target_samples, subset_reads, subset_samples, subset_fraction).
"""
import hashlib
import sys
from pathlib import Path

import pod5


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    key, target, src, out = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
    rows = {}
    with pod5.Reader(src) as r:
        for rec in r.reads():
            rows[str(rec.read_id)] = (rec.num_samples, rec.run_info.sample_rate)
    order = sorted(rows, key=lambda x: hashlib.md5(x.encode()).hexdigest())
    tgt = float("inf") if target == "all" else int(float(target))
    n = acc = 0
    while n < len(order) and acc < tgt:
        acc += rows[order[n]][0]
        n += 1
    total = sum(v[0] for v in rows.values())
    rates = sorted({v[1] for v in rows.values()})
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{key}.cpu_subset.ids.txt").write_text("".join(f"{a}\n" for a in order[:n]))
    (out / f"{key}.cpu_subset.tsv").write_text(
        "key\tset_reads\tset_samples\tsample_rate\ttarget_samples\tsubset_reads\tsubset_samples\tsubset_fraction\n"
        f"{key}\t{len(rows)}\t{total}\t{','.join(map(str, rates))}\t{target}\t{n}\t{acc}\t{acc / total:.4f}\n")
    print(f"{key}: {len(rows):,} reads, {total:,} samples; subset {n:,} reads, {acc:,} samples ({acc / total:.1%})")


if __name__ == "__main__":
    main()
