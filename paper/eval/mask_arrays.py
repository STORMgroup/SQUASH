#!/usr/bin/env python3
"""The A. thaliana reference with Col-CC's tandem arrays masked (the supplement's masked QUAST of D8; evaluation only).

For each chromosome named by --mask, from its own sequence: per 10-kb window, the fraction of its 21-mers (canonical,
without N) that occur >= 20 times in the chromosome, written to three decimals; windows of >= 0.3 are array windows,
merged across gaps < 100 kb; clusters spanning >= 100 kb are the arrays (CEN180 centromeres, the 45S rDNA arrays, further
satellite blocks). The arrays are written as N; the other records (the organelles) are copied. Writes OUT.fa (80 bases a
line, the input's header lines) and OUT.bed (record, start, end, `array`; 0-based, end exclusive); checks that the Ns
added equal the arrays' length minus the Ns already there.
Usage: mask_arrays.py REF.fa OUT.fa OUT.bed --mask NAME[,NAME...]
"""
import argparse
import sys

import numpy as np

K, MINCOPY, WINDOW, FRAC, JOIN, MIN_SPAN = 21, 20, 10_000, 0.3, 100_000, 100_000


def read_fa(path):
    name, seq, recs = None, [], []
    for line in open(path):
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].rstrip("\n"), []
        else:
            seq.append(line.strip())
    if name:
        recs.append((name, "".join(seq)))
    return recs


def high_copy_fractions(seq: str) -> list[str]:
    """Per 10-kb window: the fraction of its valid 21-mers occurring >= MINCOPY times (both strands), as '%.3f'."""
    s = seq.upper()
    code = np.frombuffer(s.encode(), dtype=np.uint8)
    lut = np.full(256, 4, np.uint8)
    for i, c in enumerate(b"ACGT"):
        lut[c] = i
    x = lut[code].astype(np.uint64)
    n = len(x) - K + 1
    valid = np.convolve((x == 4).astype(np.int32), np.ones(K, np.int32), "valid") == 0
    fw = np.zeros(n, np.uint64)
    rc = np.zeros(n, np.uint64)
    xc = np.where(x == 4, 0, x)
    for j in range(K):
        fw = (fw << np.uint64(2)) | xc[j:j + n]
        rc = rc | ((np.uint64(3) - xc[j:j + n]) << np.uint64(2 * j))
    can = np.minimum(fw, rc)[valid]
    pos = np.flatnonzero(valid)
    _, inv, cnt = np.unique(can, return_inverse=True, return_counts=True)
    high = cnt[inv] >= MINCOPY
    nw = -(-len(s) // WINDOW)
    tot = np.bincount(pos // WINDOW, minlength=nw)
    hi = np.bincount(pos // WINDOW, weights=high, minlength=nw)
    return [f"{hi[w] / max(1, tot[w]):.3f}" for w in range(nw)]


def arrays(seq: str) -> list[list[int]]:
    """The array intervals [start, end) of one chromosome."""
    fr = high_copy_fractions(seq)
    hits = [(w * WINDOW, min((w + 1) * WINDOW, len(seq))) for w, f in enumerate(fr) if float(f) >= FRAC]
    clusters = []
    for a, b in hits:
        if clusters and a - clusters[-1][1] < JOIN:
            clusters[-1][1] = b
        else:
            clusters.append([a, b])
    return [c for c in clusters if c[1] - c[0] >= MIN_SPAN]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("out")
    ap.add_argument("bed")
    ap.add_argument("--mask", required=True, help="the records to mask, comma-separated (the chromosomes)")
    a = ap.parse_args()
    names = set(a.mask.split(","))
    total = total_new = 0
    with open(a.out + ".part", "w") as fh, open(a.bed + ".part", "w") as bh:
        for header, s in read_fa(a.ref):
            cid = header.split()[0]
            if cid in names:
                names.discard(cid)
                ivs = arrays(s)
                merged = []
                for x, y in sorted(ivs):
                    if merged and x <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], y)
                    else:
                        merged.append([x, y])
                arr = bytearray(s, "ascii")
                before = arr.count(b"N") + arr.count(b"n")
                for x, y in merged:
                    arr[x:y] = b"N" * (y - x)
                    bh.write(f"{cid}\t{x}\t{y}\tarray\n")
                after = arr.count(b"N") + arr.count(b"n")
                covered = sum(y - x for x, y in merged)
                already = sum(s[x:y].upper().count("N") for x, y in merged)
                assert after - before == covered - already, (cid, after - before, covered, already)
                total, total_new = total + covered, total_new + after - before
                print(f"{cid}\t{len(s)}\tmasked {covered}\tnew N {after - before}: "
                      + ", ".join(f"{x / 1e6:.2f}-{y / 1e6:.2f}" for x, y in merged), file=sys.stderr)
                s = arr.decode()
            else:
                print(f"{cid}\t{len(s)}\tunchanged", file=sys.stderr)
            fh.write(f">{header}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    if names:
        sys.exit(f"{a.ref}: no record {sorted(names)}")
    import os
    os.replace(a.bed + ".part", a.bed)
    os.replace(a.out + ".part", a.out)
    print(f"masked {total} bp ({total_new} new N)", file=sys.stderr)


if __name__ == "__main__":
    main()
