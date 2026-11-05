#!/usr/bin/env python3
"""The D. melanogaster reference with its three repeat-rich regions masked (the supplement's masked QUAST of D7;
evaluation only): the rDNA locus, chrY and the pericentromeric heterochromatin, from UCSC's dm6 annotation tables:
- rdna: within 20 kb of an rmsk element of class rRNA (merged);
- chrY: the whole chromosome;
- het: on the five arms (chr2L, chr2R, chr3L, chr3R, chrX), the unbanded proximal (pericentromeric) end, outside the
  polytene-banded part (cytoBandIdeo: from the first named band's start to the last one's end); the centromere is at
  the high end of chr2L, chr3L and chrX and at coordinate 0 of chr2R and chr3R.
Those intervals are written as N into the release 6 nuclear reference (the seven sequences, named by accession); the
other parts (euchromatin, telomeric ends, chr4 but its rDNA) are copied. Writes OUT.fa (80 bases a line, the input's
header lines) and OUT.bed (accession, start, end, chromosome:region; 0-based, end exclusive); checks that the Ns added
equal the intervals minus the Ns already there.
Usage: mask_dmel.py REF.fa cytoBandIdeo.txt.gz rmsk.txt.gz OUT.fa OUT.bed
"""
import gzip
import os
import sys
from collections import defaultdict

ACC = {"NC_004354.4": "chrX", "NT_033779.5": "chr2L", "NT_033778.4": "chr2R", "NT_037436.4": "chr3L",
       "NT_033777.3": "chr3R", "NC_004353.4": "chr4", "NC_024512.1": "chrY"}
ARMS = ("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
DISTAL_LOW = ("chr2L", "chr3L", "chrX")   # arms whose coordinate 0 is the telomere (the centromere at the high end)
RDNA_BP = 20000
MASK = ("rdna", "chrY", "het")


def banded(path):
    """per arm: (first banded start, last banded end) of the bands with a name"""
    lo, hi = {}, {}
    with gzip.open(path, "rt") as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")[:4]
            if c in ARMS and name:
                lo[c] = min(lo.get(c, 1 << 40), int(s))
                hi[c] = max(hi.get(c, 0), int(e))
    return {c: (lo[c], hi[c]) for c in ARMS}


def rrna(path):
    rr = defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[11] == "rRNA":
                rr[x[5]].append((int(x[6]), int(x[7])))
    return rr


def merge(iv):
    iv.sort()
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def intersect(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append([s, e])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def subtract(a, b):          # a minus b, both merged
    out, j = [], 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur:
                out.append([cur, b[k][0]])
            cur = max(cur, b[k][1])
            k += 1
        if cur < e:
            out.append([cur, e])
    return out


def regions(length, band, rr):
    """chromosome -> sorted (start, end, region) covering it"""
    reg = {}
    for c, n in length.items():
        rd = merge([[max(0, s - RDNA_BP), min(n, e + RDNA_BP)] for s, e in rr.get(c, [])])
        if c in ("chr4", "chrY"):
            base = [(0, n, c)]
        else:
            lo, hi = band[c]
            near0 = "tel" if c in DISTAL_LOW else "het"
            far = "het" if c in DISTAL_LOW else "tel"
            base = [(0, lo, near0), (lo, hi, "eu"), (hi, n, far)]
        out = []
        for s, e, r in base:
            for a, b in subtract([[s, e]], rd):
                out.append((a, b, r))
            for a, b in intersect([[s, e]], rd):
                out.append((a, b, "rdna"))
        reg[c] = sorted(out)
    return reg


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


def main():
    ref, cyto, rmsk, out, bed = sys.argv[1:6]
    recs = read_fa(ref)
    acc_of = {v: k for k, v in ACC.items()}
    length = {ACC[h.split()[0]]: len(s) for h, s in recs}
    reg = regions(length, banded(cyto), rrna(rmsk))
    nonN = sum(len(s) - s.upper().count("N") for _, s in recs)
    total = total_new = 0
    with open(out + ".part", "w") as fh, open(bed + ".part", "w") as bh:
        for header, s in recs:
            c = ACC[header.split()[0]]
            ivs = merge([[a, b] for a, b, r in reg[c] if r in MASK])
            arr = bytearray(s, "ascii")
            before = arr.count(b"N") + arr.count(b"n")
            for a, b in ivs:
                assert 0 <= a < b <= len(arr), (c, a, b, len(arr))
                arr[a:b] = b"N" * (b - a)
                bh.write(f"{acc_of[c]}\t{a}\t{b}\t{c}:{next(r for x, y, r in reg[c] if x <= a < y)}\n")
            after = arr.count(b"N") + arr.count(b"n")
            covered = sum(b - a for a, b in ivs)
            already = sum(s[a:b].upper().count("N") for a, b in ivs)
            assert after - before == covered - already, (c, after - before, covered, already)
            total, total_new = total + covered, total_new + after - before
            fh.write(f">{header}\n")
            for i in range(0, len(arr), 80):
                fh.write(arr[i:i + 80].decode() + "\n")
            print(f"{c}\t{len(s)}\tmasked {covered}\t(new N {after - before})", file=sys.stderr)
    os.replace(bed + ".part", bed)
    os.replace(out + ".part", out)
    print(f"masked {total} bp ({total_new} new N) of {sum(length.values())} ({total_new / nonN:.4f} of the non-N genome "
          f"{nonN})", file=sys.stderr)


if __name__ == "__main__":
    main()
