#!/usr/bin/env python3
"""The steps of the two pool lists made with basecalls or a reference (sdys_pool.sh, hs13q_pool.sh); plain Python, no
dependencies. Every list is written sorted by read id.

Usage:
  pools.py sample READ_BC.tsv N MIN_BASES OUT.txt
      READ_BC.tsv: read id, barcode, basecalled length (one read per line, in the order of the demultiplexed BAM).
      The first N reads of each barcode with at least MIN_BASES basecalled bases, in that order -> OUT.txt (one id per
      line, in table order).
  pools.py barcodes READ_BC.tsv SAMPLE.paf N MIN_BASES OUT.tsv
      barcode -> strain: each sampled read (as `sample` picks them) is placed by its primary record (tp:A:P, MAPQ >= 1)
      with the longest aligned query span; a barcode is a strain's (the PAF's target names) when >= 50% of its sampled
      reads are placed and >= 80% of those on that strain, '-' otherwise -> OUT.tsv (barcode, sampled, placed, one count
      per strain, strain).
  pools.py strain READ_BC.tsv BARCODES.tsv STRAIN OUT.txt
      every read whose barcode is STRAIN's -> OUT.txt.
  pools.py cohort CANDIDATES.paf CONTIG START MIN_MAPQ OUT.txt
      per read the primary record (tp:A:P) with the highest chaining score s1; the read is in the cohort when that
      record has MAPQ >= MIN_MAPQ on CONTIG and starts at >= START (0-based target start) -> OUT.txt.
"""
import sys
from collections import Counter, defaultdict


def read_bc(path):
    with open(path) as f:
        for line in f:
            rid, bc, n = line.rstrip("\n").split("\t")
            yield rid, bc, int(n)


def sampled(path, n, min_bases):
    """[(read id, barcode)] of the first n reads per barcode with >= min_bases bases, in table order"""
    cnt, out = Counter(), []
    for rid, bc, length in read_bc(path):
        if length >= min_bases and cnt[bc] < n:
            cnt[bc] += 1
            out.append((rid, bc))
    return out


def write_ids(path, ids):
    with open(path + ".tmp", "w") as f:
        f.writelines(f"{r}\n" for r in ids)
    import os
    os.replace(path + ".tmp", path)


def cmd_sample(rbc, n, min_bases, out):
    write_ids(out, [r for r, _ in sampled(rbc, int(n), int(min_bases))])


def cmd_barcodes(rbc, paf, n, min_bases, out):
    picked = sampled(rbc, int(n), int(min_bases))
    bc_of = dict(picked)
    best = {}
    with open(paf) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if "tp:A:P" not in c[12:] or int(c[11]) < 1:
                continue
            span = int(c[3]) - int(c[2])
            if c[0] not in best or span > best[c[0]][0]:
                best[c[0]] = (span, c[5])
    strains = sorted({t for _, t in best.values()})
    count, placed = Counter(b for _, b in picked), defaultdict(Counter)
    for r, (_, t) in best.items():
        if r in bc_of:
            placed[bc_of[r]][t] += 1
    with open(out + ".tmp", "w") as f:
        f.write("\t".join(["barcode", "sampled", "placed"] + strains + ["strain"]) + "\n")
        for b in sorted(count):
            m = sum(placed[b].values())
            call = "-"
            for s in strains:
                if m >= 0.5 * count[b] and placed[b][s] >= 0.8 * m:
                    call = s
                    break
            f.write("\t".join([b, str(count[b]), str(m)] + [str(placed[b][s]) for s in strains] + [call]) + "\n")
    import os
    os.replace(out + ".tmp", out)
    sys.stdout.write(open(out).read())


def cmd_strain(rbc, barcodes, strain, out):
    with open(barcodes) as f:
        head = f.readline().rstrip("\n").split("\t")
        ib, js = head.index("barcode"), head.index("strain")
        mine = {c[ib] for c in (line.rstrip("\n").split("\t") for line in f) if c[js] == strain}
    if not mine:
        sys.exit(f"pools.py: no barcode of strain {strain} in {barcodes}")
    write_ids(out, sorted(r for r, b, _ in read_bc(rbc) if b in mine))


def cmd_cohort(paf, contig, start, min_mapq, out):
    best = {}
    with open(paf) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if "tp:A:P" not in c[12:]:
                continue
            s1 = next((int(x[5:]) for x in c[12:] if x.startswith("s1:i:")), 0)
            if c[0] not in best or s1 > best[c[0]][0]:
                best[c[0]] = (s1, c[5], int(c[7]), int(c[11]))
    keep = sorted(r for r, (_, t, t0, q) in best.items() if t == contig and q >= int(min_mapq) and t0 >= int(start))
    write_ids(out, keep)
    print(f"{len(keep):,} cohort reads of {len(best):,} placed", file=sys.stderr)


if __name__ == "__main__":
    cmds = {"sample": (cmd_sample, 4), "barcodes": (cmd_barcodes, 5), "strain": (cmd_strain, 4),
            "cohort": (cmd_cohort, 5)}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds or len(sys.argv) - 2 != cmds[sys.argv[1]][1]:
        sys.exit(__doc__)
    cmds[sys.argv[1]][0](*sys.argv[2:])
