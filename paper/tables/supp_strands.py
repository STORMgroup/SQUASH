#!/usr/bin/env python3
"""paper/tables/supp_strands.py: the strand of each SQUASH unitig on D1-D8 -> OUT_DIR/tables/supp_strands.tex.
The evaluation's spelled unitigs (OUT_DIR/eval/KEY/squash/asm/squash_spelled.fa) and its split by strand
(strands/squash_spelled_fwd.fa, _rev.fa: each unitig to the strand where more of it aligns); the unitigs on neither
strand are left out of the per-strand QUAST. Lengths are spelled bases."""
import lib

V = lib.Values("supp_strands")
COLS = ["unitigs", "forward unitigs", "forward Mb", "reverse unitigs", "reverse Mb", "neither unitigs", "neither Mb",
        "neither % of bases"]
lines, prev = [], None
for s in lib.main_sets():
    key, d = s["key"], s["D"]
    if prev is not None and s["bacterium"] != prev:
        lines[-1] += "\n  \\midrule"
    prev = s["bacterium"]
    files = [lib.squash_fasta(key, "asm", "NAME_spelled.fa"), lib.squash_fasta(key, "strands", "NAME_spelled_fwd.fa"),
             lib.squash_fasta(key, "strands", "NAME_spelled_rev.fa")]
    if all(files):
        allu, fw, rv = (lib.fasta_lengths(f) for f in files)
        ne = {u: x for u, x in allu.items() if u not in fw and u not in rv}
        tot = sum(allu.values())
        cells = [f"{len(allu):,}", f"{len(fw):,}", f"{sum(fw.values()) / 1e6:.2f}", f"{len(rv):,}",
                 f"{sum(rv.values()) / 1e6:.2f}", f"{len(ne):,}", f"{sum(ne.values()) / 1e6:.2f}",
                 f"{100 * sum(ne.values()) / tot:.2f}"]
    else:
        cells = [lib.NA] * len(COLS)
    src = "; ".join(lib.rel(f) for f in files if f)
    for c, v in zip(COLS, cells):
        V.add(key, "SQUASH", c, v, src)
    lines.append(f"  {d} {lib.short(d):<28} & " + " & ".join(cells) + " \\\\")
head = ("  & & \\multicolumn{2}{c}{\\textbf{Forward strand}} & \\multicolumn{2}{c}{\\textbf{Reverse strand}} & "
        "\\multicolumn{3}{c}{\\textbf{Neither strand}} \\\\\n"
        "  \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(l){7-9}\n"
        "  \\textbf{Dataset} & \\hd[r]{Unitigs} & \\hd[r]{Unitigs} & \\hd[r]{Mb} & \\hd[r]{Unitigs} & \\hd[r]{Mb} & "
        "\\hd[r]{Unitigs} & \\hd[r]{Mb} & \\hd[r]{\\% of\\\\bases} \\\\\n  \\midrule\n")
lib.write_table("supp_strands", (
    "% supp_strands.tex (paper/tables/supp_strands.py): the strand of each SQUASH unitig (its spelled sequence aligned to\n"
    "% the reference with minimap2, the strand with more aligned bases) and the unitigs aligning to neither strand, left\n"
    "% out of the per-strand QUAST. Spelled lengths (bases).\n" + lib.HD + "\\setlength\\tabcolsep{3pt}\n" +
    lib.tabular("@{}lrrrrrrrr@{}", head, lines)), V)
