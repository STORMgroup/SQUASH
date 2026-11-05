#!/usr/bin/env python3
"""paper/tables/datasets.py: the datasets table (D1-D8) -> OUT_DIR/tables/datasets.tex.
Genome length G, reads and signal-estimated bases (samples x 400 / sample rate, summed over the set's reads) from
sets.tsv; average read length = bases / reads, coverage = bases / G (both strands)."""
import lib

V = lib.Values("datasets")
if "signal_bases" not in lib.main_sets()[0]:
    raise SystemExit(f"datasets.py: {lib.SETS_TSV} has no column signal_bases (each set's signal-estimated bases)")
lines, prev = [], None
for s in lib.main_sets():
    key, d = s["key"], s["D"]
    g, reads, bases = int(s["G"]), int(s["reads"]), int(s["signal_bases"])
    if prev is not None and s["bacterium"] != prev:
        lines.append("  \\midrule")
    prev = s["bacterium"]
    cells = [f"{g / 1e6:.2f}", f"{reads:,}", f"{bases / reads / 1e3:.1f}", f"{bases / g:.1f}"]
    for c, v in zip(["genome length (Mb)", "reads", "average read length (kb)", "coverage (x)"], cells):
        V.add(key, d, c, v, lib.rel(lib.SETS_TSV))
    lines.append(f"  {d} {lib.NAMES[d][1]} & " + " & ".join(cells) + " \\\\")
lib.write_table("datasets", (
    "% datasets.tex (paper/tables/datasets.py): genome length G = the reference the assemblies are scored against;\n"
    "% reads and signal-estimated bases of each set (sets.tsv); coverage counts both strands.\n"
    "\\providecommand{\\hd}[2][l]{\\textbf{\\begin{tabular}[b]{@{}#1@{}}#2\\end{tabular}}}  % two-line header cell\n" +
    lib.tabular("@{}lrrrr@{}",
                "  \\textbf{Dataset} & \\hd[r]{Genome\\\\length (Mb)} & \\hd[r]{Reads\\\\(\\#)} & "
                "\\hd[r]{Average read\\\\length (kb)} & \\hd[r]{Coverage\\\\($\\times$)} \\\\\n  \\midrule\n", lines)), V)
