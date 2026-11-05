#!/usr/bin/env python3
"""paper/tables/supp_presets.py: the values each SQUASH preset (-x small / large) sets -> OUT_DIR/tables/supp_presets.tex.
Read from docs/USAGE.md ("Every tunable value": one row per option, its small and large values). The values both
presets share are one row; each preset's row has the values that differ or apply to it only. Every option of USAGE.md
must be placed below (LINES) or left out on purpose (LEFT_OUT), so a new option stops the script."""
import re
import lib

USAGE = lib.PAPER.parent / "docs" / "USAGE.md"
LINES = {   # the options of each row, line by line, as the table shows them
    "shared": [["--evdet-win", "--evdet-thr", "--evdet-peak"], ["--event-merge", "-q", "-e", "-w", "--bases-per-second"],
               ["--low-complexity", "--seed-repeat-fraction"],
               ["--min-anchors", "--min-overlap", "--min-cover", "--coverage"], ["--stop-below", "--band"],
               ["--density", "--repeat-ppm", "--k-min", "--k-max"],
               ["--abundance-cutoffs", "--tip-max", "--superbubble-max"],
               ["--min-support", "--max-conflict", "--anchor-nodes"], ["--anchor-edge", "--tangle-dominance"],
               ["--join-paths", "--join-threshold", "--join-min-nodes", "--join-head"],
               ["--layout-max-gap", "--layout-own-gap"]],
    "small": [["--rounds"]],
    "large": [["--rounds", "--hops-from", "--changed-from"], ["--skip-from", "--skip-pad", "--qsym-check"]]}
LEFT_OUT = {"--max-path-rows"}   # a resource cap of the graph construction, no effect unless reached

V = lib.Values("supp_presets")
text = USAGE.read_text()
section = text[text.index("## Every tunable value"):text.index("## Fixed parts")]
value, names = {}, {}
for line in section.splitlines():
    m = re.match(r"\| (`-[^|]+) \| ([^|]+) \| ([^|]+) \|", line)
    if m:
        opts = re.findall(r"`(-[^`]+)`", m.group(1))
        for o in opts:
            value[o] = (m.group(2).strip(), m.group(3).strip())
            names[o] = opts[-1]
placed = [o for row in LINES.values() for line in row for o in line]
missing = {names[o] for o in value} - {names[o] for o in placed} - LEFT_OUT
assert not missing, f"options of USAGE.md not in the table: {sorted(missing)}"
for o in (o for line in LINES["shared"] for o in line):
    assert value[o][0] == value[o][1], f"{o}: the presets differ ({value[o]}), it cannot be in the shared row"


def cell(row, i):
    lines = [" ".join(f"{o} {value[o][i]}" for o in line) for line in LINES[row]]
    for o in (o for line in LINES[row] for o in line):
        V.add(row, o, ("small", "large")[i], value[o][i], lib.rel(USAGE))
    tt = [f"\\texttt{{{x}}}" for x in lines]
    return tt[0] if len(tt) == 1 else "\\stack{" + "\\\\".join(tt) + "}"


ranges = {}
for s in lib.main_sets():
    ranges.setdefault(s["preset"], []).append(s["D"])
SIZE = {"small": "$< 10$ Mb", "large": "$\\ge 10$ Mb"}   # USAGE.md: small under 10 Mb, large from 10 Mb
use = {p: f"\\stack{{Genomes {SIZE[p]}\\\\({d[0]}--{d[-1]})}}" for p, d in ranges.items()}
rows = [f"  \\texttt{{small}}, \\texttt{{large}} & {cell('shared', 0)} & Default parameters \\\\", "  \\midrule",
        f"  \\texttt{{small}} & {cell('small', 0)} & {use['small']} \\\\", "  \\midrule",
        f"  \\texttt{{large}} & {cell('large', 1)} & {use['large']} \\\\"]
lib.write_table("supp_presets", (
    "% supp_presets.tex (paper/tables/supp_presets.py): SQUASH's two presets (docs/USAGE.md); the shared values are one\n"
    "% row. Left out: --max-path-rows (a resource cap, no effect unless reached) and the run options (-o, -t, ...).\n"
    "\\providecommand{\\stack}[2][l]{\\begin{tabular}[t]{@{}#1@{}}#2\\end{tabular}}  % multi-line cell, as \\hd\n" +
    lib.tabular("@{}lll@{}", "  \\textbf{Preset} & \\textbf{Corresponding parameters} & \\textbf{Usage} \\\\\n  \\midrule\n",
                rows)), V)
