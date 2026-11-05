#!/usr/bin/env python3
"""paper/tables/supp_hac.py: hifiasm on the sup reads (the main tables' QUAST, OUT_DIR/eval/KEY/compare/quast/) and on
the HAC reads of the same POD5s (OUT_DIR/eval/KEY/compare/quast_hac/: its hifiasm column) of D1-D8
-> OUT_DIR/tables/supp_hac.tex."""
import lib

V = lib.Values("supp_hac")
lines, prev = [], None
for s in lib.main_sets():
    key, d = s["key"], s["D"]
    if prev is not None and s["bacterium"] != prev:
        lines[-1] += "\n  \\midrule"
    prev = s["bacterium"]
    cells = []
    for view in ("quast", "quast_hac"):
        p = lib.eval_path(key, "compare", view, "report.tsv")
        rep = lib.quast(p)
        col = next((c for c in rep or {} if c.startswith("hifiasm")), "hifiasm") if view == "quast_hac" else "hifiasm"
        row = lib.strand_cells(rep, [col])
        for m, c in zip(lib.METRICS, row):
            V.add(key, "hifiasm " + ("sup" if view == "quast" else "HAC"), m, c, lib.rel(p))
        cells += row
    lines.append(f"  {d} {lib.short(d):<28} & {' & '.join(cells)} \\\\")
head = ("  & \\multicolumn{3}{c}{\\textbf{sup reads}} & \\multicolumn{3}{c}{\\textbf{HAC reads}} \\\\\n"
        "  \\cmidrule(lr){2-4}\\cmidrule(l){5-7}\n"
        f"  \\textbf{{Dataset}} & {lib.COLS_CORR} & {lib.COLS_CORR} \\\\\n  \\midrule\n")
lib.write_table("supp_hac", (
    "% supp_hac.tex (paper/tables/supp_hac.py): hifiasm on the sup reads (the main tables' runs) and on the HAC reads of\n"
    "% the same POD5s, by QUAST as the main table.\n" + lib.HD + "\\setlength\\tabcolsep{3pt}\n" +
    lib.tabular("@{}lrrr@{\\hspace{1em}}rrr@{}", head, lines)), V)
