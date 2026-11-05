#!/usr/bin/env python3
"""paper/tables/supp_athal_masked.py: D8 (A. thaliana) by QUAST against Col-CC with its organelles as is (the main
table's QUAST) and with Col-CC's tandem arrays masked as N -> OUT_DIR/tables/supp_athal_masked.tex.
OUT_DIR/eval/KEY/{squash,compare}/quast/ and quast_masked/; SQUASH forward / reverse strand."""
import lib

V = lib.Values("supp_athal_masked")
key = next(s["key"] for s in lib.main_sets() if s["D"] == "D8")
head = ("  & \\multicolumn{3}{c}{\\textbf{Col-CC}} & \\multicolumn{3}{c}{\\textbf{Arrays masked}} \\\\\n"
        "  \\cmidrule(lr){2-4}\\cmidrule(l){5-7}\n"
        f"  \\textbf{{Tool}} & {lib.COLS_CORR} & {lib.COLS_CORR} \\\\\n  \\midrule\n")
lib.write_table("supp_athal_masked", (
    "% supp_athal_masked.tex (paper/tables/supp_athal_masked.py): D8 A. thaliana by QUAST against Col-CC + organelles as\n"
    "% is (= the main table) and with Col-CC's tandem arrays masked; SQUASH forward / reverse strand.\n" + lib.HD +
    "\\setlength\\tabcolsep{3pt}\n" + lib.tabular("@{}lrrr@{\\hspace{1em}}rrr@{}", head, lib.two_views(key, "quast_masked", V))), V)
