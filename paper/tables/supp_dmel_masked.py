#!/usr/bin/env python3
"""paper/tables/supp_dmel_masked.py: D7 (D. melanogaster) by QUAST against release 6 as is (the main table's QUAST)
and with its repeat-rich regions masked as N (the rDNA locus, chrY, the pericentromeric heterochromatin)
-> OUT_DIR/tables/supp_dmel_masked.tex. OUT_DIR/eval/KEY/{squash,compare}/quast/ and quast_masked/; SQUASH forward /
reverse strand."""
import lib

V = lib.Values("supp_dmel_masked")
key = next(s["key"] for s in lib.main_sets() if s["D"] == "D7")
head = ("  & \\multicolumn{3}{c}{\\textbf{Release 6}} & \\multicolumn{3}{c}{\\textbf{Repeat-rich regions masked}} \\\\\n"
        "  \\cmidrule(lr){2-4}\\cmidrule(l){5-7}\n"
        f"  \\textbf{{Tool}} & {lib.COLS_CORR} & {lib.COLS_CORR} \\\\\n  \\midrule\n")
lib.write_table("supp_dmel_masked", (
    "% supp_dmel_masked.tex (paper/tables/supp_dmel_masked.py): D7 D. melanogaster by QUAST against release 6 as is (=\n"
    "% the main table) and with its repeat-rich regions masked; SQUASH forward / reverse strand.\n" + lib.HD +
    "\\setlength\\tabcolsep{3pt}\n" + lib.tabular("@{}lrrr@{\\hspace{1em}}rrr@{}", head, lib.two_views(key, "quast_masked", V))), V)
