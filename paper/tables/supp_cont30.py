#!/usr/bin/env python3
"""paper/tables/supp_cont30.py: contiguity of the supplement's six sets at ~30x per strand (D1-D5, D7: the same
genomes as the main sets), in contiguity_both's format and from the same files -> OUT_DIR/tables/supp_cont30.tex.
The six sets are paired in two columns of three, in the order of their labels."""
import lib

V = lib.Values("supp_cont30")
blocks = [(s["D"], lib.contiguity_cells(s["key"], V)) for s in lib.supp_sets()]
half = (len(blocks) + 1) // 2
head = (f"  \\textbf{{Tool}} & \\textbf{{Dataset}} & {lib.COLS_CONT}\n    & \\textbf{{Dataset}} & {lib.COLS_CONT} \\\\\n"
        "  \\cmidrule(r){1-6}\\cmidrule{7-11}\n")
rows = lib.two_panel(blocks[:half], blocks[half:], ["\\sysname", "Rawsamble", "nanoMDBG", "Flye", "hifiasm"], 1, 4)
lib.write_table("supp_cont30", (
    "% supp_cont30.tex (paper/tables/supp_cont30.py): contiguity at ~30x per strand, in the main contiguity table's\n"
    "% format (contiguity_both.py).\n" + lib.HD + "\\setlength\\tabcolsep{3pt}\n" +
    lib.tabular("@{}llrrrr@{\\hspace{1em}}lrrrr", head, rows)), V)
