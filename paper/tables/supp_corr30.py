#!/usr/bin/env python3
"""paper/tables/supp_corr30.py: QUAST of the supplement's six sets at ~30x per strand (D1-D5, D7), in
correctness_both's format and from the same files -> OUT_DIR/tables/supp_corr30.tex. -- = undefined in QUAST (an
NGA50 of a strand that covers less than half of the genome)."""
import lib

V = lib.Values("supp_corr30")
blocks = [(s["D"], lib.correctness_cells(s["key"], V)) for s in lib.supp_sets()]
half = (len(blocks) + 1) // 2
head = (f"  \\textbf{{Tool}} & \\textbf{{Dataset}} & {lib.COLS_CORR}\n    & \\textbf{{Dataset}} & {lib.COLS_CORR} \\\\\n"
        "  \\cmidrule(r){1-5}\\cmidrule{6-9}\n")
rows = lib.two_panel(blocks[:half], blocks[half:],
                     ["\\sysname{}, spelled", "\\sysname{}, placed", "nanoMDBG", "Flye", "hifiasm"], 2, 3)
lib.write_table("supp_corr30", (
    "% supp_corr30.tex (paper/tables/supp_corr30.py): QUAST at ~30x per strand, in the main correctness table's format\n"
    "% (SQUASH forward / reverse strand; \\best / \\besthalf against nanoMDBG as there).\n" + lib.HD + lib.BESTHALF +
    "\\setlength\\tabcolsep{2.2pt}\n" + lib.tabular("@{}llrrr@{\\hspace{0.75em}}lrrr", head, rows)), V)
