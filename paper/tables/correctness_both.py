#!/usr/bin/env python3
"""paper/tables/correctness_both.py: QUAST of D1-D8, bacteria (D1-D4) and eukaryotes (D5-D8) side by side
-> OUT_DIR/tables/correctness_both.tex.
SQUASH spelled and placed, forward / reverse strand: OUT_DIR/eval/KEY/squash/quast/report.tsv; nanoMDBG, Flye, hifiasm
(sup reads): OUT_DIR/eval/KEY/compare/quast/report.tsv. \\best / \\besthalf: SQUASH spelled against nanoMDBG, strand
by strand; nanoMDBG where at least as good as both strands; placed, Flye and hifiasm not ranked."""
import lib

V = lib.Values("correctness_both")
sets = lib.main_sets()
bact = [(s["D"], lib.correctness_cells(s["key"], V)) for s in sets if s["bacterium"]]
euk = [(s["D"], lib.correctness_cells(s["key"], V)) for s in sets if not s["bacterium"]]
assert len(bact) == len(euk), "the two panels need as many datasets each"
head = (lib.SPANNER.format(n=4, a=5, b=6, c=9) +
        f"  \\textbf{{Tool}} & \\textbf{{Dataset}} & {lib.COLS_CORR}\n    & \\textbf{{Dataset}} & {lib.COLS_CORR} \\\\\n"
        "  \\cmidrule(r){1-5}\\cmidrule{6-9}\n")
rows = lib.two_panel(bact, euk, ["\\sysname{}, spelled", "\\sysname{}, placed", "nanoMDBG", "Flye", "hifiasm"], 2, 3)
lib.write_table("correctness_both", (
    "% correctness_both.tex (paper/tables/correctness_both.py): QUAST 5.2.0, bacteria and eukaryotes side by side.\n"
    "% Upper rows = SQUASH spelled / placed (forward / reverse strand) and nanoMDBG; lower rows = Flye and hifiasm on the\n"
    "% same sup reads (not ranked). \\best / \\besthalf: the spelled row against nanoMDBG strand by strand, nanoMDBG where\n"
    "% at least as good as both strands; placed not ranked.\n" + lib.HD + lib.BESTHALF +
    "\\setlength\\tabcolsep{2.2pt}  % set inside the body, so the float needs no change\n" +
    lib.tabular("@{}llrrr@{\\hspace{0.75em}}lrrr", head, rows)), V)
