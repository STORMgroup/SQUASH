#!/usr/bin/env python3
"""paper/tables/contiguity_both.py: contiguity of D1-D8, bacteria (D1-D4) and eukaryotes (D5-D8) side by side
-> OUT_DIR/tables/contiguity_both.tex.
SQUASH: OUT_DIR/eval/KEY/squash/tier1.tsv (lengths signal-estimated; NG50 against 2G, both strands);
Rawsamble + miniasm, nanoMDBG, Flye, hifiasm: OUT_DIR/eval/KEY/compare/tier1.tsv (Rawsamble "not run" where it has no
run folder). Largest comp. = the longest contig where the output has no links. \\best: the better raw-signal
assembler of each dataset and column (fewer unitigs, larger lengths; an undefined NG50 loses)."""
import lib

V = lib.Values("contiguity_both")
sets = lib.main_sets()
bact = [(s["D"], lib.contiguity_cells(s["key"], V)) for s in sets if s["bacterium"]]
euk = [(s["D"], lib.contiguity_cells(s["key"], V)) for s in sets if not s["bacterium"]]
assert len(bact) == len(euk), "the two panels need as many datasets each"
head = (lib.SPANNER.format(n=5, a=6, b=7, c=11) +
        f"  \\textbf{{Tool}} & \\textbf{{Dataset}} & {lib.COLS_CONT}\n    & \\textbf{{Dataset}} & {lib.COLS_CONT} \\\\\n"
        "  \\cmidrule(r){1-6}\\cmidrule{7-11}\n")
rows = lib.two_panel(bact, euk, ["\\sysname", "Rawsamble", "nanoMDBG", "Flye", "hifiasm"], 1, 4)
lib.write_table("contiguity_both", (
    "% contiguity_both.tex (paper/tables/contiguity_both.py): contiguity, bacteria and eukaryotes side by side.\n"
    "% SQUASH and Rawsamble + miniasm: signal-estimated lengths, NG50 against twice the genome length (both strands);\n"
    "% nanoMDBG, Flye, hifiasm on sup reads. Largest comp. = the longest contig where the output has no links. \\best =\n"
    "% the better raw-signal value of each dataset and column (fewer unitigs, larger lengths; an undefined NG50 loses).\n"
    + lib.HD + "\\setlength\\tabcolsep{3pt}  % set inside the body, so the float needs no change\n" +
    lib.tabular("@{}llrrrr@{\\hspace{1em}}lrrrr", head, rows)), V)
