#!/usr/bin/env python3
"""paper/tables/supp_commands.py: the parameters of each tool on D1-D8 (input and output paths left out)
-> OUT_DIR/tables/supp_commands.tex. SQUASH's row: each set's preset (sets.tsv) at the paper's 64 threads; the other
rows: the commands of the paper's runs (Rawsamble's index step also takes RawHash2's R10.4.1 pore model, -p; hifiasm
ran without -l0 on D6, the one diploid genome: the caption's dagger)."""
import itertools
import lib

V = lib.Values("supp_commands")
sets = lib.main_sets()
n = len(sets)
squash = " & ".join(f"\\multicolumn{{{len(list(g))}}}{{c}}{{\\texttt{{-x {p} -t 64}}}}"
                    for p, g in itertools.groupby(s["preset"] for s in sets))
for s in sets:
    V.add(s["key"], "SQUASH", "preset", s["preset"], lib.rel(lib.SETS_TSV))
OTHERS = [("Rawsamble", "-x r10ava -t 64"), ("nanoMDBG", "asm --in-ont --threads 64"), ("Flye", "--nano-hq -t 64"),
          ("hifiasm$^\\dagger$", "--ont -l0 -t 64"), ("Dorado GPU (sup)", "basecaller sup"),
          ("Dorado GPU (HAC)", "basecaller hac"), ("Dorado CPU (HAC)", "basecaller --device cpu hac")]
rows = [f"  \\sysname & {squash} \\\\"] + [f"  {t} & \\multicolumn{{{n}}}{{c}}{{\\texttt{{{c}}}}} \\\\" for t, c in OTHERS]
head = "  \\textbf{Tool} & " + " & ".join(f"\\textbf{{{s['D']}}}" for s in sets) + " \\\\\n  \\midrule\n"
lib.write_table("supp_commands", (
    "% supp_commands.tex (paper/tables/supp_commands.py): each tool's parameters per dataset (input and output paths\n"
    "% left out). Rawsamble's index step also passes RawHash2's R10.4.1 9-mer pore model (-p).\n" +
    lib.tabular("@{}l" + "c" * n + "@{}", head, rows)), V)
