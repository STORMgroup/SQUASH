#!/usr/bin/env python3
"""paper/tables/supp_versions.py: the version and source of each tool of the evaluation -> OUT_DIR/tables/supp_versions.tex.
Versions and links from versions.tsv; ROWS names the tools shown (Dorado once for the timed basecalling and once for
the assemblers' reads, D7's older version in brackets: its 4-kHz runs need Dorado 0.9.2). SQUASH's link is the
paper's \\sysrelease."""
import csv
import lib

VERSIONS = lib.PAPER / "versions.tsv"
ROWS = [("\\sysname", "SQUASH", ["1.0.0"]), ("Rawsamble (RawHash2)", "Rawsamble (RawHash2)", ["2.5 (commit 92366af)"]),
        ("miniasm", "miniasm", ["0.3-r179"]), ("nanoMDBG (metaMDBG)", "nanoMDBG (metaMDBG)", ["1.2"]),
        ("Flye", "Flye", ["2.9.6-b1802"]), ("hifiasm", "hifiasm", ["0.25.0-r726"]),
        ("Dorado (timed basecalling)", "Dorado", ["2.1.2", "0.9.2"]),
        ("Dorado (sup reads of the assemblers)", "Dorado", ["1.4.0", "0.9.2"]),
        ("Remora", "Remora", ["3.3.0"]), ("minimap2", "minimap2", ["2.30-r1287"]), ("QUAST", "QUAST", ["5.2.0"])]
V = lib.Values("supp_versions")
tools = {}
with open(VERSIONS) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["kind"] == "tool":
            tools.setdefault(r["name"], {})[r["version"]] = r["source"]
rows = []
for label, name, vs in ROWS:
    assert name in tools and all(v in tools[name] for v in vs), f"{name} {vs}: not in {VERSIONS}"
    src = tools[name][vs[0]]
    link = "\\sysrelease" if name == "SQUASH" else f"\\url{{{src}}}"
    version = vs[0] + (f" (D7: {vs[1]})" if len(vs) > 1 else "")
    V.add(name, label, "version", version, lib.rel(VERSIONS))
    rows.append(f"  {label} & {version} & {link} \\\\")
lib.write_table("supp_versions", (
    "% supp_versions.tex (paper/tables/supp_versions.py): versions and source links (versions.tsv).\n" +
    lib.tabular("@{}lll@{}", "  \\textbf{Tool} & \\textbf{Version} & \\textbf{Link to the Source Code} \\\\\n  \\midrule\n",
                rows)), V)
