"""paper/tables/lib.py: what the table and figure scripts share: the settings (OUT_DIR, from config.sh), the sets
(sets.tsv), the readers of the run, timing and evaluation folders under OUT_DIR (README.md, "Inputs"), the number
formats, the LaTeX of the paper's two-panel tables and values.tsv (every number with its source file).

A value whose input is missing prints as ?? and its file is listed on stderr; a tool that has no run folder for a set
prints as "not run". Python 3.10 standard library."""
import csv
import os
import re
import sys
from pathlib import Path

PAPER = Path(os.environ.get("PAPER") or Path(__file__).resolve().parents[1])
OUT = Path(os.environ.get("OUT_DIR") or PAPER / "out")
SETS_TSV = Path(os.environ.get("SETS_TSV") or PAPER / "sets.tsv")
TABLES, FIGURES = OUT / "tables", OUT / "figures"
NA, NOT_RUN = "??", "not run"
# SQUASH's rows in the evaluation's files: its QUAST columns, tier-1 rows and sequence files go by either name
SQUASH = ("SQUASH", "raw-mdbg")   # SQUASH's label in the evaluation's files (the second: the paper's own files)
# each tool's timed steps under OUT_DIR/assembly/KEY/TOOL/logs/ (STEP.time.txt, /usr/bin/time -v)
STEPS = {"nanomdbg": ["nanomdbg"], "flye": ["flye"], "hifiasm": ["hifiasm"], "hifiasm_hac": ["hifiasm"],
         "rawsamble": ["rawsamble_index", "rawsamble_overlap", "miniasm"]}
COMPARATORS = ["nanoMDBG", "Flye", "hifiasm"]
RAWSAMBLE = "Rawsamble + miniasm"
MISSING = []


def rel(p):
    """a path as values.tsv names it: under OUT_DIR relative to it, under paper/ as paper/..."""
    p = Path(p)
    for base, pre in ((OUT, ""), (PAPER, "paper/")):
        try:
            return pre + str(p.relative_to(base))
        except ValueError:
            pass
    return str(p)


def need(p):
    """p if it exists, else None (and p is listed as missing)"""
    p = Path(p)
    if p.exists():
        return p
    MISSING.append(rel(p))
    return None


# ------------------------------------------------------------------------------------------------------------ the sets
NAMES = {   # D label -> (name in the result tables, name and citation in the datasets table)
    "D1": ("\\textit{E. coli}", "\\textit{E. coli} \\cite{hall_nanovarbench_2024}"),
    "D2": ("\\textit{C. jejuni}", "\\textit{C. jejuni} \\cite{hall_nanovarbench_2024}"),
    "D3": ("\\textit{V. parahaemolyticus}", "\\textit{V. parahaemolyticus} \\cite{hall_nanovarbench_2024}"),
    "D4": ("\\textit{S. dysgalactiae}", "\\textit{S. dysgalactiae} \\cite{schababerle_methylome_2025}"),
    "D5": ("\\textit{A. rabiei}", "\\textit{A. rabiei} \\cite{debler_boscalid_2025, shah_arabiei_2020}"),
    "D6": ("Human chr13q", "Human chromosome $13$ q-arm \\cite{ont_giab_2025, hansen_q100_2025}"),
    "D7": ("\\textit{D. melanogaster}", "\\textit{D. melanogaster} \\cite{ont_dmel_2023, hoskins_dmel_2015}"),
    "D8": ("\\textit{A. thaliana}", "\\textit{A. thaliana} \\cite{earlham_athal_2026, colcc_2023}")}
# names broken over two lines in the two-panel tables
TWO_LINES = {"D3": ("\\textit{V. parahae-}", "\\textit{molyticus}"), "D4": ("\\textit{S. dys-}", "\\textit{galactiae}"),
             "D7": ("\\textit{D. melano-}", "\\textit{gaster}")}


def sets():
    """sets.tsv's rows, each with 'D': its paper label, a supplement set the label of the main set of its genome"""
    with open(SETS_TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    label = {r["genome"]: r["paper"] for r in rows if r["paper"] != "supp"}
    for r in rows:
        r["D"] = label[r["genome"]]
        r["bacterium"] = r["preset"] == "small"
    return rows


def main_sets():
    return sorted((r for r in sets() if r["paper"] != "supp"), key=lambda r: int(r["D"][1:]))


def supp_sets():
    return sorted((r for r in sets() if r["paper"] == "supp"), key=lambda r: int(r["D"][1:]))


def short(d):
    return NAMES[d][0]


# ------------------------------------------------------------------------------------------------------------ readers
def time_v(p):
    """(wall s, CPU s = user + system, peak RSS kB) of a /usr/bin/time -v file of a finished step, else None"""
    p = need(p)
    if p is None:
        return None
    t = p.read_text()
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", t)
    ok = re.search(r"Exit status: (\d+)", t)
    if not m or not ok or ok.group(1) != "0":
        MISSING.append(rel(p) + " (not a finished run)")
        return None
    wall = sum(float(x) * 60 ** i for i, x in enumerate(reversed(m.group(1).split(":"))))
    cpu = sum(float(re.search(rf"{k} time \(seconds\): (\S+)", t).group(1)) for k in ("User", "System"))
    return wall, cpu, int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1))


def cpu_share(p):
    """'Percent of CPU this job got' of a /usr/bin/time -v file, as a fraction"""
    return float(re.search(r"Percent of CPU this job got: (\d+)%", Path(p).read_text()).group(1)) / 100


def total(parts):
    """the sum of several steps: wall and CPU summed, the larger peak; None unless every step is known"""
    if not parts or any(p is None for p in parts):
        return None
    return sum(p[0] for p in parts), sum(p[1] for p in parts), max(p[2] for p in parts)


def squash_time(key):
    """SQUASH's whole run on a set (OUT_DIR/KEY.time.txt, squash.sh's /usr/bin/time -v)"""
    return time_v(OUT / f"{key}.time.txt")


def ran(key, tool):
    """whether the tool has a run folder for the set (OUT_DIR/assembly/KEY/TOOL)"""
    return (OUT / "assembly" / key / tool).is_dir()


def tool_time(key, tool):
    """a tool's timed run on a set: its timed steps summed (staging not counted)"""
    d = OUT / "assembly" / key / tool / "logs"
    return total([time_v(d / f"{s}.time.txt") for s in STEPS[tool]])


def tool_sources(key, tool):
    return "; ".join(rel(OUT / "assembly" / key / tool / "logs" / f"{s}.time.txt") for s in STEPS[tool])


def basecall_file(key, model, device):
    return OUT / "timing" / key / f"basecall_{model}_{device}.time.txt"


def basecall_time(key, model, device):
    """dorado's timed basecalling of a set (OUT_DIR/timing/KEY/basecall_MODEL_DEVICE.time.txt)"""
    return time_v(basecall_file(key, model, device))


def hac_cpu(key):
    """HAC basecalling on the CPU: (wall s, CPU s, peak RSS kB) of the whole set. The wall time is set_est_s of
    OUT_DIR/timing/KEY/basecall_hac_cpu.tsv (the set's samples / the samples per second dorado reached on the timed
    reads: the whole set for the bacteria, a subset for the large sets); CPU time = that x the CPU share of the timed
    run, peak memory the timed run's (basecall_hac_cpu.time.txt)."""
    tsv, tv = need(OUT / "timing" / key / "basecall_hac_cpu.tsv"), basecall_file(key, "hac", "cpu")
    t = time_v(tv)
    if tsv is None or t is None:
        return None
    with open(tsv) as f:
        wall = float(next(csv.DictReader(f, delimiter="\t"))["set_est_s"])
    return wall, wall * cpu_share(tv), t[2]


def hac_cpu_sources(key):
    return f"{rel(OUT / 'timing' / key / 'basecall_hac_cpu.tsv')}; {rel(basecall_file(key, 'hac', 'cpu'))}"


def quast(p):
    """a QUAST report.tsv as {column: {metric: value}}, None if missing"""
    p = need(p)
    if p is None:
        return None
    rows = [line.split("\t") for line in p.read_text().splitlines()]
    return {c: {r[0]: r[1 + i] for r in rows[1:]} for i, c in enumerate(rows[0][1:])}


def squash_column(rep, kind, strand):
    """SQUASH's column of a QUAST report: kind spelled / placed, strand forward / reverse"""
    for name in SQUASH:
        c = f"{name} ({kind}; {strand}-strand unitigs)"
        if rep and c in rep:
            return c
    return f"{SQUASH[0]} ({kind}; {strand}-strand unitigs)"


def tier1(p):
    """a tier-1 table as {assembler: row}, {} if missing"""
    p = need(p)
    if p is None:
        return {}
    with open(p) as f:
        return {r["assembler"]: r for r in csv.DictReader(f, delimiter="\t")}


def squash_tier1(key):
    rows = tier1(OUT / "eval" / key / "squash" / "tier1.tsv")
    return next((rows[n] for n in SQUASH if n in rows), None)


def squash_fasta(key, sub, name):
    """one of SQUASH's sequence files of its evaluation (asm/squash_spelled.fa, strands/squash_spelled_fwd.fa, ...)"""
    d = OUT / "eval" / key / "squash" / sub
    for n in SQUASH:
        p = d / name.replace("NAME", n.lower())
        if p.is_file():
            return p
    return need(d / name.replace("NAME", SQUASH[1]))


def fasta_lengths(p):
    """{name: length} of a FASTA file"""
    d, n = {}, None
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                n = line[1:].split()[0]
                d[n] = 0
            elif n is not None:
                d[n] += len(line.strip())
    return d


def eval_path(key, who, *parts):
    """OUT_DIR/eval/KEY/squash/... (who = squash) or OUT_DIR/eval/KEY/compare/... (who = compare)"""
    return OUT.joinpath("eval", key, who, *parts)


# ------------------------------------------------------------------------------------------------------------ formats
def mb(x, d=2):
    return f"{float(x) / 1e6:.{d}f}"


def hms(s):
    s = round(s)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def gb(kb):
    return f"{kb / 1e6:.1f}"


def cpu_h(s):
    return f"{s / 3600:,.1f}"


def rank(vals, low):
    """indices of the best displayed value (ties all) among those that are numbers; [] if fewer than two"""
    nums = [(i, float(v.replace(",", ""))) for i, v in enumerate(vals) if v not in (NA, "--", "")]
    if len(nums) < 2:
        return []
    b = (min if low else max)(x for _, x in nums)
    return [i for i, x in nums if x == b]


# ------------------------------------------------------------------------------------------------------------ values.tsv
class Values:
    """the numbers of one table or figure, each with its source file: written to OUT_DIR/tables/values/NAME.tsv,
    then OUT_DIR/tables/values.tsv is remade from every such file"""

    def __init__(self, name):
        self.name, self.rows = name, []

    def add(self, key, row, column, value, source):
        """one number (as printed, without its shading) and its source; returns the value as given"""
        self.rows.append((self.name, key, row, column, re.sub(r"\\best(?:half)?\{([^{}]*)\}", r"\1", str(value)), source))
        return value

    def write(self):
        d = TABLES / "values"
        d.mkdir(parents=True, exist_ok=True)
        head = "table\tkey\trow\tcolumn\tvalue\tsource\n"
        (d / f"{self.name}.tsv").write_text(head + "".join("\t".join(map(str, r)) + "\n" for r in self.rows))
        body = "".join(line for f in sorted(d.glob("*.tsv")) for line in f.read_text().splitlines(keepends=True)[1:])
        (TABLES / "values.tsv").write_text(head + body)


def write_table(name, text, values):
    """OUT_DIR/tables/NAME.tex and the table's values; the inputs that were missing go to stderr"""
    TABLES.mkdir(parents=True, exist_ok=True)
    (TABLES / f"{name}.tex").write_text(text)
    values.write()
    report(name, sum(1 for r in values.rows if NA in str(r[4])))


def report(name, n_missing):
    seen = list(dict.fromkeys(MISSING))
    print(f"{name}: written" + (f"; {n_missing} values ?? ({len(seen)} inputs missing)" if seen or n_missing else ""),
          file=sys.stderr)
    for p in seen[:20]:
        print(f"  missing: {p}", file=sys.stderr)
    if len(seen) > 20:
        print(f"  ... and {len(seen) - 20} more", file=sys.stderr)


# ------------------------------------------------------------------------------------------------------------ LaTeX
HD = ("\\providecommand{\\hd}[2][l]{\\textbf{\\begin{tabular}[b]{@{}#1@{}}#2\\end{tabular}}}  % two-line header cell\n"
      "\\providecommand{\\twoname}[2]{\\raisebox{0pt}[\\height][0pt]{\\begin{tabular}[t]{@{}l@{}}#1\\\\#2\\end{tabular}}}"
      "  % a two-line dataset name in one cell, no extra row height\n")
BESTHALF = ("\\providecolor{bestshade}{gray}{0.9}  % = gray!20; match \\best\n"
            "\\providecommand{\\besthalf}[1]{{\\setlength{\\fboxsep}{0pt}\\colorbox{bestshade}{\\kern1pt\\strut#1\\kern1pt}}}\n")
COLS_CONT = ("\\hd[r]{Unitig\\\\Count} & \\hd[r]{Longest\\\\unitig\\\\(Mbp)} & \\hd[r]{NG50\\\\(Mbp)} & "
             "\\hd[r]{Largest\\\\comp.\\\\(Mbp)}")
COLS_CORR = "\\hd[r]{NGA50\\\\(Mbp)} & \\hd[r]{Misasm.} & \\hd[r]{Coverage\\\\(\\%)}"
COLS_TIME = "\\hd[r]{Elapsed\\\\time\\\\(h:mm:ss)} & \\hd[r]{CPU\\\\time\\\\(h)} & \\hd[r]{Peak\\\\memory\\\\(GB)}"
SPANNER = ("  & \\multicolumn{{{n}}}{{c}}{{\\textbf{{Bacteria}}}} & \\multicolumn{{{n}}}{{c}}{{\\textbf{{Eukaryotes}}}} \\\\\n"
           "  \\cmidrule(lr){{2-{a}}}\\cmidrule(l){{{b}-{c}}}\n")


def labels(d):
    """a dataset's label cells down its five rows: its number, then its name (two lines in one cell if long)"""
    name = "\\twoname{%s}{%s}" % TWO_LINES[d] if d in TWO_LINES else short(d)
    return [d, name, "", "", ""]


def two_panel(left, right, tools, cmid_after, ncols, not_run=NOT_RUN):
    """the rows of a side-by-side table: left / right = [(D label, [cells of each tool's row or None])], paired in
    order; a row None prints as 'not run' across its cells"""
    out = []
    for n, ((dl, rl), (dr, rr)) in enumerate(zip(left, right)):
        if n:
            out.append(f"  \\cmidrule(r){{1-{ncols + 2}}}\\cmidrule{{{ncols + 3}-{2 * ncols + 3}}}")
        ll, lr = labels(dl), labels(dr)
        for i, tool in enumerate(tools):
            L = " & ".join(rl[i]) if rl[i] is not None else f"\\multicolumn{{{ncols}}}{{c}}{{{not_run}}}"
            R = " & ".join(rr[i]) if rr[i] is not None else f"\\multicolumn{{{ncols}}}{{c}}{{{not_run}}}"
            out.append(f"  {tool:<22} & {ll[i]:<26} & {L}\n    & {lr[i]:<26} & {R} \\\\")
            if i == cmid_after:
                out.append(f"  \\cmidrule(lr){{3-{ncols + 2}}}\\cmidrule(l){{{ncols + 4}-{2 * ncols + 3}}}")
    return out


def tabular(spec, head, rows):
    return f"\\begin{{tabular}}{{{spec}}}\n  \\toprule\n" + head + "\n".join(rows) + "\n  \\bottomrule\n\\end{tabular}\n"


# ------------------------------------------------------------------------------------------------------------ cells
def contiguity_cells(key, values):
    """the five rows (SQUASH, Rawsamble + miniasm, nanoMDBG, Flye, hifiasm) of a set's contiguity: unitigs, longest
    unitig, NG50 and largest component (the longest contig where the output has no links) from the tier-1 tables;
    \\best marks the better raw-signal assembler of each column (fewer unitigs, larger lengths; an undefined NG50
    loses). A Rawsamble row is None (not run) when it has no run folder."""
    sq_src, cmp_src = eval_path(key, "squash", "tier1.tsv"), eval_path(key, "compare", "tier1.tsv")
    raw = squash_tier1(key)
    cmp_rows = tier1(cmp_src)

    def cells(x):
        if not x:
            return [NA] * 4
        comp = x["largest_component_bp"] if x["largest_component_bp"] not in ("-", "") else x["largest_bp"]
        ng50 = "--" if str(x["NG50"]) in ("0", "0.0") else mb(x["NG50"])
        return [f"{int(x['segments']):,}", mb(x["largest_bp"]), ng50, mb(comp)]

    r = cells(raw)
    s = cells(cmp_rows.get(RAWSAMBLE)) if ran(key, "rawsamble") else None
    if s is not None:
        for j in range(4):
            if r[j] != "--" and s[j] != "--":
                best = rank([r[j], s[j]], low=(j == 0))
            else:
                best = [0] if s[j] == "--" and r[j] != "--" else [1] if r[j] == "--" and s[j] != "--" else []
            for i in best:
                (r, s)[i][j] = f"\\best{{{(r, s)[i][j]}}}"
    rows = [r, s] + [cells(cmp_rows.get(c)) for c in COMPARATORS]
    cols = ["unitigs", "longest unitig (Mbp)", "NG50 (Mbp)", "largest component (Mbp)"]
    for tool, row, src in zip(["SQUASH", RAWSAMBLE] + COMPARATORS, rows, [sq_src] + [cmp_src] * 4):
        for c, v in zip(cols, row or [NOT_RUN] * 4):
            values.add(key, tool, c, v, rel(src))
    return rows


def quast_cell(rep, col, metric):
    """a QUAST value as shown: NGA50 in Mbp, misassemblies with thousands separators, genome fraction to 0.1;
    -- where QUAST reports it undefined (-), ?? where missing"""
    v = rep.get(col, {}).get(metric) if rep else None
    if v is None or v == "":
        return NA
    if v == "-":
        return "--"
    return {"NGA50": lambda x: mb(x), "# misassemblies": lambda x: f"{int(x):,}",
            "Genome fraction (%)": lambda x: f"{float(x):.1f}"}[metric](v)


METRICS = ["NGA50", "# misassemblies", "Genome fraction (%)"]


def correctness_cells(key, values):
    """the five rows (SQUASH spelled, SQUASH placed, nanoMDBG, Flye, hifiasm) of a set's QUAST: NGA50, misassemblies,
    genome fraction; SQUASH forward / reverse strand. \\best / \\besthalf: SQUASH spelled against nanoMDBG, strand by
    strand; nanoMDBG where at least as good as both strands; placed, Flye and hifiasm not ranked."""
    sq_src, cmp_src = eval_path(key, "squash", "quast", "report.tsv"), eval_path(key, "compare", "quast", "report.tsv")
    rr, cr = quast(sq_src), quast(cmp_src)
    SF, SR = squash_column(rr, "spelled", "forward"), squash_column(rr, "spelled", "reverse")
    PF, PR = squash_column(rr, "placed", "forward"), squash_column(rr, "placed", "reverse")
    sp, pl, na = [], [], []
    for m in METRICS:
        hi = m != "# misassemblies"
        a, b, n = quast_cell(rr, SF, m), quast_cell(rr, SR, m), quast_cell(cr, "nanoMDBG", m)
        pl.append(f"{quast_cell(rr, PF, m)}\\,/\\,{quast_cell(rr, PR, m)}")
        for tool, v, src in (("SQUASH spelled forward", a, sq_src), ("SQUASH spelled reverse", b, sq_src),
                             ("SQUASH placed forward", quast_cell(rr, PF, m), sq_src),
                             ("SQUASH placed reverse", quast_cell(rr, PR, m), sq_src), ("nanoMDBG", n, cmp_src)):
            values.add(key, tool, m, v, rel(src))
        if any(v in (NA, "--") for v in (a, b, n)):
            sp.append(f"{a}\\,/\\,{b}")
            na.append(n)
            continue
        fa, fb, fn = (float(v.replace(",", "")) for v in (a, b, n))
        ok = (lambda v: v >= fn) if hi else (lambda v: v <= fn)
        if ok(fa) and ok(fb):
            sp.append(f"\\best{{{a}\\,/\\,{b}}}")
        elif ok(fa):
            sp.append(f"\\besthalf{{{a}}}\\,/\\,{b}")
        elif ok(fb):
            sp.append(f"{a}\\,/\\,\\besthalf{{{b}}}")
        else:
            sp.append(f"{a}\\,/\\,{b}")
        na.append(f"\\best{{{n}}}" if ((fn >= fa and fn >= fb) if hi else (fn <= fa and fn <= fb)) else n)
    oth = []
    for c in ("Flye", "hifiasm"):
        oth.append([values.add(key, c, m, quast_cell(cr, c, m), rel(cmp_src)) for m in METRICS])
    return [sp, pl, na] + oth


def strand_cells(rep, cols):
    """a QUAST row of one or two columns (SQUASH forward / reverse strand), one cell per metric"""
    return ["\\,/\\,".join(quast_cell(rep, c, m) for c in cols) for m in METRICS]


def routes(key):
    """the five pipelines of a set from its raw signal: {name: (wall s, CPU s or None, peak kB) or None}; a route that
    basecalls adds dorado's time to hifiasm's and takes the larger peak. Rawsamble's is absent when it has no run."""
    sq = squash_time(key)
    hac = tool_time(key, "hifiasm_hac")
    sup = tool_time(key, "hifiasm")
    hc, hg, sg = hac_cpu(key), basecall_time(key, "hac", "gpu"), basecall_time(key, "sup", "gpu")
    r = {"SQUASH": sq,
         "HAC (CPU) + hifiasm": (hc[0] + hac[0], hc[1] + hac[1], max(hc[2], hac[2])) if hc and hac else None,
         "HAC (GPU) + hifiasm": (hg[0] + hac[0], None, max(hg[2], hac[2])) if hg and hac else None,
         "sup (GPU) + hifiasm": (sg[0] + sup[0], None, max(sg[2], sup[2])) if sg and sup else None}
    if ran(key, "rawsamble"):
        r[RAWSAMBLE] = tool_time(key, "rawsamble")
    src = {"SQUASH": rel(OUT / f"{key}.time.txt"), RAWSAMBLE: tool_sources(key, "rawsamble"),
           "HAC (CPU) + hifiasm": f"{hac_cpu_sources(key)}; {tool_sources(key, 'hifiasm_hac')}",
           "HAC (GPU) + hifiasm": f"{rel(basecall_file(key, 'hac', 'gpu'))}; {tool_sources(key, 'hifiasm_hac')}",
           "sup (GPU) + hifiasm": f"{rel(basecall_file(key, 'sup', 'gpu'))}; {tool_sources(key, 'hifiasm')}"}
    return r, src


def two_views(key, view_b, values):
    """a set's QUAST rows in two views side by side: the main QUAST (quast/) and another (view_b, e.g. quast_masked/)
    of SQUASH (spelled, placed; forward / reverse strand) and the comparators -> table rows"""
    rows, names = [], ["\\sysname{}, spelled", "\\sysname{}, placed"] + COMPARATORS
    reps = {(who, v): (quast(eval_path(key, who, v, "report.tsv")), eval_path(key, who, v, "report.tsv"))
            for who in ("squash", "compare") for v in ("quast", view_b)}
    for i, name in enumerate(names):
        who = "squash" if i < 2 else "compare"
        cells = []
        for v in ("quast", view_b):
            rep, src = reps[(who, v)]
            cols = [squash_column(rep, ("spelled", "placed")[i], s) for s in ("forward", "reverse")] if i < 2 else [name]
            row = strand_cells(rep, cols)
            for m, c in zip(METRICS, row):
                values.add(key, name.replace("\\sysname{}", "SQUASH"), f"{v} {m}", c, rel(src))
            cells += row
        rows.append(f"  {name:<22} & {' & '.join(cells)} \\\\" + ("\n  \\midrule" if i == 1 else ""))
    return rows
