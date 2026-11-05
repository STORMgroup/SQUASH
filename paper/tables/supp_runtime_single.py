#!/usr/bin/env python3
"""paper/tables/supp_runtime_single.py: elapsed time and peak memory of each tool alone on D1-D8 (the runs that the pipelines
of supp_runtime add up) -> OUT_DIR/tables/supp_runtime_single.tex.
SQUASH: OUT_DIR/KEY.time.txt; Rawsamble + miniasm (index + overlaps + miniasm), nanoMDBG, Flye, hifiasm (sup reads):
their timed steps (OUT_DIR/assembly/KEY/TOOL/logs/); dorado: OUT_DIR/timing/KEY/ (HAC on the CPU: the set's time
estimated from the timed reads, the memory of that run). Memory = peak RSS of the host process (GPU memory not
counted)."""
import lib

V = lib.Values("supp_runtime_single")
TOOLS = [("SQUASH", "\\sysname"), ("rawsamble", "\\hd[c]{Rawsamble\\\\+ miniasm}"), ("nanomdbg", "nanoMDBG"),
         ("flye", "Flye"), ("hifiasm", "hifiasm"), ("sup_gpu", "\\hd[c]{sup\\\\(GPU)}"), ("hac_gpu", "\\hd[c]{HAC\\\\(GPU)}"),
         ("hac_cpu", "\\hd[c]{HAC\\\\(CPU)$^\\dagger$}")]


def alone(key, tool):
    """(timing or the text 'not run', source)"""
    if tool == "SQUASH":
        return lib.squash_time(key), lib.rel(lib.OUT / f"{key}.time.txt")
    if tool in lib.STEPS:
        if not lib.ran(key, tool):
            return lib.NOT_RUN, lib.rel(lib.OUT / "assembly" / key / tool)
        return lib.tool_time(key, tool), lib.tool_sources(key, tool)
    if tool == "hac_cpu":
        return lib.hac_cpu(key), lib.hac_cpu_sources(key)
    model, device = tool.split("_")
    return lib.basecall_time(key, model, device), lib.rel(lib.basecall_file(key, model, device))


tl, ml, prev = [], [], None
for s in lib.main_sets():
    key, d = s["key"], s["D"]
    if prev is not None and s["bacterium"] != prev:
        tl[-1] += "\n  \\cmidrule(lr){2-9}"
        ml[-1] += "\n  \\cmidrule(lr){2-9}"
    prev = s["bacterium"]
    tc, mc = [], []
    for tool, _ in TOOLS:
        x, src = alone(key, tool)
        t = x if isinstance(x, str) else (lib.hms(x[0]) if x else lib.NA)
        m = x if isinstance(x, str) else (lib.gb(x[2]) if x else lib.NA)
        tc.append(V.add(key, tool, "elapsed (h:mm:ss)", t, src))
        mc.append(V.add(key, tool, "peak memory (GB)", m, src))
    tl.append(f"  {d} & " + " & ".join(tc) + " \\\\")
    ml.append(f"  {d} & " + " & ".join(mc) + " \\\\")
head = ("  & \\multicolumn{5}{c}{\\textbf{Assemblers}} & \\multicolumn{3}{c}{\\textbf{Basecalling (dorado)}} \\\\\n"
        "  \\cmidrule(lr){2-6}\\cmidrule(l){7-9}\n"
        "  \\textbf{Dataset} & " + " & ".join(h for _, h in TOOLS) + " \\\\\n  \\midrule\n"
        "  & \\multicolumn{8}{c}{Elapsed time (h:mm:ss)} \\\\\n  \\cmidrule(lr){2-9}\n")
rows = tl + ["  \\midrule\n  & \\multicolumn{8}{c}{Peak memory (GB)} \\\\\n  \\cmidrule(lr){2-9}"] + ml
lib.write_table("supp_runtime_single", (
    "% supp_runtime_single.tex (paper/tables/supp_runtime_single.py): each tool alone on D1-D8 (the runs of the runtime table): elapsed\n"
    "% time and peak memory. SQUASH = its whole run; Rawsamble + miniasm = index + overlaps + miniasm; assemblers on sup\n"
    "% reads at 64 threads (AMD EPYC 7502); dorado on one NVIDIA RTX A5000 and on 64 EPYC 7502 cores (HAC CPU$^\\dagger$:\n"
    "% the set's time estimated from a subset of the large sets, memory of the timed run). Memory = peak RSS of the host\n"
    "% process (GPU memory not counted).\n" + lib.HD + "\\setlength\\tabcolsep{3pt}\n" +
    lib.tabular("@{}lrrrrrrrr@{}", head, rows)), V)
