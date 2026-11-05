#!/usr/bin/env python3
"""paper/tables/supp_runtime.py: the pipelines from raw signal to an assembly on D1-D8, bacteria (D1-D4) and
eukaryotes (D5-D8) side by side -> OUT_DIR/tables/supp_runtime.tex. Each cell: elapsed time, CPU time (user + system),
peak memory (GB = 10^6 kB of /usr/bin/time -v).
SQUASH: its whole run (OUT_DIR/KEY.time.txt). Rawsamble + miniasm: index + overlaps + miniasm (OUT_DIR/assembly/KEY/
rawsamble/logs/; "not run" where it has no run folder). The basecalling routes: dorado (OUT_DIR/timing/KEY/) + hifiasm
on the reads of the same model (OUT_DIR/assembly/KEY/hifiasm_hac/ for HAC, hifiasm/ for sup), elapsed times summed,
memory the larger peak; no CPU time for GPU basecalling (N/A). HAC (CPU): dorado on 64 cores, the set's time
estimated from the timed reads (lib.hac_cpu). The three CPU-only pipelines are ranked (lowest shaded) once all three
are known."""
import lib

V = lib.Values("supp_runtime")
ROWS = ["SQUASH", lib.RAWSAMBLE, "HAC (CPU) + hifiasm", "HAC (GPU) + hifiasm", "sup (GPU) + hifiasm"]
COLS = ["elapsed (h:mm:ss)", "CPU time (h)", "peak memory (GB)"]


def cells(key):
    r, src = lib.routes(key)
    out = []
    for name in ROWS:
        if name not in r:
            out.append(None)
            continue
        x = r[name]
        out.append([lib.hms(x[0]), lib.cpu_h(x[1]) if x[1] is not None else "N/A", lib.gb(x[2])] if x else
                   [lib.NA, "N/A" if "GPU" in name else lib.NA, lib.NA])
    up = out[:3]
    if all(c is not None for c in up) and lib.NA not in sum(up, []):
        for j in range(3):
            for i in lib.rank([c[j].replace(":", "") for c in up], low=True):
                up[i][j] = f"\\best{{{up[i][j]}}}"
    for name, row in zip(ROWS, out):
        for c, v in zip(COLS, row or [lib.NOT_RUN] * 3):
            V.add(key, name, c, v, src[name])
    return out


sets = lib.main_sets()
bact = [(s["D"], cells(s["key"])) for s in sets if s["bacterium"]]
euk = [(s["D"], cells(s["key"])) for s in sets if not s["bacterium"]]
assert len(bact) == len(euk), "the two panels need as many datasets each"
head = (lib.SPANNER.format(n=4, a=5, b=6, c=9) +
        f"  \\textbf{{Pipeline}} & \\textbf{{Dataset}} & {lib.COLS_TIME}\n    & \\textbf{{Dataset}} & {lib.COLS_TIME} \\\\\n"
        "  \\cmidrule(r){1-5}\\cmidrule{6-9}\n")
tools = ["\\sysname", "Rawsamble + miniasm", "HAC (CPU) + hifiasm$^\\dagger$", "HAC (GPU) + hifiasm", "sup (GPU) + hifiasm"]
lib.write_table("supp_runtime", (
    "% supp_runtime.tex (paper/tables/supp_runtime.py): the pipelines from raw signal to an assembly, bacteria and\n"
    "% eukaryotes side by side; 64 threads on two AMD EPYC 7502 CPUs, dorado on one NVIDIA RTX A5000 (GPU) or 64 cores.\n"
    "% SQUASH: its whole run. Basecalling routes: dorado + hifiasm on its reads, elapsed times summed, memory the larger\n"
    "% peak. HAC (CPU)$^\\dagger$: the set's time estimated from dorado's rate on the timed reads (a one-hour subset for\n"
    "% the eukaryotes); its CPU time = dorado's CPU share x that time + hifiasm's. Upper rows (CPU only) ranked once all\n"
    "% three are known, the lowest value shaded; lower rows not ranked.\n" + lib.HD +
    "\\setlength\\tabcolsep{3pt}  % set inside the body, so the float needs no change\n" +
    lib.tabular("@{}llrrr@{\\hspace{1em}}lrrr", head, lib.two_panel(bact, euk, tools, 2, 3))), V)
