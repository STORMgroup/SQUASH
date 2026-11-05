#!/usr/bin/env python3
"""paper/tables/fig_runtime.py: the runtime figure -> OUT_DIR/figures/runtime.pdf. Grouped bars of the five pipelines
of supp_runtime (the same runs and sums, lib.routes): (a) elapsed time, D1-D4 in minutes and D5-D8 in hours on two
axes side by side; (b) peak memory (GB) over D1-D8; the value above each bar. Needs matplotlib (requirements.txt);
font Liberation Sans where installed."""
import lib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

PIPES = ["SQUASH", lib.RAWSAMBLE, "HAC (CPU) + hifiasm", "HAC (GPU) + hifiasm", "sup (GPU) + hifiasm"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
COL = dict(zip(PIPES, ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")))
V = lib.Values("fig_runtime")
MAIN = lib.main_sets()
val = {}
for s in MAIN:
    r, src = lib.routes(s["key"])
    # (elapsed s, peak kB); "not run" for a tool without a run folder; None where an input is missing
    val[s["D"]] = {p: (r[p][0], r[p][2]) if r.get(p) else (lib.NOT_RUN if p not in r else None) for p in PIPES}
    for p in PIPES:
        x = val[s["D"]][p]
        V.add(s["key"], p, "elapsed (s)", f"{x[0]:.0f}" if isinstance(x, tuple) else x or lib.NA, src[p])
        V.add(s["key"], p, "peak memory (GB)", f"{x[1] / 1e6:.2f}" if isinstance(x, tuple) else x or lib.NA, src[p])


def name(s):
    """the dataset's name under its number, broken as in the tables (human's too, under its bar group)"""
    if s["D"] == "D6":
        return lib.short("D6").replace(" ", "\n")
    if s["D"] in lib.TWO_LINES:
        return "\n".join(x.replace("\\textit{", "").rstrip("}") for x in lib.TWO_LINES[s["D"]])
    return lib.short(s["D"]).replace("\\textit{", "").rstrip("}")


plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7, "pdf.fonttype": 42, "axes.linewidth": 0.6,
                     "mathtext.fontset": "custom", "mathtext.it": "Liberation Sans:italic", "mathtext.rm": "Liberation Sans",
                     "mathtext.cal": "Liberation Sans",
                     "mathtext.bf": "Liberation Sans:bold"})
W = 0.155                                                   # bar width (5 bars per group of width 0.8)
OFF = [(-2 + i) * W for i in range(5)]


def draw(ax, group, kind, unit, ymax, step):
    for gi, s in enumerate(group):
        for pi, p in enumerate(PIPES):
            x, v = gi + OFF[pi], val[s["D"]][p]
            if not isinstance(v, tuple):                    # not run, or ?? where an input is missing
                ax.text(x, ymax * 0.02, v or lib.NA, rotation=90, ha="center", va="bottom", fontsize=5.6, color=MUTED)
                continue
            tot = v[0] / unit if kind == "time" else v[1] / 1e6
            ax.bar(x, tot, W, color=COL[p], edgecolor=INK, linewidth=0.35, zorder=3)
            ax.text(x, tot + ymax * 0.012, f"{tot:.1f}" if tot < 100 else f"{tot:.0f}", rotation=90, ha="center",
                    va="bottom", fontsize=5.6, color=INK)
    ax.set_xticks(range(len(group)))
    ax.set_xticklabels([s["D"] for s in group])
    for gi, s in enumerate(group):                          # the dataset's name under its number (human upright)
        ax.text(gi, -0.1, name(s), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.6,
                linespacing=1.0, color=INK2, style="italic" if "\\textit" in lib.short(s["D"]) else "normal")
    ax.set_xlim(-0.55, len(group) - 0.45)
    ax.set_ylim(0, ymax)
    ax.set_yticks(list(range(0, int(ymax) + 1, step)))
    for y in ax.get_yticks()[1:]:
        ax.axhline(y, color=GRID, lw=0.5, zorder=1)
    ax.tick_params(axis="both", length=0, colors=INK2, pad=2.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)


BG, EG = [s for s in MAIN if s["bacterium"]], [s for s in MAIN if not s["bacterium"]]
fig = plt.figure(figsize=(6.5, 4.3))
gs = fig.add_gridspec(2, 2, wspace=0.16, hspace=0.42)
a1, a1e, a2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
timed = lambda group: [v for s in group for v in val[s["D"]].values() if isinstance(v, tuple)]  # noqa: E731
tmax_b = max((v[0] / 60 for v in timed(BG)), default=1)
tmax_e = max((v[0] / 3600 for v in timed(EG)), default=1)
mmax = max((v[1] / 1e6 for v in timed(MAIN)), default=1)
draw(a1, BG, "time", 60, tmax_b * 1.2, 20)                     # bacteria in minutes
draw(a1e, EG, "time", 3600, tmax_e * 1.2, 5)                   # eukaryotes in hours
draw(a2, MAIN, "mem", None, mmax * 1.2, 50)
a1.set_ylabel("Elapsed time (min)", color=INK2)
a1e.set_ylabel("Elapsed time (h)", color=INK2)
a2.set_ylabel("Peak memory (GB)", color=INK2)
for ax, letter in ((a1, "a"), (a2, "b")):
    ax.text(-0.13 if ax is a1 else -0.065, 1.02, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left",
            va="bottom", color=INK)
fig.legend(handles=[Patch(facecolor=COL[p], edgecolor=INK, linewidth=0.35, label=p) for p in PIPES], loc="upper center",
           ncol=5, frameon=False, bbox_to_anchor=(0.53, 1.0), handlelength=1.2, handletextpad=0.4, columnspacing=1.3,
           labelcolor=INK)
fig.subplots_adjust(left=0.075, right=0.995, top=0.92, bottom=0.1)
lib.FIGURES.mkdir(parents=True, exist_ok=True)
fig.savefig(lib.FIGURES / "runtime.pdf")
V.write()
lib.report("fig_runtime", sum(1 for r in V.rows if r[4] == lib.NA))
