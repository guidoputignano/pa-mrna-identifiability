"""
Round 19: the two plain-language charts for story_so_far.html.

  story_routes : every route to the patient's enzyme capacity on one log scale,
                 against the 0.20 target. Answers "was sensor fusion the key?".
                 Values: out/t26_multiobservable.txt (sensor panel), out/t27_exvivo_kappa.txt
                 (calibration unknown), out/t22_longitudinal.txt and out/t42_serial_cv.txt
                 (change within one patient).
  story_wait   : how long you have to wait for each answer, on a 0-60 month axis.
                 Values: out/t46_timing.txt (lambda precision, 15- and 21-day steps)
                 and out/t52_knife_resolution.txt (which side of the 21-day crossing).

Palette and rcParams from plotstyle.py (validated with the dataviz skill's validator).
"""
import numpy as np
from matplotlib.patches import Patch

import plotstyle as P
from plotstyle import plt

TARGET = 0.20


def fig_routes():
    rows = [
        # (label, value, group, calibration assumed?, printed value)
        ("Monthly blood-marker record alone", 3.045, "A", False, "3.05"),
        ("+ load-turnover sensor (10 %)", 3.041, "A", False, "3.04"),
        ("+ IV tracer that measures k_e (5.4 %)", 0.454, "A", False, "0.45"),
        ("+ both of those sensors together", 0.441, "A", False, "0.44"),
        ("both of those known exactly (not obtainable)", 0.291, "A", False, "0.29"),
        ("+ ex vivo assay, κ assumed (test-tube activity)", 0.200, "A", True, "0.20"),
        ("+ all three sensors, κ assumed (test-tube activity)", 0.182, "A", True, "0.18"),
        ("ex vivo assay, κ unknown", 3.045, "A", False, "3.05"),
        ("ex vivo vs healthy controls, as whole-body capacity", 0.654, "A", False, "≈0.65"),
        ("weekly windows, model best case (5 % draw CV)", 0.133, "B", False, "0.13"),
        ("weekly windows, at the one measured CV (26.4 %)", 0.475, "B", False, "0.475"),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    # leave a gap row between the groups for the group header
    ys, y = [], 0
    for i, r in enumerate(rows):
        if i > 0 and r[2] != rows[i - 1][2]:
            y += 1.0
        ys.append(y)
        y += 1
    ys = np.array(ys)
    for (lab, v, grp, assumed, txt), yy in zip(rows, ys):
        reach = v <= TARGET + 1e-9
        col = P.C1 if reach else P.C2
        ax.barh(yy, v, left=0.0, height=0.62, color=col, edgecolor="white",
                linewidth=0.0, zorder=3,
                hatch="////" if assumed else None)
        ax.text(max(v * 1.08, TARGET * 1.08), yy, txt, va="center", ha="left", fontsize=9,
                color=P.INK, zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.1, 6.5)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(TARGET, color=P.INK, ls="--", lw=1.2, zorder=2)
    ax.text(TARGET * 0.97, ys[0] - 0.95, "target 0.20\n(about ±20 %)", ha="right",
            va="bottom", fontsize=8.5, color=P.INK)
    # group headers
    gA = ys[[r[2] == "A" for r in rows]].min() - 0.75
    gB = ys[[r[2] == "B" for r in rows]].min() - 0.75
    for yy, txt in ((gA, "CAPACITY LEVEL, ONE SNAPSHOT"),
                    (gB, "CHANGE WITHIN ONE PATIENT (a different quantity)")):
        ax.text(0.102, yy, txt, fontsize=8.2, color=P.MUTED, va="center",
                fontweight="bold")
    ax.set_xlabel("uncertainty sd(ln θ) at 5 % of normal capacity, Model A   (log scale; lower is better)")
    ax.legend(handles=[Patch(color=P.C2, label="misses the target"),
                       Patch(color=P.C1, label="reaches the target"),
                       Patch(facecolor="white", edgecolor=P.INK, hatch="////",
                             label="only with an unmeasured calibration assumed")],
              loc="upper center", bbox_to_anchor=(0.42, -0.11), ncol=3,
              fontsize=8, frameon=False)
    ax.grid(True, axis="x", color=P.GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    P.save(fig, "story_routes")


def fig_wait():
    # (label, kind, payload, note)
    rows = [
        ("15-day interval,\nevery patient (Model C)", "dot", [(0, "")], "day 0: no waiting"),
        ("Know λ to ±25 %", "dot", [(6, "")], "6 months, at any draw frequency"),
        ("Know λ to ±5 %", "dot", [(12, "weekly"), (24, "every 2 wk"), (48, "monthly")], ""),
        ("Tell a 3.6-day half-life\nfrom a 7-day one", "dot", [(6, "")], "6 months (7σ, monthly draws)"),
        ("Which side of the 21-day line,\none patient, if liver ≈ 4.3 d", "dot", [(6, "weekly"), (12, "monthly")],
         "80 % chance of a 2σ answer"),
        ("Which side of the 21-day line,\none patient, if liver ≈ 3.6 d", "never", [],
         "never: 0.6σ even after 5 years of weekly draws"),
        ("21 days for one patient\nby estimating θ (3.5-d case)", "never", [],
         "never from the passive record: margin estimate stays ±~140 µM/d"),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 4.9))
    for i, (lab, kind, pts, note) in enumerate(rows):
        y = i
        ax.axhline(y, color=P.GRID, lw=0.8, zorder=1)
        if kind == "dot":
            for x, tag in pts:
                ax.plot(x, y, "o", ms=10, color=P.C1, markeredgecolor="white",
                        markeredgewidth=1.2, zorder=4)
                if tag:
                    ax.text(x, y - 0.32, tag, ha="center", va="bottom", fontsize=7.8,
                            color=P.INK)
            if note:
                x0 = max(x for x, _ in pts)
                ax.text(x0 + 2.0, y, note, va="center", fontsize=8.6, color=P.INK)
        else:
            ax.annotate("", xy=(63.5, y), xytext=(0, y),
                        arrowprops=dict(arrowstyle="-|>", color=P.C2, lw=2.0,
                                        mutation_scale=14), zorder=3)
            ax.text(1.0, y - 0.3, note, va="bottom", fontsize=8.6, color=P.INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.8)
    ax.invert_yaxis()
    ax.set_xlim(-2, 64)
    ax.set_ylim(len(rows) - 0.4, -0.7)
    ax.set_xticks([0, 6, 12, 24, 36, 48, 60])
    ax.set_xlabel("months of on-therapy blood draws (simulated records)")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    P.save(fig, "story_wait")


if __name__ == "__main__":
    fig_routes()
    fig_wait()
