"""Supplementary figure: where the blood draws of each schedule fall in the 14-day dose cycle.

Each draw is placed on a circle by its phase in the dose cycle (days after the last infusion,
infusion at the top, time running clockwise) and on a spiral by its visit number (first visit
innermost), so a schedule that keeps hitting the same phase stacks its draws on one ray, and one
that sweeps the cycle winds round it. Draw times come from the analysis's own schedule functions
(audit/g_design_common.py), 25 draws per schedule as in main-text Figure 2b.

Run: python3 si/fig_phase_clock.py   -> si/fig/fig_S_phaseclock.pdf/.png
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "audit"))
from g_design_common import fixed, phase_spread, T  # noqa: E402

N = 25
SCHED = [("every 14 days\n(the dose interval)", fixed(N, 14.0)),
         ("every 14.5 days", fixed(N, 14.5)),
         ("monthly (30.4 days)", fixed(N, 30.4)),
         ("phase-spread", phase_spread(N, gap=30.4))]
INK, MUTED, GRID, RED, BLUE, GOOD = "#1f2933", "#5f6b7a", "#d3d9e0", "#c0392b", "#2f6fdf", "#2e9e5b"
plt.rcParams.update({"font.family": "Inter", "font.size": 7.5, "axes.edgecolor": GRID})

fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.55), subplot_kw={"projection": "polar"})
for ax, (name, t) in zip(axs, SCHED):
    t = np.asarray(t, dtype=float)
    phase = np.mod(t, T)
    theta = 2 * np.pi * phase / T
    r = 0.32 + 0.68 * np.arange(len(t)) / (len(t) - 1)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.plot(theta, r, color=GRID, lw=0.6, zorder=1)
    ax.scatter(theta, r, s=11, color=RED, edgecolor="white", linewidth=0.4, zorder=3)
    ax.scatter([theta[0]], [r[0]], s=16, color=INK, zorder=4)
    ax.set_ylim(0, 1.08)
    ax.set_yticks([])
    days = [0, 3.5, 7, 10.5]
    ax.set_xticks([2 * np.pi * d / T for d in days])
    ax.set_xticklabels(["infusion", "3.5 d", "7 d", "10.5 d"], color=MUTED, fontsize=6.3)
    ax.tick_params(axis="x", pad=-1)
    ax.grid(color=GRID, lw=0.4)
    ax.spines["polar"].set_color(GRID)
    ax.set_title(name, fontsize=7.5, color=INK, pad=13, fontweight=600)
    ax.text(0.5, -0.27, f"25 draws over {t[-1]:.0f} days", transform=ax.transAxes, ha="center",
            fontsize=6.3, color=MUTED)
fig.subplots_adjust(left=0.03, right=0.97, top=0.79, bottom=0.17, wspace=0.55)
out = HERE / "fig"
out.mkdir(exist_ok=True)
fig.savefig(out / "fig_S_phaseclock.pdf")
fig.savefig(out / "fig_S_phaseclock.png", dpi=220)
for name, t in SCHED:
    ph = np.mod(np.asarray(t, float), T)
    print(f"{name.replace(chr(10), ' '):38s} span {t[-1]:6.1f} d   phases(half-day bins) {len(np.unique(np.round(ph * 2) / 2)):2d}")
