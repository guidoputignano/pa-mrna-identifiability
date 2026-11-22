"""
Overview figure for the manuscript: what was done and what was found, at a glance.
Schematic only (no data); every statement on it is a result stated in the text.
Writes manuscript/fig/fig_overview.pdf and .png.  Run: python -W error fig_overview.py
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle, Polygon, Ellipse, FancyArrowPatch, Wedge

BLUE, ORANGE, GREEN, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, GRID, RED = "#1a1a1a", "#5a5a5a", "#d9d8d3", "#c0392b"
LIGHT = "#f4f3ef"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5})
OUT = Path(__file__).resolve().parent / "fig"
OUT.mkdir(exist_ok=True)


def panel(ax, letter, title):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.1, 0.1), 9.8, 9.8, boxstyle="round,pad=0,rounding_size=0.5",
                                fc=LIGHT, ec="none", zorder=0))
    ax.text(0.5, 9.35, letter, fontsize=8.5, fontweight="bold", color=INK, va="top")
    ax.text(1.25, 9.33, title, fontsize=7.2, fontweight="bold", color=INK, va="top")


def arrow(ax, x0, y0, x1, y1, color=INK, lw=1.0, ms=7, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=ms, color=color,
                                 lw=lw, shrinkA=0, shrinkB=0, zorder=3))


def person(ax, x, y, s=1.0):
    ax.add_patch(Circle((x, y + 1.55 * s), 0.42 * s, fc="#e9c9a8", ec=INK, lw=0.6, zorder=3))
    ax.add_patch(FancyBboxPatch((x - 0.62 * s, y - 0.2 * s), 1.24 * s, 1.35 * s,
                                boxstyle="round,pad=0,rounding_size=0.45", fc=BLUE, ec=INK, lw=0.6, zorder=3))


def iv_bag(ax, x, y):
    ax.plot([x, x], [y - 2.9, y + 0.9], color=MUTED, lw=0.9, zorder=2)            # pole
    ax.plot([x - 0.35, x + 0.35], [y + 0.9, y + 0.9], color=MUTED, lw=0.9, zorder=2)
    ax.add_patch(FancyBboxPatch((x - 0.5, y - 0.9), 1.0, 1.6, boxstyle="round,pad=0,rounding_size=0.3",
                                fc="#dbe9fb", ec=BLUE, lw=0.8, zorder=3))
    ax.add_patch(Rectangle((x - 0.5, y - 0.9), 1.0, 0.75, fc=BLUE, alpha=0.35, ec="none", zorder=3))
    ax.text(x, y + 0.05, "mRNA", ha="center", va="center", fontsize=4.8, color=BLUE, fontweight="bold", zorder=4)


def liver(ax, x, y, s=1.0):
    ax.add_patch(Ellipse((x + 0.25 * s, y), 2.1 * s, 1.05 * s, angle=-12, fc="#b5583a", ec="none", zorder=3))
    ax.add_patch(Ellipse((x - 0.75 * s, y + 0.12 * s), 1.1 * s, 0.75 * s, angle=15, fc="#b5583a", ec="none", zorder=3))


def tube(ax, x, y, h=0.9):
    ax.add_patch(FancyBboxPatch((x - 0.13, y), 0.26, h, boxstyle="round,pad=0,rounding_size=0.12",
                                fc="white", ec=INK, lw=0.5, zorder=3))
    ax.add_patch(FancyBboxPatch((x - 0.13, y), 0.26, h * 0.55, boxstyle="round,pad=0,rounding_size=0.12",
                                fc=RED, ec="none", zorder=4))


def drop(ax, x, y, s=0.32, color=BLUE):
    ax.add_patch(Circle((x, y), s, fc=color, ec="none", zorder=4))
    ax.add_patch(Polygon([[x - s * 0.93, y + s * 0.37], [x + s * 0.93, y + s * 0.37], [x, y + s * 2.1]],
                         closed=True, fc=color, ec="none", zorder=4))


def tick_icon(ax, x, y, ok=True, s=0.38):
    if ok is None:                                   # a model result, not a measurement
        ax.add_patch(Circle((x, y), s, fc="white", ec=MUTED, lw=0.9, zorder=4))
        ax.text(x, y, "M", ha="center", va="center", fontsize=5.4, color=MUTED, fontweight="bold", zorder=5)
        return
    ax.add_patch(Circle((x, y), s, fc=GREEN if ok else RED, ec="none", zorder=4))
    if ok:
        ax.plot([x - 0.18, x - 0.04, x + 0.2], [y, y - 0.15, y + 0.17], color="white", lw=1.1, zorder=5,
                solid_capstyle="round")
    else:
        ax.plot([x - 0.15, x + 0.15], [y - 0.15, y + 0.15], color="white", lw=1.1, zorder=5)
        ax.plot([x - 0.15, x + 0.15], [y + 0.15, y - 0.15], color="white", lw=1.1, zorder=5)


fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(2, 3, wspace=0.04, hspace=0.06, left=0.01, right=0.99, top=0.99, bottom=0.01)
ax = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]

# (a) the therapy and the data it generates ------------------------------------------------
a = ax[0][0]
panel(a, "a", "Therapy and routine data")
iv_bag(a, 2.0, 7.0)
a.plot([2.0, 2.6, 3.4], [6.1, 5.4, 5.6], color=BLUE, lw=0.8, zorder=2)      # drip line
person(a, 4.1, 5.0, 1.0)
liver(a, 7.3, 6.8, 1.0)
a.text(7.3, 7.75, "liver makes PCC", ha="center", fontsize=5.6, color=INK)
arrow(a, 4.9, 6.1, 6.1, 6.7, color=MUTED, lw=0.7, ms=6)
a.text(7.3, 5.75, "own capacity θ\n+ delivered enzyme Δ(t)", ha="center", va="top", fontsize=5.3, color=MUTED)
# timeline
a.plot([0.8, 9.3], [2.2, 2.2], color=INK, lw=0.7)
for k, xd in enumerate(np.linspace(1.3, 8.7, 6)):
    drop(a, xd, 2.75, s=0.22)
for xt in (2.1, 4.5, 7.0):
    tube(a, xt, 1.0, 0.8)
a.text(9.3, 3.55, "infusion every 2 weeks", fontsize=5.6, color=BLUE, ha="right")
a.text(9.3, 0.55, "blood draws, diet record", fontsize=5.6, color=RED, ha="right")

# (b) the model ---------------------------------------------------------------------------
b = ax[0][1]
panel(b, "b", "One balance, two exits")
# tank
b.add_patch(FancyBboxPatch((2.4, 2.7), 4.2, 4.3, boxstyle="round,pad=0,rounding_size=0.35", fc="white",
                           ec=INK, lw=0.8, zorder=2))
b.add_patch(Rectangle((2.45, 2.75), 4.1, 2.6, fc="#f6d9a8", ec="none", zorder=2))
b.text(4.5, 4.0, "metabolite\nlevel m", ha="center", va="center", fontsize=5.8, color=INK, zorder=3)
arrow(b, 4.5, 8.0, 4.5, 6.6, color=INK, lw=1.0)
b.text(4.85, 7.55, "load: protein, odd-chain\nfatty acids, gut bacteria", ha="left", va="center", fontsize=5.0)
# exits
arrow(b, 6.6, 3.3, 8.3, 3.3, color=GREEN, lw=2.2, ms=8)
b.text(8.35, 3.95, "enzyme\nθ + Δ(t)", ha="center", va="bottom", fontsize=5.6, color=GREEN, fontweight="bold")
arrow(b, 2.4, 3.3, 0.7, 3.3, color=ORANGE, lw=4.0, ms=9)
b.text(1.2, 4.0, "other route\n(carnitine,\nurine)", ha="center", va="bottom", fontsize=5.3, color=ORANGE,
       fontweight="bold")
# inset: three hypotheses for the other route, one mini-plot each
mm = np.linspace(0, 1, 100)
shapes = [("carnitine: rises,\nthen falls", 2.2 * mm * np.exp(-1.6 * mm) * 1.9, "-"),
          ("linear", 0.95 * mm, "--"),
          ("saturable", 1.25 * mm / (0.3 + mm), ":")]
for i, (lab, yy, ls) in enumerate(shapes):
    ins = b.inset_axes([0.08 + 0.29 * i, 0.10, 0.24, 0.12])
    ins.plot(mm, yy, color=ORANGE, lw=1.1, ls=ls)
    ins.set_xticks([]); ins.set_yticks([]); ins.set_ylim(0, 1.1)
    for sp in ins.spines.values():
        sp.set_color(MUTED); sp.set_linewidth(0.5)
    ins.set_facecolor("white")
    ins.text(0.5, -0.12, lab, fontsize=4.7, color=ORANGE, ha="center", va="top", transform=ins.transAxes)
b.text(5.0, 2.35, "three hypotheses for the other route", fontsize=4.8, color=MUTED, ha="center")

# (c) sampling: the infusions are the experiment -------------------------------------------
c = ax[0][2]
panel(c, "c", "The infusions are the experiment")
for cx, cy, lab, spread, col in ((2.75, 5.3, "draws every 14 days\n(on infusion days)", False, RED),
                                 (7.25, 5.3, "draws spread over\nthe dose cycle", True, GREEN)):
    c.add_patch(Circle((cx, cy), 1.75, fc="white", ec=INK, lw=0.7, zorder=2))
    drop(c, cx, cy + 1.95, s=0.16)
    c.text(cx - 0.3, cy + 2.05, "infusion", ha="right", va="center", fontsize=4.8, color=BLUE)
    angles = np.linspace(90, 90 - 330, 12) if spread else np.full(12, 90.0)
    for i, ang in enumerate(angles):
        r = 1.35 if spread else 1.35 - 0.0 * i
        tx, ty = cx + r * np.cos(np.radians(ang)), cy + r * np.sin(np.radians(ang))
        c.add_patch(Circle((tx, ty), 0.13, fc=RED, ec="white", lw=0.3, zorder=4))
    c.text(cx, cy - 2.3, lab, ha="center", va="top", fontsize=5.4, color=INK)
    tick_icon(c, cx + 1.45, cy + 1.45, ok=spread, s=0.33)
c.text(5.0, 8.3, "where blood draws land in the 14-day cycle", ha="center", fontsize=5.4, color=MUTED)
c.text(5.0, 0.75, "same point of the cycle every time: little information.\nBook every 14.5 days; record doses and intake.",
       ha="center", va="bottom", fontsize=5.2, color=INK)

# (d) capacity in severe PA -----------------------------------------------------------------
d = ax[1][0]
panel(d, "d", "Capacity in severe PA")
# A + B bars
d.add_patch(Rectangle((1.0, 6.2), 8.0 * 0.84, 1.0, fc=ORANGE, ec="none", zorder=3))
d.add_patch(Rectangle((1.0 + 8.0 * 0.84, 6.2), 8.0 * 0.16, 1.0, fc=GREEN, ec="none", zorder=3))
d.text(1.2, 6.7, "other route 84%", va="center", fontsize=5.6, color="white", fontweight="bold", zorder=4)
d.text(9.0, 7.45, "enzyme 16%", ha="right", fontsize=5.6, color=GREEN, fontweight="bold")
d.text(5.0, 5.55, "a 1% error in the large share\nis a 5% error in the small one", ha="center", va="top",
       fontsize=5.5, color=INK)
# profile sketch
x = np.linspace(0, 1, 200)
yprof = np.where(x < 0.72, 0.03, (x - 0.72) ** 2 * 60)
d.plot(1.0 + 7.6 * x, 1.6 + np.minimum(yprof, 2.3), color=INK, lw=1.1, zorder=3)
d.plot([1.0, 8.6], [1.6, 1.6], color=MUTED, lw=0.5)
d.axvline  # noqa
d.text(1.6, 1.85, "no lower bound", fontsize=5.3, color=MUTED)
d.text(8.6, 1.35, "candidate capacity →", fontsize=4.8, color=MUTED, ha="right", va="top")
d.text(8.0, 3.95, "upper bound", fontsize=5.3, color=MUTED, ha="right")
tick_icon(d, 8.9, 9.0, ok=False, s=0.33)
d.text(0.55, 1.6, "fit worse →", fontsize=4.6, color=MUTED, rotation=90, va="bottom")
d.text(5.0, 0.2, "only an upper bound, while the other\nroute is known only roughly", ha="center", va="bottom",
       fontsize=5.2, color=INK)

# (e) the therapeutic window ----------------------------------------------------------------
e = ax[1][1]
panel(e, "e", "The therapeutic window")
t = np.linspace(0, 42, 600)
lam = np.log(2) / 4.3
peaks = np.array([0, 14, 28])
enz = np.zeros_like(t)
for p0 in peaks:
    enz += np.where(t >= p0, np.exp(-lam * (t - p0)), 0.0)
enz = enz / enz.max()
X = 1.0 + 8.0 * t / 42.0
Y = 2.2 + 4.6 * enz
e.add_patch(Rectangle((1.0, 2.2 + 4.6 * 0.06), 8.0, 4.6 * 0.94, fc=GREEN, alpha=0.10, ec="none", zorder=1))
e.plot(X, Y, color=GREEN, lw=1.2, zorder=3)
e.plot([1.0, 9.0], [2.2 + 4.6, 2.2 + 4.6], color=RED, lw=0.8, ls="--")
e.plot([1.0, 9.0], [2.2 + 4.6 * 0.06, 2.2 + 4.6 * 0.06], color=INK, lw=0.8, ls="--")
e.text(9.0, 7.0, "ceiling: own + delivered never above normal", ha="right", va="bottom", fontsize=5.2, color=RED)
e.text(1.0, 2.25, "critical capacity", ha="left", va="top", fontsize=5.2, color=INK)
e.text(9.0, 1.75, "enzyme over three infusions", fontsize=5.2, color=GREEN, va="top", ha="right")
e.text(5.0, 0.35, "interval safe if the trough clears the line:\nneeds half-life and dose, not θ",
       ha="center", va="bottom", fontsize=5.3, color=INK)

# (f) what the trial can measure -------------------------------------------------------------
f = ax[1][2]
panel(f, "f", "What the record can identify")
rows = [(True, "the enzyme half-life,\nfrom dated infusions"),
        (True, "change within one patient, as a\ndifference, if the marker is stable"),
        (False, "absolute capacity\nin severe PA"),
        (None, "how late an infusion can be:\na model result, not a measurement")]
for i, (ok, lab) in enumerate(rows):
    yy = 7.8 - 1.65 * i
    tick_icon(f, 1.2, yy, ok=ok, s=0.36)
    f.text(1.9, yy, lab, va="center", fontsize=5.6, color=INK)
f.text(5.0, 0.35, "plus a short data specification", ha="center", fontsize=5.3, color=MUTED)

fig.savefig(OUT / "fig_overview.pdf")
fig.savefig(OUT / "fig_overview.png", dpi=220)
print("wrote", OUT / "fig_overview.pdf")
