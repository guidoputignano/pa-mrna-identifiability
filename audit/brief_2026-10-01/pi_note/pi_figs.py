"""
Figures for the technical note (pi_note.tex), drawn in the paper's vocabulary:
no in-figure titles (captions carry the reading), "local identification model"
and "parent model" rather than letters, and no kappa (the paper uses kappa for the
carnitine gain).

  pi_pencil.pdf  : the line pencil, dosed record vs one wide excursion (t51 data)
  pi_routes.pdf  : every route to theta at 5 % of normal, on one log scale
                   (values from out/t26, t27, t29, t22, t42, as in story_figs.py)
  pi_margin.pdf  : the 21-day safety margin against the PCC half-life, with the
                   range of the crossing over the parent model's admissible
                   carnitine constants (out/t53_crossing_shape.txt)

The t1 figure (fig/t1_identifiability.png) is used unchanged.
Run from this directory:  python -W error pi_figs.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import plotstyle as P  # noqa: E402
from plotstyle import plt  # noqa: E402
import t51_pencil_geometry as T  # noqa: E402

OUT = HERE / "fig"
OUT.mkdir(exist_ok=True)
TARGET = 0.20


def _save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig/{name}.pdf")


# ------------------------------------------------------------------ pencil
XLIM, YLIM = (0, 400), (0, 10)
MNORM = mpl.colors.Normalize(vmin=0, vmax=380)
CMAP = mpl.cm.viridis


def _lines(ax, Q):
    th = np.linspace(*XLIM, 200)
    for i in range(Q["a"].size):
        ke = (Q["c"][i] - Q["a"][i] * th) / Q["b"][i]
        ax.plot(th, ke, color=CMAP(MNORM(Q["m_obs"][i])), lw=1.0, alpha=0.85, zorder=2)


def _ellipse(ax, x, Cov, n=2.0, **kw):
    w, V = np.linalg.eigh(Cov)
    t = np.linspace(0, 2 * np.pi, 160)
    xy = (V @ (np.sqrt(np.maximum(w, 0))[:, None] * np.array([np.cos(t), np.sin(t)]))) * n
    ax.plot(x[0] + xy[0], x[1] + xy[1], zorder=5, **kw)


def fig_pencil():
    fig = plt.figure(figsize=(10.6, 4.2))
    fig.set_layout_engine("none")
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.2)
    ax = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    cax = fig.add_subplot(gs[2])
    for a, Q, ttl in [(ax[0], T.pencil(seed=0), "(a) dosed monthly record, 25 draws"),
                      (ax[1], T.wide_pencil(seed=1), "(b) record with one wide load excursion")]:
        _lines(a, Q)
        sig_d, sig_num = T.dist_sigma(Q)
        _, NW = T.fit_point(Q, weights=(Q["a"] ** 2 + Q["b"] ** 2) / sig_num ** 2)
        _ellipse(a, np.array([T.TH5, T.KE]), np.linalg.inv(NW), color=P.C2, lw=2.0,
                 label="2σ near-concurrence region")
        a.plot(T.TH5, T.KE, "o", color=P.C2, ms=8, zorder=6, markeredgecolor="white",
               markeredgewidth=1.0, label="true (θ, $k_e$)")
        a.set_xlim(*XLIM); a.set_ylim(*YLIM)
        a.set_title(ttl, fontsize=9.5, loc="left")
        a.set_xlabel("θ  (µM/d)")
        P.finish(a)
    ax[0].set_ylabel("$k_e$  (1/d)")
    ax[1].legend(loc="upper right", fontsize=8, framealpha=0.92)
    sm = ScalarMappable(norm=MNORM, cmap=CMAP); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("marker level $m_i$ (µM)", fontsize=8.5)
    _save(fig, "pi_pencil")


# ------------------------------------------------------------------ routes
def fig_routes():
    rows = [
        ("monthly record alone", 3.045, "A", False, "3.05"),
        ("+ load-turnover tracer (10 %)", 3.041, "A", False, "3.04"),
        ("+ IV $^{13}$C$_3$ tracer for $k_e$ (5.4 %)", 0.454, "A", False, "0.45"),
        ("+ both tracers", 0.441, "A", False, "0.44"),
        ("load and $k_e$ known exactly (idealised)", 0.291, "A", False, "0.29"),
        ("+ ex vivo PCC assay, scale factor known", 0.200, "A", True, "0.20"),
        ("+ all three, scale factor known", 0.182, "A", True, "0.18"),
        ("+ ex vivo PCC assay, scale factor unknown", 3.045, "A", False, "3.05"),
        ("ex vivo ratio to controls, read as in vivo θ", 0.654, "A", False, "≈ 0.65"),
        ("2 × 104 weekly draws, 5 % per-draw CV", 0.133, "B", False, "0.13"),
        ("2 × 104 weekly draws, 26.4 % per-draw CV", 0.475, "B", False, "0.475"),
    ]
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ys, y = [], 0
    for i, r in enumerate(rows):
        if i > 0 and r[2] != rows[i - 1][2]:
            y += 1.0
        ys.append(y)
        y += 1
    ys = np.array(ys)
    for (lab, v, grp, assumed, txt), yy in zip(rows, ys):
        col = P.C1 if v <= TARGET + 1e-9 else P.C2
        ax.barh(yy, v, left=0.0, height=0.62, color=col, edgecolor="white", linewidth=0.0,
                zorder=3, hatch="////" if assumed else None)
        ax.text(max(v * 1.08, TARGET * 1.08), yy, txt, va="center", ha="left", fontsize=8.6,
                color=P.INK, zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.6)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.1, 6.5)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(TARGET, color=P.INK, ls="--", lw=1.1, zorder=2)
    gA = ys[[r[2] == "A" for r in rows]].min() - 0.8
    gB = ys[[r[2] == "B" for r in rows]].min() - 0.8
    for yy, txt in ((gA, "absolute capacity, sd(ln θ̂)"),
                    (gB, "within-patient change, sd(ln θ̂₂ − ln θ̂₁)")):
        ax.text(0.102, yy, txt, fontsize=8.2, color=P.MUTED, va="center", style="italic")
    ax.set_xlabel("standard deviation on the log scale at 5 % of normal capacity  (target 0.20, dashed)")
    ax.legend(handles=[Patch(color=P.C2, label="misses 0.20"),
                       Patch(color=P.C1, label="reaches 0.20"),
                       Patch(facecolor="white", edgecolor=P.INK, hatch="////",
                             label="only with an unmeasured calibration supplied")],
              loc="upper center", bbox_to_anchor=(0.36, -0.12), ncol=3, fontsize=8,
              frameon=False)
    ax.grid(True, axis="x", color=P.GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _save(fig, "pi_routes")


# ------------------------------------------------------------------ margin
D14, TAU0, THC, CROSS = 6000.0, 14.0, 150.775, 3.5452
BOX = (3.294, 4.083)        # out/t53_crossing_shape.txt, eight corners of the admissible box


def _trough(tau, thalf):
    x = np.exp(-np.log(2.0) / thalf * tau)
    return D14 * (tau / TAU0) * x / (1.0 - x)


def fig_margin():
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    th = np.linspace(2.6, 6.0, 400)
    ax.axhline(0.0, color=P.INK, lw=1.1, zorder=3)
    ax.plot(th, _trough(21, th) - THC, color=P.C1, lw=2.4, zorder=4,
            label=r"$\mu(21) = \Delta_{\mathrm{trough}}(21) - \theta_c$")
    ax.plot(th, _trough(14, th) - THC, color=P.MUTED, lw=1.1, ls="--", zorder=3,
            label=r"$\mu(14)$, zero at 2.62 d")
    # range of the 21-day crossing over the admissible carnitine constants
    ax.axvspan(*BOX, color=P.GRID, alpha=0.55, zorder=0, lw=0)
    ax.text(np.mean(BOX), 640, "21-day crossing over the\nadmissible carnitine constants\n"
            f"{BOX[0]:.2f}–{BOX[1]:.2f} d", ha="center", va="top", fontsize=7.8, color=P.MUTED)
    ax.axvline(CROSS, color=P.INK, ls=":", lw=1.1, zorder=3)
    ax.annotate(f"crossing at nominal\nconstants, {CROSS:.3f} d", xy=(CROSS, -60), xytext=(2.66, -170),
                fontsize=7.8, color=P.INK, arrowprops=dict(arrowstyle="-", color=P.INK, lw=0.6))
    for x0, lab, col, off in [
            (3.5, "assumed in the design, 3.5 d\nμ = −8 µM/d", P.C2, (-92, 34)),
            (3.6, "measured, human fibroblast, 3.6 d\nμ = +10 µM/d", P.C4, (26, -46)),
            (4.3, "liver estimate, 4.3 d\n(plausible 3–6 d)", P.C3, (14, -50))]:
        y0 = _trough(21, x0) - THC
        ax.plot(x0, y0, "o", color=col, ms=7.5, zorder=6, markeredgecolor="white",
                markeredgewidth=1.0)
        ax.annotate(lab, (x0, y0), textcoords="offset points", xytext=off, fontsize=7.8,
                    color=P.INK, arrowprops=dict(arrowstyle="-", color=col, lw=0.8))
    ax.set_xlim(2.6, 6.0); ax.set_ylim(-220, 660)
    ax.set_xlabel("PCC protein half-life (d)")
    ax.set_ylabel("trough safety margin μ (µM/d)")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)
    P.finish(ax)
    _save(fig, "pi_margin")


if __name__ == "__main__":
    fig_pencil()
    fig_routes()
    fig_margin()
