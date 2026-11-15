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

  pi_t1.pdf      : t1's U-curve, correlations and profile, recomputed with t1's own
                   functions and drawn at print size
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
plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
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
    fig = plt.figure(figsize=(6.7, 2.75))
    fig.set_layout_engine("none")
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.2)
    ax = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    cax = fig.add_subplot(gs[2])
    for a, Q, ttl in [(ax[0], T.pencil(seed=0), "(a) the paper's design: 25 monthly draws"),
                      (ax[1], T.wide_pencil(seed=1), "(b) imagined: 25 loads, 0.45× to 3.2× the usual")]:
        _lines(a, Q)
        sig_d, sig_num = T.dist_sigma(Q)
        _, NW = T.fit_point(Q, weights=(Q["a"] ** 2 + Q["b"] ** 2) / sig_num ** 2)
        _ellipse(a, np.array([T.TH5, T.KE]), np.linalg.inv(NW), color=P.C2, lw=2.0,
                 label="2σ near-concurrence region")
        a.plot(T.TH5, T.KE, "o", color=P.C2, ms=8, zorder=6, markeredgecolor="white",
               markeredgewidth=1.0, label="true ($\\theta$, $k_e$)")
        a.set_xlim(*XLIM); a.set_ylim(*YLIM)
        a.set_title(ttl, loc="left")
        a.set_xlabel("$\\theta$  (µM/d)")
        P.finish(a)
    ax[0].set_ylabel("$k_e$  (1/d)")
    ax[1].legend(loc="upper right", framealpha=0.92)
    sm = ScalarMappable(norm=MNORM, cmap=CMAP); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("marker level $m_i$ (µM)")
    _save(fig, "pi_pencil")


# ------------------------------------------------------------------ routes
def fig_routes():
    rows = [
        ("monthly record alone", 3.045, "A", False, "3.05"),
        ("+ load-turnover tracer (10 %)", 3.041, "A", False, "3.04"),
        ("+ IV $^{13}$C$_3$ tracer for $k_e$ (5.4 %)", 0.454, "A", False, "0.45"),
        ("+ both tracers", 0.441, "A", False, "0.44"),
        ("unrecorded load and $k_e$ exact (idealised)", 0.291, "A", False, "0.29"),
        ("+ ex vivo PCC assay, scale factor known", 0.200, "A", True, "0.20"),
        ("+ both tracers and ex vivo assay, factor known", 0.182, "A", True, "0.18"),
        ("+ ex vivo PCC assay, scale factor unknown", 3.045, "A", False, "3.05"),
        ("ex vivo ratio to controls, read as in vivo θ", 0.654, "A", False, "≈ 0.65"),
        ("2 × 104 weekly draws, 5 % per-draw CV", 0.133, "B", True, "0.13"),
        ("2 × 104 weekly draws, 26.4 % per-draw CV", 0.475, "B", False, "0.475"),
    ]
    fig, ax = plt.subplots(figsize=(5.9, 3.9))
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
        ax.text(max(v * 1.08, TARGET * 1.08), yy, txt, va="center", ha="left",
                color=P.INK, zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.1, 6.5)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(TARGET, color=P.INK, ls="--", lw=1.1, zorder=2)
    gA = ys[[r[2] == "A" for r in rows]].min() - 0.8
    gB = ys[[r[2] == "B" for r in rows]].min() - 0.8
    for yy, txt in ((gA, r"absolute capacity, sd(ln $\hat\theta$)"),
                    (gB, r"within-patient change, sd(ln $\hat\theta_2$ − ln $\hat\theta_1$)")):
        ax.text(-0.01, yy, txt, transform=ax.get_yaxis_transform(), ha="right", color=P.INK,
                va="center", style="italic")
    ax.set_xlabel("sd on the log scale at 5 % of normal capacity (target 0.20, dashed)")
    ax.legend(handles=[Patch(color=P.C2, label="misses 0.20"),
                       Patch(color=P.C1, label="reaches 0.20"),
                       Patch(facecolor="white", edgecolor=P.INK, hatch="////",
                             label="only if an unmeasured constant is supplied")],
              loc="upper center", bbox_to_anchor=(0.30, -0.14), ncol=3, frameon=False)
    ax.grid(True, axis="x", color=P.GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _save(fig, "pi_routes")


# ------------------------------------------------------------------ margin
D14, TAU0, THC, CROSS = 6000.0, 14.0, 150.775, 3.5452
BOX = (3.29, 3.95)          # out/t53_crossing_shape.txt, corners of the paper's ranges (tab_C_shape)


def _trough(tau, thalf, D=D14):
    x = np.exp(-np.log(2.0) / thalf * tau)
    return D * (tau / TAU0) * x / (1.0 - x)


D_ADM = 2672.0             # paper, scope section: largest admissible increment at 5 % of normal


def fig_margin():
    fig, ax = plt.subplots(figsize=(5.7, 3.5))
    th = np.linspace(2.6, 6.0, 400)
    ax.axhline(0.0, color=P.INK, lw=1.1, zorder=3)
    ax.plot(th, _trough(21, th) - THC, color=P.C1, lw=2.4, zorder=4,
            label=r"$\mu(21) = \Delta_{\mathrm{tr}}(21) - \theta_{\mathrm{crit}}^{\mathrm{int}}$, $D = 6000$")
    ax.plot(th, 2850.0 * np.exp(-np.log(2.0) / th * 21) - THC, color=P.C1, lw=1.2, ls=":", zorder=4,
            label=r"$\mu(21)$, largest admissible 21-day infusion")
    ax.plot(th, _trough(14, th) - THC, color=P.MUTED, lw=1.1, ls="--", zorder=3,
            label=r"$\mu(14)$, $D = 6000$, zero at 2.62 d")
    # range of the 21-day crossing over the admissible carnitine constants
    ax.axvspan(*BOX, color=P.GRID, alpha=0.55, zorder=0, lw=0)
    ax.text(2.64, 640, "grey band: 21-day\ncrossing over the\npaper's carnitine\nranges, "
            f"{BOX[0]:.2f}–{BOX[1]:.2f} d", ha="left", va="top", fontsize=7, color=P.MUTED)
    ax.axvline(CROSS, color=P.INK, ls=":", lw=1.1, zorder=3)
    ax.annotate(f"crossing at nominal\nconstants, {CROSS:.3f} d", va="top", xy=(CROSS, -60), xytext=(2.64, -148),
                fontsize=7, color=P.INK, arrowprops=dict(arrowstyle="-", color=P.INK, lw=0.6))
    for x0, lab, col, off in [
            (3.5, "assumed in the design, 3.5 d, μ = −8 µM/d", P.C2, (4.45, -135)),
            (3.6, "human fibroblast, 3.6 d, μ = +10 µM/d", P.C4, (4.45, -75)),
            (4.3, "likely liver value, 4.3 d\n(plausible 3–6 d)", P.C3, (4.75, 235))]:
        y0 = _trough(21, x0) - THC
        ax.plot(x0, y0, "o", color=col, ms=6.5, zorder=6, markeredgecolor="white",
                markeredgewidth=1.0)
        ax.annotate(lab, (x0, y0), xytext=off, fontsize=7, va="top", color=P.INK,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.8))
    ax.set_xlim(2.6, 6.0); ax.set_ylim(-220, 660)
    ax.set_xlabel("PCC protein half-life (d)")
    ax.set_ylabel("trough safety margin μ (µM/d)")
    leg = ax.legend(loc="upper left", bbox_to_anchor=(0.33, 1.0), frameon=True, framealpha=1.0,
                    edgecolor="none")
    leg.set_zorder(10)
    P.finish(ax)
    _save(fig, "pi_margin")


# ------------------------------------------------------------------ t1
def fig_t1():
    import core as C
    import t1_identifiability as T1
    curve = T1.scan_curve()
    pmin, sdmin = T1.u_minimum()
    lo, hi = T1.band_edges()
    mech = T1.mechanism()
    pr = T1.profiles()
    crb = C.sd_theta(C.TH_5)[0]
    fig, ax = plt.subplots(1, 3, figsize=(6.9, 2.5))
    a = ax[0]
    a.fill_between(curve[:, 0], curve[:, 2], curve[:, 3], color=P.C1, alpha=0.18, lw=0)
    a.plot(curve[:, 0], curve[:, 1], color=P.C1, lw=1.6)
    a.axhline(0.20, color=P.INK, lw=0.9, ls="--")
    a.text(1.1, 0.225, "target 0.20", color=P.INK, fontsize=6.5)
    a.axvspan(100 * lo / C.TH0, 100 * hi / C.TH0, color=P.C3, alpha=0.10, lw=0)
    a.plot([5], [crb], "o", color=P.C2, ms=4.5, zorder=5)
    a.annotate(f"5%: {crb:.2f}", (5, crb), (9, 5.0), fontsize=6.5, color=P.INK,
               arrowprops=dict(arrowstyle="-", color=P.MUTED, lw=0.6))
    a.plot([pmin], [sdmin], "o", color=P.C1, ms=4.5, zorder=5)
    a.annotate(f"min {sdmin:.3f} at {pmin:.0f}%", (pmin, sdmin), (1.3, 0.128), fontsize=6.5,
               arrowprops=dict(arrowstyle="-", color=P.MUTED, lw=0.6))
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("capacity, % of normal"); a.set_ylabel(r"sd(ln $\hat\theta$)")
    a.set_title("(a) precision, 25 monthly draws", loc="left")
    P.finish(a)
    b = ax[1]
    pc = np.array([d["pct"] for d in mech])
    b.plot(pc, [d["c_ke"] for d in mech], "-o", color=P.C1, ms=2.5, lw=1.4)
    b.plot(pc, [d["c_km"] for d in mech], "-s", color=P.C2, ms=2.5, lw=1.4)
    b.plot(pc, [d["c_b"] for d in mech], "-^", color=P.C3, ms=2.5, lw=1.4)
    b.text(1.2, -0.86, r"with ln$k_e$", color=P.C1, fontsize=6.5)
    b.text(150, 0.80, r"with ln$K_m$", color=P.C2, fontsize=6.5, ha="right")
    b.text(1.2, 0.02, r"with $b$", color=P.C3, fontsize=6.5)
    b.axhline(-1, color=P.MUTED, lw=0.6, ls=":")
    b.set_xscale("log"); b.set_ylim(-1.05, 1.05)
    b.set_xlabel("capacity, % of normal"); b.set_ylabel(r"corr. with ln $\theta$")
    b.set_title(r"(b) correlation with ln $\theta$", loc="left")
    P.finish(b)
    c = ax[2]
    x = np.exp(pr["lt"]) / C.TH_5
    for k, cv in enumerate(pr["curves"]):
        c.plot(x, cv, color=P.C1, lw=0.5, alpha=0.35, label="24 records" if k == 0 else None)
    c.plot(x, pr["expected"], color=P.INK, lw=1.5, label="noiseless")
    c.plot(x, 0.5 * (np.log(x) / crb) ** 2, color=P.C2, lw=1.1, ls="--", label="Cramér–Rao")
    c.axhline(1.92, color=P.MUTED, lw=0.7, ls=":")
    c.text(0.006, 2.3, "95% cut (1.92)", fontsize=6.5, color=P.MUTED)
    c.axvline(np.median(pr["his"]) / C.TH_5, color=P.C3, lw=1.0)
    c.text(np.median(pr["his"]) / C.TH_5 * 0.7, 6.2,
           f"95% upper end\n{100*np.median(pr['his'])/C.TH0:.0f}% of normal", fontsize=6.5,
           color=P.INK, ha="right")
    c.set_xscale("log"); c.set_ylim(-0.3, 14)
    c.set_xlabel(r"$\theta$ / true $\theta$")
    c.set_ylabel(r"$-\Delta\log$ profile lik.")
    c.set_title("(c) profile, 5% of normal", loc="left")
    c.legend(fontsize=6.5, loc="upper left", frameon=False)
    P.finish(c)
    fig.tight_layout(w_pad=1.0)
    _save(fig, "pi_t1")


if __name__ == "__main__":
    fig_t1()
    fig_pencil()
    fig_routes()
    fig_margin()
