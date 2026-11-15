"""
Round-18 figures (reworked). Three separate figures:

  r18_pencil : the explanatory figure -- panels A (narrow dosed record) and B
               (wide excursion) on SHARED, POSITIVE-QUADRANT axes, lines coloured
               by marker level m_i, with an inset of the slope identity
               -1/(K_m+m). The triangulation picture: each draw is a bearing; the
               narrow record takes every bearing from nearly one direction, so the
               fix smears along the line of sight; the wide record is a second
               bearing from a different angle.
  r18_fibre  : panel C on its own -- the pencil is a family indexed by K_m, and
               the concurrence point slides along Remark 2's fibre.
  r18_knife  : the CONCLUSION figure -- Delta_trough(21) against PCC protein
               half-life, with theta_crit as a horizontal line, the 3.5452 d
               crossing, the old assumed 3.5 d (below) and measured human 3.6 d
               (above), and the 3-6 d bridging band shaded.

Palette/rcParams from plotstyle.py; sequential colour is viridis (CVD-safe).
"""
import numpy as np

import plotstyle as P
from plotstyle import plt
import t51_pencil_geometry as T

import matplotlib as mpl
from matplotlib.cm import ScalarMappable

XLIM, YLIM = (0, 400), (0, 10)
MNORM = mpl.colors.Normalize(vmin=0, vmax=380)      # shared m_i colour scale
CMAP = mpl.cm.viridis


def _lines_by_m(ax, Q, lw=1.0, alpha=0.85):
    th = np.linspace(*XLIM, 200)
    m = Q["m_obs"]
    for i in range(Q["a"].size):
        ke = (Q["c"][i] - Q["a"][i] * th) / Q["b"][i]
        ax.plot(th, ke, color=CMAP(MNORM(m[i])), lw=lw, alpha=alpha, zorder=2)


def _ellipse(ax, x, Cov, color, n=2.0, **kw):
    w, V = np.linalg.eigh(Cov)
    t = np.linspace(0, 2 * np.pi, 160)
    xy = (V @ (np.sqrt(np.maximum(w, 0))[:, None] * np.array([np.cos(t), np.sin(t)]))) * n
    ax.plot(x[0] + xy[0], x[1] + xy[1], color=color, zorder=5, **kw)


def _wls(Q):
    sig_d, sig_num = T.dist_sigma(Q)
    w = (Q["a"] ** 2 + Q["b"] ** 2) / sig_num ** 2
    return T.fit_point(Q, weights=w)


# =============================================================== figure 1
def fig_pencil():
    fig = plt.figure(figsize=(11.4, 4.7))
    fig.set_layout_engine("none")       # the colourbar column is laid out by the gridspec
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.22)
    ax = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    cax = fig.add_subplot(gs[2])
    Pn = T.pencil(seed=0)
    Pw = T.wide_pencil(seed=1)

    for a, Q, ttl in [(ax[0], Pn, "A — narrow dosed record: bearings from one direction"),
                      (ax[1], Pw, "B — wide excursion: a second bearing")]:
        _lines_by_m(a, Q)
        _, NW = _wls(Q)
        _ellipse(a, np.array([T.TH5, T.KE]), np.linalg.inv(NW), P.C2, n=2.0,
                 lw=2.2, label="2σ near-concurrence region")
        a.plot(T.TH5, T.KE, "o", color=P.C2, ms=9, zorder=6,
               markeredgecolor="white", markeredgewidth=1.0, label="true (θ, k_e)")
        a.set_xlim(*XLIM); a.set_ylim(*YLIM)
        a.set_title(ttl, fontsize=9.5, loc="left")
        a.set_xlabel("θ  (enzyme capacity, V_max)")
        P.finish(a)
    ax[0].set_ylabel("k_e  (linear clearance rate)")
    ax[1].legend(loc="upper right", fontsize=8, framealpha=0.92)

    # shared colourbar for m_i, in its own column
    sm = ScalarMappable(norm=MNORM, cmap=CMAP); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("marker level  m_i  (µM)", fontsize=8.5)
    fig.savefig(P.FIG / "r18_pencil.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig/r18_pencil.png")


# =============================================================== figure 2
def fig_fibre():
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    base = T.wide_pencil(seed=1)
    m, L = base["m_true"], base["loads"]              # noise-free lines
    Kgrid = np.linspace(20, 600, 160)
    fib = np.array([T.fit_point(dict(a=m / (Kk + m), b=m, c=L), weights=None)[0]
                    for Kk in Kgrid])
    ax.plot(fib[:, 0], fib[:, 1], color=P.GRID, lw=3.2, zorder=1, label="fibre (K_m swept)")
    for Kk, col in zip([50.0, 100.0, 200.0], [P.C1, P.C2, P.C3]):
        x, _ = T.fit_point(dict(a=m / (Kk + m), b=m, c=L), weights=None)
        ax.plot(x[0], x[1], "o", color=col, ms=10, zorder=5,
                markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(f"K_m = {Kk:.0f}", (x[0], x[1]), textcoords="offset points",
                    xytext=(8, 6), fontsize=8.5, color=col)
    ax.plot(T.TH5, T.KE, "*", color=P.INK, ms=17, zorder=6,
            markeredgecolor="white", markeredgewidth=0.8, label="truth (K_m = 100)")
    ax.set_xlim(0, 500); ax.set_ylim(3.6, 5.6)
    ax.set_xlabel("θ  (enzyme capacity, V_max)")
    ax.set_ylabel("k_e  (linear clearance rate)")
    ax.set_title("The pencil is a family in K_m: the concurrence\npoint slides along the fibre the intercept leaves open",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.92)
    P.finish(ax)
    P.save(fig, "r18_fibre")


# =============================================================== figure 3
D14, TAU0, THC = 6000.0, 14.0, 150.775
CROSS = 3.5452


def _trough(tau, thalf):
    lam = np.log(2.0) / thalf
    x = np.exp(-lam * tau)
    return D14 * (tau / TAU0) * x / (1.0 - x)


def fig_knife():
    # The decision as a MARGIN that crosses zero: mu(tau) = Delta_trough(tau) - theta_crit.
    # A zero crossing reads more cleanly than a line intersection, and it keeps the
    # y-axis a real quantity (safety margin, flux units) so the only percentage on the
    # figure -- the 1.5% -- is unambiguously a HALF-LIFE (x-axis) difference.
    fig, ax = plt.subplots(figsize=(7.7, 5.0))
    th = np.linspace(3.0, 6.0, 400)
    ax.axvspan(3.0, CROSS, color=P.C2, alpha=0.09, zorder=0)
    ax.axvspan(CROSS, 6.0, color=P.C3, alpha=0.07, zorder=0)
    ax.axhline(0.0, color=P.INK, lw=1.3, zorder=3)
    ax.text(5.97, 10, "margin = 0  (θ_crit = 151 µM/d)", color=P.INK, fontsize=8,
            ha="right", va="bottom")
    ax.plot(th, _trough(21, th) - THC, color=P.C1, lw=2.6, zorder=4,
            label="μ(21 d) = Δ_trough(21) − θ_crit")
    ax.plot(th, _trough(14, th) - THC, color=P.MUTED, lw=1.1, ls="--", zorder=3,
            label="μ(14 d), always positive (14-day step never needs θ)")
    ax.axvline(CROSS, color=P.INK, ls=":", lw=1.2, zorder=3)
    for x0, lab, col, dy in [
            (3.5, f"repo's assumed 3.5 d\nmargin {_trough(21,3.5)-THC:+.0f} µM/d (below)", P.C2, -120),
            (3.6, f"measured human 3.6 d\n(Chandler & Ballard 1985)\nmargin {_trough(21,3.6)-THC:+.0f} µM/d (above)", P.C4, 120)]:
        y0 = _trough(21, x0) - THC
        ax.plot(x0, y0, "o", color=col, ms=9, zorder=6,
                markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(lab, (x0, y0), textcoords="offset points",
                    xytext=(16, dy), fontsize=7.6, color=col,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.8))
    # annotate the x-gap explicitly (crossing -> measured), the source of the 1.5%
    ax.annotate("", xy=(3.6, 235), xytext=(CROSS, 235),
                arrowprops=dict(arrowstyle="<->", color=P.INK, lw=1.1))
    ax.annotate("crossing 3.5452 d → measured 3.6 d:\n0.055 d = 1.5 % in HALF-LIFE",
                xy=(3.573, 235), xytext=(3.95, 330), fontsize=7.8, color=P.INK,
                arrowprops=dict(arrowstyle="-", color=P.INK, lw=0.7))
    ax.text(3.03, -150, "21-day step NEEDS θ\n(θ unmeasurable → undeliverable)",
            fontsize=8, color=P.C2, style="italic", va="top")
    ax.text(4.8, 560, "21-day step FREE-FOR-ALL\n(no θ needed)",
            fontsize=8, color=P.C3, style="italic")
    ax.set_xlim(3.0, 6.0); ax.set_ylim(-220, 760)
    ax.set_xlabel("PCC protein half-life  (days)")
    ax.set_ylabel("safety margin  μ = Δ_trough(21) − θ_crit  (µM/d)")
    ax.set_title("The 21-day step turns on the PCC protein half-life; the measured value is\n"
                 "1.5 % (in half-life) past the crossing, +10 µM/d (Model C at its repo constants, K_c = 100)",
                 fontsize=9.5)
    ax.legend(loc="upper left", fontsize=7.8, framealpha=0.92)
    P.finish(ax)
    P.save(fig, "r18_knife")


if __name__ == "__main__":
    fig_pencil()
    fig_fibre()
    fig_knife()
