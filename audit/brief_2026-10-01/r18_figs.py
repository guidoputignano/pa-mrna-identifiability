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

    # inset: the slope identity -1/(K_m+m), dosed band marked (bottom-left of A)
    axi = ax[0].inset_axes([0.10, 0.10, 0.40, 0.30])
    axi.set_facecolor("white")
    mm = np.linspace(1, 600, 400)
    axi.plot(mm, -1.0 / (T.KM + mm), color=P.INK, lw=1.6)
    axi.axvspan(Pn["m_obs"].min(), Pn["m_obs"].max(), color=P.C1, alpha=0.22)
    axi.axhline(-1.0 / T.KM, color=P.MUTED, ls=":", lw=0.8)
    axi.text(0.04, 0.90, "slope −1/(K_m+m)", transform=axi.transAxes,
             fontsize=7.2, va="top")
    axi.text(Pn["m_obs"].max() + 20, -1 / T.KM + 0.0013, "dosed band:\nslope ≈ −1/100",
             fontsize=6.4, color=P.C1, va="center")
    axi.set_xlabel("m  (µM)", fontsize=7, labelpad=1); axi.tick_params(labelsize=6.5)
    axi.set_ylim(-0.011, 0); axi.set_xlim(0, 600)
    for s in ("top", "right"):
        axi.spines[s].set_visible(False)
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
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    th = np.linspace(2.8, 6.2, 400)
    # 3-6 d bridging band
    ax.axvspan(3.0, 6.0, color=P.C1, alpha=0.08, zorder=0,
               label="3–6 d bridging band (tissue-ratio work)")
    # regime shading split at the crossing
    ax.axvspan(2.8, CROSS, color=P.C2, alpha=0.08, zorder=0)
    # the trough curve for the 21-day step
    ax.plot(th, _trough(21, th), color=P.INK, lw=2.4, zorder=4,
            label="Δ_trough(21 d)")
    # faint reference: the 14-day step is always far above theta_crit
    ax.plot(th, _trough(14, th), color=P.MUTED, lw=1.1, ls="--", zorder=3,
            label="Δ_trough(14 d), for reference")
    # theta_crit horizontal
    ax.axhline(THC, color=P.C3, lw=1.8, zorder=4)
    ax.text(2.84, THC + 12, f"θ_crit = {THC:.0f} µM/d", color=P.C3, fontsize=8.5, va="bottom")
    # crossing
    ax.axvline(CROSS, color=P.INK, ls=":", lw=1.2, zorder=3)
    ax.text(CROSS, 690, f"crossing\n{CROSS:.3f} d", ha="center", fontsize=8, color=P.INK)
    # the two anchors on the curve
    for x0, lab, col, dy in [(3.5, "repo's assumed\n3.5 d (below)", P.C2, -70),
                             (3.6, "measured human\n3.6 d — Chandler &\nBallard 1985 (above)", P.C4, 60)]:
        y0 = _trough(21, x0)
        ax.plot(x0, y0, "o", color=col, ms=9, zorder=6,
                markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(lab, (x0, y0), textcoords="offset points",
                    xytext=(-96 if x0 < CROSS else 12, dy), fontsize=7.6, color=col,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.8))
    # regime labels
    ax.text(3.02, 560, "21-day step NEEDS θ\n(θ unmeasurable → undeliverable)",
            fontsize=8, color=P.C2, style="italic")
    ax.text(4.5, 560, "21-day step FREE-FOR-ALL\n(no θ needed)",
            fontsize=8, color=P.C3, style="italic")
    ax.set_xlim(2.8, 6.2); ax.set_ylim(0, 760)
    ax.set_xlabel("PCC protein half-life  (days)")
    ax.set_ylabel("trough therapeutic increment  Δ_trough  (µM/d)")
    ax.set_title("The decision turns on a protein half-life, and the measured value sits ~1.5 % above the crossing",
                 fontsize=9.5)
    ax.legend(loc="upper left", fontsize=7.8, framealpha=0.92)
    P.finish(ax)
    P.save(fig, "r18_knife")


if __name__ == "__main__":
    fig_pencil()
    fig_fibre()
    fig_knife()
