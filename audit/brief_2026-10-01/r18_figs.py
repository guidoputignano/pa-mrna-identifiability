"""
Round-18 figure: the line pencil in the (theta, k_e) plane (the PI's formulation).

Three panels:
  A  the NARROW dosed record -- 25 near-parallel lines, a long near-concurrence
     smear along the fibre; the identifiability failure, drawn.
  B  the WIDE excursion record (Task-7 forcing) -- the slopes fan out and the
     pencil concurs crisply at the true point; width buys identifiability.
  C  the pencil is a FAMILY indexed by K_m (Remark 2's fibre): reading the same
     noise-free lines with different assumed K_m slides the concurrence point
     along the fibre.

Palette and rcParams from plotstyle.py (the repo's validated palette).
"""
import numpy as np

import plotstyle as P
from plotstyle import plt
import t51_pencil_geometry as T


def _lines(ax, Q, th_lo, th_hi, color, lw=0.9, alpha=0.7, label=None):
    th = np.linspace(th_lo, th_hi, 200)
    for i in range(Q["a"].size):
        ke = (Q["c"][i] - Q["a"][i] * th) / Q["b"][i]
        ax.plot(th, ke, color=color, lw=lw, alpha=alpha, zorder=2,
                label=label if i == 0 else None)


def _ellipse(ax, x, Cov, color, n=2.0, **kw):
    w, V = np.linalg.eigh(Cov)
    t = np.linspace(0, 2 * np.pi, 120)
    xy = (V @ (np.sqrt(np.maximum(w, 0))[:, None] * np.array([np.cos(t), np.sin(t)]))) * n
    ax.plot(x[0] + xy[0], x[1] + xy[1], color=color, zorder=5, **kw)


def fig_pencil():
    fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.3))

    # ---- A: narrow dosed record --------------------------------------
    Pn = T.pencil(seed=0)
    _lines(ax[0], Pn, -200, 900, P.MUTED, lw=0.8, alpha=0.55)
    sig_d, sig_num = T.dist_sigma(Pn)
    wls = (Pn["a"] ** 2 + Pn["b"] ** 2) / sig_num ** 2
    xW, NW = T.fit_point(Pn, weights=wls)
    _ellipse(ax[0], np.array([T.TH5, T.KE]), np.linalg.inv(NW), P.C1,
             n=2.0, lw=2.0, label="2σ concurrence region")
    ax[0].plot(T.TH5, T.KE, "o", color=P.C2, ms=9, zorder=6, label="true (θ, k_e)")
    ax[0].set_xlim(-200, 900); ax[0].set_ylim(-2, 12)
    ax[0].set_title("A  narrow dosed record\n25 near-parallel lines, smear along the fibre",
                    fontsize=10)
    ax[0].legend(loc="upper right", fontsize=8, framealpha=0.9)

    # ---- B: wide excursion record ------------------------------------
    Pw = T.wide_pencil(seed=1)
    _lines(ax[1], Pw, 0, 350, P.MUTED, lw=0.8, alpha=0.6)
    sig_d2, sig_num2 = T.dist_sigma(Pw)
    wls2 = (Pw["a"] ** 2 + Pw["b"] ** 2) / sig_num2 ** 2
    xW2, NW2 = T.fit_point(Pw, weights=wls2)
    _ellipse(ax[1], np.array([T.TH5, T.KE]), np.linalg.inv(NW2), P.C1,
             n=2.0, lw=2.0, label="2σ concurrence region")
    ax[1].plot(T.TH5, T.KE, "o", color=P.C2, ms=9, zorder=6, label="true (θ, k_e)")
    ax[1].set_xlim(0, 350); ax[1].set_ylim(0, 10)
    ax[1].set_title("B  wide excursion (drug forcing)\nslopes fan out, the pencil concurs",
                    fontsize=10)
    ax[1].legend(loc="upper right", fontsize=8, framealpha=0.9)

    # ---- C: the K_m family (Remark 2's fibre) ------------------------
    # read the SAME noise-free wide lines with different assumed K_m; the
    # concurrence point slides along the fibre.
    base = T.wide_pencil(seed=1)
    m = base["m_true"]; L = base["loads"]            # noise-free
    Kms = [50.0, 100.0, 200.0]
    cols = [P.C1, P.C2, P.C3]
    fib_th, fib_ke = [], []
    Kgrid = np.linspace(20, 600, 120)
    for Kk in Kgrid:
        a = m / (Kk + m)
        Q = dict(a=a, b=m, c=L)
        x, _ = T.fit_point(Q, weights=None)
        fib_th.append(x[0]); fib_ke.append(x[1])
    ax[2].plot(fib_th, fib_ke, color=P.GRID, lw=3.0, zorder=1,
               label="fibre (K_m swept)")
    for Kk, col in zip(Kms, cols):
        a = m / (Kk + m)
        x, _ = T.fit_point(dict(a=a, b=m, c=L), weights=None)
        ax[2].plot(x[0], x[1], "o", color=col, ms=9, zorder=5)
        ax[2].annotate(f"K_m={Kk:.0f}", (x[0], x[1]), textcoords="offset points",
                       xytext=(6, 5), fontsize=8, color=col)
    ax[2].plot(T.TH5, T.KE, "*", color=P.INK, ms=15, zorder=6,
               label="truth (K_m=100)")
    ax[2].set_title("C  the pencil is a family in K_m\nconcurrence slides along the fibre (Remark 2)",
                    fontsize=10)
    ax[2].legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax[2].set_xlim(0, 500)

    for a in ax:
        a.set_xlabel("θ  (enzyme capacity, V_max)")
        a.set_ylabel("k_e  (linear clearance rate)")
        P.finish(a)
    P.save(fig, "r18_pencil")


if __name__ == "__main__":
    fig_pencil()
