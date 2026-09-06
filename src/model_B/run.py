"""
Model B: both clearance channels saturating.

    dm/dt = s(t) - theta*g(m) - ke*h(m)
    g(m) = m/(Km+m) ,  h(m) = m/(Ke+m)

Question answered: under what condition is a tolerable steady state lost?

Total clearance is now bounded by theta + ke. If the load exceeds that
ceiling no equilibrium exists and the metabolite diverges. Clinically the
tolerable steady state is lost earlier, when the equilibrium migrates past
the tolerable level, so the operative threshold is s_peak > C(m_tol).

NOTE the limit of that question. The model returns a CONDITION, not a
FREQUENCY: nothing here says how often the load crosses the threshold,
because the arrival of illness is exogenous to the vector field. Frequency
is the one thing the records hold reliably (admissions are dated), while
s_peak and m_tol are never measured. Removing the mismatch needs a hazard
formulation, not a deterministic threshold. See the open items.

Produces
    fig_B1_ceiling      clearance curve, ceiling, and three load lines
    fig_B2_bifurcation  equilibrium against load, showing the escape
    fig_B3_thetacrit    theta_crit as a function of two unmeasured quantities
    fig_B4_cost         identifiability cost of the correction, A against B

Run:  python src/model_B/run.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt

from common.models import (ParamsA, ParamsB, g, h, clearance_A, clearance_B,
                           steady_state, sensitivity_A, sensitivity_B,
                           condition_number)
from common.style import (save, figsize, GREEN, RED, BLUE, AMBER, GREY, MUTED,
                          corner, over_x, at_point)


P = ParamsB()
CEIL = P.theta + P.ke
M_MAX = 3000.0     # right edge of the Model B panels, named once


# =====================================================================
def fig_ceiling():
    """Total clearance is bounded; loads above the ceiling admit no equilibrium."""
    m = np.linspace(0, M_MAX, 1500)
    C = clearance_B(m, P)

    fig, ax = plt.subplots(figsize=figsize(0.49, 0.86))
    ax.plot(m, C, color=GREEN, lw=2.2, label=r"total clearance $C(m)$")
    ax.axhline(CEIL, color=GREEN, ls=":", lw=1.1)
    # The band between the ceiling and the peak-load line is only about 7pt
    # tall on the page, and the clearance curve crowds it from below, so the
    # label is cut to the symbol. The title already supplies the word.
    at_point(ax, (M_MAX, CEIL), r"$\theta + k_e$", dx=-3, dy=3,
             color=GREEN, fontsize=7.5, ha="right", va="bottom")

    # Each label clears its own line by a fixed number of POINTS, so the gap
    # is the same on the page whatever the axis limits become. x is data,
    # because the clearance curve crosses the middle of the panel and a single
    # x for all three would put one of them on it.
    for s, col, lab, ls, x, side in [
            (1500, BLUE, "baseline load", "-", 1200, +1),
            (2900, AMBER, "moderate illness", "-", 1500, -1),
            (3600, RED, "peak load: no equilibrium", "--", 300, +1)]:
        ax.axhline(s, color=col, lw=1.4, ls=ls)
        me = steady_state(clearance_B, P, s)
        if np.isfinite(me):
            ax.plot([me], [s], "o", ms=5.5, color=col, zorder=4)
        at_point(ax, (x, s), lab, dy=4 * side, color=col, fontsize=7.5,
                 va="bottom" if side > 0 else "top")

    ax.set_xlabel(r"metabolite $m$")
    ax.set_ylabel("flux")
    ax.set_xlim(0, M_MAX)
    ax.set_ylim(0, 4200)
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("bounded clearance: above the ceiling\nthere is no steady state")
    save(fig, "fig_B1_ceiling")


# =====================================================================
def fig_bifurcation():
    """Equilibrium metabolite level against load. It escapes, it does not collide."""
    loads = np.linspace(50, CEIL * 0.999, 700)
    eq = np.array([steady_state(clearance_B, P, s) for s in loads])
    m_tol = 400.0

    fig, ax = plt.subplots(figsize=figsize(0.49, 0.86))
    ax.plot(loads, eq, color=GREEN, lw=2.0, label="stable equilibrium")
    ax.axvline(CEIL, color=RED, ls="--", lw=1.3)
    at_point(ax, (CEIL, 0), r"$s = \theta + k_e$", dx=-4, dy=0, color=RED,
             fontsize=8, rotation=90, va="bottom", ha="right",
             annotation_clip=False)
    ax.axhline(m_tol, color=BLUE, ls=":", lw=1.2)
    at_point(ax, (0, m_tol), r"tolerable level $m_{\mathrm{tol}}$", dx=8, dy=4,
             color=BLUE, fontsize=8, va="bottom")

    s_clin = clearance_B(m_tol, P)
    ax.plot([s_clin], [m_tol], "o", ms=6, color=BLUE, zorder=4)
    at_point(ax, (s_clin, m_tol), "clinical threshold\narrives first",
             dx=-14, dy=48, arrow=True, color=BLUE, fontsize=8,
             ha="right", va="bottom")

    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"equilibrium $m^*$")
    ax.set_ylim(0, 2200)
    ax.set_title("the equilibrium escapes to infinity,\nit does not annihilate")
    ax.legend(loc="upper left")
    save(fig, "fig_B2_bifurcation")


# =====================================================================
def fig_thetacrit():
    """theta_crit depends on two quantities that are never measured."""
    s_peaks = np.linspace(1500, 4500, 160)
    m_tols = np.linspace(150, 900, 160)
    S, M = np.meshgrid(s_peaks, m_tols)
    TH = (S - P.ke * h(M, P.Ke)) / g(M, P.Km)

    fig, ax = plt.subplots(figsize=figsize(0.47, 0.95))
    cs = ax.contourf(S, M, TH, levels=18, cmap="RdYlGn_r", alpha=0.85)
    cl = ax.contour(S, M, TH, levels=8, colors="white", linewidths=0.6)
    ax.clabel(cl, inline=True, fontsize=6.5, fmt="%.0f")
    fig.colorbar(cs, ax=ax, label=r"$\theta_{\mathrm{crit}}$")
    # Both axis labels used to carry "(not measured)" and were clipped at
    # this panel width. The caption says it once instead.
    ax.set_xlabel(r"assumed peak load $s_{\mathrm{peak}}$")
    ax.set_ylabel(r"assumed tolerable level $m_{\mathrm{tol}}$")
    ax.set_title(r"$\theta_{\mathrm{crit}}$ is an operational definition,"
                 "\nnot a physiological constant")
    save(fig, "fig_B3_thetacrit")


# =====================================================================
def fig_cost():
    """Making the model biologically faithful costs identifiability."""
    m0 = 10.0
    spreads = np.logspace(-2, 1.3, 50)
    PA = ParamsA()
    cA, cB = [], []
    for sp in spreads:
        lv = m0 * np.array([1.0, 1 + sp, 1 + 2 * sp, 1 + 3 * sp])
        cA.append(condition_number(sensitivity_A, PA, lv[:3]))
        cB.append(condition_number(sensitivity_B, P, lv))
    fig, ax = plt.subplots(figsize=figsize(0.50, 0.86))
    ax.loglog(spreads * 100, cA, color=GREEN, label="Model A  (3 parameters)")
    ax.loglog(spreads * 100, cB, color=RED, ls="--",
              label="Model B  (4 parameters)")
    ax.set_xlabel("relative excursion of $m$  [%]")
    ax.set_ylabel("condition number")
    ax.set_title("the more faithful model is\nthe harder one to identify")
    ax.grid(True, which="both", alpha=0.4)
    ax.legend()
    save(fig, "fig_B4_cost")
    print(f"    at 10% excursion:  A = {cA[np.argmin(abs(spreads-0.1))]:.2e}, "
          f"B = {cB[np.argmin(abs(spreads-0.1))]:.2e}")


if __name__ == "__main__":
    print("MODEL B")
    fig_ceiling()
    fig_bifurcation()
    fig_thetacrit()
    fig_cost()
