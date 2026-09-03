"""
Model B: both clearance channels saturating.

    dm/dt = s(t) - theta*g(m) - ke*h(m)
    g(m) = m/(Km+m) ,  h(m) = m/(Ke+m)

Question answered: when does a crisis occur?

Total clearance is now bounded by theta + ke. If the load exceeds that
ceiling no equilibrium exists and the metabolite diverges. Clinically the
crisis begins earlier, when the equilibrium migrates past the tolerable
level, so the operative threshold is s_peak > C(m_bar).

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
from common.style import save, GREEN, RED, BLUE, AMBER, GREY, MUTED


P = ParamsB()
CEIL = P.theta + P.ke


# =====================================================================
def fig_ceiling():
    """Total clearance is bounded; loads above the ceiling admit no equilibrium."""
    m = np.linspace(0, 3000, 1500)
    C = clearance_B(m, P)

    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.plot(m, C, color=GREEN, lw=2.2, label=r"total clearance $C(m)$")
    ax.axhline(CEIL, color=GREEN, ls=":", lw=1.1)
    ax.text(2450, CEIL + 60, r"ceiling $\theta + k_e$", color=GREEN, fontsize=8)

    for s, col, lab, ls in [(1500, BLUE, "baseline load", "-"),
                            (2900, AMBER, "moderate illness", "-"),
                            (3600, RED, "peak load: no equilibrium", "--")]:
        ax.axhline(s, color=col, lw=1.4, ls=ls)
        me = steady_state(clearance_B, P, s)
        if np.isfinite(me):
            ax.plot([me], [s], "o", ms=5.5, color=col, zorder=4)
        ax.text(2450, s + 60, lab, color=col, fontsize=8)

    ax.set_xlabel(r"metabolite $m$")
    ax.set_ylabel("flux")
    ax.set_ylim(0, 4000)
    ax.set_title("bounded clearance: above the ceiling there is no steady state")
    save(fig, "fig_B1_ceiling")


# =====================================================================
def fig_bifurcation():
    """Equilibrium metabolite level against load. It escapes, it does not collide."""
    loads = np.linspace(50, CEIL * 0.999, 700)
    eq = np.array([steady_state(clearance_B, P, s) for s in loads])
    m_bar = 400.0

    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.plot(loads, eq, color=GREEN, lw=2.0, label="stable equilibrium")
    ax.axvline(CEIL, color=RED, ls="--", lw=1.3)
    ax.text(CEIL - 60, 1700, r"$s = \theta + k_e$", color=RED,
            fontsize=8, rotation=90, va="top", ha="right")
    ax.axhline(m_bar, color=BLUE, ls=":", lw=1.2)
    ax.text(120, m_bar + 60, r"tolerable level $\bar m$", color=BLUE, fontsize=8)

    s_clin = clearance_B(m_bar, P)
    ax.plot([s_clin], [m_bar], "o", ms=6, color=BLUE, zorder=4)
    ax.annotate("clinical threshold\narrives first", xy=(s_clin, m_bar),
                xytext=(s_clin - 1400, m_bar + 900), color=BLUE, fontsize=8,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))

    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"equilibrium $m^*$")
    ax.set_ylim(0, 2200)
    ax.set_title("the equilibrium escapes to infinity, it does not annihilate")
    ax.legend(loc="upper left")
    save(fig, "fig_B2_bifurcation")


# =====================================================================
def fig_thetacrit():
    """theta_crit depends on two quantities that are never measured."""
    s_peaks = np.linspace(1500, 4500, 160)
    m_bars = np.linspace(150, 900, 160)
    S, M = np.meshgrid(s_peaks, m_bars)
    TH = (S - P.ke * h(M, P.Ke)) / g(M, P.Km)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    cs = ax.contourf(S, M, TH, levels=18, cmap="RdYlGn_r", alpha=0.85)
    cl = ax.contour(S, M, TH, levels=8, colors="white", linewidths=0.6)
    ax.clabel(cl, inline=True, fontsize=6.5, fmt="%.0f")
    fig.colorbar(cs, ax=ax, label=r"$\theta_{\mathrm{crit}}$")
    ax.set_xlabel(r"assumed peak load $s_{\mathrm{peak}}$  (not measured)")
    ax.set_ylabel(r"assumed tolerable level $\bar m$  (not measured)")
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
    fig, ax = plt.subplots(figsize=(4.8, 3.3))
    ax.loglog(spreads * 100, cA, color=GREEN, label="Model A  (3 parameters)")
    ax.loglog(spreads * 100, cB, color=RED, ls="--",
              label="Model B  (4 parameters)")
    ax.set_xlabel("relative excursion of $m$  [%]")
    ax.set_ylabel("condition number")
    ax.set_title("the more faithful model is the harder one to identify")
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
