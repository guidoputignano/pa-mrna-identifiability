"""
Model C: carnitine depletion.

    dm/dt = s(t) - theta*g(m) - k*c*m
    dc/dt = R - k*c*m

Each molecule of propionate eliminated consumes one molecule of carnitine,
which leaves with the product. Carnitine is a consumed reagent, not a
catalyst. If depletion outpaces replenishment, effective clearance
DECREASES with m rather than merely saturating, the clearance curve becomes
non-monotone, and the system admits two equilibria: a stable one and an
unstable one which acts as a separatrix.

Question answered: why does recovery require active reduction of the load
below baseline, rather than mere removal of the trigger?

Removing the trigger returns s to baseline, and at baseline the patient who
has crossed the separatrix is still on the wrong side of it. Recovery needs
an overshoot: either s driven BELOW baseline, which raises the separatrix
above the patient, or m lowered directly. Under B, by contrast, returning s
to baseline returns m to baseline.

Produces
    fig_C1_phaseline    dm/dt against m, showing the two equilibria
    fig_C2_fold         bifurcation diagram with the saddle-node fold
    fig_C3_timeseries   B and C responses to an identical transient illness
    fig_C4_basin        the basin of attraction and the two routes to recovery

Run:  python src/model_C/run.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from common.models import (ParamsB, ParamsC, g, clearance_B,
                           clearance_C_quasistatic, equilibria_nonmonotone)
from common.style import save, GREEN, RED, BLUE, AMBER, GREY, MUTED


PC = ParamsC()
PB = ParamsB()
# Model B matched to Model C at the same operating point, for a fair
# comparison: same theta, alternative route saturating instead of declining.
PB2 = ParamsB(theta=PC.theta, Km=PC.Km, ke=520.0, Ke=50.0)


# =====================================================================
def fig_phaseline():
    """The phase line. Two equilibria at baseline load, none at peak load."""
    m = np.linspace(1e-3, 1100, 2500)
    s_base, s_peak = 600.0, 850.0

    fig, ax = plt.subplots(figsize=(5.6, 3.5))
    for s, col, lab, ls in [(s_base, GREEN, "baseline load", "-"),
                            (s_peak, RED, "peak load", "--")]:
        dm = s - np.array([clearance_C_quasistatic(x, PC) for x in m])
        ax.plot(m, dm, color=col, ls=ls, label=lab)

    ax.axhline(0, color=MUTED, lw=0.8)
    roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, PC, s_base,
                                         grid=np.linspace(1e-3, 1100, 8000))
    for r, st in zip(roots, stab):
        ax.plot([r], [0], "o", ms=7, zorder=5,
                color=GREEN if st else "white",
                markeredgecolor=GREEN, markeredgewidth=1.8)
        ax.annotate("stable\nattractor" if st else "unstable:\nthe separatrix",
                    xy=(r, 0), xytext=(r, -220 if st else 300),
                    ha="center", fontsize=7.5, color=GREEN if st else RED,
                    arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7))

    ax.set_xlabel(r"metabolite $m$")
    ax.set_ylabel(r"$dm/dt$")
    ax.set_ylim(-320, 480)
    ax.set_title("depletion makes clearance non-monotone, creating a separatrix")
    ax.legend(loc="upper left")
    save(fig, "fig_C1_phaseline")
    print(f"    equilibria at baseline load: {np.round(roots,1)}  stable={stab}")


# =====================================================================
def fig_fold():
    """Bifurcation diagram: the two branches collide at a fold."""
    loads = np.linspace(250, 800, 420)
    stable_s, stable_m, unst_s, unst_m = [], [], [], []
    for s in loads:
        roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, PC, s,
                                             grid=np.linspace(1e-3, 2000, 6000))
        for r, st in zip(roots, stab):
            (stable_s if st else unst_s).append(s)
            (stable_m if st else unst_m).append(r)

    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    ax.plot(stable_s, stable_m, ".", ms=2.0, color=GREEN, label="stable branch")
    ax.plot(unst_s, unst_m, ".", ms=2.0, color=RED, label="unstable branch")
    if unst_s:
        s_fold = max(stable_s)
        ax.axvline(s_fold, color=GREY, ls=":", lw=1.0)
        ax.text(s_fold - 30, max(stable_m) * 0.9, "saddle-node", rotation=90,
                fontsize=7.5, color=MUTED, ha="right", va="top")
    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"equilibrium $m^*$")
    ax.set_title("two branches that collide, unlike Model B")
    ax.legend(loc="upper left")
    save(fig, "fig_C2_fold")


# =====================================================================
def _simulate(rhs, m0, T=160.0, illness=(20.0, 34.0), s_base=600.0, s_ill=900.0):
    def s(t):
        return s_ill if illness[0] <= t <= illness[1] else s_base
    sol = solve_ivp(lambda t, y: [s(t) - rhs(y[0])], (0, T), [m0],
                    dense_output=True, max_step=0.2, rtol=1e-8, atol=1e-8)
    tt = np.linspace(0, T, 1400)
    return tt, sol.sol(tt)[0], np.array([s(t) for t in tt])


def fig_timeseries():
    """The discriminating prediction, as a simulation."""
    tt, mB, sB = _simulate(lambda x: clearance_B(x, PB2), m0=47.0)
    tt, mC, sC = _simulate(lambda x: clearance_C_quasistatic(x, PC), m0=47.0)

    fig, ax = plt.subplots(2, 1, figsize=(5.8, 4.4), sharex=True,
                           gridspec_kw={"height_ratios": [1, 2.4]})
    ax[0].fill_between(tt, 0, sC, color=AMBER, alpha=0.25, step="mid")
    ax[0].plot(tt, sC, color=AMBER, lw=1.2)
    ax[0].set_ylabel(r"load $s(t)$")
    ax[0].set_title("identical transient illness, two clearance hypotheses")

    ax[1].plot(tt, mB, color=GREEN, label="Model B: returns to baseline")
    ax[1].plot(tt, mC, color=RED, ls="--", label="Model C: does not return")
    ax[1].axvspan(20, 34, color=AMBER, alpha=0.12)
    ax[1].text(27, ax[1].get_ylim()[1] * 0.05, "illness", ha="center",
               fontsize=7.5, color=MUTED)
    ax[1].set_xlabel("time  [days]")
    ax[1].set_ylabel(r"metabolite $m$")
    ax[1].set_yscale("log")
    ax[1].legend(loc="upper left")
    save(fig, "fig_C3_timeseries")


# =====================================================================
def fig_basin():
    """
    There is no upper stable state in Model C: beyond the separatrix the
    metabolite grows without bound. So the clinically meaningful object is
    not a hysteresis loop but the BASIN OF ATTRACTION.

    The consequence is the question this model answers. Removing the trigger
    puts s back at baseline, and the separatrix at baseline load is still
    below the patient: nothing has been recovered. Recovery needs an
    overshoot, and the figure shows the two that standard care uses. Driving
    s BELOW baseline raises the separatrix above the patient. Lowering m
    directly, by haemofiltration, moves the patient below the separatrix
    where it stands.
    """
    loads = np.linspace(250, 781, 400)
    stab_b, unst_b, ok = [], [], []
    for s in loads:
        roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, s,
                                           grid=np.linspace(1e-3, 3000, 6000))
        if st.sum() and (~st).sum():
            stab_b.append(roots[st][0]); unst_b.append(roots[~st][0]); ok.append(s)
    ok = np.array(ok); stab_b = np.array(stab_b); unst_b = np.array(unst_b)

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.fill_between(ok, unst_b, 3000, color=RED, alpha=0.10)
    ax.fill_between(ok, 0, unst_b, color=GREEN, alpha=0.08)
    ax.plot(ok, stab_b, color=GREEN, lw=2.0, label="attractor")
    ax.plot(ok, unst_b, color=RED, lw=2.0, ls="--", label="separatrix")
    ax.text(300, 1600, "escape:\nmetabolite diverges", color=RED, fontsize=8)
    ax.text(300, 120, "recovery basin", color=GREEN, fontsize=8)

    # a patient who has crossed the separatrix, at BASELINE load
    s_base, m_crisis = 600.0, 1200.0
    ax.axvline(s_base, color=GREY, ls=":", lw=1.0)
    ax.text(s_base + 6, 2400, "baseline load", color=MUTED, fontsize=7.5, rotation=90,
            va="top")
    ax.plot([s_base], [m_crisis], "o", ms=6, color=MUTED, zorder=5)
    ax.annotate("patient above the separatrix\nat baseline load", xy=(s_base, m_crisis),
                xytext=(s_base + 30, m_crisis + 800), fontsize=8, color=MUTED,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))

    # route 1: push the load BELOW baseline (glucose, protein withdrawal)
    ax.annotate("", xy=(430, m_crisis), xytext=(s_base, m_crisis),
                arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.8))
    ax.text(408, m_crisis + 140, "suppress catabolism:\nraises the separatrix\n"
            "above the patient", color=AMBER, fontsize=7.5)

    # route 2: lower the metabolite directly
    ax.annotate("", xy=(s_base, 60), xytext=(s_base, m_crisis),
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.8))
    ax.text(s_base + 10, 420, "haemofiltration:\nlowers $m$ directly", color=BLUE,
            fontsize=7.5)

    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"metabolite $m$")
    ax.set_ylim(0, 2600)
    ax.set_title("recovery requires pushing the load below baseline,\n"
                 "not merely removing the trigger")
    ax.legend(loc="upper right")
    save(fig, "fig_C4_basin")
    print("    The model reproduces standard care: intravenous glucose lowers s")
    print("    below baseline, which raises the separatrix above the patient and")
    print("    recaptures them; haemofiltration is the alternative route, lowering")
    print("    m directly when the load cannot be reduced enough.")
    print("    TESTABLE: in Model B, cessation of the illness suffices. In Model C,")
    print("    the load must be driven BELOW baseline. Records should show whether")
    print("    decompensations resolve on return to baseline management alone.")


if __name__ == "__main__":
    print("MODEL C")
    fig_phaseline()
    fig_fold()
    fig_timeseries()
    fig_basin()
