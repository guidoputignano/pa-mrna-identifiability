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
    fig_C5_margin       attractor, separatrix and the margin between them
                        against theta, with the intrinsic theta_crit

Run:  python src/model_C/run.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dataclasses import replace

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from common.models import (ParamsB, ParamsC, g, clearance_B,
                           clearance_C_quasistatic, equilibria_nonmonotone,
                           steady_state)
from common.style import save, figsize, GREEN, RED, BLUE, AMBER, GREY, MUTED


PC = ParamsC()
PB = ParamsB()
S_BASE = 600.0          # baseline substrate load, illustrative

# Model B matched to Model C at the same operating point, for a fair
# comparison: same theta and Km, alternative route saturating instead of
# declining, and ke chosen so that the two models put the patient at the SAME
# baseline attractor. Matching the operating point is the whole purpose of the
# comparison, so ke is solved for rather than guessed.
_M_BASE_C = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE,
                                   np.linspace(1e-3, 3000.0, 12000))[0][0]
_KE_B = brentq(lambda ke: clearance_B(_M_BASE_C,
                                      ParamsB(PC.theta, PC.Km, ke, 50.0)) - S_BASE,
               1.0, 1e5)
PB2 = ParamsB(theta=PC.theta, Km=PC.Km, ke=_KE_B, Ke=50.0)


# =====================================================================
def clearance_peak(p: ParamsC, theta=None):
    """
    The maximum of the Model C clearance curve, and where it occurs.

    C(m) = theta*g(m) + alt(m) rises from zero, peaks, and then decays back
    towards theta, so the maximum is interior. Its value is the largest load
    the patient can balance at all: the saddle-node in the load sits exactly
    there.
    """
    q = p if theta is None else replace(p, theta=float(theta))
    grid = np.linspace(1e-6, 6000.0, 30000)
    C = clearance_C_quasistatic(grid, q)
    i = int(np.argmax(C))
    lo = grid[max(i - 1, 0)]
    hi = grid[min(i + 1, len(grid) - 1)]
    r = minimize_scalar(lambda m: -clearance_C_quasistatic(m, q),
                        bounds=(lo, hi), method="bounded")
    return float(-r.fun), float(r.x)


def theta_crit_intrinsic(p: ParamsC, s_base, theta_hi=None):
    """
    The theta at which the fold occurs at the given load. Defined by the
    dynamics, not by a clinical convention.

    At load s_base the equilibria of Model C are the solutions of
    C(m; theta) = s_base. Lowering theta lowers C everywhere, so it lowers the
    peak of C. theta_crit_intrinsic is the theta at which that peak has fallen
    exactly to s_base: attractor and separatrix merge at the saddle-node, and
    below it no equilibrium exists at all and the patient cannot be maintained
    at baseline load.

    The point of the construction is what it does NOT contain. The definition
    in Model B,

        theta_crit = (s_peak - ke*m_bar) / g(m_bar),

    imports two quantities that are never measured: m_bar, a clinical
    convention, and s_peak, the amplitude of an illness. The definition here
    imports neither m_bar nor any tolerability convention: the threshold is
    the point at which the model's own equilibrium structure changes. It does
    still take the load as an argument, so s has not been eliminated, only
    m_bar. See the remark in the paper.

    Returns nan when even a null patient (theta = 0) can balance s_base, in
    which case the fold lies at negative theta and the quantity is vacuous.
    """
    f = lambda th: clearance_peak(p, th)[0] - s_base
    if f(0.0) >= 0.0:
        return float("nan")
    hi = theta_hi if theta_hi is not None else max(4.0 * s_base, 1.0)
    while f(hi) < 0.0 and hi < 1e9:
        hi *= 2.0
    return float(brentq(f, 0.0, hi, xtol=1e-10, rtol=1e-12))


# =====================================================================
def fig_phaseline():
    """The phase line. Two equilibria at baseline load, none at peak load."""
    m = np.linspace(1e-3, 2400, 4000)
    s_base, s_peak = S_BASE, 850.0

    fig, ax = plt.subplots(figsize=figsize(0.49, 0.86))
    for s, col, lab, ls in [(s_base, GREEN, "baseline load", "-"),
                            (s_peak, RED, "peak load", "--")]:
        dm = s - np.array([clearance_C_quasistatic(x, PC) for x in m])
        ax.plot(m, dm, color=col, ls=ls, label=lab)

    ax.axhline(0, color=MUTED, lw=0.8)
    roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, PC, s_base,
                                         grid=np.linspace(1e-3, 2400, 12000))
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
    ax.set_title("depletion makes clearance non-monotone,\ncreating a separatrix")
    ax.legend(loc="lower right", fontsize=7.5)
    save(fig, "fig_C1_phaseline")
    print(f"    equilibria at baseline load: {np.round(roots,1)}  stable={stab}")


# =====================================================================
def fig_fold():
    """Bifurcation diagram: the two branches collide at a fold."""
    # Windowed on the fold. The unstable branch climbs steeply as the load
    # falls, so starting far below the saddle-node squashes the collision,
    # which is the whole point of the figure, into the bottom of the axes.
    s_fold = clearance_peak(PC)[0]
    loads = np.linspace(0.72 * s_fold, 1.004 * s_fold, 480)
    stable_s, stable_m, unst_s, unst_m = [], [], [], []
    for s in loads:
        roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, PC, s,
                                             grid=np.linspace(1e-3, 8000, 24000))
        for r, st in zip(roots, stab):
            (stable_s if st else unst_s).append(s)
            (stable_m if st else unst_m).append(r)

    fig, ax = plt.subplots(figsize=figsize(0.49, 0.86))
    ax.plot(stable_s, stable_m, ".", ms=2.0, color=GREEN, label="stable branch")
    ax.plot(unst_s, unst_m, ".", ms=2.0, color=RED, label="unstable branch")
    ax.axvline(s_fold, color=GREY, ls=":", lw=1.0)
    ax.text(s_fold - 4, 0.97 * max(unst_m), "saddle-node", rotation=90,
            fontsize=7.5, color=MUTED, ha="right", va="top")
    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"equilibrium $m^*$")
    ax.set_title("two branches that collide,\nunlike Model B")
    # between the branches: upper left is where the unstable one climbs
    ax.legend(loc="center left", markerscale=4, framealpha=0.95)
    save(fig, "fig_C2_fold")


# =====================================================================
def _simulate(rhs, m0, T=160.0, illness=(20.0, 34.0), s_base=S_BASE, s_ill=900.0):
    def s(t):
        return s_ill if illness[0] <= t <= illness[1] else s_base
    sol = solve_ivp(lambda t, y: [s(t) - rhs(y[0])], (0, T), [m0],
                    dense_output=True, max_step=0.2, rtol=1e-8, atol=1e-8)
    tt = np.linspace(0, T, 1400)
    return tt, sol.sol(tt)[0], np.array([s(t) for t in tt])


def fig_timeseries():
    """The discriminating prediction, as a simulation."""
    m0 = _M_BASE_C          # the shared baseline attractor, by construction
    tt, mB, sB = _simulate(lambda x: clearance_B(x, PB2), m0=m0)
    tt, mC, sC = _simulate(lambda x: clearance_C_quasistatic(x, PC), m0=m0)

    fig, ax = plt.subplots(2, 1, figsize=figsize(0.48, 1.05), sharex=True,
                           gridspec_kw={"height_ratios": [1, 2.4]})
    ax[0].fill_between(tt, 0, sC, color=AMBER, alpha=0.25, step="mid")
    ax[0].plot(tt, sC, color=AMBER, lw=1.2)
    ax[0].set_ylabel(r"load $s(t)$")
    # The label belongs over the pulse it names. In the lower panel it had to
    # share a log axis with two curves and a legend, and touched whichever
    # moved.
    ax[0].text(27, sC.max() * 1.10, "illness", ha="center", va="bottom",
               fontsize=7.5, color=MUTED)
    ax[0].set_ylim(top=sC.max() * 1.45)
    ax[0].set_title("identical illness,\ntwo clearance hypotheses")

    ax[1].plot(tt, mB, color=GREEN, label="Model B: returns to baseline")
    ax[1].plot(tt, mC, color=RED, ls="--", label="Model C: does not return")
    ax[1].axvspan(20, 34, color=AMBER, alpha=0.12)
    ax[1].set_xlabel("time  [days]")
    ax[1].set_ylabel(r"metabolite $m$")
    ax[1].set_yscale("log")
    # Headroom before the legend, or the box lands on the Model C curve where
    # it climbs. On a log axis the factor is decades, not units.
    ax[1].set_ylim(top=max(mB.max(), mC.max()) * 30)
    ax[1].legend(loc="upper left", framealpha=0.95)
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
    # The separatrix exists only for theta < s < max_m C(m), since total
    # clearance tends to theta from above: below theta it never falls back
    # through s, and above the peak there is no equilibrium at all.
    s_fold = clearance_peak(PC)[0]          # the saddle-node in the load
    loads = np.linspace(0.60 * s_fold, 0.999 * s_fold, 400)
    stab_b, unst_b, ok = [], [], []
    for s in loads:
        roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, s,
                                           grid=np.linspace(1e-3, 20000, 60000))
        if st.sum() and (~st).sum():
            stab_b.append(roots[st][0]); unst_b.append(roots[~st][0]); ok.append(s)
    ok = np.array(ok); stab_b = np.array(stab_b); unst_b = np.array(unst_b)

    fig, ax = plt.subplots(figsize=figsize(0.50, 0.88))
    ax.fill_between(ok, unst_b, 3000, color=RED, alpha=0.10)
    ax.fill_between(ok, 0, unst_b, color=GREEN, alpha=0.08)
    ax.plot(ok, stab_b, color=GREEN, lw=2.0, label="attractor")
    ax.plot(ok, unst_b, color=RED, lw=2.0, ls="--", label="separatrix")
    # Every label below is placed against one of four free regions: the red
    # wedge above the separatrix on the right, the green basin along the
    # bottom, the band between the horizontal arrow and the separatrix, and
    # the strip left of the baseline line. The earlier version put three of
    # them at the same height and they printed on top of one another. The
    # captions carry the explanations, so the labels here name the route and
    # stop.
    # low in the red wedge: the upper right is now the patient's corner
    ax.text(0.997 * s_fold, 1000, "escape:\nmetabolite\ndiverges", color=RED,
            fontsize=8, ha="right")
    ax.text(0.66 * s_fold, 300, "recovery basin", color=GREEN, fontsize=8)

    # a patient who has crossed the separatrix, at BASELINE load
    # The separatrix at baseline load now sits at about 1036, so the
    # illustrative crisis point has to clear it visibly rather than graze it.
    s_base, m_crisis = S_BASE, 1700.0
    ax.axvline(s_base, color=GREY, ls=":", lw=1.0)
    # Horizontal, in the empty strip above the separatrix. Rotated against
    # the line it ran into whichever label shared that column.
    ax.text(s_base, 2560, "baseline load", color=MUTED, fontsize=7.5,
            ha="center", va="top")
    ax.plot([s_base], [m_crisis], "o", ms=6, color=MUTED, zorder=5)
    ax.annotate("patient above\nthe separatrix", xy=(s_base, m_crisis),
                xytext=(s_base + 12, m_crisis + 170), fontsize=8, color=MUTED,
                ha="left", va="bottom",
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))

    # route 1: push the load BELOW baseline (glucose, protein withdrawal)
    s_rescue = 0.62 * s_fold
    ax.annotate("", xy=(s_rescue, m_crisis), xytext=(s_base, m_crisis),
                arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.8))
    ax.text(0.5 * (s_rescue + s_base), m_crisis - 130, "suppress catabolism",
            color=AMBER, fontsize=7.5, ha="center", va="top")

    # route 2: lower the metabolite directly
    ax.annotate("", xy=(s_base, 60), xytext=(s_base, m_crisis),
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.8))
    ax.text(s_base - 12, 640, "haemofiltration", color=BLUE,
            fontsize=7.5, ha="right", va="center")

    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"metabolite $m$")
    ax.set_ylim(0, 2600)
    ax.set_title("recovery requires pushing the load below baseline,\n"
                 "not merely removing the trigger")
    # Upper left is the one corner no curve and no label reaches.
    ax.legend(loc="upper left", framealpha=0.95)
    save(fig, "fig_C4_basin")
    print("    The model reproduces standard care: intravenous glucose lowers s")
    print("    below baseline, which raises the separatrix above the patient and")
    print("    recaptures them; haemofiltration is the alternative route, lowering")
    print("    m directly when the load cannot be reduced enough.")
    print("    TESTABLE: in Model B, cessation of the illness suffices. In Model C,")
    print("    the load must be driven BELOW baseline. Records should show whether")
    print("    decompensations resolve on return to baseline management alone.")


# =====================================================================
def fig_margin():
    """
    The safety margin against theta, at baseline load.

    Everything here is a property of the model at a fixed load. As theta
    falls, the attractor rises and the separatrix falls, and the distance
    between them, which is how far the patient can be pushed before the
    crisis becomes self-sustaining, closes to zero at theta_crit_intrinsic.
    That distance is the clinically decisive quantity and it is what the
    threshold of Model B never displays, because Model B has no separatrix.

    The upper end of the range is equally structural. As theta approaches the
    baseline load the separatrix escapes to infinity and the patient becomes
    unconditionally safe at that load: total clearance tends to theta from
    above, so it never falls back through s.
    """
    thc = theta_crit_intrinsic(PC, S_BASE)
    thetas = np.linspace(thc * 1.002, 0.86 * S_BASE, 260)
    grid = np.linspace(1e-3, 30000.0, 90000)

    th_ok, attr, sep = [], [], []
    for th in thetas:
        q = replace(PC, theta=float(th))
        roots, st = equilibria_nonmonotone(clearance_C_quasistatic, q, S_BASE, grid)
        if st.sum() and (~st).sum():
            th_ok.append(th); attr.append(roots[st][0]); sep.append(roots[~st][0])
    th_ok = np.array(th_ok); attr = np.array(attr); sep = np.array(sep)
    margin = sep - attr

    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.40))

    ax[0].fill_between(th_ok, attr, sep, color=GREEN, alpha=0.10)
    ax[0].plot(th_ok, attr, color=GREEN, lw=2.0, label="attractor")
    ax[0].plot(th_ok, sep, color=RED, lw=2.0, ls="--", label="separatrix")
    ax[0].axvline(thc, color=MUTED, ls=":", lw=1.2)
    ax[0].annotate(r"$\theta_{\mathrm{crit}}$" + f" = {thc:.0f}:"
                   "\nthe two curves merge",
                   xy=(th_ok[0], 0.5 * (attr[0] + sep[0])),
                   xytext=(thc + 55, 620), fontsize=7.5, color=MUTED,
                   arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))
    ax[0].plot([PC.theta], [attr[np.argmin(abs(th_ok - PC.theta))]], "o", ms=5,
               color=MUTED, zorder=5)
    ax[0].plot([PC.theta], [sep[np.argmin(abs(th_ok - PC.theta))]], "o", ms=5,
               color=MUTED, zorder=5)
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"residual capacity $\theta$")
    ax[0].set_ylabel(r"metabolite $m$")
    ax[0].set_title("the two equilibria collide at the fold")
    ax[0].legend(loc="upper left")

    ax[1].plot(th_ok, margin, color=BLUE, lw=2.0)
    ax[1].axvline(thc, color=MUTED, ls=":", lw=1.2)
    ax[1].fill_between(th_ok, 0, margin, color=BLUE, alpha=0.08)
    i = np.argmin(abs(th_ok - PC.theta))
    ax[1].plot([PC.theta], [margin[i]], "o", ms=5.5, color=MUTED, zorder=5)
    ax[1].annotate(f"nominal patient\nmargin = {margin[i]:.0f}",
                   xy=(PC.theta, margin[i]), xytext=(PC.theta - 20, margin[i] + 380),
                   fontsize=7.5, color=MUTED,
                   arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.9))
    ax[1].text(thc + 8, margin.max() * 0.82,
               "margin vanishes:\nsaddle-node", color=MUTED, fontsize=7.5)
    ax[1].set_xlabel(r"residual capacity $\theta$")
    ax[1].set_ylabel("separatrix $-$ attractor")
    ax[1].set_title("the safety margin, at baseline load")

    save(fig, "fig_C5_margin")
    print(f"    theta_crit_intrinsic at s = {S_BASE:.0f}: {thc:.2f}")
    print(f"    at the nominal theta = {PC.theta:.0f}: attractor "
          f"{attr[i]:.1f}, separatrix {sep[i]:.1f}, margin {margin[i]:.1f}")
    print("    m_bar does NOT appear anywhere in this construction. s does.")
    return thc


if __name__ == "__main__":
    print("MODEL C")
    fig_phaseline()
    fig_fold()
    fig_timeseries()
    fig_basin()
    fig_margin()
