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

from common.numbers import NumberBook
from common.models import (ParamsB, ParamsC, g, alt_flux_C, clearance_B,
                           clearance_C_quasistatic, equilibria_nonmonotone,
                           rhs_C, steady_state)
from common.style import (save, figsize, GREEN, RED, BLUE, AMBER, GREY, MUTED,
                          corner, over_x, at_point)


PC = ParamsC()
PB = ParamsB()
S_BASE = 600.0          # baseline substrate load, illustrative
#: the registered illness amplitude, as a fold of baseline load. It is
#: the `s_pk/s_base` row of the constants registry, and it is what the
#: reduction-cost bound below perturbs with, so the bound names a
#: perturbation the model already carries rather than one chosen for it.
from common.params import ILLNESS_FOLD  # noqa: E402
Y_MAX = 2600.0          # top of the basin panel, named once

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
#: every interpolated number this module draws on a figure, classified for
#: tools/check_artefacts.py:check_figure_numbers. All four are claims about the
#: model, and fig_margin() is their generator: it writes them into
#: results/numbers/margin.tex and both panels draw out of that same dict.
DRAWN = {
    "MARG['theta_crit']": "registered: theta_crit, where attractor and "
                          "separatrix merge",
    "MARG['margin_nominal']": "registered: the nominal patient's safety margin",
    "MARG['margin_band_lo']": "registered: the low end of that margin over the "
                              "M_i range",
    "MARG['margin_band_hi']": "registered: the high end of it",
}


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


def alt_peak(p: ParamsC):
    """The maximum of the alternative route alone, and where it occurs."""
    from common.models import alt_flux_C
    grid = np.linspace(1e-6, 6e4, 200000)
    v = alt_flux_C(grid, p)
    i = int(np.argmax(v))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    r = minimize_scalar(lambda m: -alt_flux_C(m, p), bounds=(lo, hi),
                        method="bounded")
    return float(-r.fun), float(r.x)


def R0_for_peak(p: ParamsC, target=500.0):
    """
    The replenishment that puts the alternative route's peak at `target`.

    R0 is NOT a free constant: it is fixed by requiring that the alternative
    route alone cannot carry the baseline load, so that a patient with no
    residual enzyme is not maintainable. alpha is exactly proportional to R0,
    so this is one evaluation and not a search.

    Every sensitivity analysis in this module re-derives R0 through this
    function after moving Kc, Mi or rho0. Holding R0 fixed while moving one of
    those three moves the peak as well, which breaks the constraint that
    defines R0, and reports the sensitivity of theta_crit to the calibration
    as though it were sensitivity to the shape constant. That conflation is
    what made M_i look nine times more expensive than K_c when, on the same
    footing, it is cheaper.
    """
    pk, _ = alt_peak(replace(p, R0=1.0))
    return float(target / pk)


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
        # Offset in POINTS from the equilibrium, and nudged inward from the
        # left spine: the attractor sits near m = 0, so a centred label there
        # runs off the axis.
        at_point(ax, (r, 0),
                 "stable\nattractor" if st else "unstable:\nthe separatrix",
                 dx=6 if st else 0, dy=-34 if st else 42,
                 ha="left" if st else "center", va="top" if st else "bottom",
                 fontsize=7.5, color=GREEN if st else RED,
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
    over_x(ax, s_fold, "saddle-node", frac=0.97, rotation=90, fontsize=7.5,
           color=MUTED, ha="right", va="top")
    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"equilibrium $m^*$")
    ax.set_title("two branches that collide,\nunlike Model B")
    # The free wedge is below the unstable branch and above the stable one.
    # Anchored in axes fractions, so it stays there when the fold moves.
    ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.14), markerscale=4,
              fontsize=7.5, frameon=True, facecolor="white", edgecolor="none",
              framealpha=0.95)
    save(fig, "fig_C2_fold")


# =====================================================================
ILLNESS = (20.0, 34.0)   # the transient illness window, shared by the
                         # simulation, the shading and the label


def _simulate(rhs, m0, T=160.0, illness=ILLNESS, s_base=S_BASE, s_ill=900.0):
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
    at_point(ax[0], (0.5 * (ILLNESS[0] + ILLNESS[1]), sC.max()), "illness",
             dy=5, ha="center", va="bottom", fontsize=7.5, color=MUTED)
    ax[0].set_ylim(top=sC.max() * 1.45)
    ax[0].set_title("identical illness,\ntwo clearance hypotheses")

    ax[1].plot(tt, mB, color=GREEN, label="Model B: returns to baseline")
    ax[1].plot(tt, mC, color=RED, ls="--", label="Model C: does not return")
    ax[1].axvspan(*ILLNESS, color=AMBER, alpha=0.12)
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
    # Four labelled annotations plus a legend never fitted this panel at
    # print size: they collided with each other, with the separatrix and with
    # the region shading. The plot now carries only pointers, two numbered
    # arrows and one marker, and the caption says what they are. The two
    # regions move into the legend rather than floating as text.
    import matplotlib.patches as mpatches
    ax.fill_between(ok, unst_b, Y_MAX, color=RED, alpha=0.10)
    ax.fill_between(ok, 0, unst_b, color=GREEN, alpha=0.08)
    ax.plot(ok, stab_b, color=GREEN, lw=2.0, label="attractor")
    ax.plot(ok, unst_b, color=RED, lw=2.0, ls="--", label="separatrix")

    # a patient who has crossed the separatrix, at BASELINE load
    s_base, m_crisis = S_BASE, 1700.0
    ax.axvline(s_base, color=GREY, ls=":", lw=1.0)
    over_x(ax, s_base, "baseline load", frac=0.985, color=MUTED, fontsize=7.5,
           ha="center")
    ax.plot([s_base], [m_crisis], "o", ms=6, color=MUTED, zorder=5)

    # route 1: push the load BELOW baseline (glucose, protein withdrawal)
    s_rescue = 0.62 * s_fold
    ax.annotate("", xy=(s_rescue, m_crisis), xytext=(s_base, m_crisis),
                arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.8))
    at_point(ax, (s_rescue, m_crisis), "1", dx=14, dy=-7, color=AMBER,
             fontsize=9, fontweight="bold", ha="left", va="top")

    # route 2: lower the metabolite directly
    ax.annotate("", xy=(s_base, 60), xytext=(s_base, m_crisis),
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.8))
    at_point(ax, (s_base, 0.5 * m_crisis), "2", dx=-6, color=BLUE,
             fontsize=9, fontweight="bold", ha="right", va="center")

    handles, labels = ax.get_legend_handles_labels()
    handles += [mpatches.Patch(facecolor=GREEN, alpha=0.25),
                mpatches.Patch(facecolor=RED, alpha=0.25)]
    labels += ["recovery basin", "escape"]
    ax.legend(handles, labels, loc="upper left", fontsize=7,
              handlelength=1.6, labelspacing=0.35,
              frameon=True, facecolor="white", edgecolor="none", framealpha=0.95)

    ax.set_xlabel(r"substrate load $s$")
    ax.set_ylabel(r"metabolite $m$")
    ax.set_ylim(0, Y_MAX)
    ax.set_title("recovery requires pushing the load below baseline,\n"
                 "not merely removing the trigger")
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
    #  The four numbers this figure ASSERTS about the model, registered so they
    #  are diffed against HEAD like a table cell and can be quoted with \num{}
    #  rather than retyped. Both panels draw out of this dict, so the drawn
    #  value and the registered value are one object and cannot drift apart.
    #  It is a NumberBook, so a second writer giving a key a different value
    #  fails at accumulation instead of last-write-wins.
    MARG = NumberBook({"theta_crit": float(thc)})
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
    # Anchored to the merge point and offset in POINTS, so it keeps its
    # clearance whatever theta_crit turns out to be.
    at_point(ax[0], (th_ok[0], 0.5 * (attr[0] + sep[0])),
             r"$\theta_{\mathrm{crit}}$"
             + f" = {MARG['theta_crit']:.0f}:" "\nthe curves merge",
             dx=16, dy=0, arrow=True, fontsize=7.5, color=MUTED,
             ha="left", va="center")
    ax[0].plot([PC.theta], [attr[np.argmin(abs(th_ok - PC.theta))]], "o", ms=5,
               color=MUTED, zorder=5)
    ax[0].plot([PC.theta], [sep[np.argmin(abs(th_ok - PC.theta))]], "o", ms=5,
               color=MUTED, zorder=5)
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"residual capacity $\theta$")
    ax[0].set_ylabel(r"metabolite $m$")
    ax[0].set_title("the two equilibria collide at the fold")
    ax[0].legend(loc="upper left")

    # --- the margin as a BAND, not a curve -----------------------------
    #
    # The separatrix sits where alpha has decayed, and alpha decays on the
    # scale of M_i, so the margin inherits M_i almost linearly where
    # theta_crit_intrinsic does not. Drawing one curve invites the number to be
    # read as a value. It is an order of magnitude, and the figure should say
    # so itself rather than leaving it to a remark.
    #
    # The sweep is the admissible one: R_0 is recalibrated at each M_i, so the
    # constraint that fixes it is never broken.
    MI_SWEEP = (150.0, 300.0, 600.0, 1200.0)     # the eight-fold range
    coarse = np.linspace(1e-3, 30000.0, 30000)
    band = []
    for mi in MI_SWEEP:
        q0 = replace(PC, Mi=float(mi))
        q0 = replace(q0, R0=R0_for_peak(q0))
        mg = []
        for th in th_ok:
            q = replace(q0, theta=float(th))
            roots, st = equilibria_nonmonotone(clearance_C_quasistatic, q,
                                               S_BASE, coarse)
            mg.append(roots[~st][0] - roots[st][0]
                      if st.sum() and (~st).sum() else np.nan)
        band.append(np.array(mg))
    band = np.array(band)
    lo = np.nanmin(band, axis=0)
    hi = np.nanmax(band, axis=0)

    ax[1].fill_between(th_ok, lo, hi, color=BLUE, alpha=0.16,
                       label=r"$M_i$ over an eight-fold range")
    ax[1].plot(th_ok, margin, color=BLUE, lw=2.0, label="nominal $M_i = 300$")
    ax[1].axvline(thc, color=MUTED, ls=":", lw=1.2)
    i = np.argmin(abs(th_ok - PC.theta))
    #  The four numbers this figure ASSERTS about the model, registered so
    #  they are diffed against HEAD like a table cell and can be quoted with
    #  \num{} rather than retyped. Both panels draw out of this dict, so the
    #  drawn value and the registered value are one object.
    MARG.update({"margin_nominal": float(margin[i]),
                 "margin_band_lo": float(lo[i]),
                 "margin_band_hi": float(hi[i])})
    ax[1].plot([PC.theta], [MARG["margin_nominal"]], "o", ms=5.5, color=MUTED,
               zorder=5)
    at_point(ax[1], (PC.theta, MARG["margin_nominal"]),
             "nominal patient\n"
             f"margin = {MARG['margin_nominal']:.0f}\n"
             f"band {MARG['margin_band_lo']:.0f}--{MARG['margin_band_hi']:.0f}",
             dx=10, dy=34, arrow=True, fontsize=7.5, color=MUTED,
             ha="left", va="bottom")
    over_x(ax[1], thc, "margin vanishes:\nsaddle-node", frac=0.97,
           color=MUTED, fontsize=7.5, ha="left")
    ax[1].set_xlabel(r"residual capacity $\theta$")
    ax[1].set_ylabel("separatrix $-$ attractor")
    ax[1].set_title("the safety margin: sign robust, magnitude to an "
                    "order of magnitude", fontsize=8)
    ax[1].legend(loc="upper left", fontsize=6.5)

    save(fig, "fig_C5_margin")
    print(f"    theta_crit_intrinsic at s = {S_BASE:.0f}: {thc:.2f}")
    print(f"    at the nominal theta = {PC.theta:.0f}: attractor "
          f"{attr[i]:.1f}, separatrix {sep[i]:.1f}, margin {margin[i]:.1f}")
    pos = np.nanmin(band)
    print(f"    THE SIGN: the margin stays positive at every theta and every "
          f"M_i in the sweep (smallest {pos:.1f})"
          if pos > 0 else
          f"    THE SIGN FLIPS somewhere in the sweep (smallest {pos:.1f})")
    print(f"    THE MAGNITUDE: at the nominal patient the band runs "
          f"{lo[i]:.0f} to {hi[i]:.0f}, a factor of {hi[i]/lo[i]:.1f}. "
          f"Quote it as an order of magnitude.")
    print("    m_bar does NOT appear anywhere in this construction. s does.")
    from common.numbers import write_numbers
    write_numbers("margin", MARG,
                  {"theta_crit": "{:.0f}", "margin_nominal": "{:.0f}",
                   "margin_band_lo": "{:.0f}", "margin_band_hi": "{:.0f}"})
    return thc


# =====================================================================
def shape_sensitivity(alt_target=500.0):
    """
    What the three unmeasured shape constants cost, on one footing.

    Remark 12 bounded K_c and stopped there. M_i and rho_0 have exactly the
    same epistemic status -- no measurement, chosen -- and M_i was not even
    written down, so Figure 16 could not be reproduced from the paper.

    Two columns, and the difference between them is the point.

      R0 re-derived   the constraint that fixes R0 is held, so only the SHAPE
                      of the alternative route moves. This is the admissible
                      sweep: a parameter set that violates the constraint is
                      not a parameter set.
      R0 held         the constraint is broken, so the peak of alpha moves
                      with the shape constant and carries theta_crit with it.

    On the admissible footing all three constants are cheap and M_i is cheaper
    than K_c. On the inadmissible one rho_0 is catastrophic and M_i is second.
    The paper's own 1.56 for K_c is an R0-re-derived number, so the two were
    never comparable.

    The margin is reported alongside, because it is a different story: the
    separatrix sits where alpha has decayed, and alpha decays on the scale of
    M_i, so the margin inherits M_i almost linearly where theta_crit does not.
    """
    out = {}
    for name, field, vals in (("K_c", "Kc", (25, 50, 100, 200, 400, 800)),
                              ("M_i", "Mi", (150, 212, 300, 424, 600, 849, 1200)),
                              ("rho_0", "rho0", (0.35, 0.4, 0.5, 0.7, 0.9, 1.0))):
        rec = []
        for v in vals:
            q = replace(PC, **{field: float(v)})
            q_cal = replace(q, R0=R0_for_peak(q, alt_target))
            tc_cal = theta_crit_intrinsic(q_cal, S_BASE)
            tc_raw = theta_crit_intrinsic(q, S_BASE)
            mg = margin_at(q_cal, PC.theta)
            rec.append((v, tc_cal, tc_raw, mg))
        out[name] = rec
    return out


def margin_at(p: ParamsC, theta, s=S_BASE):
    """Separatrix minus attractor at this theta and load, or nan if none."""
    q = replace(p, theta=float(theta))
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, q, s,
                                       np.linspace(1e-3, 3e5, 400000))
    st = np.asarray(st, dtype=bool)
    if st.sum() and (~st).sum():
        return float(roots[~st][0] - roots[st][0])
    return float("nan")


def table_shape():
    """Render the sensitivity of theta_crit and the margin to the shape constants."""
    from common.style import write_table
    d = shape_sensitivity()
    rows = []
    for name, rec in d.items():
        cal = [r[1] for r in rec if np.isfinite(r[1])]
        raw = [r[2] for r in rec if np.isfinite(r[2])]
        mgs = [r[3] for r in rec if np.isfinite(r[3])]
        n_vac = sum(1 for r in rec if not np.isfinite(r[2]))
        span = lambda xs: f"{min(xs):.0f}--{max(xs):.0f}" if xs else "---"
        fac = lambda xs: f"{max(xs) / min(xs):.2f}" if xs else "---"
        rows.append([
            f"${name.replace('_', '_{') + '}' if '_' in name else name}$".replace("$K_{c}$", "$K_c$"),
            f"{rec[0][0]:g}--{rec[-1][0]:g}",
            f"{rec[-1][0] / rec[0][0]:.0f}$\\times$",
            fac(cal), fac(raw),
            f"{n_vac:d}", fac(mgs),
        ])
    write_table("tab_C_shape",
                ["constant", "range swept", "range", "$\\thc$, $R_0$ re-derived",
                 "$\\thc$, $R_0$ held", "vacuous", "margin, $R_0$ re-derived"],
                rows,
                note="factors are max/min over the swept range; 'vacuous' counts "
                     "sweep points at which a null-genotype patient becomes "
                     "maintainable and the threshold ceases to exist")
    print("    shape-constant sensitivity of theta_crit and of the margin:")
    for name, rec in d.items():
        cal = [r[1] for r in rec if np.isfinite(r[1])]
        mgs = [r[3] for r in rec if np.isfinite(r[3])]
        print(f"      {name:>6}: theta_crit {min(cal):7.1f}..{max(cal):7.1f} "
              f"(factor {max(cal)/min(cal):5.2f})   margin {min(mgs):8.1f}..{max(mgs):8.1f} "
              f"(factor {max(mgs)/min(mgs):5.2f})")
    return rows


def linearise_at_attractor(p: ParamsC = None, s=None):
    """
    The tangent (a, b) of the parent model's total clearance at its attractor.

    Reported alongside what Model A's intercept identity would give, because
    the two are not the same object and the note previously wrote as though
    they were. a = theta g(m_bar)^2 holds only when the alternative channel is
    linear through the origin. alpha is not, and its own intercept term
    alpha(m) - m alpha'(m) is the larger part of a here by a factor of six.
    """
    p = PC if p is None else p
    s = S_BASE if s is None else s
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, p, s,
                                       np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    mbar = float(roots[st][0])
    d = 1e-4 * mbar
    C = lambda m: p.theta * g(m, p.Km) + alt_flux_C(m, p)
    b = (C(mbar + d) - C(mbar - d)) / (2 * d)
    a = C(mbar) - mbar * b
    ap = (alt_flux_C(mbar + d, p) - alt_flux_C(mbar - d, p)) / (2 * d)
    return dict(m_bar=mbar, a=a, b=b,
                a_enzymatic=p.theta * g(mbar, p.Km) ** 2,
                a_alternative=alt_flux_C(mbar, p) - mbar * ap,
                alpha=alt_flux_C(mbar, p), alpha_prime=ap)


def recover_from_ab(a, b, m_bar, p: ParamsC):
    """
    (theta, K_m) from a fitted tangent, GIVEN alpha. Unique when it exists.

    This is why the parent model has no fibre: alpha is curved, so it
    contributes to both coefficients, and subtracting its contribution leaves
    an enzymatic (a_e, b_e) from which a_e/b_e = m_bar^2/K_m fixes K_m
    outright. Inadmissibility is informative rather than a failure: it says
    the assumed alpha already accounts for more of the intercept, or more of
    the slope, than the record contains.
    """
    d = 1e-4 * m_bar
    ap = (alt_flux_C(m_bar + d, p) - alt_flux_C(m_bar - d, p)) / (2 * d)
    a_e = a - (alt_flux_C(m_bar, p) - m_bar * ap)
    b_e = b - ap
    if a_e <= 0.0:
        return None, None, "alpha already owns the intercept"
    Km = m_bar ** 2 * b_e / a_e
    if Km <= 0.0:
        return None, None, "alpha already owns the slope"
    return a_e * (Km + m_bar) ** 2 / m_bar ** 2, Km, "ok"


def identification_cost_of_Kc(values=(25, 50, 100, 150, 200, 300, 400, 800),
                              rederive=True):
    """
    What K_c costs in the theta recovered from one fitted tangent.

    The published figures for this -- 300 to 457, a factor of 1.52, and no
    admissible solution above K_c = 200 -- do not reproduce under either
    convention. With the calibration held the admissible window is two-sided
    and tight, K_c between about 100 and 300, and inside it theta moves by
    1.41. The bound being two-sided is the useful part: the record excludes
    most of the range Remark 12 sweeps.
    """
    lin = linearise_at_attractor()
    out = []
    for Kc in values:
        q = replace(PC, Kc=float(Kc))
        if rederive:
            q = replace(q, R0=R0_for_peak(q))
        th, Km, why = recover_from_ab(lin["a"], lin["b"], lin["m_bar"], q)
        out.append((Kc, th, Km, why))
    ok = [r[1] for r in out if r[1] is not None]
    print(f"    K_c and the recovered theta ({'calibration held' if rederive else 'R0 held'}):")
    for Kc, th, Km, why in out:
        print(f"      K_c={Kc:5.0f}: " + ("inadmissible -- " + why if th is None
                                          else f"theta={th:8.2f}  Km={Km:8.2f}"))
    if ok:
        print(f"      -> admissible theta spans {min(ok):.1f}..{max(ok):.1f}, "
              f"factor {max(ok)/min(ok):.2f}")
    return out


#: the perturbations the reduction's cost is quoted against, as a fold of the
#: baseline load. The cost is a FUNCTION of the perturbation and not a
#: constant, which is the whole reason the row names one: a 1.10x step returns
#: 9.9% where a 1.50x step returns 31.8%, and an earlier version of the row
#: gave four percentages and named neither.
#:
#: The third entry is the one item 11 of the September review needs. An
#: isotope-dilution measurement of whole-body propionate turnover is itself a
#: load perturbation, and a tracer is administered in tracer quantities by
#: construction, so it is a SMALL one: the reduction has to be re-read for
#: that step rather than for either of the two the registry already carries.
#: If measuring the load opens anything, this is the number that says whether
#: the measurement is readable off the reduced model at all.
PERTURBATIONS = {"": None,               # filled from ILLNESS_FOLD below
                 "_step10": 1.10,
                 "_tracer": 1.01}


def boundary_numbers():
    """
    Section 8's validity boundary, generated rather than typed.

    Every number in that paragraph was a literal, in the one paragraph where
    the paper states the limits of its own model, and one of them described
    the wrong event: "loses its stable equilibrium at 12.6 times it" names the
    SEPARATRIX at baseline load, while the stable equilibrium is lost at the
    FOLD, where attractor and separatrix collide, which Section 2.3 defines
    and which sits at the clearance maximum. Raising the load walks the two
    roots together and they merge at m = 275.7 at a load of 706.4, i.e. at
    3.4 times the attractor and not 12.6.

    THE TWO MANAGED LEVELS, which is the trap underneath it. The parent's
    reference level is its OWN attractor, 81.93 at its own theta and load; the
    local model's is m_bar = 10. The paragraph used one phrase for both and
    therefore meant two different bands by "eight times the managed level":
    655 in the 42% figure and 80 in the factor of 4.6. Both are computed here
    so the referents are visible.
    """
    from common.numbers import write_numbers
    from design.run import P as PA, M_BAR, load
    from common.models import clearance_A
    L = float(load(0))
    peak, m_peak = clearance_peak(PC)
    eq = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE,
                                np.linspace(1e-3, 3000.0, 120000))[0]
    mbar, sep = float(eq[0]), float(eq[1])
    #  the fold: raise the load until the two roots merge, which is exactly
    #  the clearance maximum. Reported as the load and the level.
    N = {
        "parent_mbar": mbar,
        "peak_over_mbar": m_peak / mbar,
        "sep_over_mbar": sep / mbar,
        "fold_m": m_peak,
        "fold_load": peak,
        "fold_load_ratio": peak / S_BASE,
        "band_parent": 8.0 * mbar,
        "band_local": 8.0 * M_BAR,
        "band_beyond_pct": 100.0 * max(
            0.0, (np.log(8.0) - np.log(m_peak / mbar)) / np.log(8.0)),
        "hold8_ratio": clearance_A(8 * M_BAR, PA) / clearance_A(M_BAR, PA),
        "clear_ratio": (clearance_A(8 * M_BAR, PA) / clearance_A(M_BAR, PA))
        / (peak / S_BASE),
        #  the structural margin: no local operating point can lie past the
        #  parent's peak, because the local marker is capped and the peak has
        #  a floor that no capacity can lower
        "local_ceiling": L / PA.ke,
        "peak_floor": float(np.sqrt(PC.Kc * PC.Mi)),
        "peak_margin": float(np.sqrt(PC.Kc * PC.Mi)) / (L / PA.ke),
    }
    #  and how often the margin reverses over this module's own prior box.
    #  Analytic first, because sqrt(Kc Mi) = 10^((u+v)/2) with u, v uniform
    #  makes P(floor < ceiling) a triangle, and a closed form is not a seed.
    a = np.log10(L / PA.ke)
    slack = 2 * a - (1.4 + 2.0)
    N["prior_floor_pct"] = 100.0 * (0.5 * max(0.0, slack) ** 2) / (1.5 * 1.1)
    #  A Monte-Carlo count of the cells that actually reverse is NOT written
    #  here, deliberately. The event is rare: at 600 draws it reads 0.14% and
    #  at 4000 it reads 0.23%, which is sampling noise on five or so events,
    #  and a generated number the artefact guard diffs against HEAD must not
    #  wander with the seed. The closed form above is an upper bound on it
    #  (the floor being below the ceiling is necessary for the peak to be),
    #  it is exact, and it is the one the paper quotes. The swept figure is
    #  recorded in REVIEW_LOG.md instead, where a measurement belongs.
    write_numbers("boundary", N,
                  {"parent_mbar": "{:.1f}", "peak_over_mbar": "{:.2f}",
                   "sep_over_mbar": "{:.1f}", "fold_m": "{:.0f}",
                   "fold_load": "{:.1f}", "fold_load_ratio": "{:.3f}",
                   "band_parent": "{:.0f}", "band_local": "{:.0f}",
                   "band_beyond_pct": "{:.0f}", "hold8_ratio": "{:.2f}",
                   "clear_ratio": "{:.1f}", "local_ceiling": "{:.1f}",
                   "peak_floor": "{:.0f}", "peak_margin": "{:.2f}",
                   "prior_floor_pct": "{:.1f}"})
    for k in sorted(N):
        print(f"    {k:22} {N[k]}")
    return N


def reduction_cost_numbers(ks=(0.1, 1.0, 10.0, 100.0)):
    """
    What the quasi-steady reduction costs OFF equilibrium, at four k.

    WHY THIS IS NOW LOAD-BEARING. Section 2.3 justifies eliminating the
    carnitine state algebraically: `c-dot = 0` holds exactly at any fixed
    point, so at equilibrium the two formulations are one equation rearranged
    and the cost is zero. Every steady-state number in the paper is computed
    at equilibrium, and while that was the whole of the paper the k row could
    bound the off-equilibrium cost loosely. It cannot now: the moment a
    transient enters the specification, the reduction's validity OFF
    equilibrium is what the transient rests on, and this row is the only thing
    bounding it.

    THE PERTURBATION, STATED, which is the point of this function. The row
    used to quote 21.6% at k = 1 without saying what produced it, and the
    figure was not reproducible: a 10% load step gives 9.91% and the model's
    own illness gives 31.79%, a factor of three apart, both with the same
    asymptotic 1/k. A bound whose two readings differ by three is not a bound.

    So the perturbation is the one the model already carries: a step in the
    substrate load to `ILLNESS_FOLD` times baseline, which is the registered
    illness amplitude, with carnitine starting at its pre-perturbation value.
    That is the largest defensible reading, and a cost bound should be the
    largest defensible reading rather than the median of several.
    """
    from common.numbers import write_numbers
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE)
    m0 = float(roots[st][0])
    folds = dict(PERTURBATIONS)
    folds[""] = ILLNESS_FOLD

    def cost(fold, k):
        s_new = S_BASE * fold
        c0 = PC.R0 / (k * (m0 / (PC.Kc + m0)
                           + PC.rho0 * (1.0 + m0 / PC.Mi)))
        red = solve_ivp(
            lambda t, y: [s_new - clearance_C_quasistatic(y[0], PC)],
            (0.0, 160.0), [m0], rtol=1e-10, atol=1e-12,
            dense_output=True, method="Radau", max_step=0.5)
        two = solve_ivp(
            lambda t, y: rhs_C(t, y, PC, lambda _: s_new, k=k),
            (0.0, 160.0), [m0, c0], rtol=1e-10, atol=1e-12,
            dense_output=True, method="Radau", max_step=0.5)
        tt = np.linspace(0.0, 160.0, 8000)
        a, b = two.sol(tt)[0], red.sol(tt)[0]
        return 100.0 * float(np.max(np.abs(a - b) / b))

    N = {"illness_fold": ILLNESS_FOLD,
         "tracer_fold": PERTURBATIONS["_tracer"],
         "step10_fold": PERTURBATIONS["_step10"]}
    suffix_k = {1.0: "", 10.0: "_k10", 100.0: "_k100", 0.1: "_k0p1"}
    for tag, fold in folds.items():
        for k in ks:
            N[f"reduction_cost{tag}{suffix_k[k]}"] = cost(fold, k)
        print(f"    a step to {fold:.2f}x baseline:")
        for k in ks:
            print(f"      k = {k:6.1f}: peak relative error in m "
                  f"{N[f'reduction_cost{tag}{suffix_k[k]}']:6.2f}%")
    fmt = {k: "{:.1f}" for k in N}
    fmt.update({"illness_fold": "{:.2f}", "tracer_fold": "{:.2f}",
                "step10_fold": "{:.2f}",
                "reduction_cost_tracer": "{:.2f}",
                "reduction_cost_tracer_k10": "{:.2f}",
                "reduction_cost_tracer_k100": "{:.3f}",
                "reduction_cost_tracer_k0p1": "{:.2f}"})
    write_numbers("reduction", N, fmt)
    return N


def collapse_check(n=12, seed=0):
    """
    Is theta_crit a function of the peak of alpha alone?

    Sample the three shape constants over their whole admissible ranges, hold
    the calibration, and see how much of theta_crit is left moving. The answer
    is a factor of about 1.4 over three constants spanning a 32-fold, a
    13-fold and a 3-fold range between them, which is the case for saying that
    the threshold's exposure is concentrated in the one quantity the clinic
    can argue about -- whether a null genotype is maintainable -- and not in
    three constants nobody has measured.
    """
    rng = np.random.default_rng(seed)
    tcs = []
    for _ in range(n):
        q = replace(PC, Kc=float(10 ** rng.uniform(1.4, 2.9)),
                    Mi=float(10 ** rng.uniform(2.0, 3.1)),
                    rho0=float(rng.uniform(0.35, 1.0)))
        q = replace(q, R0=R0_for_peak(q))
        tcs.append(theta_crit_intrinsic(q, S_BASE))
    tcs = [t for t in tcs if np.isfinite(t)]
    print(f"    collapse: {len(tcs)} random shape sets at a held calibration give "
          f"theta_crit {min(tcs):.1f}..{max(tcs):.1f}, factor {max(tcs)/min(tcs):.2f}")
    return min(tcs), max(tcs)


if __name__ == "__main__":
    print("MODEL C")
    fig_phaseline()
    fig_fold()
    fig_timeseries()
    fig_basin()
    fig_margin()
    table_shape()
    collapse_check()
    identification_cost_of_Kc()
    reduction_cost_numbers()
    boundary_numbers()
