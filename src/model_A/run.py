"""
Model A: Michaelis-Menten enzyme, linear alternative clearance.

    dm/dt = s(t) - theta*g(m) - ke*m ,   g(m) = m/(Km+m)

Question answered: can theta be identified from clinical data?

Produces
    fig_A1_channels          the two clearance channels and the linear regime
    fig_A2_conditioning      condition number against excursion of m
    fig_A3_indistinguishable two patients with identical measurements, and the
                             response of each to a selective perturbation
    fig_A4_population        deficiency of the population system against N
    fig_A5_load_error        amplification of misspecification in s

Run:  python src/model_A/run.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares

from common.models import (ParamsA, g, clearance_A, steady_state,
                           sensitivity_A, condition_number)
from common.style import save, write_table, figsize, GREEN, RED, BLUE, AMBER, GREY, MUTED, corner, over_x, at_point


P = ParamsA()
M_OP = 10.0        # the operating point the illustrative patients are matched at
M_TOL = 40.0       # tolerable metabolite level, a clinical convention
# Illness amplitude as a MULTIPLE of the patient's own baseline load, not as an
# absolute figure. Stated absolutely it disagreed across the note: 1500 here is
# 4.65x this patient's baseline, where Model C's Figures 11 and 12 apply 1.42x
# and 1.50x and its fold sits at 1.18x. The ratio is the physiological claim,
# and each model applies it to its own parameters.
ILLNESS_FOLD = 1.50
ASSAY_CV = 5.0     # per cent, never flux units. A sigma comes from sigma_from_load().


# =====================================================================
#  The fibre of observationally identical patients.
#
#  A band of data fixes the affine fit (a, b) of Proposition 1 and nothing
#  more. Since a = theta g(m_bar)^2 and b = theta g'(m_bar) + ke both carry
#  Km, a known (a, b) leaves a one-parameter family, not a point. These are
#  its parametrisation and its endpoints.
# =====================================================================
#: every interpolated number this module draws on a figure, classified for
#: tools/check_artefacts.py:check_figure_numbers.
DRAWN = {
    "thB": "input: patient B's capacity, a setting of the comparison",
    "keB": "input: patient B's k_e, likewise",
    "theta": "input: the capacity this panel was drawn at",
    "slope": "UNREGISTERED: the load-error amplification factor, which "
             "Section 2.5 discusses and no generator owns",
    "th.max() / th.min()": "UNREGISTERED: the fold of theta consistent with "
        "(a, b) across the K_m band under noiseless data. NOT scope.tex's "
        "band_fold, which is the working band in theta",
    "WIN['mtol_window_pct']": "registered: the same object write_numbers "
                              "takes as mtol_window_pct",
    "M_TOL": "input: the tolerable-level convention this panel is drawn "
             "against",
    "100 * sf": "input: the load fraction this curve was drawn at",
    "nn": "input: the visit count this curve was drawn at",
}


def fit_ab(p: ParamsA, m_bar=M_OP):
    """The tangent intercept and slope of total clearance at m_bar."""
    gm = g(m_bar, p.Km)
    gp = p.Km / (p.Km + m_bar) ** 2
    return p.theta * gm ** 2, p.theta * gp + p.ke


def fibre(a, b, Km, m_bar=M_OP):
    """The (theta, ke) on the fibre at this Km. Inverts fit_ab exactly."""
    return a * (Km + m_bar) ** 2 / m_bar ** 2, b - a * Km / m_bar ** 2


def fibre_Km_max(a, b, m_bar=M_OP):
    """Where ke reaches zero, which bounds the fibre from above."""
    return b * m_bar ** 2 / a


def fibre_span(a, b, m_bar=M_OP):
    """
    How many fold theta spans along the fibre. Closed form, and observable.

    theta(K_m) = a(K_m + m_bar)^2/m_bar^2 on 0 < K_m <= K_max = b m_bar^2/a, so

        span = theta(K_max)/theta(0) = (1 + b m_bar / a)^2 = (C(m_bar)/a)^2 ,

    the second equality because the tangent satisfies a + b m_bar = C(m_bar).
    At steady state C(m_bar) is the RECORDED LOAD and a is the FITTED
    INTERCEPT, so both are observable and the span is computable per patient
    without knowing a single parameter.

    That turns the illustrative 169 into a per-patient diagnostic. It is the
    square of the ratio of the load the patient carries to the part of the
    clearance curve that extrapolates back to the origin, so it is large
    exactly when the operating point sits far down the saturation curve, which
    is where management puts it.
    """
    return (1.0 + b * m_bar / a) ** 2


def baseline_load(p: ParamsA = None, m_bar=M_OP):
    """The load this patient balances at the operating point."""
    return clearance_A(m_bar, P if p is None else p)


def s_peak_default(p: ParamsA = None, m_bar=M_OP, fold=ILLNESS_FOLD):
    return fold * baseline_load(p, m_bar)


def margin(theta, ke, m_tol=M_TOL, s_peak=None, Km=None):
    """theta - theta_crit, the safety margin of equation (17)."""
    Km = P.Km if Km is None else Km
    s_peak = s_peak_default() if s_peak is None else s_peak
    return theta - (s_peak - ke * m_tol) * (Km + m_tol) / m_tol


def margin_on_fibre(a, b, Km, m_tol=M_TOL, s_peak=None, m_bar=M_OP):
    """The same margin, reached through the fibre. Linear in Km."""
    return margin(*fibre(a, b, Km, m_bar), m_tol=m_tol, s_peak=s_peak, Km=Km)


def margin_window(a, b, m_tol=M_TOL, m_bar=M_OP):
    """
    The interval of s_peak on which the SIGN of the margin is undetermined.

    mu is linear in the fibre coordinate, so its sign is constant along the
    fibre unless the two endpoints disagree. They change sign at

        mu(0) = 0        =>  s_peak = a + b*m_tol
        mu(Kmax) = 0     =>  s_peak = (Kmax*S + a)/(1 + Kmax/m_tol) + b*m_tol

    with S = 2a/m_bar - a*m_tol/m_bar^2 and Kmax = b*m_bar^2/a. Between them
    the decision is not determined by the data, however wide the band. Outside
    them it is, despite theta being unknown to a factor of a hundred or more.
    """
    Kmax = fibre_Km_max(a, b, m_bar)
    S = 2 * a / m_bar - a * m_tol / m_bar ** 2
    hi = a + b * m_tol
    lo = (Kmax * S + a) / (1.0 + Kmax / m_tol) + b * m_tol
    return (lo, hi) if lo <= hi else (hi, lo)


# =====================================================================
def fig_channels():
    """The two clearance channels. In the low-m band only their sum is seen."""
    m = np.linspace(0, 600, 600)
    # Band edges named once, so the shading and its label cannot drift apart.
    band_hi, op_hi = 40.0, 25.0
    enz = P.theta * g(m, P.Km)
    alt = P.ke * m
    tot = enz + alt

    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.40))

    ax[0].axvspan(0, band_hi, color=BLUE, alpha=0.06)
    ax[0].plot(m, enz, color=GREEN, label=r"enzymatic  $\theta\,g(m)$")
    ax[0].plot(m, alt, color=AMBER, label=r"alternative  $k_e m$")
    ax[0].axhline(P.theta, color=GREEN, ls=":", lw=0.9)
    at_point(ax[0], (m[-1], P.theta), r"$\theta$", dx=-4, dy=5,
             color=GREEN, fontsize=9, ha="right", va="bottom")
    # The band is the range a stable patient's marker actually visits. Labelling
    # it "m << Km" would compare a plasma marker with an intracellular enzyme
    # constant, which is the comparison the paper declines to make.
    over_x(ax[0], band_hi * 1.4, "observed\nrange", frac=0.80, color=BLUE)
    ax[0].set_xlabel(r"metabolite $m$")
    ax[0].set_ylabel("elimination flux")
    ax[0].set_title("the two channels")
    ax[0].legend(loc="lower right")

    # A TANGENT-matched partner, not merely a point-matched one. Matching the
    # value alone leaves the two patients with different intercepts, and
    # Proposition 1 says the intercept is the discriminating quantity, so the
    # old pair illustrated a weaker degeneracy than the one being claimed.
    # This partner sits on the fibre of Remark 4: same a, same b, different Km.
    m_op = M_OP
    a0, b0 = fit_ab(P)
    thB, keB = fibre(a0, b0, 50.0)
    P2 = ParamsA(thB, 50.0, keB)

    ax[1].axvspan(0, op_hi, color=BLUE, alpha=0.06)
    ax[1].plot(m, clearance_A(m, P), color=GREEN,
               label=r"A: $\theta=3000,\ k_e=5.0$")
    ax[1].plot(m, clearance_A(m, P2), color=RED, ls="--",
               label=rf"B: $\theta={thB:.0f},\ K_m=50,\ k_e={keB:.1f}$")
    ax[1].plot([m_op], [clearance_A(m_op, P)], "o", ms=5, color=MUTED, zorder=4)
    #  the label went BELOW the point, where it landed on the x-axis label
    #  and printed over it; a small offset above put it on the curves
    #  instead. The clear air is the upper middle of the panel, so the label
    #  goes there in axes coordinates and reaches the point with an arrow.
    ax[1].annotate("same value\nand slope",
                   xy=(m_op, clearance_A(m_op, P)), xycoords="data",
                   xytext=(0.30, 0.78), textcoords="axes fraction",
                   fontsize=7.5, color=MUTED, ha="left", va="center",
                   arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7,
                                   shrinkA=2, shrinkB=3))
    over_x(ax[1], op_hi * 1.4, "operating\nregime", frac=0.97, color=BLUE,
           fontsize=7.5)
    ax[1].set_xlim(0, 300)
    ax[1].set_ylim(0, 4200)
    ax[1].set_xlabel(r"metabolite $m$")
    ax[1].set_title("identical to first order\nwhere they live")
    ax[1].legend(loc="lower right")

    save(fig, "fig_A1_channels")


# =====================================================================
def fig_conditioning():
    """Condition number against the relative excursion of m."""
    spreads = np.logspace(-2, 1.3, 60)
    conds = []
    m0 = 10.0
    for sp in spreads:
        levels = m0 * np.array([1.0, 1 + sp, 1 + 2 * sp])
        conds.append(condition_number(sensitivity_A, P, levels))
    conds = np.array(conds)

    fig, ax = plt.subplots(figsize=figsize(0.48, 0.88))
    ax.loglog(spreads * 100, conds, color=BLUE)
    ref = conds[0] * (spreads / spreads[0]) ** -2
    ax.loglog(spreads * 100, ref, color=GREY, ls="--", lw=1.0,
              label=r"slope $-2$")
    ax.axvspan(1, 10, color=RED, alpha=0.07)
    corner(ax, "clinically\nstable patient", "lower left", color=RED, fontsize=7.5)
    ax.set_xlabel("relative excursion of $m$  [%]")
    ax.set_ylabel("condition number")
    ax.set_title("conditioning degrades quadratically\nas the patient is stabilised")
    ax.legend()
    ax.grid(True, which="both", alpha=0.4)
    save(fig, "fig_A2_conditioning")


# =====================================================================
def fig_indistinguishable():
    """
    Two patients with different biology and identical measurements, and the
    response of each to a selective perturbation.

    The perturbation used is CARNITINE, which acts on the alternative channel
    k_e and is present in every patient's chart. The enzymatic channel, which
    is what the mRNA therapy would perturb, is shown alongside it only to make
    the structural point: the 2x2 system solved for (theta, k_e) has the same
    coefficient matrix whichever channel is perturbed, so its determinant
    depends on the two metabolite levels alone and not on the channel. The
    perturbation enters the right-hand side, not the matrix.
    """
    thA, keA = 3000.0, 5.0
    pA = ParamsA(thA, P.Km, keA)
    s = clearance_A(10.0, pA)
    thB = 1000.0
    keB = brentq(lambda k: clearance_A(10.0, ParamsA(thB, P.Km, k)) - s, 1e-6, 1e4)
    pB = ParamsA(thB, P.Km, keB)
    m0 = 10.0

    # accessible input: a carnitine increment on the alternative channel
    dks = np.linspace(0, 60, 200)
    mA_c = np.array([steady_state(clearance_A, ParamsA(thA, P.Km, keA + d), s)
                     for d in dks])
    mB_c = np.array([steady_state(clearance_A, ParamsA(thB, P.Km, keB + d), s)
                     for d in dks])

    # unavailable input: a therapeutic increment on the enzymatic channel
    dths = np.linspace(0, 12000, 200)
    mA_t = np.array([steady_state(clearance_A, ParamsA(thA + d, P.Km, keA), s)
                     for d in dths])
    mB_t = np.array([steady_state(clearance_A, ParamsA(thB + d, P.Km, keB), s)
                     for d in dths])

    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.40))

    ax[0].plot(dks, mA_c, color=GREEN, label=r"A: $\theta=3000,\ k_e=5.0$")
    ax[0].plot(dks, mB_c, color=RED, ls="--",
               label=rf"B: $\theta=1000,\ k_e={keB:.1f}$")
    ax[0].plot([0], [m0], "o", ms=5, color=MUTED, zorder=4)
    at_point(ax[0], (0.0, m0), "identical untreated", dx=9, dy=-5,
             color=MUTED, fontsize=7.5, ha="left", va="top")
    ax[0].set_xlabel(r"carnitine increment $\kappa u_c$ on $k_e$")
    ax[0].set_ylabel(r"steady-state $m$")
    ax[0].set_title("response to carnitine")
    ax[0].legend()

    # separation against the excursion actually achieved, both channels
    exc_c = 100 * (m0 - mA_c) / m0
    rel_c = 100 * np.abs(mA_c - mB_c) / mA_c
    exc_t = 100 * (m0 - mA_t) / m0
    rel_t = 100 * np.abs(mA_t - mB_t) / mA_t

    ax[1].plot(exc_c, rel_c, color=BLUE, label="carnitine, acts on $k_e$")
    ax[1].plot(exc_t, rel_t, color=AMBER, ls="--",
               label=r"therapy, acts on $\theta$")
    ax[1].axhline(ASSAY_CV, color=GREY, ls=":", lw=1.0)
    at_point(ax[1], (0.0, ASSAY_CV), "typical assay variability", dx=8, dy=4,
             color=GREY, fontsize=7.5, ha="left", va="bottom")
    ax[1].set_xlabel("excursion of $m$ achieved  [%]")
    ax[1].set_ylabel("separation between the patients  [%]")
    ax[1].set_xlim(0, 70)
    ax[1].set_ylim(0, 8)
    ax[1].set_title("the excursion decides")
    ax[1].legend(loc="upper left")

    save(fig, "fig_A3_indistinguishable")
    for e in (10.0, 30.0, 50.0):
        rc = np.interp(e, exc_c, rel_c)
        rt = np.interp(e, exc_t, rel_t)
        print(f"    at {e:.0f}% excursion: carnitine {rc:.3f}%, "
              f"enzymatic {rt:.3f}% separation")
    return pA, pB, s


# =====================================================================
def fig_population(N_max=30, seed=0):
    """Deficiency of the steady-state population system, and anchor recovery."""
    Ns = np.arange(1, N_max + 1)
    unknowns = Ns + 2
    equations = Ns

    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.40))

    ax[0].plot(Ns, unknowns, color=RED, label="unknowns  $N+2$")
    ax[0].plot(Ns, equations, color=GREEN, label="equations  $N$")
    ax[0].fill_between(Ns, equations, unknowns, color=RED, alpha=0.08)
    at_point(ax[0], (0.7 * N_max, 0.7 * N_max + 1), "deficiency 2\nfor every $N$",
             dx=-6, dy=-40, arrow=True, color=RED, fontsize=8,
             ha="center", va="top")
    ax[0].set_xlabel("cohort size $N$")
    ax[0].set_ylabel("count")
    ax[0].set_title("cohort size does not help")
    ax[0].legend(loc="upper left")

    # anchored recovery
    rng = np.random.default_rng(seed)
    N = 12
    thetas = np.concatenate([[0.0], rng.uniform(300, 4000, N - 1)])
    loads = rng.uniform(200, 500, N)
    ms = np.array([steady_state(clearance_A, ParamsA(t, P.Km, P.ke), s)
                   for t, s in zip(thetas, loads)])

    ke_hat = loads[0] / ms[0]                      # anchor (i): null genotype
    i = 1                                          # anchor (ii): two loads
    s2 = loads[i] * 1.6
    m2 = steady_state(clearance_A, ParamsA(thetas[i], P.Km, P.ke), s2)

    def resid(v):
        th, km = v
        return [loads[i] - ke_hat * ms[i] - th * ms[i] / (km + ms[i]),
                s2 - ke_hat * m2 - th * m2 / (km + m2)]

    th_hat, Km_hat = least_squares(resid, [1000.0, 50.0],
                                   bounds=([0, 1e-3], [1e6, 1e5])).x
    rec = [(loads[j] - ke_hat * ms[j]) * (Km_hat + ms[j]) / ms[j]
           for j in range(1, N)]

    ax[1].scatter(thetas[1:], rec, s=26, color=GREEN, zorder=3,
                  label="recovered")
    lim = [0, 4200]
    ax[1].plot(lim, lim, color=GREY, ls="--", lw=1.0, label="identity")
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel(r"true $\theta_i$")
    ax[1].set_ylabel(r"recovered $\theta_i$")
    ax[1].set_title("with the two anchors")
    ax[1].legend(loc="upper left")

    save(fig, "fig_A4_population")
    print(f"    ke recovered {ke_hat:.4f} (true {P.ke}), "
          f"Km recovered {Km_hat:.3f} (true {P.Km})")


# =====================================================================
def fig_load_error():
    """How an error in the assumed substrate load propagates into theta."""
    theta, Km, ke = 1500.0, P.Km, P.ke
    p = ParamsA(theta, Km, ke)
    m = 20.0
    s = clearance_A(m, p)
    eps = np.linspace(-0.3, 0.3, 200)
    th_hat = np.array([(s * (1 + e) - ke * m) * (Km + m) / m for e in eps])
    rel = 100 * (th_hat - theta) / theta

    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.40))
    ax[0].plot(eps * 100, rel, color=BLUE)
    ax[0].plot(eps * 100, eps * 100, color=GREY, ls="--", lw=1.0,
               label="no amplification")
    ax[0].axhline(0, color=GREY, lw=0.6)
    ax[0].axvline(0, color=GREY, lw=0.6)
    slope = np.polyfit(eps * 100, rel, 1)[0]
    ax[0].set_xlabel(r"error in assumed load $s$  [%]")
    ax[0].set_ylabel(r"resulting error in $\theta$  [%]")
    ax[0].set_title(f"at $\\theta = {theta:.0f}$: factor {slope:.2f}")
    ax[0].legend()

    # The factor is not a property of the model. It is C(m)/(theta g(m)),
    # total clearance over enzymatic clearance, so it is 1.4 only because
    # these illustrative values put most of the flux through the enzyme. It
    # diverges as theta -> 0, which is the severe end of the disease.
    ths = np.logspace(np.log10(30.0), np.log10(6000.0), 300)
    amp = np.array([clearance_A(m, ParamsA(T, Km, ke)) / (T * g(m, Km))
                    for T in ths])
    ax[1].loglog(ths, amp, color=RED)
    ax[1].plot([theta], [clearance_A(m, p) / (theta * g(m, Km))], "o",
               ms=5.5, color=MUTED, zorder=5)
    at_point(ax[1], (theta, clearance_A(m, p) / (theta * g(m, Km))),
             "left panel", dx=10, dy=10, arrow=True, fontsize=7.5, color=MUTED,
             ha="left", va="bottom")
    ax[1].axhline(1.0, color=GREY, ls=":", lw=0.9)
    ax[1].set_xlabel(r"residual capacity $\theta$")
    ax[1].set_ylabel(r"$C(m)\,/\,\theta g(m)$")
    ax[1].set_title("the factor is a flux ratio,\nand diverges as $\\theta$ falls")
    save(fig, "fig_A5_load_error")
    print(f"    amplification C(m)/(theta g(m)) at theta = "
          f"{theta:.0f}: {clearance_A(m, p) / (theta * g(m, Km)):.3f}; "
          f"at theta = 100: "
          f"{clearance_A(m, ParamsA(100.0, Km, ke)) / (100.0 * g(m, Km)):.3f}")
    return slope


# =====================================================================
def fig_fibre():
    """
    What a perfectly resolved band still leaves open.

    (a, b) known exactly is a curve in parameter space, not a point, because
    both coefficients carry Km. This is that curve.
    """
    a0, b0 = fit_ab(P)
    Kmax = fibre_Km_max(a0, b0)
    Ks = np.linspace(1e-3, Kmax, 400)
    th = np.array([fibre(a0, b0, K)[0] for K in Ks])
    kes = np.array([fibre(a0, b0, K)[1] for K in Ks])

    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.40))
    ax[0].semilogy(Ks, th, color=GREEN, lw=2.0)
    ax[0].plot([P.Km], [P.theta], "o", ms=6, color=MUTED, zorder=5)
    at_point(ax[0], (P.Km, P.theta), "the true patient", dx=-12, dy=-12,
             arrow=True, fontsize=7.5, color=MUTED, ha="right", va="top")
    ax[0].set_xlabel(r"$K_m$")
    ax[0].set_ylabel(r"$\theta$ consistent with $(a,b)$")
    ax[0].set_title(f"one band, {th.max() / th.min():.0f}-fold in "
                    r"$\theta$" "\nwith noiseless data")

    ax[1].plot(Ks, kes, color=AMBER, lw=2.0)
    ax[1].plot([P.Km], [P.ke], "o", ms=6, color=MUTED, zorder=5)
    ax[1].axhline(0, color=GREY, lw=0.8)
    at_point(ax[1], (Kmax, 0.0), r"$k_e = 0$ bounds the fibre", dx=-6, dy=8,
             arrow=True, fontsize=7.5, color=MUTED, ha="right", va="bottom")
    ax[1].plot([P.Km], [P.ke], "o", ms=6, color=MUTED, zorder=5)
    corner(ax[1], "a null genotype fixes $k_e$,\nand one point on this line\n"
                  "is the whole parameter set",
           "upper right", color=MUTED, fontsize=7.5)
    ax[1].set_xlabel(r"$K_m$")
    ax[1].set_ylabel(r"$k_e$ consistent with $(a,b)$")
    ax[1].set_title("and the anchor that closes it")
    save(fig, "fig_A6_fibre")
    print(f"    fibre: Km in (0, {Kmax:.3f}], theta in "
          f"[{th.min():.2f}, {th.max():.1f}], span {th.max() / th.min():.1f}x")
    print(f"    closed form (1 + b m_bar/a)^2 = (C(m_bar)/a)^2 = "
          f"{fibre_span(a0, b0):.2f}x   [load {baseline_load():.2f}, intercept {a0:.3f}]")
    Km_rec = (b0 - P.ke) * M_OP ** 2 / a0
    print(f"    null-genotype closure: Km = {Km_rec:.6f}, "
          f"theta = {fibre(a0, b0, Km_rec)[0]:.3f}")


# =====================================================================
def fig_margin_fibre():
    """
    The decision quantity is better conditioned than the parameter, and where
    that stops being enough.

    theta is quadratic in Km along the fibre; the margin is linear, the
    quadratic terms cancelling identically. So the SIGN of the margin is
    constant along the fibre unless its two endpoints disagree, which happens
    only for s_peak inside a computable window.
    """
    a0, b0 = fit_ab(P)
    Kmax = fibre_Km_max(a0, b0)
    base = baseline_load()
    sp0 = s_peak_default()
    Ks = np.linspace(1e-3, Kmax, 400)
    mu = np.array([margin_on_fibre(a0, b0, K, s_peak=sp0) for K in Ks])

    fig, axg = plt.subplots(2, 2, figsize=figsize(1.0, 0.76))
    ax = axg.ravel()

    # (a) the cancellation itself
    ax[0].axhspan(mu.min(), mu.max(), color=GREEN, alpha=0.07)
    ax[0].plot(Ks, mu, color=GREEN, lw=2.0)
    ax[0].plot([P.Km], [margin_on_fibre(a0, b0, P.Km, s_peak=sp0)], "o", ms=6,
               color=MUTED, zorder=5)
    ax[0].axhline(0, color=RED, ls="--", lw=1.0)
    at_point(ax[0], (P.Km, margin_on_fibre(a0, b0, P.Km, s_peak=sp0)),
             "the true patient", dx=-10, dy=-10, arrow=True, fontsize=7.5,
             color=MUTED, ha="right", va="top")
    ax[0].set_xlabel(r"$K_m$")
    ax[0].set_ylabel(r"margin $\mu$")
    ax[0].set_title("linear along the fibre,\nwhere $\\theta$ is quadratic",
                    fontsize=8)

    # (b) the sign at both endpoints, swept over illness amplitude
    lo, hi = margin_window(a0, b0)
    folds = np.linspace(1.0, 6.0, 500)
    e0 = np.array([margin_on_fibre(a0, b0, 1e-9, s_peak=f * base) for f in folds])
    e1 = np.array([margin_on_fibre(a0, b0, Kmax, s_peak=f * base) for f in folds])
    ax[1].axvspan(lo / base, hi / base, color=RED, alpha=0.12)
    ax[1].plot(folds, e0, color=GREEN, lw=1.8, label=r"$\mu$ at $K_m \to 0$")
    ax[1].plot(folds, e1, color=AMBER, lw=1.8, ls="--",
               label=r"$\mu$ at $K_m = K_{\max}$")
    ax[1].axhline(0, color=RED, ls=":", lw=1.0)
    ax[1].axvline(ILLNESS_FOLD, color=GREY, ls=":", lw=1.0)
    over_x(ax[1], ILLNESS_FOLD, "  Model C's illness", frac=0.34,
           color=MUTED, fontsize=7)
    over_x(ax[1], 0.5 * (lo + hi) / base, "sign\nundetermined", frac=0.97,
           color=RED, fontsize=7, ha="center")
    ax[1].set_xlabel("illness amplitude  [$\\times$ baseline load]")
    ax[1].set_ylabel(r"margin at the fibre ends")
    ax[1].set_title("the decision is determined\noutside one window", fontsize=8)
    ax[1].legend(loc="lower left", fontsize=6.5, frameon=True,
                 facecolor="white", edgecolor="none", framealpha=0.95)

    # (c) how wide that window is, against the convention that sets it
    mts = np.linspace(12.0, 320.0, 400)
    wid = np.array([(lambda w: 100 * (w[1] - w[0]) / (0.5 * (w[0] + w[1])))
                    (margin_window(a0, b0, x)) for x in mts])
    #  The three widths the registry's m_tol row quotes, computed HERE so the
    #  number this panel draws and the number results/numbers/local.tex holds
    #  are one object. They were two: the panel read its dot off the sweep
    #  grid with np.interp while the registered value solved margin_window at
    #  m_tol exactly, so the figure and the registry were two computations of
    #  one quantity with nothing comparing them.
    WIN = {x: 100 * (lambda w: (w[1] - w[0]) / (0.5 * (w[0] + w[1])))
                 (margin_window(a0, b0, x)) for x in (15.0, M_TOL, 160.0)}
    WIN = {"mtol_window_lo_pct": WIN[15.0], "mtol_window_pct": WIN[M_TOL],
           "mtol_window_hi_pct": WIN[160.0]}
    ax[2].semilogx(mts, wid, color=BLUE, lw=2.0)
    ax[2].plot([M_TOL], [WIN["mtol_window_pct"]], "o", ms=6, color=MUTED,
               zorder=5)
    ax[2].axvline(M_OP, color=GREY, ls=":", lw=1.0)
    over_x(ax[2], M_OP, r"  $\bar m$", frac=0.94, color=MUTED, fontsize=7)
    at_point(ax[2], (M_TOL, WIN["mtol_window_pct"]),
             f"{WIN['mtol_window_pct']:.0f}% at "
             f"$m_{{\\mathrm{{tol}}}} = {M_TOL:.0f}$",
             dx=10, dy=-4, arrow=True, fontsize=7, color=MUTED,
             ha="left", va="top")
    ax[2].set_xlabel(r"tolerable level $m_{\mathrm{tol}}$")
    ax[2].set_ylabel("window width  [% of centre]")
    ax[2].set_title("and it narrows as $m_{\\mathrm{tol}}$\napproaches $\\bar m$",
                    fontsize=8)


    # (d) two patients matched at the operating point, separating exactly
    thB = 1000.0
    keB = brentq(lambda k: clearance_A(M_OP, ParamsA(thB, P.Km, k))
                 - clearance_A(M_OP, P), 1e-6, 1e4)
    mtol = np.linspace(M_OP, 60.0, 300)
    muA = np.array([margin(P.theta, P.ke, m_tol=x) for x in mtol])
    muB = np.array([margin(thB, keB, m_tol=x) for x in mtol])
    ax[3].plot(mtol, muA, color=GREEN, lw=1.8, label=r"A: $\theta = 3000$")
    ax[3].plot(mtol, muB, color=RED, ls="--", lw=1.8,
               label=rf"B: $\theta = {thB:.0f}$")
    ax[3].axvline(M_OP, color=GREY, ls=":", lw=0.9)
    ax[3].set_xlabel(r"tolerable level $m_{\mathrm{tol}}$")
    ax[3].set_ylabel(r"margin $\mu$")
    ax[3].set_title(r"matched patients separate as "
                    r"$\Delta k_e (m_{\mathrm{tol}} - \bar m)$", fontsize=8)
    ax2 = ax[3].twinx()
    ax2.plot(mtol, muA - muB, color=MUTED, lw=1.2,
             label=r"$\mu_A - \mu_B$ (right)")
    ax2.set_ylabel(r"$\mu_A - \mu_B$", color=MUTED, fontsize=8)
    ax2.tick_params(axis="y", colors=MUTED, labelsize=7)
    h1, l1 = ax[3].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax[3].legend(h1 + h2, l1 + l2, loc="lower left", fontsize=6.5,
                 frameon=True, facecolor="white", edgecolor="none",
                 framealpha=0.95)

    save(fig, "fig_A7_margin")
    print(f"    baseline load C(m_bar) = {base:.1f}; illness at "
          f"{ILLNESS_FOLD:.2f}x = {sp0:.1f}")
    print(f"    mu over the fibre at that illness: [{mu.min():.1f}, {mu.max():.1f}], "
          f"sign constant: {bool(np.all(mu > 0) or np.all(mu < 0))}")
    print(f"    undetermined window: s_peak in ({lo:.1f}, {hi:.1f}) = "
          f"({lo / base:.2f}x, {hi / base:.2f}x) baseline, "
          f"width {100 * (hi - lo) / (0.5 * (lo + hi)):.1f}% of centre")
    #  the width at the two ends of the m_tol sweep is what the registry's
    #  m_tol row quotes as the cost of that convention, so it is generated
    #  here rather than read off this function's stdout and typed in.
    for k in ("mtol_window_lo_pct", "mtol_window_pct", "mtol_window_hi_pct"):
        print(f"      {k}: width {WIN[k]:6.2f}% of centre")
    from common.numbers import write_numbers
    write_numbers("local", WIN,
                  {"mtol_window_lo_pct": "{:.1f}", "mtol_window_pct": "{:.1f}",
                   "mtol_window_hi_pct": "{:.0f}"})


# =====================================================================
#  What the design would have to deliver.
#
#  Equation (7) is derived in the note and never evaluated. These are the
#  quantities it needs, all of them computable from fitted coefficients.
# =====================================================================
def rel_err_theta(sigma, n, Delta, a=None, m_bar=M_OP):
    """
    Relative precision of theta_hat, equation (7).

    theta_hat = a_hat/g(m_bar)^2 is a_hat times a constant, so
    sd(theta_hat)/theta = sd(a_hat)/a exactly. The g(m_bar)^-2 that appears in
    (7) is not a separate unknown: a = theta g(m_bar)^2, so
    (sigma/theta)*g(m_bar)^-2 IS sigma/a. The RELATIVE precision is therefore
    identified along the whole fibre even though theta is not.
    """
    a = fit_ab(P)[0] if a is None else a
    return (sigma / a) * np.sqrt(12.0 / n) * (m_bar / Delta)


def delta_for_precision(target, sigma, n, a=None, m_bar=M_OP):
    """Band width needed for a given relative precision on theta."""
    a = fit_ab(P)[0] if a is None else a
    return (sigma / a) * np.sqrt(12.0 / n) * m_bar / target


# For a quadratic least-squares fit on a band of full width D carrying n
# points, Var(c2_hat) = 180 sigma^2/(n D^4) in the uniform-design limit, so
# the design constant is sqrt(180) = 13.42 rather than an order of magnitude.
CURVATURE_CONST = np.sqrt(180.0)


def delta_for_curvature(sigma, n, z=2.0, p: ParamsA = None, m_bar=M_OP):
    """
    Band width at which the curvature of C is resolvable at z standard
    errors, which is what makes condition (3) checkable against a record.
    The true quadratic coefficient is |C''(m_bar)/2| = theta*Km/(Km+m_bar)^3.
    """
    p = P if p is None else p
    c2 = p.theta * p.Km / (p.Km + m_bar) ** 3
    return np.sqrt(z * CURVATURE_CONST * sigma / (np.sqrt(n) * c2))


def sigma_from_load(frac, p: ParamsA = None, m_bar=M_OP):
    """
    sigma in flux units, from a relative error in the recorded substrate load.

    The ordinate of the fit is the load, so sigma belongs on that scale. An
    earlier version passed ASSAY_CV, documented as a PERCENTAGE, as if it were
    absolute flux units: 5 units against a baseline of 322.7 is a load known to
    1.55%, which the note says throughout is not available. Nothing below
    quotes a single sigma.
    """
    return frac * baseline_load(p, m_bar)


# Relative errors in the recorded load worth reporting the criteria over. The
# note says s is "partially known" and never bounds it, so the conclusion is
# presented as a function of this rather than at a point.
SIGMA_FRACS = (0.02, 0.05, 0.10, 0.20)


# =====================================================================
#  The exact variance of the intercept, and the floor it carries.
#
#  rel_err_theta above keeps only the extrapolation term of Var(a_hat). The
#  term it drops is 1/n:
#
#      Var(a_hat) = sigma^2 [ 1/n + m_bar^2 / S ] ,   S = sum_i (m_i - m_bar)^2
#
#  and that term is the BINDING one in exactly the wide-band regime the note
#  recommends buying. Dropping it is why the design table could quote band
#  widths of 460%, 920% and 1841% of m_bar as though they were achievable.
#
#  They are not achievable, for a reason that needs no arithmetic. A band is
#  [m_bar - Delta/2, m_bar + Delta/2] and the metabolite cannot be negative,
#  so Delta <= 2 m_bar, i.e. m_bar/Delta >= 1/2. Every entry above 200% of
#  m_bar describes a band whose lower edge is below zero. Substituting the
#  bound into (7) gives a floor on the precision that no width can beat:
#
#      any design at all   sd(theta)/theta >=        (sigma/a) / sqrt(n)
#      n evenly spaced     sd(theta)/theta >= c(n) * (sigma/a) / sqrt(n)
#
#  with c(n) = sqrt(1 + 3(n-1)/(n+1)), which is 1.77 at n = 6 and rises to 2.
#  Both are attained only in the limit where the band's foot reaches m = 0,
#  the widest band the positivity of the metabolite allows. (The note's own
#  form, evaluated at its ceiling m_bar/Delta = 1/2, gives sqrt(3) = 1.73,
#  which is the size of the error it makes by dropping 1/n.)
#
#  The second follows from a closed form worth recording, since it says what
#  the best possible use of a given excursion is. For n samples split between
#  the two ends of [L, H] in the optimal proportion,
#
#      sd(a_hat) = sigma (H + L) / ( sqrt(n) (H - L) ) ,
#
#  which is (7) with sqrt(12) replaced by 2 and m_bar read as the band's
#  midpoint. The uniform design is therefore sqrt(3) = 1.73 times worse than
#  the best use of the same excursion, and neither can escape the floor.
#
#  The consequence for the paper's question is a change of kind. The design
#  requirement splits into two pieces that trade against nothing:
#
#    * a floor on n, the number of instrument events, set by sigma alone;
#    * given n above that floor, a requirement on the excursion.
#
#  Below the floor on n, no excursion whatever delivers the precision, and the
#  right answer to "how wide a band is needed" is that there is no such band.
# =====================================================================
def rel_err_exact(levels, sigma, a=None):
    """
    sd(theta_hat)/theta for a given design, from the exact Var(a_hat).

    `levels` are the metabolite levels the n observations are taken at. No
    narrow-band approximation and no assumption that they are evenly spread.
    """
    a = fit_ab(P)[0] if a is None else a
    m = np.asarray(levels, dtype=float)
    n = m.size
    mb = m.mean()
    S = float(((m - mb) ** 2).sum())
    if S <= 0.0:
        return np.inf
    return sigma * np.sqrt(1.0 / n + mb ** 2 / S) / a


def sd_a_two_point(L, H, n, sigma):
    """
    The closed form above: the best any design supported on [L, H] can do.

    Derived by minimising 1/n + m_bar^2/S over the fraction placed at H. With
    r = L/(H-L) the minimum of m_bar^2 n/S is 4r(1+r), so the bracket is
    (1+2r)^2 and the square root is exact.
    """
    return sigma * (H + L) / (np.sqrt(n) * (H - L))


def floor_constant(n, design="uniform"):
    """
    The constant c in sd(theta)/theta >= c*sigma/(a sqrt n), by design.

    "optimal"  c = 1. The 1/n term of Var(a_hat) on its own, which no design
               of any kind escapes.
    "uniform"  n evenly spaced samples with the band's foot at m = 0. Then
               S/n = D^2 (n+1)/(12(n-1)) with D = 2 m_bar, so m_bar^2 n/S is
               3(n-1)/(n+1) and c = sqrt(1 + 3(n-1)/(n+1)). This rises to 2 as
               n grows and is 1.77 at n = 6: a finite evenly spaced design
               beats its own continuous limit, because the endpoints carry
               more leverage than a uniform density does.
    """
    if design == "optimal":
        return 1.0
    return float(np.sqrt(1.0 + 3.0 * (n - 1) / (n + 1)))


def rel_err_floor(sigma, n, design="uniform", a=None):
    """
    The precision no width of band can beat, at this sigma and this n.

    Attained only in the limit where the band's foot reaches m = 0, which is
    the widest band the positivity of the metabolite allows.
    """
    a = fit_ab(P)[0] if a is None else a
    return floor_constant(n, design) * sigma / (a * np.sqrt(n))


def n_for_precision(target, sigma, design="uniform", a=None, n_max=100000):
    """
    Instrument events needed before `target` is reachable at any width.

    The uniform constant depends on n, so this is a search rather than a
    formula. Below the value returned there is no band that delivers the
    precision, and the question "how wide a band" has no answer.
    """
    a = fit_ab(P)[0] if a is None else a
    n = 1
    while n < n_max and rel_err_floor(sigma, n, design, a) > target:
        n += 1
    return n


def band_for_precision(target, sigma, n, m_hi=M_OP, a=None):
    """
    How far down a uniform band from m_hi has to reach for `target`.

    Returns None when the target is below the floor at this n, which is the
    case the old table reported as a band width instead.
    """
    f = lambda lo: rel_err_exact(np.linspace(lo, m_hi, int(n)), sigma, a) - target
    if f(1e-12) > 0.0:
        return None
    return float(brentq(f, 1e-12, m_hi * (1.0 - 1e-9)))


def top_for_curvature(sigma, n, z=2.0, m_lo=None, p: ParamsA = None,
                      m_bar=M_OP, hi_cap=1e4):
    """
    How far UP the band has to reach for the curvature of C to be resolvable.

    The curvature criterion is the one that makes condition (3) checkable
    against a record, and unlike the precision criterion it is a statement
    about the top of the excursion: the bending is at high m. Stated as a
    multiple of the operating point rather than as a band width, because the
    top of the excursion is what the safety constraint bounds and what a
    clinician can be asked for.
    """
    p = P if p is None else p
    m_lo = 0.5 * m_bar if m_lo is None else m_lo
    c2 = p.theta * p.Km / (p.Km + m_bar) ** 3

    def sd_c2(hi):
        m = np.linspace(m_lo, hi, int(n))
        X = np.vstack([np.ones_like(m), m, m ** 2]).T
        return sigma * np.sqrt(np.linalg.inv(X.T @ X)[2, 2])

    f = lambda hi: c2 / sd_c2(hi) - z
    if f(hi_cap) < 0.0:
        return None
    return float(brentq(f, m_bar * 1.0001, hi_cap))


# =====================================================================
#  The noise budget: which term binds, and what a marker panel buys.
#
#  At steady state the ordinate of the fit is the recorded load and the
#  abscissa is the assayed metabolite, so an abscissa error of sigma_m acts
#  like an ordinate error of b*sigma_m. Section 3.7 asserted that assay
#  variability is the smaller of the two and never checked it. It is not:
#  b*m_bar and the baseline load differ by under 8% at these values, so the
#  crossover sits at an assay CV of about 1.08 times the relative error in the
#  recorded load, and the assay term dominates at every load recorded better
#  than roughly its own CV.
# =====================================================================
def sigma_total(load_frac, assay_cv=None, channels=1, p: ParamsA = None,
                m_bar=M_OP):
    """
    The ordinate noise of the steady-state fit, both terms.

    `channels` is the number of markers observed as scaled copies of one
    latent pool. Independent assay noise averages as sqrt(channels); the load
    error is COMMON to every channel and does not average at all. That
    asymmetry is the whole answer to whether a marker panel is affordable.
    """
    p = P if p is None else p
    cv = (ASSAY_CV / 100.0) if assay_cv is None else assay_cv
    b = fit_ab(p, m_bar)[1]
    return float(np.hypot(load_frac * baseline_load(p, m_bar),
                          b * cv * m_bar / np.sqrt(channels)))


def assay_load_crossover(p: ParamsA = None, m_bar=M_OP):
    """
    The assay CV at which the two terms are equal, per unit load error.

    Returns k such that the terms balance at CV = k * load_frac.
    """
    p = P if p is None else p
    b = fit_ab(p, m_bar)[1]
    return float(baseline_load(p, m_bar) / (b * m_bar))


def anchor_factor(theta=None, ke=None, Km=None, m_bar=M_OP):
    """
    How a relative error in the null patient's load reaches Km and theta.

    |dln Km / dln ke|    = ke (Km + m_bar)^2 / (theta Km)
    |dln theta / dln ke| = 2 ke (Km + m_bar) / theta

    Both scale as 1/theta, so the anchor ATTENUATES for a mild patient and
    AMPLIFIES for a severe one. The crossings are at theta = ke(Km+m_bar)^2/Km
    and theta = 2 ke (Km + m_bar). This is the same shape as the load-error
    amplification of Figure 6: a ratio, not a constant.
    """
    theta = P.theta if theta is None else theta
    ke = P.ke if ke is None else ke
    Km = P.Km if Km is None else Km
    return (ke * (Km + m_bar) ** 2 / (theta * Km),
            2.0 * ke * (Km + m_bar) / theta)


def anchor_crossover(ke=None, Km=None, m_bar=M_OP):
    """The theta at which each anchor factor passes 1."""
    ke = P.ke if ke is None else ke
    Km = P.Km if Km is None else Km
    return ke * (Km + m_bar) ** 2 / Km, 2.0 * ke * (Km + m_bar)


def eiv_fit(sigma_m, Delta, n=12, reps=1500, seed=0, p: ParamsA = None,
            m_bar=M_OP, m_tol=M_TOL, constrain=True):
    """
    Errors in the abscissa, with admissibility imposed or not.

    Unconstrained least squares drives ke_hat = b_hat - a_hat*Km/m_bar^2
    negative, which is not a physical parameter set. An earlier version read a
    sign flip in the margin off such fits; it is an artefact. What survives is
    that the fit ANNOUNCES its own failure: ke_hat < 0 is checkable without
    knowing the truth.

    Returns a dict of means, plus the fraction of fits that are inadmissible
    and the fraction whose implied fibre excludes the true theta.
    """
    from scipy.optimize import lsq_linear
    p = P if p is None else p
    rng = np.random.default_rng(seed)
    sp = s_peak_default(p, m_bar)
    mt = np.linspace(m_bar - Delta / 2, m_bar + Delta / 2, n)
    s = clearance_A(mt, p)
    mu, th_r, neg, out = [], [], [], []
    for _ in range(reps):
        mo = mt + rng.normal(0.0, sigma_m, n)
        X = np.vstack([np.ones(n), mo]).T
        ah, bh = np.linalg.lstsq(X, s, rcond=None)[0]
        ke = bh - ah * p.Km / m_bar ** 2
        neg.append(ke < 0.0)
        hi = (ah + bh * m_bar) ** 2 / ah if ah > 0 else np.inf
        out.append(p.theta > hi or p.theta < ah)
        if constrain:
            # b = ke + a*Km/m_bar^2 with ke >= 0, so the design matrix is
            # rebuilt in (a, ke) and ke is bounded below at zero.
            Z = np.vstack([np.ones(n) + mo * p.Km / m_bar ** 2, mo]).T
            ah, ke = lsq_linear(Z, s, bounds=([-np.inf, 0.0],
                                              [np.inf, np.inf])).x
        th = ah / g(m_bar, p.Km) ** 2
        th_r.append(th / p.theta)
        mu.append(margin(th, ke, m_tol=m_tol, s_peak=sp, Km=p.Km))
    return dict(mu=float(np.mean(mu)),
                mu_true=margin(p.theta, p.ke, m_tol=m_tol, s_peak=sp, Km=p.Km),
                theta_ratio=float(np.mean(th_r)),
                frac_inadmissible=float(np.mean(neg)),
                frac_theta_outside=float(np.mean(out)))


def eiv_margin_bias(sigma_m, Delta, n=12, reps=2000, seed=0, p: ParamsA = None,
                    m_bar=M_OP, m_tol=M_TOL):
    """
    Errors in the ABSCISSA, which is what the steady-state design has.

    At steady state s = C(m), with s the recorded load and m the assayed
    metabolite, so the noise sits on x. Classical attenuation pulls b_hat down
    and pushes a_hat up. The two do not cancel in the margin: ke_hat =
    b_hat - a_hat*Km/m_bar^2 collapses, theta_crit rises, and mu comes out
    UNDERestimated. Returns (mean mu_hat, true mu, mean theta_hat/theta).
    """
    p = P if p is None else p
    rng = np.random.default_rng(seed)
    sp = s_peak_default(p, m_bar)
    mt = np.linspace(m_bar - Delta / 2, m_bar + Delta / 2, n)
    s = clearance_A(mt, p)
    mus, ths = [], []
    for _ in range(reps):
        mo = mt + rng.normal(0.0, sigma_m, n)
        X = np.vstack([np.ones(n), mo]).T
        ah, bh = np.linalg.lstsq(X, s, rcond=None)[0]
        th = ah / g(m_bar, p.Km) ** 2
        ke = bh - ah * p.Km / m_bar ** 2
        mus.append(margin(th, ke, m_tol=m_tol, s_peak=sp, Km=p.Km))
        ths.append(th / p.theta)
    return (float(np.mean(mus)),
            margin(p.theta, p.ke, m_tol=m_tol, s_peak=sp, Km=p.Km),
            float(np.mean(ths)))


def anchor_propagation(rel_s, reps=20000, seed=1, a=None, b=None, m_bar=M_OP):
    """
    The strengthened Proposition 5, with a variance.

    A null patient has theta = 0, so s_j = ke*m_j and a relative error in the
    recorded load enters ke_hat 1:1. Through Km = (b-ke)*m_bar^2/a and
    theta = a*(Km+m_bar)^2/m_bar^2 it reaches theta ATTENUATED, not amplified.
    """
    if a is None:
        a, b = fit_ab(P)
    rng = np.random.default_rng(seed)
    ke = P.ke * (1.0 + rng.normal(0.0, rel_s, reps))
    Km = (b - ke) * m_bar ** 2 / a
    ok = Km > 0
    th = a * (Km[ok] + m_bar) ** 2 / m_bar ** 2
    return float(np.std(Km[ok]) / P.Km), float(np.std(th) / P.theta)


# =====================================================================
def fig_design():
    """What the design would have to deliver, and what it does deliver."""
    a0, b0 = fit_ab(P)
    n_ss = 6.0                      # quasi-exogenous events, not clinic visits
    fig, axg = plt.subplots(2, 2, figsize=figsize(1.0, 0.78))
    ax = axg.ravel()

    # (a) the precision against the band, with the FLOOR restored and the
    #     region past Delta = 2 m_bar marked as what it is: not a band.
    #     Drawn on the exact Var(a_hat), not on (7)'s narrow-band form, since
    #     the two part company exactly where the note recommends going.
    fr = np.linspace(0.02, 1.999, 300)          # Delta/m_bar, up to the ceiling
    cols = [GREY, GREEN, AMBER, RED]
    for sf, col in zip(SIGMA_FRACS, cols):
        sg = sigma_total(sf)
        ys = [100 * rel_err_exact(
                  np.linspace(M_OP * (1 - f / 2), M_OP * (1 + f / 2), int(n_ss)), sg)
              for f in fr]
        ax[0].loglog(100 * fr, ys, color=col, lw=1.7, label=f"{100 * sf:.0f}% of load")
        ax[0].axhline(100 * rel_err_floor(sg, n_ss), color=col, ls=":", lw=0.9)
    ax[0].axhline(20.0, color=RED, ls="--", lw=1.0)
    ax[0].axvspan(200.0, 400.0, color=GREY, alpha=0.18, hatch="//", lw=0)
    ax[0].set_xlim(2.0, 400.0)
    corner(ax[0], "past $200\\%$ the band's\nfoot is below zero", "lower right",
           color=MUTED, fontsize=6.5)
    ax[0].set_xlabel(r"band width $\Delta$  [\% of $\bar m$]")
    ax[0].set_ylabel(r"rel. error of $\hat\theta$  [\%]")
    ax[0].set_title(r"the floor $(\sigma/a)\sqrt{3/n}$, dotted,""\nwhich no width beats",
                    fontsize=8)
    ax[0].legend(loc="lower left", fontsize=6, title=r"$\sigma$ as", title_fontsize=6)

    # (b) the curvature criterion, restated as the TOP of the excursion, which
    #     is the quantity a safety constraint bounds and a clinician can be
    #     asked for. A band width was the wrong unit: it hid the fact that the
    #     bending of C is at high m, so the band has to run up, not down.
    sfr = np.linspace(0.005, 0.25, 60)
    for nn, col, ls in [(3, GREY, ":"), (6, GREEN, "-"), (20, BLUE, "--")]:
        ys = [top_for_curvature(sigma_total(x), nn) for x in sfr]
        ax[1].plot(100 * sfr, [np.nan if y is None else y / M_OP for y in ys],
                   color=col, ls=ls, lw=1.7, label=f"$n = {nn}$")
    ax[1].axhline(1.0, color=BLUE, alpha=0.6, lw=1.2)
    corner(ax[1], "the instruments never go above $\\bar m$", "lower right",
           color=BLUE, fontsize=6.5)
    ax[1].set_xlabel(r"error in the recorded load  [\%]")
    ax[1].set_ylabel(r"top of the excursion  [$\times\,\bar m$]")
    ax[1].set_title("condition (3) is checkable\nonly above this", fontsize=8)
    ax[1].legend(loc="upper left", fontsize=6.5)

    # (c) the abscissa-noise fit announces its own failure
    ratios = np.linspace(0.02, 0.5, 18)
    D = 0.5 * M_OP
    free = [eiv_fit(r * D, D, reps=400, seed=7, constrain=False) for r in ratios]
    cons = [eiv_fit(r * D, D, reps=400, seed=7, constrain=True) for r in ratios]
    ax[2].plot(100 * ratios, [q["mu"] for q in free], color=RED, ls=":", lw=1.6,
               label=r"$\hat\mu$, fit unconstrained")
    ax[2].plot(100 * ratios, [q["mu"] for q in cons], color=RED, lw=2.0,
               label=r"$\hat\mu$, $k_e \geq 0$ imposed")
    ax[2].axhline(cons[0]["mu_true"], color=GREEN, ls="--", lw=1.5,
                  label=r"true $\mu$")
    ax[2].axhline(0.0, color=GREY, lw=0.8)
    ax2 = ax[2].twinx()
    ax2.plot(100 * ratios, [100 * q["frac_inadmissible"] for q in free],
             color=MUTED, lw=1.2)
    ax2.set_ylabel(r"fits with $\hat k_e < 0$  [%]", color=MUTED, fontsize=7)
    ax2.tick_params(axis="y", colors=MUTED, labelsize=7)
    ax[2].set_xlabel(r"$\sigma_m / \Delta$  [%]")
    ax[2].set_ylabel(r"estimated margin $\hat\mu$")
    ax[2].set_title("the sign flip was an inadmissible\nfit, and it is detectable",
                    fontsize=8)
    ax[2].legend(loc="lower left", fontsize=6)

    # (d) the anchor factor is a ratio, and it crosses 1
    ths = np.logspace(np.log10(80.0), np.log10(6000.0), 300)
    fk = [anchor_factor(theta=x)[0] for x in ths]
    ft = [anchor_factor(theta=x)[1] for x in ths]
    ax[3].loglog(ths, ft, color=GREEN, lw=2.0, label=r"into $\theta$")
    ax[3].loglog(ths, fk, color=AMBER, ls="--", lw=1.7, label=r"into $K_m$")
    ax[3].axhline(1.0, color=GREY, ls=":", lw=1.0)
    x_km, x_th = anchor_crossover()
    for xv, col in ((x_km, AMBER), (x_th, GREEN)):
        ax[3].axvline(xv, color=col, ls=":", lw=0.9)
    ax[3].plot([P.theta], [anchor_factor()[1]], "o", ms=5.5, color=MUTED, zorder=5)
    ax[3].plot([300.0], [anchor_factor(theta=300.0)[1]], "o", ms=5.5,
               color=RED, zorder=5)
    at_point(ax[3], (300.0, anchor_factor(theta=300.0)[1]),
             "Model C's $\\theta$", dx=8, dy=6, arrow=True, fontsize=7,
             color=RED, ha="left", va="bottom")
    ax[3].set_xlabel(r"residual capacity $\theta$")
    ax[3].set_ylabel("amplification of load error")
    ax[3].set_title("the anchor amplifies for\nthe severe patient", fontsize=8)
    ax[3].legend(loc="lower left", fontsize=6.5)

    save(fig, "fig_A8_design")
    print("    (7) evaluated, and the floor the old evaluation dropped:")
    for sf in SIGMA_FRACS:
        sg = sigma_total(sf)
        d_old = delta_for_precision(0.20, sigma_from_load(sf), n_ss) / M_OP
        print(f"      sigma = {100 * sf:4.0f}% of load: old table said "
              f"Delta = {100 * d_old:6.0f}% of m_bar"
              f"{'  -- NOT A BAND' if d_old > 2.0 else ''}; "
              f"floor at n=6 is {100 * rel_err_floor(sg, n_ss):6.0f}%, "
              f"so 20% needs n >= {n_for_precision(0.20, sg)}")
    q = eiv_fit(0.20 * D, D, reps=1500, constrain=False)
    r = eiv_fit(0.20 * D, D, reps=1500, constrain=True)
    print(f"    EIV at sigma_m/Delta = 20%: unconstrained mu {q['mu']:.0f}, "
          f"admissible mu {r['mu']:.0f}, true {r['mu_true']:.0f}; "
          f"{100 * q['frac_inadmissible']:.0f}% of free fits have ke < 0")
    print(f"    anchor factor into theta: {anchor_factor()[1]:.2f} at theta=3000, "
          f"{anchor_factor(theta=300.0)[1]:.2f} at theta=300; "
          f"crosses 1 at theta = {x_th:.0f}")


# =====================================================================
#  table_rank() used to write results/tables/tab_A_rank.tex here, on every
#  build, and no document ever showed it. Its claim -- the sensitivity matrix
#  reaches full rank at three distinct marker levels and not before -- is
#  Proposition 1's, and it is now asserted in
#  tests/test_three_distinct_levels_identify_the_local_model rather than
#  printed into a table nobody reads.


def table_design():
    """
    What the design has to deliver, with the floor restored.

    This replaces the table that was typed into section 3.7. Two things
    changed and both change the answer. The 1/n term of Var(a_hat) is back, so
    there is a precision no width of band can beat; and the assay term
    b*sigma_m is in the noise budget alongside the load term, where the old
    table had only the load.
    """
    n0 = 6
    cv = ASSAY_CV / 100.0
    rows = []
    for f in SIGMA_FRACS:
        s_load = f * baseline_load()
        s_assay = fit_ab(P)[1] * cv * M_OP
        s_tot = sigma_total(f)
        dom = "assay" if s_assay > s_load else "load"
        floor = rel_err_floor(s_tot, n0)
        n20 = n_for_precision(0.20, s_tot)
        n_use = 2 * n20
        lo = band_for_precision(0.20, s_tot, n_use)
        hi = top_for_curvature(s_tot, n0)
        rows.append([
            f"{100 * f:.0f}\\%",
            f"{s_load:.1f}",
            f"{s_assay:.1f}",
            dom,
            f"{100 * floor:.0f}\\%",
            f"{n20:d}",
            "---" if lo is None else f"{lo / M_OP:.2f}$\\,\\bar m$ at $n={n_use}$",
            "none" if hi is None else f"{hi / M_OP:.1f}$\\,\\bar m$",
        ])
    write_table(
        "tab_A_design",
        ["load error", "$\\sigma_{\\mathrm{load}}$", "$b\\sigma_m$", "dominant",
         f"floor at $n={n0}$", "$n$ for $20\\%$", "band foot", "excursion top"],
        rows,
        note=("floor is 2 (sigma/a)/sqrt(n), which no band width beats; "
              "band foot is where a uniform band from m_bar must reach for 20% "
              "at twice the n in the previous column; excursion top is where it "
              f"must reach for the curvature of C at two standard errors, n = {n0}"))
    print("    the design, with the floor restored:")
    for f, r in zip(SIGMA_FRACS, rows):
        print(f"      load to {100 * f:4.0f}%: sigma_load {r[1]:>7}, b*sigma_m {r[2]:>7}"
              f" ({r[3]} dominates); floor {r[4]:>5} at n={n0}; n for 20% = {r[5]:>4};"
              f" band foot {r[6]:>14}; curvature needs m up to {r[7]}")
    print("    crossover: the assay term dominates above CV = "
          f"{assay_load_crossover():.2f} x the load error")
    for J in (1, 3, 4):
        print(f"      marker panel J={J}: sigma_total at a 2% load = "
              f"{sigma_total(0.02, channels=J):.2f}, at 10% = "
              f"{sigma_total(0.10, channels=J):.2f}")
    return rows


if __name__ == "__main__":
    print("MODEL A")
    fig_channels()
    fig_conditioning()
    fig_indistinguishable()
    fig_population()
    slope = fig_load_error()
    fig_fibre()
    fig_margin_fibre()
    fig_design()
    table_design()
    from common.params import write_params_table
    write_params_table()                 # full, for the supplement
    write_params_table(compact=True)     # main text, cost bounds only
    print(f"  load-error amplification factor: {slope:.3f}")
