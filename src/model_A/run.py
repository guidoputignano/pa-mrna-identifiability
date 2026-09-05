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
    tab_A_rank               rank against number of distinct levels

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
ASSAY_CV = 5.0     # typical analytical and biological variability, per cent


# =====================================================================
#  The fibre of observationally identical patients.
#
#  A band of data fixes the affine fit (a, b) of Proposition 1 and nothing
#  more. Since a = theta g(m_bar)^2 and b = theta g'(m_bar) + ke both carry
#  Km, a known (a, b) leaves a one-parameter family, not a point. These are
#  its parametrisation and its endpoints.
# =====================================================================
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
    at_point(ax[1], (m_op, clearance_A(m_op, P)), "same value\nand slope",
             dx=26, dy=-16, arrow=True, fontsize=7.5, color=MUTED,
             ha="left", va="top")
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
    ax[2].semilogx(mts, wid, color=BLUE, lw=2.0)
    ax[2].plot([M_TOL], [np.interp(M_TOL, mts, wid)], "o", ms=6, color=MUTED,
               zorder=5)
    ax[2].axvline(M_OP, color=GREY, ls=":", lw=1.0)
    over_x(ax[2], M_OP, r"  $\bar m$", frac=0.94, color=MUTED, fontsize=7)
    at_point(ax[2], (M_TOL, np.interp(M_TOL, mts, wid)),
             f"{np.interp(M_TOL, mts, wid):.0f}% at "
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
    for x in (15.0, 40.0, 160.0):
        w = margin_window(a0, b0, x)
        print(f"      m_tol = {x:5.1f}: width "
              f"{100 * (w[1] - w[0]) / (0.5 * (w[0] + w[1])):6.2f}% of centre")


# =====================================================================
def table_rank():
    rows = []
    for levels in ([10.0], [10.0, 12.0], [10.0, 12.0, 14.0],
                   [10.0, 50.0, 200.0], [10.0, 12.0, 14.0, 16.0]):
        S = np.array([sensitivity_A(m, P) for m in levels])
        rows.append([len(levels), np.linalg.matrix_rank(S, tol=1e-10),
                     "yes" if np.linalg.matrix_rank(S, tol=1e-10) == 3 else "no"])
    write_table("tab_A_rank",
                ["distinct levels of $m$", "rank of $S$", "all three identifiable"],
                rows,
                note="three distinct levels are necessary and generically sufficient")


if __name__ == "__main__":
    print("MODEL A")
    fig_channels()
    fig_conditioning()
    fig_indistinguishable()
    fig_population()
    slope = fig_load_error()
    fig_fibre()
    fig_margin_fibre()
    table_rank()
    print(f"  load-error amplification factor: {slope:.3f}")
