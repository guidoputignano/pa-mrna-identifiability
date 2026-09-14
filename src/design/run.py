"""
The design: what the excitation condition asks for, and where it can be bought.

This module produces the part of the note that terminates in a data
specification. Nothing here is a result about the disease; everything is a
statement about what a record would have to contain.

    fig_D1_excitation   the excitation functional, its optimum, and the two
                        degenerate asymptotes
    fig_D2_collapse     an augmented-state estimator run under sparse and
                        dense sampling, with P_theta_theta over time and the
                        (theta, K_m) covariance ellipse
    fig_D3_decision     the operating characteristic of the decision: the two
                        asymmetric sign errors on mu, against assay CV,
                        sampling interval and excursion
    tab_D_excursion     excursion available from each source against the
                        excursion the condition requires
    tab_D_spec          the one-page data specification

Run:  python src/design/run.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from common.models import ParamsA, g, clearance_A
from common.excitation import (phi, excitation, excitation_max, fold_range_for,
                               optimal_band, sd_log_theta, band_levels,
                               fold_range_for_precision,
                               top_of_excursion_for_precision)
from common.style import (save, write_table, figsize, GREEN, RED, BLUE, AMBER,
                          GREY, MUTED, INK, corner, over_x, at_point)

P = ParamsA()
M_BAR = 10.0
M_TOL = 40.0
ASSAY_CV = 0.05          # a fraction, never flux units
LOAD_FRAC = 0.05         # relative error in the recorded load


# --------------------------------------------------------------------------
#  The therapeutic schedule. These are the two constants the note names and
#  never assigns, and they are assigned here because 1.8 cannot be run without
#  them.
# --------------------------------------------------------------------------
DOSE_INTERVAL = 14.0     # days. The biweekly schedule selected for mRNA-3927
                         # by model-informed dose selection. ASSUMED here only
                         # in the sense that the note adopts the trial's
                         # schedule; the schedule itself is published.
LAMBDA = np.log(2.0) / 3.5   # /day. ASSUMED. A 3.5-day half-life on the
                         # therapeutic increment, so that the increment has
                         # decayed sixteen-fold by the next dose and the
                         # inter-dose recovery is a real excursion rather than
                         # a ripple. Nothing measures it. Its role is to set
                         # how much excursion a fixed schedule creates, and
                         # fig_D2 reports the answer as a function of it.
D_AMP = 6000.0           # ASSUMED. Chosen so the marker falls by about half at
                         # peak effect, which is the order of the reductions
                         # reported in the phase 1/2 interim analysis.


def delta(t):
    """The therapeutic increment from a biweekly schedule of decaying doses."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.zeros_like(t)
    for k in range(0, 200):
        tk = k * DOSE_INTERVAL
        on = t >= tk
        out[on] += D_AMP * np.exp(-LAMBDA * (t[on] - tk))
    return out if out.size > 1 else float(out[0])


def load(t, base=None):
    """The recorded component of the load. Constant here; L_0 carries the rest."""
    base = clearance_A(M_BAR, P) if base is None else base
    return base


#: day-to-day variability of ordinary intake, as a coefficient of variation.
#: By Proposition prop:rank a PERFECTLY STEADY load is the exactly singular
#: case, so this is not noise to be averaged away -- it is what closes the
#: fibre, and a dietary record smoothed to a weekly mean destroys it.
LOAD_CV = 0.10


def design_record(times, load_cv=None, seed=3):
    """
    Levels and increments for a design, with intake varying day to day.

    Returns (levels, increments). The levels solve the balance at the
    per-sample recorded load, so a varying load moves the marker as it
    physically would rather than being added as measurement error.
    """
    load_cv = LOAD_CV if load_cv is None else load_cv
    times = np.asarray(times, dtype=float)
    inc = np.asarray(delta(times), dtype=float)
    Lb = load(0.0)
    Ls = (Lb * (1.0 + load_cv * np.random.default_rng(seed).standard_normal(times.size))
          if load_cv else np.full(times.size, Lb))
    B = P.ke * P.Km + P.theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
    return lv, inc


# ==========================================================================
def fig_excitation():
    """
    The excitation functional, the two degeneracies, and where the optimum is.

    The left panel is the statement; the right panel is the price of being in
    the wrong place, which is the one that decides a design.
    """
    Km = P.Km
    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.42))

    # (left) phi(m) = K_m/(K_m+m), and the excitation as the drop across a band
    m = np.geomspace(0.1, 1e4, 800)
    ax[0].semilogx(m, phi(m, Km), color=BLUE, lw=2.0)
    for lo, hi, col, lab in [(M_BAR / np.sqrt(2), M_BAR * np.sqrt(2), RED,
                              "where the patient is held"),
                             (Km / np.sqrt(2), Km * np.sqrt(2), GREEN,
                              "where the information is")]:
        ax[0].axvspan(lo, hi, color=col, alpha=0.14)
        ax[0].plot([lo, hi], [phi(lo, Km), phi(hi, Km)], "o-", color=col, ms=4,
                   lw=1.4, label=f"{lab}: $E$ = {excitation(lo, hi, Km):.3f}")
    ax[0].axvline(Km, color=MUTED, ls=":", lw=1.0)
    over_x(ax[0], Km, r"$K_m$", frac=0.94, color=MUTED, fontsize=8, ha="left")
    ax[0].set_xlabel(r"metabolite $m$")
    ax[0].set_ylabel(r"$\varphi(m) = K_m/(K_m+m)$")
    ax[0].set_title("the direction that separates\n"
                    r"$\theta$ from $K_m$ is flat at both ends", fontsize=8)
    ax[0].legend(loc="lower left", fontsize=6.5)

    # (right) sd(ln theta) against where the band sits, load known and not
    cen = np.geomspace(0.02, 100.0, 120) * Km
    for lb, col, ls, lab in [(True, RED, "-", "load partly unrecorded"),
                             (False, GREY, "--", "load known exactly")]:
        ys = [sd_log_theta(band_levels(c, 2.0, 9), P.theta, Km, P.ke,
                           load_bias=lb) for c in cen]
        ys = np.array(ys) / min(y for y in ys if np.isfinite(y))
        ax[1].loglog(cen / Km, ys, color=col, ls=ls, lw=1.9, label=lab)
    ax[1].axvline(M_BAR / Km, color=MUTED, ls=":", lw=1.0)
    over_x(ax[1], M_BAR / Km, "the patient\nis here", frac=0.90, color=MUTED,
           fontsize=7, ha="left")
    ax[1].set_xlabel(r"$\sqrt{m_{\min} m_{\max}} \;/\; K_m$")
    ax[1].set_ylabel(r"$\mathrm{sd}(\ln\hat\theta)$, relative to its best")
    ax[1].set_title("with the load known, saturation is best;\n"
                    "with it unrecorded, the optimum is interior", fontsize=8)
    ax[1].legend(loc="upper right", fontsize=6.5)

    save(fig, "fig_D1_excitation")
    print(f"    E at the operating band [m_bar/sqrt2, m_bar*sqrt2]: "
          f"{excitation(M_BAR / np.sqrt(2), M_BAR * np.sqrt(2), Km):.4f}")
    print(f"    E at the same fold range straddling K_m:            "
          f"{excitation_max(2.0):.4f}  ({excitation_max(2.0) / excitation(M_BAR / np.sqrt(2), M_BAR * np.sqrt(2), Km):.1f}x)")
    best = min(sd_log_theta(band_levels(c, 2.0, 9), P.theta, Km, P.ke)
               for c in np.geomspace(0.02, 100.0, 400) * Km)
    here = sd_log_theta(band_levels(M_BAR, 2.0, 9), P.theta, Km, P.ke)
    print(f"    sd(ln theta) where the patient is held is {here / best:.0f}x its best")


# ==========================================================================
#  1.8  The covariance collapse.
#
#  The estimator is the exact information filter for CONSTANT parameters,
#
#      P_k^{-1} = P_0^{-1} + sum_{i<=k} q(m_i) q(m_i)^T / sigma^2 ,
#
#  on the augmented coordinates (ln theta, ln K_m, load bias), none of which
#  carries dynamics. "No process noise on theta" is not a setting here, it is
#  the construction: a random walk on theta manufactures a finite stationary
#  covariance in a direction the data may never see, which is the exact error
#  this note exists to warn about. What does not collapse, stays.
#
#  An extended Kalman filter on the joint state (m, theta, K_m, L_0) was
#  written first and is not used, for a reason worth recording rather than
#  hiding. The marker's time constant is 1/(theta g'(m_bar) + k_e), under an
#  hour at these values, while the therapeutic increment decays with a
#  half-life of days. The marker is therefore in quasi-steady state at every
#  sampling instant a clinic could use, a design is fully described by the
#  LEVELS it visits, and the filter above is exact rather than approximate.
#  quasi_steady_check() measures the separation instead of asserting it.
# ==========================================================================
def quasi_steady_level(t, Lb=None):
    """
    The level that balances the load at the instantaneous increment.

    Closed form, because the balance is a quadratic in m:

        k_e m^2 + (k_e K_m + theta + Delta(t) - L) m - L K_m = 0 .
    """
    Lb = load(0.0) if Lb is None else Lb
    d = np.asarray(delta(t), dtype=float)
    B = P.ke * P.Km + P.theta + d - Lb
    return (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Lb * P.Km)) / (2.0 * P.ke)


def quasi_steady_check(tol=0.05, cycles=3):
    """
    How long after a dose the marker is still catching up, and how closely it
    tracks the rest of the time.

    The separation of time scales is large -- the marker's time constant is
    well under an hour and the increment decays over days -- so between doses
    the marker sits on quasi_steady_level to a fraction of a per cent. The
    dose itself is a step, and immediately after one the marker is a long way
    from the level the new increment implies, so quasi-steadiness fails in a
    short window and holds through the rest of the cycle.

    That window is a design constraint rather than a nuisance, and the data
    specification carries it: a sample drawn inside it is not a steady state
    and does not belong in a steady-state fit.

    Returns (settling time in hours, worst relative departure after it).
    """
    Lb = load(0.0)
    tau = 1.0 / (P.theta * P.Km / (P.Km + M_BAR) ** 2 + P.ke)

    # the pre-dose level, reached by running a cycle out from the post-dose one
    m0 = float(quasi_steady_level(0.0))
    for c in range(int(cycles)):
        t0 = c * DOSE_INTERVAL
        m0 = float(solve_ivp(
            lambda t, y: [Lb - (P.theta + float(delta(t))) * y[0] / (P.Km + y[0])
                          - P.ke * y[0]],
            (t0 + 1e-9, t0 + DOSE_INTERVAL - 1e-9), [m0],
            rtol=1e-9, atol=1e-12, method="Radau").y[0, -1])

    # now the dose lands: integrate the cycle that follows it from that level
    t0 = cycles * DOSE_INTERVAL
    sol = solve_ivp(
        lambda t, y: [Lb - (P.theta + float(delta(t))) * y[0] / (P.Km + y[0])
                      - P.ke * y[0]],
        (t0, t0 + DOSE_INTERVAL - 1e-9), [m0], dense_output=True,
        rtol=1e-10, atol=1e-13, method="Radau")
    tt = t0 + np.geomspace(1e-4, DOSE_INTERVAL - 1e-3, 4000)
    qs = quasi_steady_level(tt)
    rel = np.abs(sol.sol(tt)[0] - qs) / qs
    over = np.where(rel > tol)[0]
    settle = float(tt[over[-1]] - t0) if over.size else 0.0
    after = rel[(tt - t0) > settle]
    worst = float(after.max()) if after.size else 0.0
    print(f"    marker time constant {tau*24:.2f} h against a dose half-life of "
          f"{np.log(2)/LAMBDA:.1f} d, a ratio of {np.log(2)/LAMBDA/tau:.0f}")
    print(f"    after a dose the marker reaches {100*tol:.0f}% of quasi-steady in "
          f"{settle*24:.2f} h, and tracks to within {100*worst:.3f}% thereafter")
    return settle, worst


def _filter(levels, sigma, increments=None, prior_rel=(None, None),
            load_frac=LOAD_FRAC, ke_prior=0.8):
    """
    Sequential covariance of (ln theta, ln K_m, load bias) as samples arrive.

    `increments` is Delta(t_i) at each sample, and carrying it is not optional.
    The model is

        s_i = (theta + Delta_i) g(m_i; K_m) + k_e m_i + b + eps_i ,

    so the ordinate's derivative in ln K_m is -(theta + Delta_i) K_m m / (K_m+m)^2,
    NOT -theta K_m m / (K_m+m)^2. The levels a design visits are generated by
    quasi_steady_level(t), which already responds to Delta; evaluating the
    sensitivity as though Delta were absent linearises the likelihood about a
    model the record was not generated from. That is not a duplicate term
    removed, it is the wrong tangent. See audit/c1_sensitivity.py, which
    finite-differences the clearance and agrees with the corrected column to
    1.1e-10 while the uncorrected one is wrong by up to 68% of itself.

    `prior_rel` are the log-priors on (ln theta, ln K_m), and BOTH DEFAULT TO
    None, meaning flat. Stripping the filter's priors one at a time shows that
    only the one on ln k_e does any work: on the accessible design, dropping
    the ln K_m prior changes sd(ln theta) by nothing at all, and a K_m prior
    twenty times looser still changes nothing, so nothing here depends on
    knowing K_m -- which Remark rem:no_km forbids. Dropping the ln k_e prior
    takes it from 0.15 to 0.65. The paper therefore claims one prior and lets
    theta and K_m float.

    `ke_prior` is the log-prior on ln k_e. The default 0.8 FITS k_e, which is
    the footing this paper claims: the non-enzymatic channel is not known, and
    a design evaluated as though it were flatters itself (audit/PHASE2A.md, A1).
    Passing ke_prior=None holds k_e exactly known and reproduces the three-column
    covariance the earlier drafts used; it is kept for the audit, not for a
    headline.

    Passing increments=None reproduces the published (uncorrected) column and
    is kept only so the audit's numbers stay reproducible.

    A warning the caller has to own: with Delta carried and the recorded load
    CONSTANT, the three columns are exactly rank deficient -- the load-bias
    column lies in the span of the other two at

        c_theta = (1 + K_m k_e / L)/theta ,   c_Km = -1/L

    (audit/c2b_degeneracy.py proves it symbolically). The degeneracy is an
    artefact of a constant recorded load and is broken by any real day-to-day
    variation in intake, but where it holds the covariance below is filled in
    by the load-record prior rather than by the data, and the direction is
    identified only at second order.
    """
    Lb = load(0.0)
    lv = np.asarray(levels, dtype=float)
    inc = (np.zeros_like(lv) if increments is None
           else np.broadcast_to(np.asarray(increments, dtype=float), lv.shape))
    pri = [np.inf if prior_rel[0] is None else prior_rel[0] ** 2,
           np.inf if prior_rel[1] is None else prior_rel[1] ** 2,
           (load_frac * Lb) ** 2]
    if ke_prior is not None:
        pri.append(ke_prior ** 2)
    Finv = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])

    def _cov(A):
        w = np.linalg.eigvalsh(A)
        return (np.full_like(A, np.inf) if w[0] <= 1e-14 * max(w[-1], 1e-300)
                else np.linalg.inv(A))

    out = [_cov(Finv)]
    for m, d in zip(lv, inc):
        q = [P.theta * m / (P.Km + m),
             -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
             1.0]
        if ke_prior is not None:
            q.append(P.ke * m)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / sigma ** 2
        out.append(_cov(Finv))
    return np.array(out)


def fig_collapse(T=730.0):
    """
    Sparse clinic sampling against sampling dense enough to catch the
    inter-dose recovery. The gate: does the dense design collapse P_theta?
    """
    import model_A.run as LOC
    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
    quasi_steady_check()

    designs = {
        "six-monthly clinic visits": np.arange(180.0, T + 1.0, 180.0),
        "one dose cycle, every 2 days": np.arange(0.0, DOSE_INTERVAL, 2.0),
        "every 2 days for two years": np.arange(0.0, T, 2.0),
    }
    fig, ax = plt.subplots(1, 2, figsize=figsize(1.0, 0.44))
    cols = dict(zip(designs, (RED, AMBER, GREEN)))
    finals = {}
    for name, times in designs.items():
        lv = quasi_steady_level(times)
        inc = np.asarray(delta(times), dtype=float)
        Ps = _filter(lv, sigma, inc)
        finals[name] = (times, lv, Ps)
        ax[0].semilogy(np.concatenate([[0.0], times]), np.sqrt(Ps[:, 0, 0]),
                       color=cols[name], lw=1.9,
                       label=f"{name} ($n={len(times)}$)")
        rho = Ps[-1][0, 1] / np.sqrt(Ps[-1][0, 0] * Ps[-1][1, 1])
        stale = np.sqrt(_filter(lv, sigma, ke_prior=None)[-1][0, 0])
        kek = np.sqrt(_filter(lv, sigma, inc, ke_prior=None)[-1][0, 0])
        print(f"    {name:>30}: n={len(times):4d}  levels {lv.min():6.3f}-{lv.max():6.3f}"
              f"  sd(ln theta) = {np.sqrt(Ps[-1][0,0]):8.4f}  corr = {rho:+.4f}"
              f"   [k_e held known: {kek:.4f};  Delta dropped too: {stale:.4f}]")
    ax[0].axhline(0.20, color=MUTED, ls="--", lw=1.0)
    corner(ax[0], r"$20\%$ on $\theta$", "lower left", color=MUTED, fontsize=7)
    ax[0].set_xlabel("time  [days]")
    ax[0].set_ylabel(r"$\mathrm{sd}(\ln\hat\theta)$")
    ax[0].set_title(r"no process noise on $\theta$:" "\n"
                    "what does not collapse, stays", fontsize=8)
    ax[0].legend(loc="upper right", fontsize=5.8)

    th = np.linspace(0, 2 * np.pi, 240)
    for name, (times, lv, Ps) in finals.items():
        w, V = np.linalg.eigh(Ps[-1][:2, :2])
        E = V @ (np.sqrt(np.maximum(w, 1e-30))[:, None]
                 * np.array([np.cos(th), np.sin(th)]))
        ax[1].plot(E[0], E[1], color=cols[name], lw=1.7,
                   label=name.split(",")[0])
    ax[1].plot([0], [0], "o", ms=5, color=MUTED, zorder=5)
    ax[1].set_xlabel(r"$\ln\theta$ error")
    ax[1].set_ylabel(r"$\ln K_m$ error")
    ax[1].set_title("one-sigma ellipse, after two years", fontsize=8)
    ax[1].legend(loc="upper left", fontsize=6)
    save(fig, "fig_D2_collapse")

    best = min(np.sqrt(f[2][-1][0, 0]) for f in finals.values())
    print(f"    GATE: the best of these reaches sd(ln theta) = {best:.3f}. "
          f"{'It collapses.' if best < 0.2 else 'NONE of them collapses to the 20% target.'}")
    print("    The gate is decided by whether the therapeutic increment is carried in")
    print("    the sensitivity. It is a known regressor, so it is carried; see the")
    print("    docstring of _filter and audit/c1_sensitivity.py.")
    times = designs["every 2 days for two years"]
    lv = quasi_steady_level(times)
    inc = np.asarray(delta(times), dtype=float)
    sig2 = LOC.sigma_total(LOAD_FRAC / 10.0, ASSAY_CV)
    print(f"    the same dense design with the load recorded to "
          f"{100*LOAD_FRAC/10:.1f}%: sd(ln theta) = "
          f"{np.sqrt(_filter(lv, sig2, inc)[-1][0,0]):.3f}")
    lv3 = np.geomspace(P.Km / 2.0, P.Km * 2.0, len(lv))
    print(f"    a design of the same size straddling K_m: sd(ln theta) = "
          f"{np.sqrt(_filter(lv3, sigma, inc)[-1][0,0]):.4f}")
    return finals


# ==========================================================================
#  1.9  The estimand is a decision.
# ==========================================================================

def fig_protocol(T=112.0):
    """
    What the schedule looks like, and what the phase actually buys.

    Section 5 states the sampling rule in words and Figure D2 shows what each
    choice costs. Neither shows the schedule itself, so a protocol writer has
    to reconstruct it from a standard deviation. This figure is the rule.

    The right panel plots the increment rather than the phase. Phase is only
    the proxy: what identifies theta is the RANGE of Delta/theta a record
    contains, and a schedule locked to the dose interval sees one value of it.
    Faceting the three schedules keeps that comparison a spatial one, so the
    reader never has to separate two hues to make it.
    """
    SCHEDULES = [
        ("monthly, every 30.4 d", 30.4, GREEN, "o", "0.170"),
        ("every 14.5 d, half a day off", 14.5, BLUE, "^", "0.161"),
        ("fortnightly, on the dose", DOSE_INTERVAL, RED, "X", "3.51"),
    ]

    fig, ax = plt.subplot_mosaic(
        [["traj", "p0"], ["traj", "p1"], ["rails", "p2"]],
        figsize=figsize(1.0, 0.62),
        gridspec_kw={"width_ratios": [1.45, 1.0], "height_ratios": [1.0, 1.0, 1.0],
                     "hspace": 0.20, "wspace": 0.34}, layout="constrained")

    # ---------------- the marker under the dose train ----------------
    tt = np.linspace(0.0, T, 2000)
    a = ax["traj"]
    a.plot(tt, quasi_steady_level(tt), color=INK, lw=1.5, zorder=3)
    a.set_ylabel(r"marker $m(t)$")
    a.set_xlim(0, T)
    a.set_xticklabels([])
    doses = np.arange(0.0, T + 1e-9, DOSE_INTERVAL)
    a.plot(doses, quasi_steady_level(doses), linestyle="none", marker="v",
           ms=3.8, color=AMBER, mec="white", mew=0.5, zorder=4, clip_on=False)
    a.annotate("dose, every 14 d", xy=(doses[2], quasi_steady_level(doses[2])),
               xytext=(0, -13), textcoords="offset points",
               fontsize=6.3, color=AMBER, ha="center")
    a.set_title("the schedule in time", fontsize=7.4)

    # ---------------- the three schedules as rails ----------------
    r = ax["rails"]
    r.set_xlim(0, T)
    r.set_xlabel("time  [days]")
    r.set_ylim(-0.7, len(SCHEDULES) - 0.3)
    r.set_yticks(range(len(SCHEDULES)))
    r.set_yticklabels([s[0].split(",")[0] for s in reversed(SCHEDULES)], fontsize=6.4)
    for i, (lab, dt, col, mk, sd) in enumerate(SCHEDULES):
        y = len(SCHEDULES) - 1 - i
        ts = np.arange(0.0, T + 1e-9, dt)
        r.axhline(y, color=GREY, lw=0.5, zorder=1)
        r.plot(ts, np.full_like(ts, y), linestyle="none", marker=mk, ms=4.4,
               color=col, mec="white", mew=0.6, zorder=5)
    r.spines["left"].set_visible(False)
    r.tick_params(axis="y", length=0)

    # ---------------- what each schedule sees of the increment ----------------
    # one steady cycle, well past the first dose
    ph = np.linspace(0.0, DOSE_INTERVAL, 400)
    curve = np.asarray(delta(20 * DOSE_INTERVAL + ph)) / P.theta
    ymax = float(curve.max()) * 1.18

    for i, (lab, dt, col, mk, sd) in enumerate(SCHEDULES):
        b = ax[f"p{i}"]
        ts = np.arange(0.0, 730.0, dt)
        dph = np.mod(ts, DOSE_INTERVAL)
        dval = np.asarray(delta(ts)) / P.theta
        rng = float(np.ptp(dval))

        b.plot(ph, curve, color=GREY, lw=1.0, zorder=2)
        b.axvspan(0.0, 0.22, color=RED, alpha=0.18, lw=0, zorder=1)
        b.axvspan(DOSE_INTERVAL - 0.22, DOSE_INTERVAL, color=RED, alpha=0.18,
                  lw=0, zorder=1)
        b.plot(dph, dval, linestyle="none", marker=mk, ms=4.0, color=col,
               mec="white", mew=0.6, alpha=0.9, zorder=5)

        # the range each schedule actually spans, drawn as a bar at the right
        xr = DOSE_INTERVAL * 1.045
        b.plot([xr, xr], [dval.min(), dval.max()], color=col, lw=2.4,
               solid_capstyle="butt", zorder=6, clip_on=False)
        b.text(DOSE_INTERVAL * 1.10, (dval.min() + dval.max()) / 2.0,
               f"range\n{rng:.2f}", fontsize=5.9, color=col, va="center",
               ha="left")

        b.set_xlim(0, DOSE_INTERVAL)
        b.set_ylim(0, ymax)
        b.set_yticks([0, 1, 2])
        b.tick_params(labelsize=6.2)
        b.text(0.42, ymax * 0.93, lab, fontsize=6.4, color=col, va="top")
        if i < len(SCHEDULES) - 1:
            b.set_xticklabels([])
        else:
            b.set_xlabel("phase within the 14-day dose cycle  [days]")
        if i == 1:
            b.set_ylabel(r"increment $\Delta/\theta$")

    ax["p0"].set_title("what each schedule sees of the increment", fontsize=7.4)
    ax["p2"].annotate("every draw at one phase,\nso one value of $\Delta$",
                      xy=(0.35, curve.max() * 0.99), xytext=(2.4, ymax * 0.42),
                      fontsize=6.0, color=RED,
                      arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))

    save(fig, "fig_D4_protocol")

def theta_crit_local(ke, Km, m_tol=M_TOL, s_peak=None):
    """The operational threshold, in the local identification model."""
    s_peak = 1.5 * clearance_A(M_BAR, P) if s_peak is None else s_peak
    return (s_peak - ke * m_tol) * (Km + m_tol) / m_tol


def margin_and_gradient(theta, Km, ke, m_tol=M_TOL, s_peak=None):
    """
    mu and its gradient in LOG parameters, which is what the design delivers.

        mu = theta - (s_peak - ke m_tol)(K_m + m_tol)/m_tol
        d mu/d ln theta = theta
        d mu/d ln K_m   = -K_m (s_peak - ke m_tol)/m_tol
        d mu/d ln k_e   = +k_e (K_m + m_tol)

    The two parameters enter with opposite signs in log coordinates, which is
    the cancellation of the margin theorem seen from the estimator's side.
    """
    s_peak = 1.5 * clearance_A(M_BAR, P) if s_peak is None else s_peak
    mu = theta - (s_peak - ke * m_tol) * (Km + m_tol) / m_tol
    grad = np.array([theta,
                     -Km * (s_peak - ke * m_tol) / m_tol,
                     ke * (Km + m_tol),
                     0.0])
    return float(mu), grad


def sign_error_linear(levels, theta, Km, ke, sigma, m_tol=M_TOL, s_peak=None):
    """
    The Gaussian tail a delta-method calculation would give. NOT used for the
    figure, and kept because the reason it is wrong is itself a result.

    It propagates the estimator's covariance to mu and reports
    Phi(-|mu|/sd(mu_hat)). That is only correct if mu_hat is symmetric about
    mu, and it is not: the design's near-null direction is very nearly the
    direction along which mu is monotone and of one sign, so the sampling
    distribution of mu_hat is heavily skewed and does not cross zero anything
    like as often as a Gaussian of the same width would. This is the margin
    theorem showing up in the estimator, and taking the Gaussian tail at face
    value would throw the paper's own central result away.
    """
    from common.excitation import information
    from scipy.stats import norm
    F = information(levels, theta, Km, ke, sigma=sigma, load_bias=True)
    mu, grad = margin_and_gradient(theta, Km, ke, m_tol, s_peak)
    try:
        var = float(grad @ np.linalg.inv(F) @ grad)
    except np.linalg.LinAlgError:
        return dict(mu=mu, sd=np.inf, p=0.5, kind="n/a")
    if not np.isfinite(var) or var <= 0:
        return dict(mu=mu, sd=np.inf, p=0.5, kind="n/a")
    sd = np.sqrt(var)
    return dict(mu=mu, sd=sd, p=float(norm.cdf(-abs(mu) / sd)),
                kind="fail to treat" if mu < 0 else "treat unneeded")


_ENSEMBLE_CACHE = {}


def fit_ensemble(levels, theta, Km, ke, sigma_load=None, cv=None, reps=400,
                 seed=0, m_bar=M_BAR, increments=None, rho=1.0, kappa_cv=0.0):
    """
    An ensemble of fitted parameter sets from synthetic records of one design.

    The fit does not depend on the illness amplitude -- only the margin
    evaluated afterwards does -- so the ensemble is computed once per design
    and reused across the whole sweep. Each replicate fits
    (theta, K_m, k_e) together with an unknown constant in the load, started
    away from the truth so that a flat likelihood shows up as a spread rather
    than as a fixed point at the starting value.

    The assay error is placed on the ORDINATE as b*sigma_m, where the analytic
    Fisher matrix puts it, rather than on the abscissa where it physically
    sits. Perturbing m and fitting against the perturbed m would make the fit
    look better than it is, because the model would then be evaluated at the
    level the assay reported rather than the level the patient was at. That
    errors-in-variables gap has its own remark; it is not what this is for.
    """
    sigma_load = LOAD_FRAC if sigma_load is None else sigma_load
    cv = ASSAY_CV if cv is None else cv
    key = (tuple(np.round(np.asarray(levels, float), 9)), theta, Km, ke,
           sigma_load, cv, int(reps), int(seed), float(rho), float(kappa_cv),
           None if increments is None
           else tuple(np.round(np.atleast_1d(increments).astype(float), 9)))
    if key in _ENSEMBLE_CACHE:
        return _ENSEMBLE_CACHE[key]
    rng = np.random.default_rng(seed)
    lv = np.asarray(levels, dtype=float)
    inc = (np.zeros_like(lv) if increments is None
           else np.broadcast_to(np.asarray(increments, dtype=float),
                                lv.shape).astype(float))
    b_loc = theta * Km / (Km + m_bar) ** 2 + ke
    s_true = (theta + inc) * lv / (Km + lv) + ke * lv
    L_base = clearance_A(m_bar, P)
    s_rec = L_base + rho * (s_true - L_base)
    out = []
    for _ in range(int(reps)):
        s_obs = (s_rec * (1.0 + rng.normal(0.0, sigma_load, lv.size))
                 + rng.normal(0.0, b_loc * cv * m_bar, lv.size))

        inc_fit = inc * (1.0 + rng.normal(0.0, kappa_cv)) if kappa_cv else inc

        def res(q):
            t, K, k = np.exp(np.clip(q[:3], -25.0, 25.0))
            return (t + inc_fit) * lv / (K + lv) + k * lv + q[3] - s_obs

        start = np.concatenate([np.log([theta, Km, ke])
                                + rng.normal(0.0, 0.5, size=3), [0.0]])
        try:
            r = least_squares(res, start, method="lm", max_nfev=3000)
        except Exception:
            continue
        out.append(np.exp(np.clip(r.x[:3], -25.0, 25.0)))
    arr = np.array(out)
    _ENSEMBLE_CACHE[key] = arr
    return arr


def sign_error(levels, theta, Km, ke, sigma_load=None, cv=None, reps=400,
               seed=0, m_tol=M_TOL, s_peak=None, m_bar=M_BAR, increments=None,
               rho=1.0, kappa_cv=0.0):
    """
    The operating characteristic of the decision, by simulation.

    mu = theta - theta_crit. mu > 0 is a patient with margin; mu < 0 is one
    who needs the therapy intensified. The two errors are not interchangeable
    and the `kind` field says which side of the boundary the patient is on:

        treat unneeded   true mu > 0, estimated mu < 0
        fail to treat    true mu < 0, estimated mu > 0
    """
    ens = fit_ensemble(levels, theta, Km, ke, sigma_load, cv, reps, seed, m_bar,
                       increments, rho, kappa_cv)
    mu_true = margin_and_gradient(theta, Km, ke, m_tol, s_peak)[0]
    if ens.size == 0:
        return dict(mu=mu_true, p=np.nan, n=0, kind="n/a")
    mu_hat = np.array([margin_and_gradient(t, K, k, m_tol, s_peak)[0]
                       for t, K, k in ens])
    return dict(mu=mu_true, p=float(np.mean(np.sign(mu_hat) != np.sign(mu_true))),
                n=len(ens),
                kind="fail to treat" if mu_true < 0 else "treat unneeded")


def excursion_for_decision(target=0.05, mu_frac=0.30, reps=300, seed=31, n=9,
                           tops=(2.0, 5.0, 10.0, 20.0, 40.0, 80.0)):
    """
    How far the marker has to go for the DECISION to be reliable, and at what
    distance from the boundary that is asked.

    The note argues that a confidence interval on theta is the wrong output, so
    deriving the excursion requirement from a precision target on theta would
    be inconsistent. This derives it from the decision instead.

    `mu_frac` is where the question is put: at a patient whose true margin is
    -mu_frac * theta, that is who sits that far BELOW the threshold. The choice
    matters and it is not cosmetic. Asking the question arbitrarily close to
    the boundary is asking for the sign of a quantity that is nearly zero, and
    no design answers that: at |mu| = 0.11 theta the best excursion tried here
    still errs 30% of the time. That is a property of taking a sign near a
    boundary, not a property of the design, and quoting it as a design
    requirement would be a category error.

    Returns (top of the excursion as a multiple of m_bar, error there), or
    (None, best error reached).
    """
    base = clearance_A(M_BAR, P)
    # the amplitude at which the true margin is -mu_frac * theta
    grid = np.linspace(1.0, 20.0, 8001)
    mus = np.array([margin_and_gradient(P.theta, P.Km, P.ke, s_peak=f * base)[0]
                    for f in grid])
    # mus decreases with the amplitude, so -mus increases: interpolate on that.
    # The target is mu = -mu_frac*theta, i.e. -mu = +mu_frac*theta.
    f_at = float(np.interp(mu_frac * P.theta, -mus, grid))
    assert margin_and_gradient(P.theta, P.Km, P.ke,
                               s_peak=f_at * base)[0] < 0, \
        'the question must be asked on the unsafe side of the boundary'
    s_peak = f_at * base
    curve = []
    for top in tops:
        lv = np.geomspace(M_BAR, top * M_BAR, n)
        curve.append((top, sign_error(lv, P.theta, P.Km, P.ke, reps=reps,
                                      seed=seed, s_peak=s_peak)["p"]))
    hit = [t for t, e in curve if e <= target]
    return (hit[0] if hit else None), min(e for _, e in curve), f_at, curve


def cohort_size(p_err, delta=0.10, frac=0.30, z=1.96):
    """
    Patients needed to estimate the fraction with mu < 0 to +/- delta.

    Per-patient misclassification attenuates the observed fraction towards
    one half, and correcting for it inflates the variance by
    (1 - p_0 - p_1)^{-2}. Taking the two error rates equal at p_err,

        N = z^2 f_obs (1 - f_obs) / (delta^2 (1 - 2 p_err)^2) ,

    which diverges as p_err approaches one half: at that point the record
    carries no information about the sign and no cohort is large enough. That
    divergence is the point of quoting a cohort size at all.
    """
    if not np.isfinite(p_err) or p_err >= 0.5:
        return np.inf
    f_obs = frac * (1 - p_err) + (1 - frac) * p_err
    return float(np.ceil(z ** 2 * f_obs * (1 - f_obs)
                         / (delta ** 2 * (1 - 2 * p_err) ** 2)))


def fig_decision():
    """
    The sign error on mu, against the three things a protocol controls.

    A patient of fixed biology is swept across the decision boundary by moving
    the illness amplitude, so both error types are exercised on one axis and
    the asymmetry is visible. What a trial has to live with is the WIDTH of
    the region where the error rate is appreciable, not the rate at one point.
    """
    import model_A.run as LOC
    base = clearance_A(M_BAR, P)
    folds = np.linspace(1.5, 5.5, 21)
    reps = 400
    fig, ax = plt.subplots(1, 3, figsize=figsize(1.0, 0.34))

    def curve(levels, **kw):
        return [100 * sign_error(levels, P.theta, P.Km, P.ke, reps=reps,
                                 seed=11, s_peak=f * base, **kw)["p"]
                for f in folds]

    lv9 = np.geomspace(M_BAR, 5 * M_BAR, 9)
    lv_km = np.geomspace(P.Km / 3.0, P.Km * 3.0, 9)
    for cv, col, ls in [(0.02, GREY, ":"), (0.05, GREEN, "-"), (0.128, RED, "--")]:
        ax[0].plot(folds, curve(lv9, cv=cv), color=col, ls=ls, lw=1.7,
                   label=rf"{100*cv:.0f}\%")
    ax[0].set_ylabel(r"sign error on $\mu$  [\%]")
    ax[0].set_title("against assay variability", fontsize=8)
    ax[0].legend(fontsize=6, title="assay CV", title_fontsize=6)

    for nn, col, ls in [(4, RED, "--"), (9, GREEN, "-"), (40, BLUE, ":")]:
        ax[1].plot(folds, curve(np.geomspace(M_BAR, 5 * M_BAR, nn)),
                   color=col, ls=ls, lw=1.7, label=f"$n = {nn}$")
    ax[1].set_title("against the number of samples", fontsize=8)
    ax[1].legend(fontsize=6)

    for top, col, ls in [(1.5, RED, "--"), (5.0, GREEN, "-"), (20.0, BLUE, ":")]:
        ax[2].plot(folds, curve(np.geomspace(M_BAR, top * M_BAR, 9)),
                   color=col, ls=ls, lw=1.7, label=rf"top $={top:g}\,\bar m$")
    ax[2].plot(folds, curve(lv_km), color=AMBER, lw=1.7,
               label=r"straddling $K_m$")
    ax[2].set_title("against where the band sits", fontsize=8)
    ax[2].legend(fontsize=6)

    mu0 = np.array([margin_and_gradient(P.theta, P.Km, P.ke, s_peak=f * base)[0]
                    for f in folds])
    cross = float(np.interp(0.0, -mu0, folds))
    for a in ax:
        a.axvline(cross, color=MUTED, ls=":", lw=0.9)
        a.set_xlabel(r"illness amplitude  [$\times$ baseline]")
        a.set_ylim(0, 60)
    over_x(ax[0], cross, r"$\mu = 0$", frac=0.95, color=MUTED, fontsize=7,
           ha="left")
    save(fig, "fig_D3_decision")

    print("    the decision, not a confidence interval "
          f"(the boundary is at {cross:.2f}x baseline):")
    for f in (1.5, 2.5, cross - 0.3, cross + 0.3, 4.0, 5.5):
        d = sign_error(lv9, P.theta, P.Km, P.ke, reps=400, seed=12,
                       s_peak=f * base)
        dl = sign_error_linear(lv9, P.theta, P.Km, P.ke,
                               LOC.sigma_total(LOAD_FRAC, ASSAY_CV),
                               s_peak=f * base)
        print(f"      illness {f:4.2f}x: mu = {d['mu']:9.1f} ({d['kind']:>14})"
              f"   simulated {100*d['p']:5.1f}%   delta method would say "
              f"{100*dl['p']:5.1f}%")
    bad = [f for f in folds
           if sign_error(lv9, P.theta, P.Km, P.ke, reps=250, seed=13,
                         s_peak=f * base)["p"] > 0.10]
    if bad:
        print(f"      simulated sign error exceeds 10% for amplitudes "
              f"{min(bad):.2f}x to {max(bad):.2f}x baseline, a window "
              f"{100*(max(bad)-min(bad))/cross:.0f}% of the boundary wide")
    else:
        print("      simulated sign error never exceeds 10% on this sweep")


# ==========================================================================
#  1.7  The excursion table, and the specification it terminates in.
# ==========================================================================
def excursion_sources():
    """
    What each accessible source actually moves the marker by, and in which
    direction. Directions matter: the curvature the identification needs is at
    HIGH m, and most of what a chart contains moves the patient down.
    """
    # inter-dose recovery on the trial schedule, measured off the simulation
    ts = np.linspace(3 * DOSE_INTERVAL, 5 * DOSE_INTERVAL, 4000)
    sol = solve_ivp(
        lambda t, y: [load(t) - (P.theta + float(delta(t))) * y[0] / (P.Km + y[0])
                      - P.ke * y[0]],
        (0.0, 6 * DOSE_INTERVAL), [M_BAR], dense_output=True, rtol=1e-10,
        atol=1e-12, method="LSODA")
    mm = sol.sol(ts)[0]
    inter = (float(mm.min()), float(mm.max()))
    return [
        ("inter-dose recovery, biweekly schedule", inter[0], inter[1], "down"),
        ("carnitine dose step (quasi-exogenous)", 0.3 * M_BAR, M_BAR, "down"),
        ("cyclical intestinal antibiotics", 0.6 * M_BAR, M_BAR, "down"),
        ("carnitine supply interruption", M_BAR, 1.6 * M_BAR, "up"),
        ("decompensation and recovery", M_BAR, 8.0 * M_BAR, "up"),
        ("labelled-propionate bolus (not administered)", M_BAR, 2.0 * M_BAR, "up"),
    ]


def table_designs(T=730.0):
    """
    The design table of Section 7, generated rather than transcribed.

    Every row is on the footing this paper claims: the therapeutic increment
    carried, theta and K_m FREE, and priors only on ln k_e and the load record.
    The last column is the abandoned footing, where the increment was dropped
    from the sensitivity while still being used to generate the levels.
    """
    import model_A.run as LOC
    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
    dense_t = np.arange(0.0, T, 2.0)
    designs = [
        ("six-monthly clinic visits", np.arange(180.0, T + 1.0, 180.0), None),
        ("one dose cycle, every two days", np.arange(0.0, DOSE_INTERVAL, 2.0), None),
        (r"\textbf{monthly, across dose cycles}", np.arange(0.0, T, 30.4), None),
        ("two-monthly, across dose cycles", np.arange(0.0, T, 60.0), None),
        ("every two days for two years", dense_t, None),
        (r"\emph{fortnightly: the dose interval itself}",
         np.arange(0.0, T, DOSE_INTERVAL), None),
        (r"\emph{the same $365$ straddling $K_m$}", dense_t,
         np.geomspace(P.Km / 2.0, P.Km * 2.0, len(dense_t))),
    ]
    rows, raw = [], []
    for name, times, override in designs:
        lv, inc = design_record(times)
        if override is not None:
            lv = override
        Pc = _filter(lv, sigma, inc)[-1]
        sd = float(np.sqrt(Pc[0, 0]))
        rho = float(Pc[0, 1] / np.sqrt(Pc[0, 0] * Pc[1, 1]))
        dropped = float(np.sqrt(_filter(lv, sigma)[-1][0, 0]))
        drng = float((inc.max() - inc.min()) / P.theta)
        raw.append((name, len(times), lv.min(), lv.max(), drng, sd, rho, dropped))
        cell = (r"\textbf{%.4f}" % sd) if sd < 0.20 else ("%.4f" % sd)
        rows.append([name, len(times), f"${lv.min():.2f}$--${lv.max():.2f}$",
                     f"${drng:.2f}$", cell, f"${rho:+.4f}$", f"${dropped:.3f}$"])
    write_table("tab_D_designs",
                ["design", "$n$", "levels visited",
                 r"range of $\Delta/\theta$",
                 r"$\operatorname{sd}(\ln\hat\theta)$",
                 r"$\operatorname{corr}(\ln\hat\theta, \ln \hat K_m)$",
                 r"\emph{$\Delta$ dropped}"],
                rows,
                note="theta and K_m free; priors only on ln k_e and the load record")
    print("  wrote tab_D_designs.tex")
    for nm, n, lo, hi, drng, sd, rho, drop in raw:
        print(f"    {nm[:40]:>40}: n={n:4d}  range(D/th)={drng:5.2f}  sd={sd:.4f}"
              f"  corr={rho:+.4f}  [Delta dropped {drop:.3f}]")
    return raw


def table_excursion():
    rows = []
    for name, lo, hi, direction in excursion_sources():
        r = hi / lo
        E = excitation(lo, hi, P.Km)
        gm = np.sqrt(lo * hi) / P.Km
        rows.append([name,
                     f"{lo:.2f}--{hi:.1f}",
                     f"{r:.1f}$\\times$",
                     direction,
                     f"{gm:.3f}",
                     f"{E:.4f}",
                     f"{excitation_max(r):.4f}",
                     f"{E / excitation_max(r):.2f}"])
    write_table("tab_D_excursion",
                ["source", "$m$ range", "fold", "direction",
                 "$\\sqrt{m_{\\min}m_{\\max}}/K_m$", "$E$", "$E_{\\max}(r)$",
                 "frac.\\ of $E_{\\max}$"],
                rows,
                colspec=r"p{0.29\textwidth}lllllll",
                note="E is the excitation the band delivers; E_max(r) is the most "
                     "a band of that fold range could deliver if it straddled K_m")
    print("    excursion available from each source:")
    for name, lo, hi, d in excursion_sources():
        print(f"      {name:<44} m {lo:7.2f}-{hi:7.2f}  E = {excitation(lo,hi,P.Km):.4f}"
              f"  ({100*excitation(lo,hi,P.Km)/excitation_max(hi/lo):5.1f}% of the best at that fold)")
    return rows


def spec_numbers():
    """
    Every number the data specification quotes, computed in one place, so the
    table and the prose cannot drift apart.
    """
    import model_A.run as LOC
    from scipy.optimize import brentq
    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
    Lb = load(0.0)

    def sd_of(times, kap=None, load_cv=None, ke_prior=0.8, load_prior=True):
        lv, inc = design_record(times, load_cv=load_cv)
        pri = [np.inf, np.inf,
               (LOAD_FRAC * Lb) ** 2 if load_prior else np.inf,
               ke_prior ** 2]
        if kap is not None:
            pri.append(kap ** 2)
        F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
        for m, d in zip(lv, inc):
            g = m / (P.Km + m)
            q = [P.theta * g, -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                 1.0, P.ke * m]
            if kap is not None:
                q.append(d * g)
            q = np.array(q)
            F = F + np.outer(q, q) / sigma ** 2
        w = np.linalg.eigvalsh(F)
        if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
            return np.inf
        v = np.linalg.inv(F)[0, 0]
        return float(np.sqrt(v)) if v > 0 else np.inf

    monthly = np.arange(0.0, 730.0, 30.4)
    two_monthly = np.arange(0.0, 730.0, 60.0)
    N = {"sigma": sigma,
         "sd_monthly": sd_of(monthly),
         "sd_two_monthly": sd_of(two_monthly),
         "sd_aliased": sd_of(np.arange(0.0, 730.0, DOSE_INTERVAL)),
         "sd_aliased_off": sd_of(np.arange(0.0, 730.0, DOSE_INTERVAL + 0.5)),
         "sd_load_1pct": sd_of(monthly, load_cv=0.01),
         "sd_load_20pct": sd_of(monthly, load_cv=0.20),
         "sd_load_steady": sd_of(monthly, load_cv=0.0),
         # the same, with NO prior on the load bias, which is where the rank
         # deficiency of prop:rank actually shows: the prior fills that
         # direction, and only real variation lets the DATA fill it instead
         "free_steady": sd_of(monthly, load_cv=0.0, load_prior=False),
         "free_cv2": sd_of(monthly, load_cv=0.02, load_prior=False),
         "free_cv10": sd_of(monthly, load_cv=0.10, load_prior=False),
         "free_cv20": sd_of(monthly, load_cv=0.20, load_prior=False)}
    for tag, tt in (("monthly", monthly), ("two_monthly", two_monthly)):
        try:
            N["kappa_" + tag] = brentq(lambda x: sd_of(tt, kap=x) - 0.20,
                                       1e-3, 3.0, xtol=1e-5)
        except ValueError:
            N["kappa_" + tag] = None
    try:
        N["ke_needed"] = brentq(lambda x: sd_of(monthly, ke_prior=x) - 0.20,
                                1e-3, 3.0, xtol=1e-5)
    except ValueError:
        N["ke_needed"] = None
    return N


def table_spec():
    """
    The one-page data specification, in two tiers.

    The order is not cosmetic: Tier 1 measures the only quantity the method's
    feasibility turns on, and until it reports there is no defensible number of
    breath tests to ask for.
    """
    N = spec_numbers()
    kap = N["kappa_monthly"]
    kefac = np.exp(N["ke_needed"]) if N["ke_needed"] else None
    B = "\\textbf{"
    rows = [
        [B + "TIER 1} --- run first", "", ""],
        ["what to measure",
         "same-day repeatability of $1$-$^{13}$C-propionate oxidation: two "
         "breath tests bracketing one dose, on one day",
         "the only quantity the design turns on, and no published value for "
         "it exists"],
        ["what it must reach",
         ("a coefficient of variation on $\\Delta$ below $%.1f\\%%$" % (100 * kap)
          if kap else "a CV on $\\Delta$ inside the gate"),
         "derived from the gate, not chosen: at that precision monthly "
         "sampling reaches $20\\%$ on $\\theta$"],
        ["what is known now",
         "between-visit CV, mean $14.0\\%$, range $0.4$--$50.8\\%$, $n=15$, "
         "intervals $2$ months to $4$ years",
         "an upper bound on same-day error, already close to the requirement, "
         "but it does not establish it; and it is MMA, not PA"],
        [B + "TIER 2} --- conditional on Tier 1", "", ""],
        ["marker", "plasma 3-hydroxypropionate",
         "endpoint association (RR $0.7$, $95\\%$ CI $0.6$--$0.9$) and "
         "immunity to the carnitine sign reversal that disqualifies C3"],
        ["sampling schedule",
         "monthly across dose cycles, $n \\approx 25$ over two years "
         "($%.3f$); two-monthly also passes ($%.3f$)"
         % (N["sd_monthly"], N["sd_two_monthly"]),
         "ordinary outpatient density; what matters is covering the phase of "
         "the dose cycle, not the sampling rate"],
        ["sampling, the hard constraint",
         B + "never at the dose interval or a multiple of it}: fortnightly "
         "gives $%.2f$, $14.5$-daily gives $%.3f$"
         % (N["sd_aliased"], N["sd_aliased_off"]),
         "every draw lands at the same phase, the increment stops varying, "
         "and the therapy channel closes"],
        ["sampling, the other exclusion",
         "no sample within one hour of a dose",
         "the marker is not at quasi-steady state there and a steady-state "
         "fit does not apply to it"],
        ["dosing record",
         "date, dose and product lot for every infusion",
         "$\\Delta(t_i)$ is a regressor, not a nuisance: it carries the "
         "excitation"],
        ["load record, precision",
         "ordinary outpatient precision suffices ($%.3f$ at a $1\\%%$ record "
         "against $%.3f$ at $20\\%%$)" % (N["sd_load_1pct"], N["sd_load_20pct"]),
         "the load-bias column is nearly orthogonal to $\\theta$ once the "
         "increment is carried"],
        ["load record, variability",
         B + "record intake as eaten; do not smooth, average or impute}",
         "a steady load is the exactly singular case of "
         "Proposition~\\ref{prop:rank}. With the $5\\%%$ load prior in place "
         "that direction is filled by the prior and the answer barely moves; "
         "with no prior it is filled only by real variation "
         "($%s$ steady, $%.2f$ at a $2\\%%$ CV, $%.2f$ at $20\\%%$), so the "
         "variability is what makes the result independent of the prior"
         % ("\\infty" if not np.isfinite(N["free_steady"])
            else "%.2f" % N["free_steady"], N["free_cv2"], N["free_cv20"])],
        ["$k_e$",
         ("a literature prior within a factor of $%.1f$" % kefac if kefac
          else "a literature prior"),
         "the only prior the result depends on; no measurement programme is "
         "proposed for it"],
        ["$K_m$", "nothing at all",
         "every number here is unchanged by a $K_m$ prior twenty times "
         "looser, which is what Remark~\\ref{rem:no_km} demanded"],
    ]
    write_table("tab_D_spec", ["", "what the record must contain", "why"], rows,
                note="two tiers; Tier 2 is conditional on Tier 1 clearing the gate",
                colspec=(r"p{0.17\textwidth}@{\hskip 1em}"
                         r"p{0.40\textwidth}@{\hskip 1em}"
                         r"p{0.40\textwidth}"))
    print("  wrote tab_D_spec.tex")
    for k, v in N.items():
        print("    %18s = %s" % (k, v))
    return N


if __name__ == "__main__":
    print("DESIGN")
    fig_excitation()
    fig_collapse()
    fig_decision()
    fig_protocol()
    table_excursion()
    table_spec()
