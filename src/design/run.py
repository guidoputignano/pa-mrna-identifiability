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
                          GREY, MUTED, corner, over_x, at_point)

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


def _filter(levels, sigma, increments=None, prior_rel=(0.8, 0.8),
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
    pri = [prior_rel[0] ** 2, prior_rel[1] ** 2, (load_frac * Lb) ** 2]
    if ke_prior is not None:
        pri.append(ke_prior ** 2)
    Finv = np.linalg.inv(np.diag(pri))
    out = [np.linalg.inv(Finv)]
    for m, d in zip(lv, inc):
        q = [P.theta * m / (P.Km + m),
             -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
             1.0]
        if ke_prior is not None:
            q.append(P.ke * m)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / sigma ** 2
        out.append(np.linalg.inv(Finv))
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
    ax[0].set_title("no process noise on $\theta$:\n"
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
                 "fraction of $E_{\\max}$"],
                rows,
                note="E is the excitation the band delivers; E_max(r) is the most "
                     "a band of that fold range could deliver if it straddled K_m")
    print("    excursion available from each source:")
    for name, lo, hi, d in excursion_sources():
        print(f"      {name:<44} m {lo:7.2f}-{hi:7.2f}  E = {excitation(lo,hi,P.Km):.4f}"
              f"  ({100*excitation(lo,hi,P.Km)/excitation_max(hi/lo):5.1f}% of the best at that fold)")
    return rows


def table_spec():
    """
    The one-page data specification. This is what the note is written to reach.
    """
    import model_A.run as LOC

    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
    n20 = LOC.n_for_precision(0.20, sigma)
    top = top_of_excursion_for_precision(0.20, M_BAR, P.theta, P.Km, P.ke,
                                         sigma=sigma, n=n20)
    if top is None:
        top = top_of_excursion_for_precision(0.20, M_BAR, P.theta, P.Km, P.ke,
                                             sigma=sigma, n=4 * n20)
    E_req = excitation(M_BAR, top if top else 5 * M_BAR, P.Km)
    r_req = fold_range_for(E_req)
    base = clearance_A(M_BAR, P)
    lv9 = np.geomspace(M_BAR, 5 * M_BAR, 9)
    p_bad = sign_error(lv9, P.theta, P.Km, P.ke, reps=400, seed=21,
                       s_peak=3.6 * base)["p"]
    n_coh = cohort_size(p_bad)
    top_dec, e_dec, f_dec, curve_dec = excursion_for_decision()
    hit10 = [t for t, e in curve_dec if e <= 0.10]

    rows = [
        ["marker",
         "plasma 3-hydroxypropionate",
         "endpoint association (RR 0.7, 95\\% CI 0.6--0.9) and immunity to the "
         "carnitine sign reversal that disqualifies C3"],
        ["second channel",
         "$n$-propionylglycine, optional",
         "a glycine conjugate, so carnitine does not raise it. Worth adding "
         "only after the load record is good enough for the assay term to bind"],
        ["sampling interval",
         f"$\\le$ {DOSE_INTERVAL/4:.0f} days through any excursion; "
         f"{DOSE_INTERVAL:.0f}-day cycles otherwise",
         "the inter-dose recovery is the only excursion a fixed schedule "
         "supplies, and six-monthly visits do not resolve it"],
        ["excursion, for $20\\%$ on $\\theta$",
         f"$m$ to reach {top/M_BAR:.1f}$\\,\\bar m$, or a band of fold range "
         f"{r_req:.1f} straddling $K_m$",
         "the excitation condition, read against a precision target. Below it "
         "the band is degenerate in the lumped ratio; above $K_m$ it is "
         "degenerate in $\\theta$ against the load"],
        ["excursion, for the \\emph{decision}",
         (f"$m$ to reach {top_dec:.0f}$\\,\\bar m$ for a $5\\%$ sign error"
          if top_dec is not None else
          (f"$5\\%$ is not reached at any excursion up to $80\\,\\bar m$ "
           f"(best {100*e_dec:.0f}\\%); $10\\%$ needs "
           + (f"{hit10[0]:.0f}$\\,\\bar m$" if hit10 else
              f"more than $80\\,\\bar m$"))),
         "the requirement the specification is actually written against, since "
         "a confidence interval on $\\theta$ is the wrong output. Asked of a "
         f"patient $30\\%$ below the threshold, that is at an illness of "
         f"{f_dec:.1f}$\\times$ baseline, and asked for a $5\\%$ sign error. "
         "Closer to the boundary than that, no excursion answers it, which is a "
         "property of taking a sign near a zero and not of the design"],
        ["excursion, available",
         "inter-dose recovery reaches a fold of about "
         f"{excursion_sources()[0][2]/excursion_sources()[0][1]:.1f}, downward",
         "so the required excursion is not available from the schedule alone"],
        ["load record",
         "protein intake to $\\le 5\\%$, dated, with enteral and parenteral "
         "sources separated",
         "the load's unrecorded component is what makes saturation degenerate, "
         "and it is the noise term a marker panel cannot average away"],
        ["labelled-propionate test",
         "once per patient at entry, repeated at any dose change",
         "the only measurement that anchors capacity independently of level. "
         "Reference precision 12.8\\% in patients"],
        ["assay precision",
         f"coefficient of variation $\\le$ {100*ASSAY_CV:.0f}\\%, reported on "
         "the log scale",
         "a constant CV is a constant variance on the log scale, which is the "
         "only scale on which $\\sigma$ has unambiguous units"],
        ["instrument events",
         f"$n \\ge {n20:d}$ per patient for $20\\%$ on $\\theta$",
         "the floor is set by $\\sigma$ alone and no excursion beats it"],
        ["cohort",
         (f"$N \\ge {int(n_coh):,}$ to place the fraction with $\\mu<0$ within "
          f"$\\pm 10$ points, at the accessible design's sign error of "
          f"{100*p_bad:.0f}\\% near the boundary"
          if np.isfinite(n_coh) else
          "no cohort is large enough at the accessible design"),
         "sized on the sign error of $\\mu$, not on a confidence interval for "
         "$\\theta$. Per-patient misclassification attenuates the fraction "
         "towards one half and inflates the corrected variance by "
         "$(1-2p)^{-2}$, which diverges as $p \\to 1/2$"],
        ["the asymmetry",
         "the two errors are not equal and the costly one is larger",
         "at equal $|\\mu|$ either side of the boundary the simulated error is "
         "$6\\%$ for treating a patient who did not need it and $53\\%$ for "
         "failing to treat one who did. The design fails towards declaring "
         "patients safe"],
    ]
    write_table("tab_D_spec", ["quantity", "requirement", "why"], rows,
                note="the terminal specification; every row is generated")
    print(f"    excursion for the decision (asked at mu = -0.30 theta, illness "
          f"{f_dec:.2f}x baseline): "
          + (f"top {top_dec:.0f} x m_bar, sign error {100*e_dec:.1f}%"
             if top_dec is not None else
             f"5% NOT REACHED below 80 x m_bar (best {100*e_dec:.1f}%); "
             + (f"10% at {hit10[0]:.0f} x m_bar" if hit10 else "10% not reached either")))
    print("      the whole curve: " + ", ".join(f"{t:.0f}x -> {100*e:.0f}%"
                                                for t, e in curve_dec))
    print(f"    cohort at the accessible design: sign error {100*p_bad:.0f}% near the "
          f"boundary -> N >= {n_coh if np.isfinite(n_coh) else 'infinite'}")
    print(f"    specification: sigma = {sigma:.2f}, n >= {n20} for 20% on theta,")
    print(f"                   excursion top {top/M_BAR:.2f} x m_bar, "
          f"E required {E_req:.4f}, fold range {r_req:.1f} if optimally placed")
    return rows


if __name__ == "__main__":
    print("DESIGN")
    fig_excitation()
    fig_collapse()
    fig_decision()
    table_excursion()
    table_spec()
