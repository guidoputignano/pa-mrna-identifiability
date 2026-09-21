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
from common.numbers import NumberBook
from common.excitation import (phi, excitation, excitation_max, fold_range_for,
                               optimal_band, sd_log_theta, band_levels,
                               fold_range_for_precision,
                               top_of_excursion_for_precision)
from common.style import (save, write_table, figsize, GREEN, RED, BLUE, AMBER,
                          GREY, MUTED, INK, corner, over_x, beside_y, at_point)

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


#: every interpolated number this module draws on a figure, classified for
#: tools/check_artefacts.py:check_figure_numbers. `registered:` is checked;
#: the rest are declarations, and `UNREGISTERED:` is debt the guard prints on
#: every run until a generator owns the number.
DRAWN = {
    "lab": "text: which of the two bands the marker sits on",
    "PAPER_NUMBERS[key]": "registered: E for whichever of the two bands this "
        "marker is on, as E_held and E_info",
    "trough": "registered: the same object write_numbers takes as "
              "train_trough_frac",
    "DOSE_INTERVAL": "input: a model constant from common.params, which the "
                     "constants table registers",
    "np.log(2) / LAMBDA": "input: the dose half-life, likewise a constant of "
                          "the parameter set this figure was drawn at",
    "name": "text: the schedule's name",
    "len(times)": "input: the schedule's own visit count",
    "PAPER_NUMBERS[rkey]": "registered: the Delta/theta range this schedule "
        "spans, under the same delta_range_* key tab_D_designs uses, so the "
        "NumberBook itself is the comparison between figure and table",
    "100 * cv": "input: the load CV this curve was drawn at",
    "nn": "input: the visit count this curve was drawn at",
    "top": "input: the band top this curve was drawn at",
}


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

#: every quantity this module computes that the PAPER QUOTES IN PROSE, keyed
#: by the name the text cites it under. Written to results/numbers/design.tex
#: at the end of a run and quoted as \num{<key>}, so a literal in a sentence
#: is checked the same way a table is.
#: NumberBook, not a plain dict: this module's generators are called in
#: different orders by `make figures` and by the artefact guard, and a key two
#: of them disagree on then takes whichever value ran last. See NumberBook.
PAPER_NUMBERS = NumberBook()

#: the precision each quoted number is printed at, chosen to match the
#: sentence that quotes it. A key with no entry is printed to four decimals.
PAPER_FORMATS = {
    "alias_cost": "{:.0f}", "alias_cost_lo": "{:.0f}", "alias_cost_hi": "{:.0f}",
    "excitation_ratio": "{:.0f}", "load_accuracy_gate": "{:.0f}",
    "ke_needed_factor": "{:.2f}", "ke_needed_factor_lo": "{:.2f}",
    "ke_needed_factor_hi": "{:.2f}",
    #  four decimals throughout, because Table 3 prints four and the same
    #  quantity must not reach the page as 0.1699 in a table and 0.170 in a
    #  sentence about that table
    "sd_monthly": "{:.4f}", "sd_monthly_lo": "{:.4f}", "sd_monthly_hi": "{:.4f}",
    "sd_two_monthly": "{:.4f}", "sd_two_monthly_lo": "{:.4f}",
    "sd_two_monthly_hi": "{:.4f}",
    "sd_aliased": "{:.4f}", "sd_aliased_lo": "{:.4f}", "sd_aliased_hi": "{:.4f}",
    "sd_aliased_off": "{:.4f}", "free_cv2": "{:.3f}",
    "kappa_monthly_pct": "{:.1f}", "kappa_monthly_pct_lo": "{:.1f}",
    "kappa_monthly_pct_hi": "{:.1f}",
    "settle_hours": "{:.2f}", "settle_hours_1pct": "{:.2f}",
    "track_pct": "{:.4f}", "track_pct_1pct": "{:.4f}",
    "sigma": "{:.1f}",
    "timescale_ratio": "{:.0f}",
    "dropped_lo": "{:.2f}", "dropped_hi": "{:.1f}",
    #  two decimals, matching the "range of Delta/theta" column of Table 3
    "delta_range_monthly": "{:.2f}", "delta_range_two_monthly": "{:.2f}",
    "delta_range_dense": "{:.2f}", "delta_range_aliased": "{:.2f}",
    "delta_range_six_monthly": "{:.2f}", "delta_range_one_cycle": "{:.2f}",
    "delta_range_off": "{:.2f}",
    "sd_dense": "{:.4f}", "sd_dense_lo": "{:.4f}", "sd_dense_hi": "{:.4f}",
    "sd_straddle": "{:.4f}", "sd_straddle_lo": "{:.4f}",
    "sd_straddle_hi": "{:.4f}", "straddle_gain": "{:.1f}",
    "sd_one_cycle": "{:.4f}", "sd_one_cycle_lo": "{:.4f}",
    "sd_one_cycle_hi": "{:.4f}",
    "sd_six_monthly": "{:.4f}", "sd_six_monthly_lo": "{:.4f}",
    "sd_six_monthly_hi": "{:.4f}",
    "sd_dense_flat": "{:.4f}", "sd_dense_flat_lo": "{:.4f}",
    "sd_dense_flat_hi": "{:.4f}", "ke_prior_worth": "{:.1f}",
    #  the dose train's own shape, quoted in four places from memory
    "train_peak_frac": "{:.2f}", "train_trough_frac": "{:.3f}",
    "train_ratio": "{:.2f}", "train_drop": "{:.2f}", "phi_drop": "{:.4f}",
    #  the quasi-steady separation, quoted as 0.8 h and 0.81 h and 10^2 and 104
    "tau_hours": "{:.2f}", "settle_hours_1pct": "{:.2f}",
    #  the K_m sweep and the excursion table's decompensation row
    "decomp_E": "{:.2f}", "decomp_frac_pct": "{:.0f}",
    #  the two bands fig_excitation draws, at the precision its legend prints
    "E_held": "{:.3f}", "E_info": "{:.3f}",
}

#: the seed `design_record` has always drawn its intake record from. Kept as
#: the default so nothing computed before September 2026 changes silently.
BASE_SEED = 3

#: every design number this module reports is a MEDIAN over these records with
#: its range quoted beside it. Until September 2026 each was a single draw at
#: BASE_SEED, because no caller ever passed a seed. Mostly that did not matter
#: -- the monthly headline moves 0.11% from seed 3 to the median -- but on the
#: weak designs it does: six-monthly moves 8.6%, and 68.6% with the k_e prior
#: removed. See audit/n3_ke_prior_worth.py, which found it.
SEEDS = tuple(range(24))


def median_spread(values):
    """(median, min, max) over a seed sweep, infinities dropped."""
    v = np.asarray([x for x in values], dtype=float)
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


def tex_spread(med, lo, hi, gate=0.20, digits=4):
    """A median with its seed range, bolded when it clears the gate."""
    if not np.isfinite(med):
        return r"$\infty$"
    body = ("\\textbf{%.*f}" % (digits, med)) if med < gate else ("%.*f" % (digits, med))
    return r"%s\,\scriptsize$[%.*f, %.*f]$" % (body, digits, lo, digits, hi)


def design_record(times, load_cv=None, seed=BASE_SEED):
    """
    Levels and increments for a design, with intake varying day to day.

    Returns (levels, increments). The levels solve the balance at the
    per-sample recorded load, so a varying load moves the marker as it
    physically would rather than being added as measurement error.

    `seed` selects the intake record. It has a default so that existing callers
    keep their numbers, but a number that reaches the paper is swept over
    `SEEDS` and reported as a median with its range; see `median_spread`.
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
    #  E for the two bands the panel draws. decomp_E is a THIRD band, the
    #  decompensation excursion of tab_D_excursion, so neither of these was
    #  registered anywhere and the legend was asserting two numbers about the
    #  design that nothing compared against anything.
    for lo, hi, col, lab, key in [(M_BAR / np.sqrt(2), M_BAR * np.sqrt(2), RED,
                                   "where the patient is held", "E_held"),
                                  (Km / np.sqrt(2), Km * np.sqrt(2), GREEN,
                                   "where the information is", "E_info")]:
        PAPER_NUMBERS[key] = float(excitation(lo, hi, Km))
        ax[0].axvspan(lo, hi, color=col, alpha=0.14)
        ax[0].plot([lo, hi], [phi(lo, Km), phi(hi, Km)], "o-", color=col, ms=4,
                   lw=1.4, label=f"{lab}: $E$ = {PAPER_NUMBERS[key]:.3f}")
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
    #  the window is a TOLERANCE, not a constant, and the paper quotes it at
    #  two: 5% in Section 4.1 and 1% in Section 4.4 and Figure 6's caption.
    #  `track_pct` is the departure AFTER that window, so it moves with the
    #  tolerance too and needs the same treatment -- it was keyed as though it
    #  did not, so whichever call ran last decided its value, and the guard's
    #  fast path and `make figures` called this function a different number of
    #  times in a different order and wrote 0.9861 and 4.9745 for the same
    #  key. tau and the timescale ratio do not depend on the tolerance at all.
    suffix = "" if abs(tol - 0.05) < 1e-12 else "_%dpct" % round(100 * tol)
    PAPER_NUMBERS.update({
        "settle_hours" + suffix: settle * 24.0,
        "track_pct" + suffix: 100.0 * worst,
        "tau_hours": tau * 24.0,
        "timescale_ratio": np.log(2) / LAMBDA / tau,
    })
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
    takes it from 0.1500 to 1.1096, a factor of 7.4 taken record by record
    (6.9 to 7.9 over the 24 records of SEEDS).

    This docstring previously said "from 0.15 to 0.65", which was the rounded
    form of a pair, 0.1474 and 0.6451, quoted in Section 7 and in the
    supplement and REPRODUCIBLE FROM NOWHERE: audit/n3_ke_prior_worth.py
    searched every design Section 7 lists at 24 records each, and 60
    combinations of horizon, sampling interval and phase offset at the seed
    the paper used, and found neither figure. It is withdrawn here too, so the
    code does not assert what the paper has retracted (REVIEW_LOG.md R19).

    The paper therefore claims one prior and lets theta and K_m float.

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


def table_factors(T=730.0):
    """
    The two factors of the sensitivity ratio, on the same record.

    Section 3 typed this table into the .tex. The numbers are computed on the
    two-daily grid over two years, which is the densest design the paper
    prices, and the point of the table is the asymmetry between the factors
    rather than either one.
    """
    from common.excitation import phi
    t = np.arange(0.0, T, 2.0)
    lv = quasi_steady_level(t)
    inc = np.asarray(delta(t), dtype=float)
    ph = phi(lv, P.Km)
    fac = 1.0 + inc / P.theta
    pr = ph * fac
    def row(name, v, ratio=True):
        return [name, f"${v.min():.4f}$--${v.max():.4f}$",
                (f"${v.max() / v.min():.2f}\\times$" if ratio else ""),
                f"${v.max() - v.min():.4f}$"]
    rows = [row(r"$\varphi(m)$, the marker band", ph),
            row(r"$1 + \Delta/\theta$, the dose train", fac)]
    rows.append([r"the product", f"${pr.min():.3f}$--${pr.max():.3f}$", "",
                 r"$\mathbf{%.3f}$" % (pr.max() - pr.min())])
    write_table("tab_D_factors",
                ["factor", "range over the record", "as a ratio", "drop"],
                rows,
                colspec="lrrr",
                note=("the dose train supplies %.1f times the drop the marker "
                      "band does" % ((pr.max() - pr.min()) / (ph.max() - ph.min()))))
    PAPER_NUMBERS.update({
        "excitation_ratio": (pr.max() - pr.min()) / (ph.max() - ph.min()),
        "phi_drop": ph.max() - ph.min(),
        "train_drop": fac.max() - fac.min(),
        "train_ratio": fac.max() / fac.min(),
    })
    print("  wrote tab_D_factors.tex")
    print(f"    ratio of drops: {PAPER_NUMBERS['excitation_ratio']:.1f}")


def table_load_precision(T=730.0):
    """
    What recording the substrate load more accurately is worth.

    The accuracy enters twice and both have to move together: through
    `sigma`, which is the ordinate noise of the fit, and through the
    load-bias prior, which is the same fraction of the load. Holding either
    one fixed while sweeping the other prices something the specification
    does not ask for.
    """
    import model_A.run as LOC
    monthly = np.arange(0.0, T, 30.4)
    rows, raw = [], []
    for frac in (0.01, 0.05, 0.10, 0.20, 0.40):
        sig = LOC.sigma_total(frac, ASSAY_CV)
        sds = []
        for s in SEEDS:
            lv, inc = design_record(monthly, seed=s)
            sds.append(float(np.sqrt(
                _filter(lv, sig, inc, load_frac=frac)[-1][0, 0])))
        med, lo, hi = median_spread(sds)
        raw.append((frac, sig, med, lo, hi))
        rows.append([f"${100 * frac:.0f}\\%$", f"${sig:.1f}$",
                     tex_spread(med, lo, hi)])
    write_table("tab_D_loadprec",
                ["substrate load recorded to", r"$\sigma$",
                 r"$\operatorname{sd}(\ln\hat\theta)$"],
                rows, colspec="lcc",
                note=("medians over %d intake records at a %.0f%% day-to-day "
                      "variability, range in brackets" % (len(SEEDS), 100 * LOAD_CV)))
    print("  wrote tab_D_loadprec.tex")
    for frac, sig, med, lo, hi in raw:
        print(f"    load to {100 * frac:5.1f}%  sigma {sig:7.1f}  "
              f"sd {med:.4f} [{lo:.4f}, {hi:.4f}]")
    #  where the design crosses the gate, by interpolation in log-log, which
    #  is what the text quotes as the recording accuracy the gate requires
    import numpy as _np
    xs = _np.log([r[0] for r in raw])
    ys = _np.log([r[2] for r in raw])
    cross = float(_np.exp(_np.interp(_np.log(0.20), ys, xs)))
    #  the seed spreads as well as the medians, because the registry row for
    #  sigma quotes both and a spread typed beside a generated median is the
    #  same hazard one column over
    PAPER_NUMBERS.update({
        "load_accuracy_gate": 100 * cross,
        "sd_load_acc_1pct": raw[0][2],
        "sd_load_acc_1pct_lo": raw[0][3], "sd_load_acc_1pct_hi": raw[0][4],
        "sd_load_acc_20pct": raw[3][2],
        "sd_load_acc_20pct_lo": raw[3][3], "sd_load_acc_20pct_hi": raw[3][4],
        "sd_load_acc_40pct": raw[4][2],
        "sd_load_acc_40pct_lo": raw[4][3], "sd_load_acc_40pct_hi": raw[4][4],
    })
    print(f"    reaches 0.20 at a recording accuracy of {100 * cross:.1f}%")
    return raw, cross


def table_instruments(T=730.0):
    """
    Three sources of excitation on one metric.

    Same footing throughout: theta and K_m free, a log-prior of 0.8 on ln k_e,
    the load recorded to 5%, twenty-five monthly draws over two years. The
    first two rows carry no therapy at all, so their levels solve the balance
    without the increment; the carnitine row adds recorded steps in k_e.
    """
    import model_A.run as LOC
    monthly = np.arange(0.0, T, 30.4)
    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
    Lb = load(0.0)

    def med_of(fn):
        return median_spread([fn(s) for s in SEEDS])

    def sd_untreated(load_cv, ke_step=0.0):
        """
        No therapy at all, so the levels solve the balance without the
        increment. A carnitine step is a RECORDED step in k_e on alternate
        visits with its gain treated as known, which is the favourable case
        for carnitine, and the level is re-solved at the stepped k_e: the
        step moves the patient as well as the column, and moving the patient
        is most of what it does (audit/e7_excitation_budget.py).
        """
        def one(s):
            rng = np.random.default_rng(s)
            Ls = Lb * (1.0 + load_cv * rng.standard_normal(monthly.size))
            ke = P.ke * (1.0 + ke_step * (np.arange(monthly.size) % 2))
            B = ke * P.Km + P.theta - Ls
            lv = (-B + np.sqrt(B ** 2 + 4.0 * ke * Ls * P.Km)) / (2.0 * ke)
            F = np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * Lb) ** 2, 1.0 / 0.8 ** 2])
            for m, k in zip(lv, ke):
                g = m / (P.Km + m)
                q = np.array([P.theta * g,
                              -P.theta * P.Km * m / (P.Km + m) ** 2, 1.0, k * m])
                F = F + np.outer(q, q) / sigma ** 2
            w = np.linalg.eigvalsh(F)
            if w[0] <= 1e-14 * max(w[-1], 1e-300):
                return np.inf
            return float(np.sqrt(np.linalg.inv(F)[0, 0]))
        return med_of(one)

    def sd_dose(s):
        lv, inc = design_record(monthly, seed=s)
        return float(np.sqrt(_filter(lv, sigma, inc)[-1][0, 0]))

    gate = 0.20
    spont10 = sd_untreated(0.10)
    spont40 = sd_untreated(0.40)
    carn_lo = sd_untreated(0.10, ke_step=0.25)
    carn_hi = sd_untreated(0.10, ke_step=1.00)
    train = med_of(sd_dose)
    rows = [
        [r"spontaneous load variability alone, at $10\%$ day-to-day CV",
         f"${spont10[0]:.2f}$", r"$%.1f\times$ worse" % (spont10[0] / gate)],
        [r"\quad the same at an implausible $40\%$ CV",
         f"${spont40[0]:.2f}$", r"$%.1f\times$ worse" % (spont40[0] / gate)],
        [r"plus recorded carnitine steps of $25$--$100\%$ in $k_e$",
         f"${min(carn_lo[0], carn_hi[0]):.2f}$--${max(carn_lo[0], carn_hi[0]):.2f}$",
         r"\emph{worse still}"],
        [r"\textbf{the dose train}", r"$\mathbf{%.4f}$" % train[0],
         r"\textbf{passes}"],
    ]
    write_table("tab_D_instruments",
                ["source of excitation",
                 r"$\operatorname{sd}(\ln\hat\theta)$",
                 r"against $0.20$"],
                rows,
                #  the first column holds sentences, so it has to wrap
                colspec=(r">{\raggedright\arraybackslash}p{0.56\linewidth}cc"),
                note=("medians over %d intake records, one footing throughout"
                      % len(SEEDS)))
    PAPER_NUMBERS.update({
        "sd_spont_10": spont10[0], "sd_spont_40": spont40[0],
        "sd_carn_lo": min(carn_lo[0], carn_hi[0]),
        "sd_carn_hi": max(carn_lo[0], carn_hi[0]),
    })
    print("  wrote tab_D_instruments.tex")
    for lab, v in (("load variability 10%", spont10), ("load variability 40%", spont40),
                   ("carnitine 25%", carn_lo), ("carnitine 100%", carn_hi),
                   ("the dose train", train)):
        print(f"    {lab:>22}: {v[0]:.4f} [{v[1]:.4f}, {v[2]:.4f}]")


def fig_dose_train(T=91.0):
    """
    The experiment, drawn: the dose train, what it does to the increment and
    to the marker, and where two schedules land in it.

    The paper's thesis is that the therapy is already the identification
    experiment, and a reader has to be able to see the experiment. The top
    panel is the increment the schedule delivers, the bottom is the marker
    it drives, and the two schedules are marked on both: monthly across dose
    cycles, which walks around the cycle, and fortnightly, which is the dose
    interval itself and therefore lands at one phase for ever.
    """
    t = np.linspace(0.0, T, 4000)
    inc = np.asarray(delta(t), dtype=float) / P.theta
    lv = quasi_steady_level(t)

    monthly = np.arange(0.0, T, 30.4)
    fortnightly = np.arange(0.0, T, DOSE_INTERVAL)
    doses = np.arange(0.0, T, DOSE_INTERVAL)

    fig, ax = plt.subplots(2, 1, figsize=figsize(1.0, 0.58), sharex=True)

    #  (top) the increment the schedule delivers
    ax[0].plot(t, inc, color=GREEN, lw=1.8)
    for d in doses:
        ax[0].plot([d, d], [0, inc.max() * 1.04], color=GREY, lw=0.7, ls=":",
                   zorder=0)
    trough = float(inc[(t > doses[-2]) & (t < doses[-1])].min())
    ax[0].axhline(trough, color=MUTED, ls="--", lw=0.9)
    ax[0].annotate(r"trough  $%.2f\,\theta$" % trough,
                   xy=(0.30 * T, trough), xytext=(0.30 * T, trough + 0.28),
                   color=MUTED, fontsize=7, ha="left", va="bottom",
                   arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    ax[0].set_ylabel(r"increment  $\Delta(t)/\theta$")
    ax[0].set_ylim(0, inc.max() * 1.12)
    corner(ax[0], "a dose every %.0f d, half-life %.1f d"
           % (DOSE_INTERVAL, np.log(2) / LAMBDA), "upper right",
           color=MUTED, fontsize=7)

    #  (bottom) the marker the increment drives, and where each schedule lands
    ax[1].plot(t, lv, color=BLUE, lw=1.8)
    for times, col, mk, lab in ((monthly, GREEN, "o", "monthly, across cycles"),
                                (fortnightly, RED, "s",
                                 "fortnightly = the dose interval")):
        y = quasi_steady_level(times)
        ax[1].plot(times, y, mk, color=col, ms=5.5, mec="white", mew=0.8,
                   label=lab, zorder=5)
        for x, yy in zip(times, y):
            ax[1].plot([x, x], [lv.min() * 0.97, yy], color=col, lw=0.8,
                       alpha=0.55, zorder=1)
    ax[1].set_xlabel("days")
    ax[1].set_ylabel(r"marker  $m$  [$\mu$M]")
    #  headroom above the peaks so the legend sits in clear air rather than
    #  over the curve it is labelling
    ax[1].set_ylim(lv.min() * 0.97, lv.max() * 1.30)
    ax[1].legend(loc="upper center", ncol=2, fontsize=7, framealpha=0.0,
                 borderpad=0.2, handletextpad=0.5, columnspacing=1.6)

    #  the phase each schedule visits over the WHOLE record, which is the
    #  point: three months is what fits on a page, two years is the design
    full = np.arange(0.0, 730.0, 30.4)
    ph_m = np.sort((full % DOSE_INTERVAL) / DOSE_INTERVAL)
    ph_f = (np.arange(0.0, 730.0, DOSE_INTERVAL) % DOSE_INTERVAL) / DOSE_INTERVAL
    PAPER_NUMBERS.update({
        "phase_gap_monthly": float(np.max(np.diff(np.r_[ph_m, ph_m[0] + 1.0]))),
        "phase_spread_fortnightly": float(ph_f.max() - ph_f.min()),
    })
    for a in ax:
        a.set_xlim(0, T)
    fig.align_ylabels(ax)
    save(fig, "fig_D0_dosetrain")
    PAPER_NUMBERS.update({
        "train_trough_frac": trough,
        "train_peak_frac": float(inc.max()),
        "marker_low": float(lv.min()),
        "marker_high": float(lv.max()),
    })
    print(f"    increment runs {inc.min():.3f} to {inc.max():.3f} of theta, "
          f"trough {trough:.3f}; marker {lv.min():.2f} to {lv.max():.2f}")
    print(f"    over two years the monthly schedule leaves no gap in cycle "
          f"phase wider than {PAPER_NUMBERS['phase_gap_monthly']:.2f} of a "
          f"cycle; the fortnightly one spans "
          f"{PAPER_NUMBERS['phase_spread_fortnightly']:.2f}")


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
    #  the target label sat under its own line, where the dense design's
    #  curve runs; above the line is clear
    ax[0].annotate(r"$20\%$ on $\theta$", xy=(0.34, 0.21),
                   xycoords=("axes fraction", "data"),
                   color=MUTED, fontsize=7, ha="left", va="bottom")
    ax[0].set_xlabel("time  [days]")
    ax[0].set_ylabel(r"$\mathrm{sd}(\ln\hat\theta)$")
    ax[0].set_title(r"no process noise on $\theta$:" "\n"
                    "what does not collapse, stays", fontsize=8)
    #  both panels carry the same three series and the left panel's legend
    #  sat on top of the six-monthly curve, so one legend goes under the
    #  figure instead of two inside it

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
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
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
    #  no sd field: it was carried here as a hand-typed string, the aliased
    #  entry still reading 3.51, which is that design's seed-3 MINIMUM and not
    #  its median of 4.5845. Nothing drew it -- the panels compute their own
    #  range -- but a stale number in a figure is invisible to the artefact
    #  guard, which compares rendered bytes, so it does not get to live here.
    #  the fifth field is the delta_range_* key each panel's range bar draws.
    #  These three schedules are exactly three rows of tab_D_designs over the
    #  same two-year grid, so the bar and the table are ONE number: writing it
    #  into PAPER_NUMBERS under that key makes the NumberBook the comparison,
    #  and a figure that drifted from the table would fail at accumulation.
    SCHEDULES = [
        ("monthly, every 30.4 d", 30.4, GREEN, "o", "monthly"),
        ("every 14.5 d, half a day off", 14.5, BLUE, "^", "off"),
        ("fortnightly, on the dose", DOSE_INTERVAL, RED, "X", "aliased"),
    ]

    # Nested gridspecs rather than one mosaic: the trajectory and its rails are
    # one object and want almost no gap, while the three increment panels want
    # a real one. A single grid cannot hold both spacings.
    fig = plt.figure(figsize=figsize(1.0, 0.62))
    fig.set_layout_engine("none")
    outer = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.34,
                             left=0.085, right=0.90, top=0.93, bottom=0.105)
    gl = outer[0, 0].subgridspec(2, 1, height_ratios=[2.3, 1.0], hspace=0.04)
    gr = outer[0, 1].subgridspec(3, 1, hspace=0.28)
    ax = {"traj": fig.add_subplot(gl[0]), "rails": fig.add_subplot(gl[1])}
    for i in range(3):
        ax[f"p{i}"] = fig.add_subplot(gr[i])

    # ---------------- the marker under the dose train ----------------
    tt = np.linspace(0.0, T, 2000)
    a = ax["traj"]
    mm = quasi_steady_level(tt)
    a.plot(tt, mm, color=INK, lw=1.5, zorder=3)
    a.set_ylabel(r"marker $m(t)$")
    a.set_xlim(0, T)
    a.set_xticklabels([])
    # a clear strip below the trough, so the dose label sits under the
    # trajectory rather than across it
    lo, hi = float(mm.min()), float(mm.max())
    span = hi - lo
    a.set_ylim(lo - span * 0.20, hi + span * 0.05)
    doses = np.arange(0.0, T + 1e-9, DOSE_INTERVAL)
    a.plot(doses, quasi_steady_level(doses), linestyle="none", marker="v",
           ms=3.8, color=AMBER, mec="white", mew=0.5, zorder=4)
    a.text(T * 0.012, lo - span * 0.15, "dose, every 14 d",
           fontsize=6.3, color=AMBER, ha="left", va="center")
    a.set_title("the schedule in time", fontsize=7.4)

    # ---------------- the three schedules as rails ----------------
    r = ax["rails"]
    r.set_xlim(0, T)
    r.set_xlabel("time  [days]")
    r.set_ylim(-0.7, len(SCHEDULES) - 0.3)
    r.set_yticks(range(len(SCHEDULES)))
    r.set_yticklabels([s[0].split(",")[0] for s in reversed(SCHEDULES)], fontsize=6.4)
    for i, (lab, dt, col, mk, _key) in enumerate(SCHEDULES):
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

    for i, (lab, dt, col, mk, key) in enumerate(SCHEDULES):
        b = ax[f"p{i}"]
        ts = np.arange(0.0, 730.0, dt)
        dph = np.mod(ts, DOSE_INTERVAL)
        dval = np.asarray(delta(ts)) / P.theta
        rkey = f"delta_range_{key}"
        PAPER_NUMBERS[rkey] = float(np.ptp(dval))

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
               f"range\n{PAPER_NUMBERS[rkey]:.2f}", fontsize=5.9, color=col,
               va="center", ha="left")

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
                   label=f"{100*cv:.0f}%")
    ax[0].set_ylabel(r"sign error on $\mu$  [%]")
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
        #  headroom above the peaks so the per-panel legends sit in clear
        #  air rather than on the curves they label. The same limit on all
        #  three, because the panels are meant to be read against each other
        a.set_ylim(0, 96)
    #  the boundary label sat at the top of its own line, where all three
    #  curves peak, and the mu was lost behind them. Low and to the left is
    #  empty in every panel, because that is where the error rate is zero.
    ax[0].annotate(r"$\mu = 0$", xy=(cross, 0.0), xycoords="data",
                   xytext=(0.06, 0.14), textcoords="axes fraction",
                   color=MUTED, fontsize=7, ha="left", va="center",
                   arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6,
                                   shrinkA=1, shrinkB=1))
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
        # one intake record per seed; every number below is a median with its
        # range, because a single draw is not a property of the design
        sds, rhos, drops = [], [], []
        for s in SEEDS:
            lv, inc = design_record(times, seed=s)
            if override is not None:
                lv = override
            Pc = _filter(lv, sigma, inc)[-1]
            sds.append(float(np.sqrt(Pc[0, 0])))
            rhos.append(float(Pc[0, 1] / np.sqrt(Pc[0, 0] * Pc[1, 1])))
            drops.append(float(np.sqrt(_filter(lv, sigma)[-1][0, 0])))
        sd, sd_lo, sd_hi = median_spread(sds)
        rho, _, _ = median_spread(rhos)
        dropped, _, _ = median_spread(drops)
        lv, inc = design_record(times, seed=BASE_SEED)
        if override is not None:
            lv = override
        drng = float((inc.max() - inc.min()) / P.theta)
        raw.append((name, len(times), lv.min(), lv.max(), drng, sd, rho,
                    dropped, sd_lo, sd_hi))
        rows.append([name, len(times), f"${lv.min():.2f}$--${lv.max():.2f}$",
                     f"${drng:.2f}$", tex_spread(sd, sd_lo, sd_hi),
                     f"${rho:+.4f}$", f"${dropped:.3f}$"])
    head = ["design", "$n$", "levels visited",
            r"range of $\Delta/\theta$",
            r"$\operatorname{sd}(\ln\hat\theta)$",
            r"$\operatorname{corr}(\ln\hat\theta, \ln \hat K_m)$"]
    note = ("theta and K_m free; priors only on ln k_e and the load record. "
            "Medians over %d intake records, range in brackets." % len(SEEDS))
    #  Seven columns do not fit the page: the last one ran off the right edge
    #  and printed as a truncated word. The main text carries six, and the
    #  column it drops -- the same design evaluated on a footing this paper
    #  does not use -- goes to the supplement whole, where the validation
    #  section reads it.
    write_table("tab_D_designs", head, [r[:6] for r in rows], note=note)
    print("  wrote tab_D_designs.tex")
    write_table("tab_D_designs_full", head + [r"\emph{$\Delta$ dropped}"],
                rows, note=note)
    print("  wrote tab_D_designs_full.tex")
    #  what the dropped-increment column spans over the designs a protocol
    #  could adopt, which is what the text quotes now that the column itself
    #  has moved to the supplement. `raw` is
    #  (name, n, lv_min, lv_max, range, sd, corr, dropped, sd_lo, sd_hi).
    acc = [r[7] for r in raw if "straddling" not in r[0]]
    PAPER_NUMBERS.update({"dropped_lo": min(acc), "dropped_hi": max(acc)})
    #  the range of Delta/theta each design sees, which the text quoted from
    #  memory in nine different values. One key per design, and the half-day
    #  offset added here because the text quotes it and the table does not.
    by_name = {"monthly": r"\textbf{monthly, across dose cycles}",
               "two_monthly": "two-monthly, across dose cycles",
               "dense": "every two days for two years",
               "aliased": r"\emph{fortnightly: the dose interval itself}",
               "six_monthly": "six-monthly clinic visits",
               "one_cycle": "one dose cycle, every two days"}
    dense = next(r for r in raw if r[0] == "every two days for two years")
    PAPER_NUMBERS.update({"sd_dense": dense[5], "sd_dense_lo": dense[8],
                          "sd_dense_hi": dense[9]})
    #  the same design with the k_e prior taken away, which is what that prior
    #  is worth and what the registry's s_lnke row quotes. A width of 1e9 is
    #  flat here and NOT the same as ke_prior=None, which holds k_e exactly
    #  known and drops the column altogether.
    flat = []
    for s in SEEDS:
        lv, inc = design_record(dense_t, seed=s)
        flat.append(float(np.sqrt(
            _filter(lv, sigma, inc, ke_prior=1e9)[-1][0, 0])))
    flat = median_spread(flat)
    PAPER_NUMBERS.update({"sd_dense_flat": flat[0], "sd_dense_flat_lo": flat[1],
                          "sd_dense_flat_hi": flat[2],
                          "ke_prior_worth": flat[0] / dense[5]})
    lut = {r[0]: r[4] for r in raw}
    PAPER_NUMBERS.update({f"delta_range_{k}": lut[v]
                          for k, v in by_name.items() if v in lut})
    #  the two rows the text quotes that spec_numbers() does not recompute:
    #  the K_m-straddling design, which is the placement argument, and the
    #  single dose cycle. Both were hand-typed in Section 4, one of them to
    #  three decimals against the table's four.
    sd_lut = {r[0]: (r[5], r[8], r[9]) for r in raw}
    for key, nm in (("straddle", r"\emph{the same $365$ straddling $K_m$}"),
                    ("one_cycle", "one dose cycle, every two days"),
                    ("six_monthly", "six-monthly clinic visits")):
        if nm in sd_lut:
            med, slo, shi = sd_lut[nm]
            PAPER_NUMBERS.update({f"sd_{key}": med, f"sd_{key}_lo": slo,
                                  f"sd_{key}_hi": shi})
    if "sd_straddle" in PAPER_NUMBERS:
        PAPER_NUMBERS["straddle_gain"] = dense[5] / PAPER_NUMBERS["sd_straddle"]
    off = design_record(np.arange(0.0, T, DOSE_INTERVAL + 0.5),
                        seed=BASE_SEED)[1]
    PAPER_NUMBERS["delta_range_off"] = float((off.max() - off.min()) / P.theta)
    for nm, n, lo, hi, drng, sd, rho, drop, slo, shi in raw:
        print(f"    {nm[:40]:>40}: n={n:4d}  range(D/th)={drng:5.2f}  "
              f"sd={sd:.4f} [{slo:.4f}, {shi:.4f}]"
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
    for name, lo, hi, _d in excursion_sources():
        if "decompensation" in name:
            PAPER_NUMBERS.update({
                "decomp_E": excitation(lo, hi, P.Km),
                "decomp_frac_pct": 100.0 * excitation(lo, hi, P.Km)
                / excitation_max(hi / lo)})
    write_table("tab_D_excursion",
                ["source", "$m$ range", "fold", "direction",
                 "$\\sqrt{m_{\\min}m_{\\max}}/K_m$", "$E$", "$E_{\\max}(r)$",
                 "frac.\\ of $E_{\\max}$"],
                rows,
                #  ragged rather than justified: at this width LaTeX was
                #  hyphenating the source names mid-word to make the lines
                #  meet, which is unreadable in a column of labels
                colspec=r">{\raggedright\arraybackslash}p{0.27\textwidth}lllllll",
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

    def sd_of(times, kap=None, load_cv=None, ke_prior=0.8, load_prior=True,
              seed=BASE_SEED):
        lv, inc = design_record(times, load_cv=load_cv, seed=seed)
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

    def over_seeds(**kw):
        """
        Median and range over SEEDS, because one intake record is not a
        property of the design. Until September 2026 every number below was a
        single draw at BASE_SEED, and on the weak designs that mattered: the
        aliased row's seed-3 value is the MINIMUM of its own 24-seed range.
        """
        return median_spread([sd_of(seed=s, **kw) for s in SEEDS])

    monthly = np.arange(0.0, 730.0, 30.4)
    two_monthly = np.arange(0.0, 730.0, 60.0)
    spread = {
        "sd_monthly": over_seeds(times=monthly),
        "sd_two_monthly": over_seeds(times=two_monthly),
        "sd_aliased": over_seeds(times=np.arange(0.0, 730.0, DOSE_INTERVAL)),
        "sd_aliased_off": over_seeds(
            times=np.arange(0.0, 730.0, DOSE_INTERVAL + 0.5)),
        "sd_load_1pct": over_seeds(times=monthly, load_cv=0.01),
        "sd_load_20pct": over_seeds(times=monthly, load_cv=0.20),
        "sd_load_steady": over_seeds(times=monthly, load_cv=0.0),
        # the same, with NO prior on the load bias, which is where the rank
        # deficiency of prop:rank actually shows: the prior fills that
        # direction, and only real variation lets the DATA fill it instead
        "free_steady": over_seeds(times=monthly, load_cv=0.0, load_prior=False),
        "free_cv2": over_seeds(times=monthly, load_cv=0.02, load_prior=False),
        "free_cv10": over_seeds(times=monthly, load_cv=0.10, load_prior=False),
        "free_cv20": over_seeds(times=monthly, load_cv=0.20, load_prior=False),
    }
    N = {"sigma": sigma}
    #: the headline keys keep their names and now carry the MEDIAN; the
    #: <key>_lo / <key>_hi keys carry the range that median came from, and the
    #: paper quotes both wherever it quotes the median.
    for k, (med, lo, hi) in spread.items():
        N[k], N[k + "_lo"], N[k + "_hi"] = med, lo, hi
    #: the ratio the aliasing claim rests on, taken per record and then
    #: summarised, rather than as a ratio of two separately-taken medians
    ratios = [sd_of(times=np.arange(0.0, 730.0, DOSE_INTERVAL), seed=s)
              / sd_of(times=monthly, seed=s) for s in SEEDS]
    (N["alias_cost"], N["alias_cost_lo"],
     N["alias_cost_hi"]) = median_spread(ratios)
    def solve_over_seeds(f):
        """
        brentq per intake record, then the median and range of the answers.

        These are the requirements Tier 2 is written against, so a single-draw
        answer here is a single-draw REQUIREMENT, which is worse than a
        single-draw precision.
        """
        out = []
        for s in SEEDS:
            try:
                out.append(brentq(lambda x: f(x, s) - 0.20, 1e-3, 3.0,
                                  xtol=1e-7))
            except ValueError:
                pass
        return median_spread(out) if out else (None, None, None)

    for tag, tt in (("monthly", monthly), ("two_monthly", two_monthly)):
        med, lo, hi = solve_over_seeds(
            lambda x, s, tt=tt: sd_of(tt, kap=x, seed=s))
        N["kappa_" + tag] = med
        N["kappa_" + tag + "_lo"], N["kappa_" + tag + "_hi"] = lo, hi
    med, lo, hi = solve_over_seeds(
        lambda x, s: sd_of(monthly, ke_prior=x, seed=s))
    N["ke_needed"], N["ke_needed_lo"], N["ke_needed_hi"] = med, lo, hi
    return N


def table_spec():
    """
    The one-page data specification, in three tiers.

    The order is not cosmetic and it is not cost: it is dependency. Tier 1 is
    one blood draw and it rules out a family of structures under which theta
    is not identifiable at all, so the tiers below it are calibrating a
    parameter that may not exist until it reports. Tier 2 then measures the
    only quantity the method's feasibility turns on, and until THAT reports
    there is no defensible number of breath tests to ask for.
    """
    N = spec_numbers()
    PAPER_NUMBERS.update(N)
    if N.get("ke_needed"):
        PAPER_NUMBERS["ke_needed_factor"] = float(np.exp(N["ke_needed"]))
        PAPER_NUMBERS["ke_needed_factor_lo"] = float(np.exp(N["ke_needed_lo"]))
        PAPER_NUMBERS["ke_needed_factor_hi"] = float(np.exp(N["ke_needed_hi"]))
    for tag in ("", "_lo", "_hi"):
        v = N.get("kappa_monthly" + tag)
        if v is not None and np.isfinite(v):
            PAPER_NUMBERS["kappa_monthly_pct" + tag] = 100.0 * v
    kap = N["kappa_monthly"]
    kefac = np.exp(N["ke_needed"]) if N["ke_needed"] else None

    def rng(key, digits=3, scale=1.0, unit=""):
        """
        The seed range beside a median, in the units the median is quoted in.

        `scale` and `unit` exist because some of these are percentages and a
        range printed in the wrong units is worse than no range at all.
        """
        lo, hi = N.get(key + "_lo"), N.get(key + "_hi")
        if lo is None or not np.isfinite(lo) or not np.isfinite(hi):
            return ""
        return (r", $%.*f$ to $%.*f$%s across %d records"
                % (digits, scale * lo, digits, scale * hi, unit, len(SEEDS)))
    B = "\\textbf{"
    rows = [
        [B + "TIER 1} --- run first, and it is one draw", "", ""],
        ["what to measure",
         r"two or three plasma draws inside the first half hour after one "
         r"scheduled infusion, at the assay precision assumed throughout",
         "rules out the published alternative and any structure whose marker "
         "has no relaxation of its own, under which $\\theta$ is not "
         "identifiable at any design; see Section~\\ref{sec:tierstruct}"],
        ["what it must reach",
         r"a draw within fifteen minutes of the dose; the separation is "
         r"$16.3\sigma$ at three minutes and $0.68\sigma$ at one hour",
         "the two structures share the steady map exactly and differ only in "
         "the travel between levels, so the whole of the test is in the "
         "first hour"],
        [B + "TIER 2} --- conditional on Tier 1", "", ""],
        ["what to measure",
         "same-day repeatability of $1$-$^{13}$C-propionate oxidation: two "
         "breath tests bracketing one dose, on one day",
         "the only quantity the design turns on, and no published value for "
         "it exists"],
        ["what it must reach",
         ("a coefficient of variation on $\\Delta$ below $%.1f\\%%$%s"
          % (100 * kap, rng("kappa_monthly", digits=1, scale=100.0,
                            unit=r"\%"))
          if kap else "a CV on $\\Delta$ inside the gate"),
         "derived from the gate, not chosen: at that precision monthly "
         "sampling reaches $20\\%$ on $\\theta$. This is a \\emph{requirement}, "
         "so it is quoted as a median over intake records rather than as one "
         "draw"],
        ["what is known now",
         "between-visit CV, mean $14.0\\%$, range $0.4$--$50.8\\%$, $n=15$, "
         "intervals $2$ months to $4$ years",
         "an upper bound on same-day error, already close to the requirement, "
         "but it does not establish it; and it is MMA, not PA"],
        [B + "TIER 3} --- conditional on Tier 2", "", ""],
        ["marker", "plasma 3-hydroxypropionate",
         "endpoint association (RR $0.7$, $95\\%$ CI $0.6$--$0.9$) and "
         "immunity to the carnitine sign reversal that disqualifies C3"],
        ["sampling schedule",
         "monthly across dose cycles, $n \\approx 25$ over two years "
         "($%.4f$%s); two-monthly also passes ($%.4f$%s). These are "
         "\\emph{sampling spreads} of $\\ln\\hat\\theta$, not interval "
         "half-widths; see Remark~\\ref{rem:curvature}"
         % (N["sd_monthly"], rng("sd_monthly", digits=4),
            N["sd_two_monthly"], rng("sd_two_monthly", digits=4)),
         "ordinary outpatient density; what matters is covering the phase of "
         "the dose cycle, not the sampling rate"],
        ["sampling, the hard constraint",
         B + "never at the dose interval or a multiple of it}: fortnightly "
         "gives $%.4f$%s against the monthly design's $%.4f$, a cost of "
         "$%.0f$-fold ($%.0f$ to $%.0f$ across records); $14.5$-daily gives "
         "$%.4f$%s, which recovers it"
         % (N["sd_aliased"], rng("sd_aliased", digits=4), N["sd_monthly"],
            N["alias_cost"], N["alias_cost_lo"], N["alias_cost_hi"],
            N["sd_aliased_off"], rng("sd_aliased_off", digits=4)),
         "every draw lands at the same phase, the increment stops varying, "
         "and the therapy channel closes"],
        #  the exclusion is on the FIT, not on the draw. A first-hour sample
        #  is good data; it is inadmissible in a steady-state fit and in no
        #  other sense, and the row used to say "no sample within one hour",
        #  which forecloses the only window any dynamic estimator could use.
        ["sampling, the other exclusion",
         r"\textbf{do not put a draw taken within one hour of a dose into "
         r"the steady-state fit}. Such a draw may be taken, and is the only "
         r"window that distinguishes the two structures of "
         r"Section~\ref{sec:structuredep}; it is excluded from THIS fit, not "
         r"from the protocol",
         "the marker is not at quasi-steady state there and a steady-state "
         "fit does not apply to it. That is a constraint on the estimator "
         "and not on the patient"],
        ["dosing record",
         "date, dose and product lot for every infusion",
         "$\\Delta(t_i)$ is a regressor, not a nuisance: it carries the "
         "excitation"],
        ["load record, precision",
         "ordinary outpatient \\emph{variability} suffices ($%.4f$ at a $1\\%%$ "
         "day-to-day CV against $%.4f$ at $20\\%%$); what binds is the "
         "\\emph{per-visit recording error}, and it is the row below. It is "
         "\\textbf{not} the mean offset: $98\\%%$ of the degradation at a $20\\%%$ "
         "recording error is per-visit scatter and $2\\%%$ is offset, so a "
         "weighed-intake calibration of the mean is not a substitute for a "
         "repeated as-eaten record"
         % (N["sd_load_1pct"], N["sd_load_20pct"]),
         "the offset is the direction the fibre runs along, and a varying "
         "recorded load closes it from the data, leaving its prior inert"],
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
         ("a prior within a factor of $%.2f$ ($%.2f$ to $%.2f$ across records)"
          % (kefac, np.exp(N["ke_needed_lo"]), np.exp(N["ke_needed_hi"]))
          if kefac else "a prior of some width"),
         B + "no published source supplies it}, so this is the one number the "
         "design asks for and cannot cite. It is also the cheapest second "
         "external measurement: tightening it is worth $2.1\\times$ at nominal "
         "capacity and $6.9\\times$ in the disease range, though even an exact "
         "$k_e$ leaves $0.388$  there"],
        ["$K_m$", "nothing at all",
         "every number here is unchanged by a $K_m$ prior twenty times "
         "looser, which is what Remark~\\ref{rem:no_km} demanded"],
    ]
    write_table("tab_D_spec", ["", "what the record must contain", "why"], rows,
                note="three tiers, ordered by dependency: Tier 2 is "
                     "conditional on Tier 1 ruling out the alternative "
                     "structure, and Tier 3 on Tier 2 clearing the gate",
                #  every column ragged. Justified at these widths LaTeX
                #  hyphenates the row labels mid-word ("preci-sion",
                #  "vari-ability") and opens gaps inside the sentences, which
                #  is what made this unreadable as a table rather than the
                #  content being wrong for one.
                colspec=(r">{\raggedright\arraybackslash}p{0.19\textwidth}"
                         r"@{\hskip 0.8em}"
                         r">{\raggedright\arraybackslash}p{0.39\textwidth}"
                         r"@{\hskip 0.8em}"
                         r">{\raggedright\arraybackslash}p{0.38\textwidth}"))
    print("  wrote tab_D_spec.tex")
    for k, v in N.items():
        print("    %18s = %s" % (k, v))
    return N


if __name__ == "__main__":
    print("DESIGN")
    fig_excitation()
    #  the settling window and the two time constants are quoted in prose,
    #  not drawn, so nothing else here would produce them and the numbers
    #  file would come out five keys short of what the paper inputs. The
    #  paper quotes the window at two tolerances, 5% and 1%, so both run.
    quasi_steady_check()
    quasi_steady_check(tol=0.01)
    fig_dose_train()
    fig_collapse()
    fig_decision()
    fig_protocol()
    table_designs()
    table_excursion()
    table_factors()
    table_load_precision()
    table_instruments()
    table_spec()
    from common.numbers import write_numbers
    write_numbers("design", PAPER_NUMBERS, PAPER_FORMATS)
