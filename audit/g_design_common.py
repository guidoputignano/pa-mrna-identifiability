"""
Shared pieces for the Phase 3 design figures, so g1, g2 and g3 cannot drift.

Everything here is stated in the paper's own coordinates. The information
matrix is built in (ln theta, ln K_m, b, ln k_e), which is what
src/design/run.py:324-329 differentiates in and what prop:rank is stated in
(06_excitation.tex:332). The k_e convention is `sd_from`'s, from
audit/e9_rank_vs_precision.py:49-58: ke_prior=None means k_e is held EXACTLY
KNOWN and the column is DROPPED, not given a flat prior.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from design.run import P, delta, load, DOSE_INTERVAL, LOAD_FRAC, ASSAY_CV
import model_A.run as LOC

TH, KM, KE = P.theta, P.Km, P.ke
LB = load(0.0)
SIG = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
T = DOSE_INTERVAL
GATE = 0.20

#: clinic visits land on a half-day grid, which is the resolution the paper's
#: own 14.5 d schedule already assumes.
GRID = 0.5


# --------------------------------------------------------------------------
def fixed(n, step, t0=0.0):
    """n visits at a fixed interval."""
    return t0 + step * np.arange(n, dtype=float)


def phase_spread(n, gap=30.4, t0=0.0, period=T, grid=GRID, phases=None):
    """
    n visits whose dose-cycle phases are as evenly spread as the visit grid
    allows, subject to a minimum gap between visits.

    The construction is the inverse map of g3: a target phase is chosen first
    and a calendar day is then found for it. The targets are visited in a
    GREEDY order, nearest-reachable first, rather than in numerical order.
    That matters: taking them in order forces a wait of up to a whole cycle
    between visits and stretches the record by about 40% at equal n, which
    would hand the comparison to the phase-spread design for a reason that has
    nothing to do with phase.

    `phases` supplies the target phases explicitly, for a caller that cannot
    use the whole dial. Left None the targets are n phases evenly spaced
    around the cycle, which is what g1 and g2 use and what this returned
    before the argument existed.
    """
    targets = (list((np.arange(n, dtype=float) * period / n) % period)
               if phases is None else [float(p) % period for p in phases])
    out, remaining = [], list(targets)
    for k in range(n):
        earliest = t0 if k == 0 else out[-1] + gap
        base = np.ceil(earliest / grid) * grid
        best, best_t = None, None
        for p in remaining:
            cyc = np.floor(base / period) * period
            cand = np.round(np.array([cyc + p, cyc + period + p]) / grid) * grid
            cand = cand[cand >= base - 1e-9]
            if not cand.size:
                continue
            if best_t is None or cand[0] < best_t:
                best, best_t = p, float(cand[0])
        if best is None:
            best, best_t = remaining[0], float(base)
        remaining.remove(best)
        out.append(best_t)
    return np.array(out)


#: far enough in that the dose train has reached steady state. A patient who
#: is ALREADY ON THERAPY when enrolment starts is at this end of the train, not
#: at t = 0, and for a fortnightly schedule that is the whole difference
#: between "nearly singular" and "exactly singular".
STEADY = 100 * T

SCHEDULES = {
    "14 d (the dose interval)": lambda n: fixed(n, 14.0),
    "14 d, already on therapy": lambda n: fixed(n, 14.0, t0=STEADY),
    "14.5 d": lambda n: fixed(n, 14.5),
    "monthly (30.4 d)": lambda n: fixed(n, 30.4),
    "phase-spread": lambda n: phase_spread(n, gap=30.4),
}


# --------------------------------------------------------------------------
def record(times, load_cv=0.0, seed=0):
    """
    Levels, increments and the per-sample recorded load.

    The load is lognormal, which is what audit/e9_rank_vs_precision.py:33-36
    settled on: a Gaussian draw at a large CV puts the load and the marker
    negative, and 14.1's 40% row was computed on a partly unphysical draw.
    """
    times = np.asarray(times, dtype=float)
    inc = np.asarray(delta(times), dtype=float)
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        Ls = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        Ls = np.full(times.size, LB)
    B = KE * KM + TH + inc - Ls
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    return m, inc, Ls


def columns(m, inc, kappa=False):
    """
    The design columns in the paper's coordinates,
    (ln theta, ln K_m, b, ln k_e) and optionally ln kappa.

    The fifth column is the DOSE SCALE, `Delta -> kappa Delta`, whose
    derivative is `Delta g(m)`. It is design/run.py:943's column exactly. Left
    out, as g1 to g4 left it out, the design is being told the dose scale is
    known to infinite precision, which no assay establishes.
    """
    g = m / (KM + m)
    cols = [TH * g,
            -(TH + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m),
            KE * m]
    if kappa:
        cols.append(inc * g)
    return np.stack(cols, axis=1)


def sd_ln_theta(m, inc, ke_prior=0.8, load_prior=True, kappa_prior=None,
                sigma=SIG):
    """
    sd(ln theta_hat) from the information filter.

    ke_prior=None holds k_e exactly known and DROPS its column, which is
    _filter's own convention and e9's. Giving it a flat prior instead would
    make "k_e known" read as "k_e free" and invert the comparison.

    load_prior=False drops the prior on the load bias `b`, NOT its column. At a
    constant recorded load that is the exactly singular case of P5's fibre, so
    it returns inf, and that is the answer rather than a failure: with the
    prior on, every load-CV comparison is being carried by the prior.

    kappa_prior=None leaves the dose scale out of the parameter vector, which
    asserts it is known exactly. A number puts it in with that log-prior.
    """
    Q = columns(m, inc, kappa=kappa_prior is not None)
    pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2 if load_prior else np.inf,
           ke_prior ** 2 if ke_prior is not None else None]
    if kappa_prior is not None:
        pri.append(kappa_prior ** 2)
    if ke_prior is None:                      # drop the column, do not free it
        keep = [i for i in range(Q.shape[1]) if i != 3]
        Q, pri = Q[:, keep], [pri[i] for i in keep]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def sd_over_seeds(times, load_cv, ke_prior=0.8, seeds=range(8),
                  load_prior=True, kappa_prior=None):
    """Median and range of sd(ln theta_hat) over load realisations.

    At load_cv = 0 the record is deterministic and one seed is the answer;
    above it the spread is large enough that a single seed misleads, which is
    the failure corrected in PHASE4_PI.md C5.
    """
    kw = dict(ke_prior=ke_prior, load_prior=load_prior,
              kappa_prior=kappa_prior)
    if not load_cv:
        v = sd_ln_theta(*record(times, 0.0)[:2], **kw)
        return v, v, v
    vals = np.array([sd_ln_theta(*record(times, load_cv, seed=s)[:2], **kw)
                     for s in seeds], float)
    fin = vals[np.isfinite(vals)]
    if fin.size == 0:
        return np.inf, np.inf, np.inf
    return float(fin.min()), float(np.median(fin)), float(fin.max())


# --------------------------------------------------------------------------
def cycle_phases(times, period=T, tol=1e-6):
    """The distinct dose-cycle phases a schedule visits."""
    ph = np.sort(np.asarray(times, float) % period)
    keep = [ph[0]]
    for x in ph[1:]:
        if x - keep[-1] > tol:
            keep.append(x)
    if len(keep) > 1 and (period - keep[-1]) + keep[0] < tol:
        keep.pop()
    return np.array(keep)


def delta_range(times):
    """The range of Delta/theta a schedule covers, which f4 found is what binds."""
    inc = np.asarray(delta(times), dtype=float)
    return float((inc.max() - inc.min()) / TH)


def settling_days(tol=0.01):
    """
    The post-dose window inside which a draw is not a steady state, at a
    stated relative tolerance.

    The window depends on the tolerance and the repository's own default is
    5%, which is not tight enough to shade as "off the locus": at that
    tolerance the worst departure AFTER settling is 4.97%, i.e. right at the
    bound. The figures use 1%, and both are reported.

    This is ONE of three quantities that get called a settling time, and it is
    the smallest. It is the departure from the MOVING quasi-steady locus,
    which is the question the exclusion rule is for. The paper's 4.2 h and
    P1f's 3.71 h close a full displacement to a FIXED target, at two different
    operating points. g6_settling.py computes all three and shows which choice
    accounts for which gap; every Phase 3 script calls this function, so they
    cannot disagree with one another.
    """
    import io
    import contextlib
    import design.run as D
    with contextlib.redirect_stdout(io.StringIO()):
        settle, worst = D.quasi_steady_check(tol=tol)
    return float(settle), float(worst)
