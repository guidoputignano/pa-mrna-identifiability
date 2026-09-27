"""
The constant-theta assumption, priced against record length.

THE TENSION. eq:input fixes theta_base over the whole record, and
Section sec:local warrants it on the kinetics: most pathogenic variants reduce
the QUANTITY of functional enzyme rather than the turnover number, so theta is
genotype-determined and constant. But Tier 1 needs 118 monthly draws, which is
9.8 years, and "constant over a decade" is a materially stronger assumption
than "constant over the afebrile range", which is what the paper actually
warrants. Growth, disease progression and any change in endogenous expression
all act on the decade and none acts on the fever.

THE ARITHMETIC. Let capacity drift log-linearly, ln theta(t) = ln theta_ref
+ r (t - t_ref), and let the estimator fit a CONSTANT theta, which is what
eq:input licenses it to do. Then

    variance falls with the number of visits,
    bias grows with the span those visits cover,

and a monthly schedule cannot buy one without the other. So there is an
interior optimum in record length and a floor that no schedule beats, and the
question is whether that floor sits above or below the 0.20 gate.

THE BIAS, exactly, inside the linearisation the design is already stated in.
Data are generated with a parameter vector that varies by visit, beta_i =
beta_ref + delta_i e_theta with delta_i = r (t_i - t_ref). The estimator that
fits a constant beta returns, in expectation,

    E[beta_hat] - beta_ref = F^{-1} (1/sigma^2) SUM_i q_i q_{i,theta} delta_i,

F being the same information matrix every other script here builds. The first
component is the bias in ln theta_hat. PART A checks that against simulation
rather than asserting it.

THE REFERENCE TIME IS A CHOICE AND IT IS REPORTED BOTH WAYS. A dosing
decision needs capacity NOW, so t_ref is the last visit and the bias is first
order in r. Taking t_ref at the record's midpoint makes the first-order term
vanish, and that is reported alongside to show the finding is about which
theta is wanted and not about a badly chosen origin.

FOOTING. theta = 3000 with the recorded load held at L = 322.727 umol/L/d,
K_m = 100, k_e = 5/d, D = 6000 on a 14 d train with a 3.5 d half-life, monthly
draws, load recorded to 5% with a day-to-day intake CV of 10%, k_e log-prior
0.8, load-bias prior 0.05 L, gate 0.20. Medians over 24 intake records with
the seed range reported. The Tier 1 arm additionally carries the increment's
scale at the published 14.02%.

Nothing here records an expected result in advance.

Run:  python audit/u7_drift.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, load)

KM, KE, THETA = 100.0, 5.0, 3000.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
YEAR = 365.25
KE_PRIOR = 0.8
SB = LOAD_FRAC * LB
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
PUBLISHED_KAPPA_CV = 0.1402
SEEDS = range(24)


def root(theta, inc, Ltot, Km=KM, ke=KE):
    theta = np.asarray(theta, dtype=float)
    inc = np.asarray(inc, dtype=float)
    Ltot = np.asarray(Ltot, dtype=float)
    q = ke * Km + theta + inc - Ltot
    return (-q + np.sqrt(q * q + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


def schedule(n, step=MONTH):
    t = step * np.arange(n, dtype=float)
    return t, PEAK * np.exp(-LAMBDA * np.mod(t - T_DOSE, DOSE_INTERVAL))


def sigma_at(theta=THETA, assay_cv=ASSAY_CV):
    mbar = float(root(theta, 0.0, LB))
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(LOAD_FRAC * LB, slope * assay_cv * mbar))


def design(n, r, t_ref, seed, step=MONTH, kappa_cv=None, fit_drift=False):
    """
    Columns, priors and the per-visit log-capacity offsets.

    `fit_drift` appends a column for the drift rate itself, d/dr, which is
    the theta column times (t_i - t_ref). Its prior is flat: the point of
    that arm is to ask what the record can do with no outside knowledge of
    the drift at all.
    """
    t, inc = schedule(n, step)
    rng = np.random.default_rng(seed)
    Ls = LB * (1.0 + LOAD_CV * rng.standard_normal(t.size))
    tref = {"end": t[-1], "mid": 0.5 * (t[0] + t[-1]),
            "start": t[0]}[t_ref]
    dl = r * (t - tref) / YEAR                  # r is per YEAR
    theta_t = THETA * np.exp(dl)
    m = root(theta_t, inc, Ls)
    g = m / (KM + m)
    cols = [theta_t * g,
            -(theta_t + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m),
            KE * m]
    pri = [np.inf, np.inf, SB ** 2, KE_PRIOR ** 2]
    if kappa_cv is not None:
        cols.append(inc * g)
        pri.append(kappa_cv ** 2)
    if fit_drift:
        cols.append(theta_t * g * (t - tref) / YEAR)
        pri.append(np.inf)
    return np.stack(cols, axis=1), np.array(pri, dtype=object), dl, t


def exact_bias(n, r, t_ref="end", seed=0, step=MONTH, kappa_cv=None,
               sigma=None):
    """
    The constant-theta estimator's offset on NOISELESS drifting data,
    computed by solving the nonlinear weighted least squares exactly.

    WHY NOT THE FIRST-ORDER PROJECTION. The linearised bias of PART A is
    exact as a projection but it is a projection of a perturbation assumed
    small, and r times the record's span is not small here: at 0.05 a year
    over sixty years it is 3.3, a factor of 27 in capacity. Used there it
    returns a bias that is not monotone in the record length and changes
    sign, which is the linearisation failing and not the estimator doing
    something interesting. This function has no such range restriction.
    """
    sigma = sigma_at() if sigma is None else sigma
    t, inc = schedule(n, step)
    rng = np.random.default_rng(seed)
    Ls = LB * (1.0 + LOAD_CV * rng.standard_normal(t.size))
    tref = {"end": t[-1], "mid": 0.5 * (t[0] + t[-1]), "start": t[0]}[t_ref]
    theta_t = THETA * np.exp(r * (t - tref) / YEAR)
    m_dat = root(theta_t, inc, Ls)

    def resid(z):
        """
        The flux residual the design actually fits, evaluated AT THE
        OBSERVED LEVEL.

        Writing it instead as a difference between the model's own root and
        the data's own root is the error this function made first, and it is
        silent: the model flux at the model's root equals the load by
        definition, so the residual loses every dependence on theta and the
        fit returns its start point. The exact bias then came out as zero in
        some cells and as nonsense in others.
        """
        lt, lK, bs, lk = z[0], z[1], z[2], z[3]
        kap = np.exp(z[4]) if kappa_cv is not None else 1.0
        theta, Km, ke = np.exp(lt), np.exp(lK), np.exp(lk)
        if Km <= 0 or ke <= 0:
            return np.full(len(m_dat) + 3, 1e6)
        fl = (theta + kap * inc) * m_dat / (Km + m_dat) + ke * m_dat
        out = [(fl - (Ls + bs * SB)) / sigma, [bs],
               [(lk - np.log(KE)) / KE_PRIOR]]
        out.append([z[4] / kappa_cv] if kappa_cv is not None else [0.0])
        return np.concatenate([np.atleast_1d(x) for x in out])

    z0 = np.array([np.log(THETA), np.log(KM), 0.0, np.log(KE), 0.0])
    best, cost = None, np.inf
    for j in (0.0, 0.3, -0.3):
        try:
            rr = least_squares(resid, z0 + j * np.array([1.0, -.5, .3, .2, 0.]),
                               method="trf", xtol=1e-14, ftol=1e-14,
                               gtol=1e-14, max_nfev=20000)
        except Exception:
            continue
        if rr.cost < cost:
            best, cost = rr, rr.cost
    if best is None:
        return np.nan
    return float(best.x[0] - np.log(THETA))


def sd_bias(n, r, t_ref="end", seed=0, step=MONTH, kappa_cv=None,
            fit_drift=False, sigma=None, exact=True):
    """(sd, bias, rmse) of ln theta_hat for a record of n draws."""
    sigma = sigma_at() if sigma is None else sigma
    Q, pri, dl, t = design(n, r, t_ref, seed, step, kappa_cv, fit_drift)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / float(x) for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf, np.inf, np.inf
    Fi = np.linalg.inv(F)
    sd = float(np.sqrt(Fi[0, 0]))
    if fit_drift:
        return sd, 0.0, sd       # the drift is in the model
    if exact:
        bias = exact_bias(n, r, t_ref, seed, step, kappa_cv, sigma)
    else:
        rhs = (Q.T @ (Q[:, 0] * dl)) / sigma ** 2
        bias = float((Fi @ rhs)[0])
    if not np.isfinite(bias):
        return sd, np.inf, np.inf
    return sd, bias, float(np.hypot(sd, bias))


def over_seeds(fn, *a, **kw):
    v = np.array([fn(*a, seed=s, **kw) for s in SEEDS], dtype=float)
    fin = v[np.all(np.isfinite(v), axis=1)] if v.ndim > 1 else v
    return fin


def med_rmse(n, r, **kw):
    out = np.array([sd_bias(n, r, seed=s, **kw) for s in SEEDS])
    out = out[np.isfinite(out).all(axis=1)]
    if not out.size:
        return np.inf, np.inf, np.inf, np.inf, np.inf
    return (float(np.median(out[:, 0])), float(np.median(out[:, 1])),
            float(np.median(out[:, 2])), float(out[:, 2].min()),
            float(out[:, 2].max()))


# ======================================================================
def part_a():
    print("=" * 78)
    print("PART A.  the bias formula, against simulation")
    print("=" * 78)
    print("\n  A record is generated with a drifting capacity and fitted with")
    print("  a constant one, by weighted least squares in the linearised")
    print("  model, and the mean of ln theta_hat is compared with the")
    print("  formula. No noise is added: the quantity under test is the")
    print("  systematic part, so noise would only add scatter to a check of")
    print("  a mean.\n")
    sigma = sigma_at()
    print(f"    {'n':>5} {'r /yr':>7} {'span yr':>8} {'r*span':>8} "
          f"{'linearised':>11} {'exact':>10} {'gap':>9}")
    print("    " + "-" * 62)
    for n in (13, 25, 60, 118, 240):
        for r in (0.01, 0.05, 0.10):
            t, _ = schedule(n)
            span = t[-1] / YEAR
            Q, pri, dl, tt = design(n, r, "end", seed=0)
            F = np.diag([0.0 if np.isinf(x) else 1.0 / float(x) for x in pri]) \
                + Q.T @ Q / sigma ** 2
            rhs = (Q.T @ (Q[:, 0] * dl)) / sigma ** 2
            lin = float((np.linalg.inv(F) @ rhs)[0])
            ex = exact_bias(n, r, "end", seed=0)
            print(f"    {n:5d} {r:7.2f} {span:8.2f} {r * span:8.3f} "
                  f"{lin:11.5f} {ex:10.5f} {abs(lin - ex):9.4f}")
    print("\n  The two agree while r times the span is small and part")
    print("  company once it is not, which is the linearisation's own")
    print("  validity condition showing up as a number. EVERY BIAS REPORTED")
    print("  BELOW IS THE EXACT COLUMN. The linearised one is kept here")
    print("  because it is the natural first thing to compute and it is")
    print("  wrong in exactly the regime the Tier 1 record sits in.")


def part_b():
    print()
    print("=" * 78)
    print("PART B.  sd, bias and their combination against record length")
    print("=" * 78)
    print("\n  monthly draws, capacity wanted at the LAST visit, medians over")
    print("  24 intake records. rmse is hypot(sd, bias).\n")
    for r in (0.0, 0.02, 0.05, 0.10, 0.20):
        print(f"  drift r = {r:.2f} per year "
              f"(a factor of {np.exp(r):.3f} a year)")
        print(f"    {'n':>5} {'years':>7} {'sd':>9} {'bias':>10} {'rmse':>9} "
              f"{'seed range of rmse':>26}")
        print("    " + "-" * 70)
        best = (np.inf, None)
        for n in (13, 25, 50, 75, 100, 118, 160, 240):
            sd, bi, rm, lo, hi = med_rmse(n, r)
            flag = ""
            if rm < best[0]:
                best = (rm, n)
            print(f"    {n:5d} {n * MONTH / YEAR:7.2f} {sd:9.4f} {bi:10.4f} "
                  f"{rm:9.4f} {f'[{lo:.4f}, {hi:.4f}]':>26}{flag}")
        print(f"    best of the grid: n = {best[1]}, rmse = {best[0]:.4f}"
              f"{'  (meets the gate)' if best[0] <= GATE else '  (misses the gate)'}")
        print()


def optimum(r, kappa_cv=None, nmax=240):
    """The record length minimising rmse, and the floor there."""
    best = (np.inf, None)
    for n in [6, 8, 10, 13, 16, 20, 25, 32, 40, 50, 63, 79, 100, 126, 158,
              200, 240]:
        if n > nmax:
            break
        rm = med_rmse(n, r, kappa_cv=kappa_cv)[2]
        if rm < best[0]:
            best = (rm, n)
    n0 = best[1]
    for n in range(max(6, int(n0 * 0.75)), min(nmax, int(n0 * 1.35)) + 1,
                   max(1, int(n0 * 0.06))):
        rm = med_rmse(n, r, kappa_cv=kappa_cv)[2]
        if rm < best[0]:
            best = (rm, n)
    return best


def part_c():
    print()
    print("=" * 78)
    print("PART C.  the floor, swept over the drift rate")
    print("=" * 78)
    print("\n  The optimum is interior because sd falls with n and bias grows")
    print("  with the span. The floor is the rmse there: no schedule of")
    print("  monthly draws beats it at that drift rate. Record lengths are")
    print("  capped at 240 monthly draws, which is 20 years and already")
    print("  longer than any study that could be run.\n")
    print(f"    {'r /yr':>8} {'halving/doubling':>17} {'n*':>6} {'years*':>8} "
          f"{'floor':>9} {'vs gate':>10}")
    print("    " + "-" * 62)
    rows = []
    for r in (0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20,
              0.35, 0.50, 0.69):
        fl, n = optimum(r)
        rows.append((r, n, fl))
        print(f"    {r:8.3f} {np.log(2.0) / r:17.1f} {n:6d} "
              f"{n * MONTH / YEAR:8.2f} {fl:9.4f} "
              f"{'PASSES' if fl <= GATE else 'FAILS':>10}")
    #  where does the floor cross the gate
    rs = np.array([x[0] for x in rows])
    fs = np.array([x[2] for x in rows])
    cross = None
    for i in range(len(rs) - 1):
        if fs[i] <= GATE < fs[i + 1]:
            cross = float(np.exp(np.log(rs[i]) + (GATE - fs[i])
                                 * (np.log(rs[i + 1]) - np.log(rs[i]))
                                 / (fs[i + 1] - fs[i])))
    print(f"\n  the floor crosses the 0.20 gate at a drift of about "
          f"{cross:.3f} per year" if cross else
          "\n  the floor does not cross the gate on this grid")
    if cross:
        print(f"  which is a halving time of {np.log(2.0) / cross:.1f} years,")
        print(f"  or {100 * (np.exp(cross) - 1):.1f}% a year.")
    return rows


def part_d():
    print()
    print("=" * 78)
    print("PART D.  fitting the drift instead of assuming it away")
    print("=" * 78)
    print("\n  A sixth column, d/dr, with a FLAT prior: the record is asked")
    print("  to find the drift itself. The bias then vanishes by")
    print("  construction and the whole cost appears in the variance.\n")
    print(f"    {'r /yr':>7} {'n':>5} {'years':>7} {'constant: rmse':>15} "
          f"{'drift fitted: sd':>18} {'ratio':>8}")
    print("    " + "-" * 64)
    for r in (0.02, 0.05, 0.10, 0.20):
        for n in (25, 118, 200):
            c = med_rmse(n, r)[2]
            d = med_rmse(n, r, fit_drift=True)[0]
            print(f"    {r:7.2f} {n:5d} {n * MONTH / YEAR:7.2f} {c:15.4f} "
                  f"{d:18.4f} {d / c:8.2f}")
    print("\n  Fitting the drift removes the bias and buys nothing if the")
    print("  variance it costs is larger, which is what the ratio column")
    print("  says at each footing.")

    print("\n  THE HONEST PRICE OF THE COLUMN is what it costs when there is")
    print("  no drift to find, because that is the premium paid for the")
    print("  insurance. At r > 0 the drift is itself excitation in theta and")
    print("  flatters the arm that models it.\n")
    print(f"    {'n':>5} {'years':>7} {'constant, r=0':>14} "
          f"{'drift column, r=0':>18} {'premium':>9}")
    print("    " + "-" * 56)
    for n in (13, 25, 50, 118, 240):
        a = med_rmse(n, 0.0)[0]
        b = med_rmse(n, 0.0, fit_drift=True)[0]
        print(f"    {n:5d} {n * MONTH / YEAR:7.2f} {a:14.4f} {b:18.4f} "
              f"{b / a:9.3f}")


def part_f():
    print()
    print("=" * 78)
    print("PART F.  the floor on the arm that actually needs the long record")
    print("=" * 78)
    print("\n  Everything above is at nominal capacity with the increment")
    print("  known exactly, where the gate is met at 13 to 25 draws and the")
    print("  drift has no span to act over. The arm that needs a decade is")
    print("  Tier 1, whose increment is the published 14.02%, and that is")
    print("  where the two effects meet.")
    print()
    print("  The increment's prior does not shrink with n, so the variance")
    print("  has a floor of its own; the bias has no floor. Below, the best")
    print("  record length at each drift rate and the rmse there.\n")
    print(f"    {'r /yr':>8} {'n*':>6} {'years*':>8} {'sd at n*':>10} "
          f"{'bias at n*':>11} {'floor':>9} {'vs gate':>9}")
    print("    " + "-" * 66)
    rows = []
    for r in (0.0, 0.002, 0.005, 0.0075, 0.010, 0.015, 0.020, 0.030, 0.050):
        fl, n = optimum(r, kappa_cv=PUBLISHED_KAPPA_CV)
        sd, bi, rm, lo, hi = med_rmse(n, r, kappa_cv=PUBLISHED_KAPPA_CV)
        rows.append((r, n, fl))
        print(f"    {r:8.4f} {n:6d} {n * MONTH / YEAR:8.2f} {sd:10.4f} "
              f"{bi:11.4f} {fl:9.4f} "
              f"{'PASSES' if fl <= GATE else 'FAILS':>9}")
    rs = [x[0] for x in rows]
    fs = [x[2] for x in rows]
    cross = None
    for i in range(len(rs) - 1):
        if fs[i] <= GATE < fs[i + 1]:
            cross = rs[i] + (GATE - fs[i]) * (rs[i + 1] - rs[i]) \
                / (fs[i + 1] - fs[i])
    if cross:
        print(f"\n  the floor crosses the gate at a drift of about "
              f"{cross:.4f} per year,")
        print(f"  which is {100 * (np.exp(cross) - 1):.2f}% a year, a halving "
              f"time of {np.log(2.0) / cross:.0f} years.")
    else:
        print("\n  the floor does not cross the gate on this grid")
    print()
    print("  and the same arm with the drift FITTED rather than assumed away:")
    print()
    print(f"    {'r /yr':>8} {'n':>6} {'constant: floor':>16} "
          f"{'drift fitted: sd':>18} {'vs gate':>9}")
    print("    " + "-" * 56)
    for r in (0.005, 0.010, 0.020, 0.050):
        fl, n = optimum(r, kappa_cv=PUBLISHED_KAPPA_CV)
        d = med_rmse(240, r, kappa_cv=PUBLISHED_KAPPA_CV, fit_drift=True)[0]
        print(f"    {r:8.4f} {240:6d} {fl:16.4f} {d:18.4f} "
              f"{'PASSES' if d <= GATE else 'FAILS':>9}")
    return rows


def part_e():
    print()
    print("=" * 78)
    print("PART E.  what this does to Tier 1 specifically")
    print("=" * 78)
    print("\n  Tier 1's 118 monthly draws span 9.82 years. Carrying the")
    print("  increment at the published 14.02% throughout, which is the")
    print("  footing that made 118 the requirement in the first place.\n")
    print(f"    {'r /yr':>8} {'sd at n=118':>13} {'bias at n=118':>15} "
          f"{'rmse':>9} {'verdict':>10}")
    print("    " + "-" * 60)
    for r in (0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10):
        sd, bi, rm, lo, hi = med_rmse(118, r, kappa_cv=PUBLISHED_KAPPA_CV)
        print(f"    {r:8.3f} {sd:13.4f} {bi:15.4f} {rm:9.4f} "
              f"{'PASSES' if rm <= GATE else 'FAILS':>10}")
    print("\n  and the same with the reference time at the record's midpoint,")
    print("  which is the most favourable choice of estimand:\n")
    print(f"    {'r /yr':>8} {'bias (end)':>12} {'bias (mid)':>12} "
          f"{'rmse (mid)':>12} {'verdict':>10}")
    print("    " + "-" * 58)
    for r in (0.01, 0.02, 0.05, 0.10):
        _, be, _, _, _ = med_rmse(118, r, kappa_cv=PUBLISHED_KAPPA_CV)
        sd, bm, rm, _, _ = med_rmse(118, r, t_ref="mid",
                                    kappa_cv=PUBLISHED_KAPPA_CV)
        print(f"    {r:8.3f} {be:12.4f} {bm:12.4f} {rm:12.4f} "
              f"{'PASSES' if rm <= GATE else 'FAILS':>10}")
    print("\n  The midpoint column is what the assumption costs when the")
    print("  quantity wanted is the record's average capacity. The 'end'")
    print("  column is what it costs when the quantity wanted is the")
    print("  patient's capacity now, which is what a dosing interval needs.")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
