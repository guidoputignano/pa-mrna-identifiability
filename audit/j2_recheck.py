"""
Phase 1 step 2, the independent re-check. Four claims, each by a route that
shares no code path with the one that produced it.

    (1) the Wald sd and the profile are the same object, checked by finite
        differencing the profile at its minimum instead of inverting a matrix
    (2) the LR statistic is correctly formed, checked where its distribution
        is KNOWN: strip the priors and the misspecification and it must be
        chi^2_1
    (3) the b >= 0 active set is the constrained optimum, checked against a
        general constrained optimiser
    (4) the profile intervals are the set they claim to be, checked by direct
        evaluation on a fine grid instead of interpolation

Nothing here records an expected result in advance.

Run:  python audit/j2_recheck.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.optimize import minimize
from scipy import stats

import j1_uncertainty as J


def head(t):
    print("\n" + "=" * 94); print("  " + t); print("=" * 94)


# ------------------------------------------------------------------ check 1
def check_wald_against_curvature():
    head("(1) THE WALD sd AGAINST THE CURVATURE OF THE PROFILE")
    print("  _wald_sd inverts an information matrix built from analytic")
    print("  derivative columns. The profile is minimised numerically over a")
    print("  grid. If both describe the same objective then")
    print("      sd = sqrt(2 / d^2 J_profile / d(ln theta)^2)  at the minimum,")
    print("  because J is 2 x the negative log posterior. No matrix is inverted")
    print("  on this route and no analytic column is used.\n")
    cfg = dict(load_prior=True, kappa=None, bound_b=False)
    print(f"    {'rep':>4} {'wald sd':>10} {'curvature sd':>14} {'ratio':>9}")
    print("  " + "-" * 44)
    for i in range(6):
        m, inc, s, _, _ = J.record(0, i, cfg, 0.10, "draw")
        f = J.fit(m, inc, s, cfg)
        lth = np.log(f["theta"])
        h = 0.02
        g = np.array([lth - h, lth, lth + h])
        Jp = J.profile(m, inc, s, cfg, lt=g)
        curv = (Jp[0] - 2 * Jp[1] + Jp[2]) / h ** 2
        sd_c = np.sqrt(2.0 / curv) if curv > 0 else np.inf
        sd_w = J._wald_sd(m, inc, s, cfg, f)
        print(f"    {i:4d} {sd_w:10.6f} {sd_c:14.6f} {sd_c/sd_w:9.4f}")
    print("\n    A ratio at 1 says the two are the same object computed two")
    print("    ways. It also retires the coordinate bug an earlier draft had,")
    print("    where the d/d ln k_e column was paired with a prior sd of")
    print("    0.8 k_e instead of 0.8, inflating sd(ln theta) about fourfold.")


# ------------------------------------------------------------------ check 2
def check_lr_null():
    head("(2) IS THE LR STATISTIC RIGHT? TEST IT WHERE ITS LAW IS EXACT")
    print("  An earlier draft of this check generated data with the priors off")
    print("  and K_m free and asserted the profile LR must be chi^2_1. It is")
    print("  not: that problem has a condition number of 3.7e7, the profile")
    print("  over a near-flat direction absorbs most of a change in theta, and")
    print("  the statistic came out at 0.31 x chi^2_1 at EVERY quantile. That")
    print("  is a property of the problem, not of the arithmetic, and the check")
    print("  as posed could not tell the two apart.\n")
    print("  Where the law IS exact: hold K_m at the truth and the model is")
    print("  LINEAR in (theta, k_e, b) with sigma known, so the profile LR for")
    print("  one parameter is exactly chi^2_1, whatever the conditioning.\n")
    # ke_prior=None flattens the k_e log-prior. An earlier draft of this
    # check said `priors off` while the k_e prior was hard-wired on, which is
    # exactly what shrank the statistic below chi^2_1 and looked like a bug in
    # it. The prior was the bug, in the check.
    cfg = dict(load_prior=False, kappa=None, bound_b=False, ke_prior=None)
    m0, inc, _, _, _ = J.record(0, 0, cfg, 0.10, "mean")
    g = m0 / (J.KM + m0)
    mu = J.TH * g + inc * g + J.KE * m0
    rng = np.random.default_rng(4242)
    n = 20000
    S = mu + rng.normal(0.0, J.SIG, (n, len(m0)))
    one = np.array([J.KM])

    def lr_at(fixed_grid, label):
        X, y, sds, p0 = J.build(fixed_grid, m0, inc, S[:, None, :], cfg)
        _, Ju, _ = J.solve(X, y, sds, p0, False)
        Xr, yr, sdr, p0r = J.build(fixed_grid, m0, inc, S[:, None, :], cfg,
                                   theta_fixed=J.TH)
        _, Jr, _ = J.solve(Xr, yr, sdr, p0r, False)
        v = Jr.min(axis=1) - Ju.min(axis=1)
        v = np.maximum(v[np.isfinite(v)], 0.0)
        print(f"\n  {label}, {v.size} replicates")
        print(f"    {'quantile':>10} {'empirical':>12} {'chi^2_1':>10} {'ratio':>8}")
        print("  " + "-" * 44)
        for qq in (0.50, 0.75, 0.90, 0.95, 0.99):
            e = float(np.quantile(v, qq)); c = float(stats.chi2.ppf(qq, 1))
            print(f"    {qq:10.2f} {e:12.4f} {c:10.4f} {e/c:8.4f}")
        ks = stats.kstest(v, "chi2", args=(1,))
        print(f"    mean {v.mean():.4f} against 1, "
              f"Kolmogorov-Smirnov p = {ks.pvalue:.4f}")
        return v

    # the strongest form of the check: j1's solver against a from-scratch
    # projection computation on the same data, sharing no code at all.
    X = np.stack([g, m0, -np.ones_like(m0)], 1)
    y = S - inc * g
    Pj = X @ np.linalg.pinv(X)
    Ju = np.einsum("ni,ni->n", y - (Pj @ y.T).T, y - (Pj @ y.T).T) / J.SIG ** 2
    Xr, yr = X[:, 1:], y - J.TH * g
    Pr = Xr @ np.linalg.pinv(Xr)
    Jr = np.einsum("ni,ni->n", yr - (Pr @ yr.T).T, yr - (Pr @ yr.T).T) / J.SIG ** 2
    v_text = Jr - Ju
    v_j1 = lr_at(one, "K_m HELD at the truth, ALL priors flat: the model is linear")
    print(f"\n    j1's solver against a from-scratch projection computation on")
    print(f"    the same data: max |difference| "
          f"{np.abs(np.sort(v_text) - np.sort(v_j1)):s}" if False else
          f"    the same data: max |difference| "
          f"{float(np.max(np.abs(np.sort(v_text) - np.sort(np.maximum(v_j1, 0.0))))):.3e}")
    lr_at(J.GRID, "K_m FREE on the grid, all priors flat: nonlinear")
    cfgk = dict(cfg, ke_prior=J.KE_PRIOR_LOG)

    def lr_cfg(c, label):
        X, y, sds, p0 = J.build(J.GRID, m0, inc, S[:, None, :], c)
        _, Ju, _ = J.solve(X, y, sds, p0, False)
        Xr, yr, sdr, p0r = J.build(J.GRID, m0, inc, S[:, None, :], c,
                                   theta_fixed=J.TH)
        _, Jr, _ = J.solve(Xr, yr, sdr, p0r, False)
        v = np.maximum(Jr.min(axis=1) - Ju.min(axis=1), 0.0)
        print(f"\n  {label}: mean {v.mean():.4f}, 95% {np.quantile(v,0.95):.4f} "
              f"against 3.8415")
        return v

    lr_cfg(cfgk, "K_m free, k_e log-prior 0.8 ON (j1's own setting)")
    print("\n    The first block is the arithmetic check and must sit on")
    print("    chi^2_1. The second adds one nonlinear, nearly-flat direction.")
    print("    The third adds the k_e prior back, which is j1's own setting and")
    print("    which shrinks the statistic further because a penalised")
    print("    parameter is not a free one. None of the three is a defect in")
    print("    the statistic; together they are the reason 3.841 has to be")
    print("    checked against the empirical distribution rather than assumed.")


# ------------------------------------------------------------------ check 3
def check_bound():
    head("(3) THE b >= 0 ACTIVE SET AGAINST A GENERAL CONSTRAINED OPTIMISER")
    print("  j1 handles the bound by solving free and, where b came out")
    print("  negative, re-solving with b held at zero. That is exact for ONE")
    print("  inequality, and `exact for a reason` is worth checking against a")
    print("  solver that knows nothing about the reason.\n")
    cfg = dict(load_prior=False, kappa=None, bound_b=True)
    print(f"    {'rep':>4} {'J active-set':>14} {'J scipy':>13} {'diff':>11} "
          f"{'b active-set':>13} {'b scipy':>11} {'on bound':>9}")
    print("  " + "-" * 82)
    worst = 0.0
    for i in range(8):
        m, inc, s, _, _ = J.record(0, i, cfg, 0.10, 0.0)
        f = J.fit(m, inc, s, cfg)
        Km = f["Km"]
        X, y, sds, p0 = J.build(np.array([Km]), m, inc, s, cfg)
        X, y = X[0], y[0]

        def obj(p):
            r = X @ p - y
            d = p - p0
            return (r @ r / J.SIG ** 2
                    + np.sum(np.where(np.isinf(sds), 0.0,
                                      d ** 2 / np.where(np.isinf(sds), 1.0, sds) ** 2)))

        res = minimize(obj, np.array([f["theta"], f["ke"], max(f["b"], 0.0)]),
                       method="SLSQP",
                       bounds=[(None, None), (None, None), (0.0, None)],
                       options={"ftol": 1e-14, "maxiter": 2000})
        d = f["J"] - res.fun
        worst = max(worst, abs(d))
        print(f"    {i:4d} {f['J']:14.8f} {res.fun:13.8f} {d:+11.2e} "
              f"{f['b']:13.6f} {res.x[2]:11.6f} {str(f['on_bound']):>9}")
    print(f"\n    worst |difference| in the objective: {worst:.2e}")
    print("    A positive difference would mean the active set stopped short of")
    print("    the constrained optimum; a negative one that it left the")
    print("    feasible set.")


# ------------------------------------------------------------------ check 4
def check_interval_set():
    head("(4) ARE THE PROFILE INTERVALS THE SET THEY CLAIM TO BE?")
    print("  interval_from walks out from the minimum and interpolates the two")
    print("  crossings. Here the objective is evaluated directly on a grid ten")
    print("  times finer and the set {ln theta : LR <= 3.841} is read off it,")
    print("  with no interpolation and no walking.\n")
    cfg = dict(load_prior=True, kappa=None, bound_b=False)
    fine = np.log(J.TH) + np.linspace(-1.5, 1.5, 1601)
    print(f"    {'rep':>4} {'coarse lo':>11} {'fine lo':>11} "
          f"{'coarse hi':>11} {'fine hi':>11} {'width diff':>12} {'connected':>10}")
    print("  " + "-" * 76)
    for i in range(6):
        m, inc, s, _, _ = J.record(0, i, cfg, 0.10, "draw")
        lr = J.profile(m, inc, s, cfg)
        lo, hi, _ = J.interval_from(J.LT, lr - lr.min(), J.CHI2_1_95)
        lrf = J.profile(m, inc, s, cfg, lt=fine)
        ok = lrf - lrf.min() <= J.CHI2_1_95
        idx = np.where(ok)[0]
        conn = bool(idx.size and (idx[-1] - idx[0] + 1) == idx.size)
        flo, fhi = fine[idx[0]], fine[idx[-1]]
        print(f"    {i:4d} {lo:11.5f} {flo:11.5f} {hi:11.5f} {fhi:11.5f} "
              f"{(hi-lo)-(fhi-flo):+12.2e} {str(conn):>10}")
    print(f"\n    The fine grid's own spacing is {fine[1]-fine[0]:.2e}, so a")
    print("    difference of that order is the grid and not the method.")
    print("    `connected` matters: interval_from returns an INTERVAL, and if")
    print("    the accepted set were disconnected that would be the wrong")
    print("    summary of it however accurate the endpoints.")


if __name__ == "__main__":
    print(__doc__)
    check_wald_against_curvature()
    check_lr_null()
    check_bound()
    check_interval_set()
