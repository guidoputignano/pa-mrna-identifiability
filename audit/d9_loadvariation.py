"""
The decisive experiment for the specification.

Proposition rank says a record at a CONSTANT recorded load is exactly rank
deficient, and that any real variation in intake closes it. Both estimators in
this note are then reporting artefacts on the constant-load design:

  - the information filter fills the null direction from the load prior and
    returns 0.1475, which the data did not supply;
  - the prior-free Monte Carlo starts half a decade from the truth, cannot
    return along a direction whose gradient is third order, and reports the
    spread it started with.

So the question is not which estimator to believe. It is whether putting the
ordinary day-to-day variation of eating back into the record fixes both.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import least_squares
import design.run as D, model_A.run as LOC

P = D.P
L0 = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def record(t, load_cv, seed=0):
    rng = np.random.default_rng(seed)
    dd = np.asarray(D.delta(t), dtype=float)
    Ls = (L0 * (1.0 + load_cv * rng.standard_normal(t.size)) if load_cv
          else np.full(t.size, L0))
    B = P.ke * P.Km + P.theta + dd - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
    return lv, dd, Ls


def fisher_sd(lv, dd, load_frac=D.LOAD_FRAC):
    F = np.linalg.inv(np.diag([0.8 ** 2, 0.8 ** 2, (load_frac * L0) ** 2, 0.8 ** 2]))
    for m, d in zip(lv, dd):
        q = np.array([P.theta * m / (P.Km + m),
                      -(P.theta + d) * P.Km * m / (P.Km + m) ** 2, 1.0, P.ke * m])
        F = F + np.outer(q, q) / SIG ** 2
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def mc(lv, dd, Ls, reps=300, seed=5, s_peak=None):
    rng = np.random.default_rng(seed)
    b_loc = P.theta * P.Km / (P.Km + D.M_BAR) ** 2 + P.ke
    out = []
    for _ in range(reps):
        s_obs = (Ls * (1.0 + rng.normal(0.0, D.LOAD_FRAC, lv.size))
                 + rng.normal(0.0, b_loc * D.ASSAY_CV * D.M_BAR, lv.size))

        def res(q):
            t_, K, k = np.exp(np.clip(q[:3], -25.0, 25.0))
            return (t_ + dd) * lv / (K + lv) + k * lv + q[3] - s_obs
        start = np.concatenate([np.log([P.theta, P.Km, P.ke])
                                + rng.normal(0.0, 0.5, 3), [0.0]])
        try:
            r = least_squares(res, start, method="lm", max_nfev=3000)
        except Exception:
            continue
        out.append(np.exp(np.clip(r.x[:3], -25.0, 25.0)))
    ens = np.array(out)
    sd = float(np.std(np.log(ens[:, 0]), ddof=1))
    mu_t = D.margin_and_gradient(P.theta, P.Km, P.ke, D.M_TOL, s_peak)[0]
    mh = np.array([D.margin_and_gradient(a, b, c, D.M_TOL, s_peak)[0] for a, b, c in ens])
    return sd, float(np.mean(np.sign(mh) != np.sign(mu_t)))


def main():
    print(__doc__)
    base = D.clearance_A(D.M_BAR, P)
    t = np.arange(0.0, 730.0, 30.4)
    print("=" * 86)
    print("MONTHLY SAMPLING ACROSS DOSE CYCLES, n = 24, two years")
    print("=" * 86)
    print(f"{'day-to-day load CV':>20} {'Fisher sd':>11} {'MC sd':>9} "
          f"{'MC sign error at mu = -0.3 theta':>34}")
    for cv in (0.0, 0.02, 0.05, 0.10, 0.20, 0.35):
        lv, dd, Ls = record(t, cv, seed=3)
        f = fisher_sd(lv, dd)
        sd, p = mc(lv, dd, Ls, s_peak=3.6 * base)
        print(f"{cv:19.0%} {f:11.4f} {sd:9.4f} {p:33.1%}")
    print()
    print("  The Fisher barely moves: the prior was already holding that direction.")
    print("  The Monte Carlo is the one that changes, because the data now pin what")
    print("  the prior was standing in for.")
    print()
    print("=" * 86)
    print("THE SAME, WITH THE INCREMENT DROPPED FROM THE FIT (the pre-2C footing)")
    print("=" * 86)
    for cv in (0.0, 0.10, 0.35):
        lv, dd, Ls = record(t, cv, seed=3)
        sd, p = mc(lv, np.zeros_like(dd), Ls, s_peak=3.6 * base)
        print(f"{cv:19.0%} {'---':>11} {sd:9.4f} {p:33.1%}")
    print()
    print("  Load variation alone does not rescue a fit that ignores the increment.")


if __name__ == "__main__":
    main()
