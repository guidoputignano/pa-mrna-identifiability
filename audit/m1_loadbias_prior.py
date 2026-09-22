"""
m1 -- how much of the headline precision is the load-bias prior, and what the
September 2026 dietary review does to it.

The design's information matrix carries a prior on the load-bias column of
sd = LOAD_FRAC * L0, with LOAD_FRAC = 5%. Section 8 of the paper states that
when the recorded load is constant the null direction of Proposition rank "is
filled by the 5% prior on the load bias", and that only real day-to-day
variation lets the data fill it instead.

Two things the dietary review supplies bear on that sentence:

  (i) European PA dietary records in the published record are PRESCRIPTIONS
      (Daly 2017, PMID 29021961, verified: a cross-sectional survey of what 186
      patients were "prescribed"), and no PA study reports repeated as-eaten
      daily intake. A prescription record is a CONSTANT recorded load.
 (ii) The best analogue for how far a prescription sits from what is eaten is
      PKU, where actual natural protein ran at 167% of prescribed, a median
      +40% (Nutrients 2023;15:4903). That is the load BIAS, and it is the
      quantity LOAD_FRAC stands for.

So this sweeps LOAD_FRAC under both record regimes. Nothing here is asserted in
advance; the table is whatever it is.
"""
import sys
import numpy as np
sys.path.insert(0, 'src')
import matplotlib
matplotlib.use('Agg')
import design.run as D, model_A.run as LOC

P = D.P
L0 = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
T = np.arange(0.0, 730.0, 30.4)          # the paper's monthly schedule, n = 24


def levels(load_cv, seed):
    """The marker levels a record would contain, as in d9_loadvariation."""
    rng = np.random.default_rng(seed)
    dd = np.asarray(D.delta(T), dtype=float)
    Ls = (L0 * (1.0 + load_cv * rng.standard_normal(T.size)) if load_cv
          else np.full(T.size, L0))
    B = P.ke * P.Km + P.theta + dd - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
    return lv, dd


def columns(lv, dd):
    """Sensitivity columns in (ln theta, ln Km, load bias, ln ke)."""
    return np.array([[P.theta * m / (P.Km + m),
                      -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                      1.0,
                      P.ke * m] for m, d in zip(lv, dd)])


def sd_inverse(Q, load_frac):
    """Path A: build F, invert it, take the (0,0) entry."""
    F = np.linalg.inv(np.diag([0.8 ** 2, 0.8 ** 2, (load_frac * L0) ** 2, 0.8 ** 2]))
    F = F + Q.T @ Q / SIG ** 2
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def sd_schur(Q, load_frac):
    """Path B, independent of A: solve the same normal equations as a weighted
    least-squares problem with the prior written as four extra pseudo-rows, and
    read sd(ln theta) off the QR of the stacked design. No matrix inverse of F,
    no diagonal extraction."""
    pri = np.array([0.8, 0.8, load_frac * L0, 0.8])
    A = np.vstack([Q / SIG, np.diag(1.0 / pri)])
    R = np.linalg.qr(A, mode='r')
    e0 = np.zeros(4); e0[0] = 1.0
    # sd = || R^{-T} e0 ||  since Cov = (R^T R)^{-1}
    v = np.linalg.solve(R.T, e0)
    return float(np.linalg.norm(v))


def main():
    print(__doc__)
    fracs = [0.01, 0.05, 0.10, 0.20, 0.40, 0.80, np.inf]
    for cv, label in ((0.10, "recorded load varies day to day at 10% CV  (what the paper assumes)"),
                      (0.00, "recorded load CONSTANT  (a prescription record)")):
        print("=" * 92)
        print(label)
        print("=" * 92)
        print(f"{'load-bias prior':>16} {'sd(ln theta) median':>21} {'over 24 seeds':>26} {'check B':>10}")
        for lf in fracs:
            vals, chk = [], []
            for s in range(24):
                lv, dd = levels(cv, seed=100 + s)
                Q = columns(lv, dd)
                vals.append(sd_inverse(Q, lf if np.isfinite(lf) else 1e12))
                chk.append(sd_schur(Q, lf if np.isfinite(lf) else 1e12))
            v = np.array(vals); c = np.array(chk)
            tag = "flat" if not np.isfinite(lf) else f"{lf:.0%}"
            print(f"{tag:>16} {np.median(v):21.4f} "
                  f"{'[' + format(v.min(), '.4f') + ', ' + format(v.max(), '.4f') + ']':>26} "
                  f"{np.max(np.abs(c - v) / v):10.2e}")
        print()
    print("check B is the largest relative disagreement between the two code paths.")


if __name__ == "__main__":
    main()


def sigma_sweep():
    """
    The second lever, and the one spec_numbers() holds fixed.

    spec_numbers() computes sigma ONCE at LOAD_FRAC = 5% and reuses it on every
    row, so the table's "load record kept to 1% or to 20%" rows vary the
    day-to-day VARIATION in the record and not the recording ERROR. The error
    enters only through sigma. If a prescription record sits a median 40% from
    what is eaten, that divergence lands in the residual.
    """
    import model_A.run as LOC
    print("=" * 92)
    print("sd(ln theta), monthly draws, as the RECORDING ERROR feeding sigma is varied")
    print("  (record still varies day to day at 10% CV; load-bias prior held at 5%)")
    print("=" * 92)
    print(f"{'recording error':>16} {'sigma':>10} {'sd(ln theta) median':>21} {'over 24 seeds':>26}")
    for lf in (0.01, 0.05, 0.10, 0.20, 0.40):
        sig = LOC.sigma_total(lf, D.ASSAY_CV)
        vals = []
        for s in range(24):
            lv, dd = levels(0.10, seed=100 + s)
            Q = columns(lv, dd)
            F = np.linalg.inv(np.diag([0.8 ** 2, 0.8 ** 2,
                                       (D.LOAD_FRAC * L0) ** 2, 0.8 ** 2]))
            F = F + Q.T @ Q / sig ** 2
            vals.append(float(np.sqrt(np.linalg.inv(F)[0, 0])))
        v = np.array(vals)
        print(f"{lf:16.0%} {sig:10.2f} {np.median(v):21.4f} "
              f"{'[' + format(v.min(), '.4f') + ', ' + format(v.max(), '.4f') + ']':>26}")


if __name__ == "__main__":
    sigma_sweep()
