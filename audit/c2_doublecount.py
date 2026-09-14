"""
C2 -- is carrying Delta in the sensitivity CONSISTENCY, or the same term twice?

THE DERIVATION, stated once and explicitly.

A record is a set of PAIRS (m_i, s_i): m_i the marker measured at time t_i,
s_i the propiogenic load recorded over the same window.  The model asserts a
balance between them,

    s_i = (theta + Delta_i) g(m_i; K_m) + k_e m_i + b + eps_i ,          (*)

with Delta_i = Delta(t_i) the therapeutic increment -- a KNOWN regressor, read
off the dosing record and the product's PK, not a parameter.  The unknowns are
(theta, K_m, k_e, b); the noise sits on the ordinate.

The Fisher column is the derivative of (*) at FIXED m_i, because m_i is the
abscissa -- it is data, not a model output:

    ds/d ln theta = theta g(m_i)
    ds/d ln K_m   = -(theta + Delta_i) K_m m_i / (K_m + m_i)^2
    ds/d ln k_e   = k_e m_i
    ds/d b        = 1

Delta_i appears exactly once, in (*).  quasi_steady_level(t) is used only to
GENERATE the synthetic m_i; it is the data-generating process, not part of the
likelihood.  Delta therefore enters in two DIFFERENT roles -- it moves the
abscissa (through the data) and it scales the K_m coefficient (through the
model) -- and neither substitutes for the other.  Dropping it from the column
does not remove a duplicate; it linearises the likelihood about a model the
data were not generated from.

A derivation can be wrong, so the operational test is below: if Delta were
double-counted the Fisher matrix would OVERSTATE the information and the Monte
Carlo spread would exceed the Fisher prediction.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import least_squares
import design.run as D, model_A.run as LOC

P = D.P
Lb = D.load(0.0)
SIGMA = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def fisher_sd(lv, dd, free_ke, use_delta, priors=True, ke_prior=0.8):
    pri = [0.8 ** 2, 0.8 ** 2, (D.LOAD_FRAC * Lb) ** 2]
    if free_ke:
        pri.append(ke_prior ** 2)
    Finv = (np.linalg.inv(np.diag(pri)) if priors
            else np.zeros((len(pri), len(pri))))
    for m, dl in zip(lv, dd):
        eff = P.theta + (dl if use_delta else 0.0)
        q = [P.theta * m / (P.Km + m), -eff * P.Km * m / (P.Km + m) ** 2, 1.0]
        if free_ke:
            q.append(P.ke * m)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / SIGMA ** 2
    C = np.linalg.inv(Finv)
    return float(np.sqrt(C[0, 0])), float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]))


def mc(lv, dd, free_ke, use_delta, reps=400, seed=7):
    """
    Monte Carlo of the SAME estimator.  Data are generated from the true model
    with the true Delta; the fit uses Delta or not, per `use_delta`.
    """
    rng = np.random.default_rng(seed)
    s_true = (P.theta + dd) * lv / (P.Km + lv) + P.ke * lv
    x0 = np.array([np.log(P.theta), np.log(P.Km), 0.0]
                  + ([np.log(P.ke)] if free_ke else []))

    def resid(x, s):
        th, km = np.exp(x[0]), np.exp(x[1])
        ke = np.exp(x[3]) if free_ke else P.ke
        eff = th + (dd if use_delta else 0.0)
        return (eff * lv / (km + lv) + ke * lv + x[2] - s) / SIGMA

    out = []
    for _ in range(reps):
        s = s_true + rng.normal(0.0, SIGMA, size=lv.size)
        r = least_squares(resid, x0, args=(s,), method="lm", xtol=1e-12, ftol=1e-12)
        out.append(r.x[0])
    out = np.array(out)
    return float(out.std(ddof=1)), float(out.mean() - np.log(P.theta))


def main():
    print(__doc__)
    T = 730.0
    cases = {
        "one dose cycle (n=7)": np.arange(0.0, D.DOSE_INTERVAL, 2.0),
        "two-daily x 2y (n=365)": np.arange(0.0, T, 2.0),
    }
    print("=" * 86)
    print("C2  Fisher prediction against Monte Carlo of the same estimator")
    print("    (prior-free Fisher, so it is the Cramer-Rao bound for the MC's own fit)")
    print("=" * 86)
    print(f"{'design':>24} {'k_e':>7} {'Delta':>7} | {'Fisher sd':>10} {'MC sd':>10} "
          f"{'MC/Fisher':>10} | {'MC bias':>9}")
    for name, t in cases.items():
        lv = D.quasi_steady_level(t)
        dd = np.asarray(D.delta(t), dtype=float)
        for free_ke in (False, True):
            for use_delta in (False, True):
                f, _ = fisher_sd(lv, dd, free_ke, use_delta, priors=False)
                m, bias = mc(lv, dd, free_ke, use_delta)
                print(f"{name:>24} {'free' if free_ke else 'known':>7} "
                      f"{'in' if use_delta else 'dropped':>7} | {f:10.4f} {m:10.4f} "
                      f"{m/f:10.3f} | {bias:+9.4f}")
    print()
    print("  MC/Fisher ~ 1 means the Fisher matrix describes the estimator exactly:")
    print("  no information is being counted twice (that would push the ratio ABOVE 1)")
    print("  and none is being thrown away (that would push it below).")
    print()
    print("=" * 86)
    print("C2b  the same, WITH the module's priors -- the published footing")
    print("=" * 86)
    for name, t in cases.items():
        lv = D.quasi_steady_level(t); dd = np.asarray(D.delta(t), dtype=float)
        for free_ke in (False, True):
            a, ra = fisher_sd(lv, dd, free_ke, False)
            b, rb = fisher_sd(lv, dd, free_ke, True)
            print(f"{name:>24} {'k_e free' if free_ke else 'k_e known':>10}:  "
                  f"published {a:.4f} (corr {ra:+.4f})   Delta reinstated {b:.4f} (corr {rb:+.4f})")


if __name__ == "__main__":
    main()
