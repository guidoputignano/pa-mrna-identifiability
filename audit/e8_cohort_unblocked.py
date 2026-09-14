"""
14.3 -- how many patients for the question that is NOT blocked.

The note's cohort arithmetic (~2,200) is for the DECISION question: placing the
fraction of patients whose margin is negative within ten points. That question is
declared blocked, for reasons that have nothing to do with cohort size -- the
threshold it is taken against is stipulated, not measured.

A different question is not blocked, and it is the one the design specification
actually makes: does a monthly outpatient record over two years identify a
patient's residual capacity to the stated precision? That is a WITHIN-patient
claim, so patients are needed for two things only:

    (a) to cover the heterogeneity, since the achieved precision depends on the
        patient's own parameters -- most of all on Delta/theta, because the dose
        is fixed in absolute terms while residual capacity is not; and
    (b) to confirm the predicted precision against the realised spread, which is
        a variance-estimation problem and needs far fewer patients than any
        question about a population fraction.

Ranges below are ASSUMPTIONS and are stated as such. Nothing measures them.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from design.run import P, load, delta, LOAD_FRAC, ASSAY_CV, LOAD_CV
import model_A.run as LOC

SIGMA = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
LB = load(0.0)
MONTHLY = np.arange(0.0, 730.0, 30.4)
TARGET = 0.20


def banner(t):
    print("\n" + "=" * 82)
    print("  " + t)
    print("=" * 82)


def sd_for(theta, Km, ke, load_mult, seed):
    """sd(ln theta) for one patient on the monthly design with the dose train."""
    inc = delta(MONTHLY)
    rng = np.random.default_rng(seed)
    Lp = LB * load_mult
    L = Lp * (1.0 + LOAD_CV * rng.standard_normal(len(MONTHLY)))
    B = ke * Km + theta + inc - L
    m = (-B + np.sqrt(B ** 2 + 4.0 * ke * L * Km)) / (2.0 * ke)
    F = np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * Lp) ** 2, 1.0 / 0.8 ** 2])
    for mi, di in zip(m, inc):
        q = np.array([theta * mi / (Km + mi),
                      -(theta + di) * Km * mi / (Km + mi) ** 2,
                      1.0, ke * mi])
        F = F + np.outer(q, q) / SIGMA ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def cohort(n=4000, seed=11):
    """A heterogeneous cohort. Every range here is assumed, not measured."""
    rng = np.random.default_rng(seed)
    lu = lambda lo, hi: np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    theta = P.theta * lu(1 / 3, 3.0)      # residual capacity, 9-fold span
    Km = P.Km * lu(1 / 3, 3.0)            # unknown, 9-fold span
    ke = P.ke * lu(1 / 2.2, 2.2)          # the prior's own factor
    lm = lu(0.5, 2.0)                     # managed load, 4-fold span
    sd = np.array([sd_for(t, k, e, l, 100 + i)
                   for i, (t, k, e, l) in enumerate(zip(theta, Km, ke, lm))])
    return theta, Km, ke, lm, sd


if __name__ == "__main__":
    print(__doc__)
    theta, Km, ke, lm, sd = cohort()
    ok = np.isfinite(sd)
    sdf = sd[ok]

    banner("THE DISTRIBUTION OF ACHIEVED PRECISION ACROSS A HETEROGENEOUS COHORT")
    for q in (5, 25, 50, 75, 95):
        print(f"    {q:2d}th percentile of sd(ln theta-hat): {np.percentile(sdf, q):.4f}")
    p_ok = float(np.mean(sdf <= TARGET))
    print(f"\n    fraction meeting the {TARGET} target: {100*p_ok:.1f}%  "
          f"({int(p_ok*sdf.size)} of {sdf.size})")

    print("\n    what drives the spread -- correlation of sd with each factor:")
    for nm, v in (("Delta/theta (dose over capacity)", delta(MONTHLY).mean() / theta),
                  ("theta", theta), ("K_m", Km), ("k_e", ke), ("load", lm)):
        r = np.corrcoef(np.log(v[ok]), np.log(sdf))[0, 1]
        print(f"      {nm:34s} r(log,log) = {r:+.3f}")

    banner("(a) PATIENTS TO BOUND THE FRACTION THAT MEETS THE TARGET")
    print("    A binomial proportion at p ~ %.2f, half-width w at 95%%:" % p_ok)
    for w in (0.20, 0.15, 0.10, 0.05):
        n = int(np.ceil(1.96 ** 2 * p_ok * (1 - p_ok) / w ** 2))
        print(f"      +/- {100*w:2.0f} points  ->  {n:4d} patients")

    banner("(b) PATIENTS TO CONFIRM THE PREDICTED PRECISION AGAINST THE REALISED SPREAD")
    print("    The design predicts a per-patient sd. Confirming it means estimating")
    print("    a standard deviation from N independent patients, whose own relative")
    print("    standard error is about 1/sqrt(2(N-1)).")
    for rel in (0.40, 0.30, 0.20, 0.15):
        n = int(np.ceil(1 + 1 / (2 * rel ** 2)))
        print(f"      confirm sd to +/- {100*rel:2.0f}% relative  ->  {n:3d} patients")

    banner("THE ANSWER, AND THE CONTRAST")
    n_frac = int(np.ceil(1.96 ** 2 * p_ok * (1 - p_ok) / 0.10 ** 2))
    n_conf = int(np.ceil(1 + 1 / (2 * 0.20 ** 2)))
    print(f"    Confirming the design achieves its stated precision: ~{n_conf} patients.")
    print(f"    Bounding the fraction of patients it works for to +/-10 points: ~{n_frac}.")
    print(f"    Both are inside a phase 1/2 trial of sixteen, or just outside it.")
    print(f"    The decision question, for contrast, needs about 2,200.")
    print()
    print("    That contrast is the point. The blocked question and the unblocked")
    print("    one differ by two orders of magnitude in cohort, and the note's")
    print("    specification should ask for the one that is affordable.")
