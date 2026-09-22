"""
THE TIER 1 MARGIN: WHAT EACH READING OF THE INCREMENT'S REPEATABILITY COSTS.

Scoping only. Nothing here reaches the paper.

The specification asks for the delivered increment to be characterised to a
coefficient of variation below 11.6%, and that requirement is derived rather
than chosen: it is the CV at which monthly sampling lands exactly on
sd(ln theta) = 0.20. The only published repeatability for the assay is
between-visit, 14.0% with a range of 0.4 to 50.8%, in fifteen methylmalonic
acidemia patients at visits two months to four years apart. Nobody has measured
the same-day figure the design actually needs.

We cannot measure it either. What we can do is price every value it might take,
so that when the measurement arrives its consequence is already known.

HOW THE INCREMENT'S PRECISION ENTERS. The increment is carried as a regressor
with an unknown scale kappa, so the design matrix gains a fifth column,
d s / d kappa = Delta(t) g(m), and kappa carries a log-prior of the stated CV.
sd(ln theta) is then the (1,1) entry of the inverse of the five-parameter
information matrix. That is the same construction spec_numbers() uses to solve
for the 11.6%, run forwards instead of inverted.

THE CAPACITY QUESTION IS ASKED SEPARATELY, because Section 6.2 reports that at
a disease-range capacity the estimate is gone for reasons that have nothing to
do with the increment. If that is right, increment precision should not matter
there at all, and a sweep is the way to find out rather than assume.

Run:  python audit/t2_tier1_margin.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

import model_A.run as LOC
from common.models import ParamsA
from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, SEEDS, load, median_spread)

KM, KE = 100.0, 5.0
LB = load(0.0)
KE_PRIOR_LOG = 0.8
GATE = 0.20
MONTH = 30.4

#: the published between-visit figure, and the requirement derived from the gate
PUBLISHED_CV = 0.1402
REQUIRED_CV = 0.1162

T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))


def op_level(theta):
    B = KE * KM + theta - LB
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)


def sigma_of(theta):
    return LOC.sigma_total(LOAD_FRAC, ASSAY_CV,
                           p=ParamsA(theta=theta, Km=KM, ke=KE),
                           m_bar=op_level(theta))


def sd_of(theta, n_visits, kappa_cv, seed):
    """
    sd(ln theta) on n monthly visits with the increment known to `kappa_cv`.

    kappa_cv None means the increment is known exactly, which is the footing
    every other number in the paper is computed on.
    """
    t = np.arange(0.0, n_visits * MONTH, MONTH)[:int(n_visits)]
    tau = np.mod(t - T_DOSE, DOSE_INTERVAL)
    inc = PEAK * np.exp(-LAMBDA * tau)
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    sig = sigma_of(theta)

    pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2, KE_PRIOR_LOG ** 2]
    if kappa_cv is not None:
        pri.append(kappa_cv ** 2)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    for m, d in zip(lv, inc):
        gm = m / (KM + m)
        q = [theta * gm, -(theta + d) * KM * m / (KM + m) ** 2, 1.0, KE * m]
        if kappa_cv is not None:
            q.append(d * gm)
        q = np.asarray(q)
        F = F + np.outer(q, q) / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over_seeds(theta, n_visits, kappa_cv):
    return median_spread([sd_of(theta, n_visits, kappa_cv, s) for s in SEEDS])


def visits_to_gate(theta, kappa_cv, nmax=600):
    """Monthly visits to reach the gate, median over records. None if never."""
    for n in range(4, nmax + 1):
        if over_seeds(theta, n, kappa_cv)[0] < GATE:
            return n
    return None


def floor_at(theta, kappa_cv, n_big=20000):
    """
    Where the design lands with the visit count unbounded.

    The kappa prior does not shrink with n, so there is a floor, and the CV at
    which that floor crosses the gate is the CV beyond which no number of
    visits is enough.
    """
    return over_seeds(theta, n_big, kappa_cv)[0]


def main():
    print(__doc__)
    bar = "=" * 78
    N25 = 25                                   # the recommended record

    print(bar)
    print("1. THE SWEEP, ON THE RECOMMENDED MONTHLY RECORD, AT NOMINAL CAPACITY")
    print(bar)
    print(f"  {N25} monthly draws over two years, medians over {len(SEEDS)} "
          f"intake records\n")
    print(f"  {'increment CV':>13} {'sd(ln theta) at n=25':>34} "
          f"{'visits to 0.20':>15}")
    exact = over_seeds(3000.0, N25, None)
    print(f"  {'exact':>13} {exact[0]:11.4f} [{exact[1]:.4f}, {exact[2]:.4f}]"
          f"      {visits_to_gate(3000.0, None) or 'never':>15}")
    rows = []
    for cv in (0.02, 0.04, 0.06, 0.08, 0.10, 0.1162, 0.14, 0.1402,
               0.16, 0.20, 0.25, 0.30):
        med, lo, hi = over_seeds(3000.0, N25, cv)
        v = visits_to_gate(3000.0, cv)
        rows.append((cv, med, lo, hi, v))
        mark = " *" if med < GATE else "  "
        tag = ""
        if abs(cv - REQUIRED_CV) < 1e-9:
            tag = "   <- the requirement"
        if abs(cv - PUBLISHED_CV) < 1e-9:
            tag = "   <- the published between-visit figure"
        print(f"  {100*cv:12.2f}% {med:11.4f} [{lo:.4f}, {hi:.4f}]{mark}"
              f"     {(v if v else 'never'):>15}{tag}")

    print("\n" + bar)
    print("2. THE CV AT WHICH NO NUMBER OF VISITS IS ENOUGH")
    print(bar)
    print("  With the visit count unbounded the kappa prior still does not")
    print("  shrink, so the design has a floor. Where that floor crosses the")
    print("  gate is where the requirement stops being a question of effort.\n")
    for cv in (0.10, 0.1162, 0.14, 0.20, 0.30, 0.50):
        print(f"    CV {100*cv:5.2f}%: floor {floor_at(3000.0, cv):.4f}"
              f"{'  <- inside the gate' if floor_at(3000.0, cv) < GATE else ''}")
    try:
        cv_star = brentq(lambda c: floor_at(3000.0, c) - GATE, 0.05, 2.0,
                         xtol=1e-4)
        print(f"\n  the floor reaches the gate at an increment CV of "
              f"{100*cv_star:.1f}%")
        print(f"  beyond that, no number of monthly visits reaches "
              f"sd(ln theta) = {GATE}")
    except ValueError:
        print("\n  the floor does not cross the gate over [5%, 200%]")

    print("\n" + bar)
    print("3. THE TWO QUESTIONS")
    print(bar)
    pub = over_seeds(3000.0, N25, PUBLISHED_CV)
    vpub = visits_to_gate(3000.0, PUBLISHED_CV)
    print(f"  a. At the published {100*PUBLISHED_CV:.1f}%, does the design pass?")
    print(f"     on the recommended {N25} draws: {pub[0]:.4f} "
          f"[{pub[1]:.4f}, {pub[2]:.4f}] against a gate of {GATE}")
    print(f"     -> {'PASSES' if pub[0] < GATE else 'FAILS'} at n = {N25}")
    print(f"     visits needed at that CV: "
          f"{vpub if vpub else 'no number up to 600'}")
    if vpub:
        print(f"     which is {vpub} monthly draws, "
              f"{vpub * MONTH / 365.25:.1f} years of follow-up")
    print(f"\n  b. How much better must a same-day pair be?")
    print(f"     the requirement is {100*REQUIRED_CV:.2f}% and the published "
          f"between-visit figure is {100*PUBLISHED_CV:.2f}%,")
    print(f"     so a same-day pair must be a factor of "
          f"{PUBLISHED_CV/REQUIRED_CV:.2f} tighter,")
    print(f"     a {100*(1-REQUIRED_CV/PUBLISHED_CV):.0f}% reduction in the "
          f"coefficient of variation.")
    #  and the same question asked the other way: how many visits buy the gap
    print(f"     Bought with visits instead: at the published CV the gate is "
          f"reached at {vpub if vpub else '>600'} draws")
    print(f"     against {visits_to_gate(3000.0, None)} with the increment known "
          f"exactly.")

    print("\n" + bar)
    print("4. DOES ANY OF THIS MATTER IN THE DISEASE RANGE?")
    print(bar)
    print("  Section 6.2 says the estimate is gone at a disease-range capacity")
    print("  for reasons that have nothing to do with the increment. If that is")
    print("  right, increment precision should change nothing here.\n")
    for theta, lab in ((3000.0, "100%"), (300.0, "10%"), (150.0, "5%")):
        line = [f"  theta = {theta:6.0f} ({lab:>4}):"]
        for cv, nm in ((None, "exact"), (0.02, "2%"), (REQUIRED_CV, "11.6%"),
                       (0.30, "30%")):
            line.append(f"{nm} {over_seeds(theta, 25, cv)[0]:8.4f}")
        print("   ".join(line))

    print("\n" + bar)
    print("R1. RE-CHECK: the solved requirement reproduces the gate")
    print(bar)
    got = over_seeds(3000.0, N25, REQUIRED_CV)[0]
    print(f"  at the requirement's own CV of {100*REQUIRED_CV:.2f}%, "
          f"sd(ln theta) = {got:.4f}")
    print(f"  which should sit on the gate of {GATE}: "
          f"off by {abs(got-GATE)/GATE:.2%}")

    print("\nR2. RE-CHECK: an exact increment reproduces the paper")
    print(f"  {exact[0]:.4f} [{exact[1]:.4f}, {exact[2]:.4f}] against Table 3's "
          f"0.1699 [0.1684, 0.1709]")

    print("\nR3. RE-CHECK: the sweep is monotone in the CV")
    seq = [r[1] for r in rows]
    print(f"  sd rises with every loosening of the prior: "
          f"{all(b >= a for a, b in zip(seq, seq[1:]))}")


if __name__ == "__main__":
    main()
