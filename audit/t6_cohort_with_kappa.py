"""
4. COHORT COVERAGE WITH THE INCREMENT CV DRAWN ALONGSIDE THE OTHER FOUR.

Scoping only. Nothing here reaches the paper.

Section 5.5 sweeps four thousand simulated patients over theta, K_m, k_e and
managed load, and reports that just over half meet the target. It holds the
delivered increment EXACTLY KNOWN, which is the footing every other number in
the paper is computed on and which measurement B has since shown is not the
patient's situation: the published per-subject repeatability of the breath
test runs from 0.39% to 50.84%, and beyond an increment CV of about 18.7% no
number of monthly visits reaches the gate at all.

So the cohort as published varies the four things that are properties of the
patient and holds fixed the one that is a property of the MEASUREMENT, and it
is the measurement that has a floor. This draws the increment CV as a fifth
factor, from the published per-subject range, and asks what the coverage is.

THE INCREMENT ENTERS AS Proposition B's FIFTH COLUMN. The delivered increment
is Delta(t) with an unknown scale kappa, so the design matrix gains a column
d s / d kappa = Delta(t) g(m) and kappa carries a log-prior of the stated CV.
That prior does not shrink with the record, which is where the floor comes
from.

FOUR ARMS, so the change is attributable rather than merely observed:

    A  the published cohort, increment exact            (the baseline)
    B  the same cohort, increment at the published MEAN of 14.02%
    C  the same cohort, increment drawn log-uniform over 0.39%-50.84%
    D  arm C with theta held at the illustrative 3000, to separate the
       increment's contribution from the capacity heterogeneity's

Every patient in every arm carries the same seed, so the arms differ in the
increment and in nothing else.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import model_A.run as LOC
from design.run import (ASSAY_CV, LOAD_CV, LOAD_FRAC, P, delta, load)

SIGMA = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
LB = load(0.0)
MONTHLY = np.arange(0.0, 730.0, 30.4)
INC = np.asarray(delta(MONTHLY), float)
TARGET = 0.20
KE_PRIOR_LOG = 0.8

#: the published per-subject range of the breath test's repeatability
CV_LO, CV_HI, CV_MEAN = 0.0039, 0.5084, 0.1402
#: the floor measurement B found: beyond this no number of visits suffices
CV_FLOOR = 0.187


def sd_for(theta, Km, ke, load_mult, seed, kappa_cv=None):
    """
    sd(ln theta) for one patient, optionally with the increment scale free.

    `kappa_cv is None` reproduces Section 5.5's own footing exactly; a value
    adds Proposition B's fifth column and its log-prior.
    """
    rng = np.random.default_rng(seed)
    Lp = LB * load_mult
    L = Lp * (1.0 + LOAD_CV * rng.standard_normal(MONTHLY.size))
    B = ke * Km + theta + INC - L
    m = (-B + np.sqrt(B ** 2 + 4.0 * ke * L * Km)) / (2.0 * ke)

    pri = [np.inf, np.inf, (LOAD_FRAC * Lp) ** 2, KE_PRIOR_LOG ** 2]
    if kappa_cv is not None:
        pri.append(kappa_cv ** 2)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    for mi, di in zip(m, INC):
        gm = mi / (Km + mi)
        q = [theta * gm, -(theta + di) * Km * mi / (Km + mi) ** 2, 1.0, ke * mi]
        if kappa_cv is not None:
            q.append(di * gm)
        q = np.asarray(q)
        F = F + np.outer(q, q) / SIGMA ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def draw(n=4000, seed=11):
    """Section 5.5's own cohort, seed and ranges unchanged."""
    rng = np.random.default_rng(seed)
    lu = lambda lo, hi: np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    theta = P.theta * lu(1 / 3, 3.0)
    Km = P.Km * lu(1 / 3, 3.0)
    ke = P.ke * lu(1 / 2.2, 2.2)
    lm = lu(0.5, 2.0)
    #: drawn AFTER the other four, from a separate generator, so the four
    #: published columns are bit-identical to Section 5.5's
    cv = np.exp(np.random.default_rng(seed + 1).uniform(
        np.log(CV_LO), np.log(CV_HI), n))
    return theta, Km, ke, lm, cv


def summarise(name, sd, extra=""):
    """
    Coverage over the WHOLE cohort, not over the identifiable part of it.

    A patient the arm cannot identify at all returns inf, and dropping those
    before taking the fraction would count a design's failures as absences and
    report coverage that rises as the design gets worse. They are misses.
    """
    ok = np.isfinite(sd)
    s = sd[ok]
    q = [np.percentile(s, p) for p in (5, 25, 50, 75, 95)]
    frac = 100.0 * np.mean(np.where(ok, sd, np.inf) <= TARGET)
    print(f"  {name:<44} {q[2]:7.4f}  [{q[1]:.3f}, {q[3]:.3f}]  "
          f"{frac:5.1f}%  {sd.size - ok.sum():5d}  {extra}")
    return frac, q[2]


def cohort_numbers():
    """
    Section 5.5's coverage, with the increment CV drawn as a fifth factor.

    The published sweep varies the four things that are properties of the
    PATIENT -- capacity, half-saturation, alternative clearance, managed load
    -- and holds fixed the one that is a property of the MEASUREMENT. It is
    the measurement that has a floor: beyond an increment CV of about 18.7%
    no number of monthly visits reaches the gate, and the published
    per-subject repeatability runs to 50.84%.

    Arm A is kept as a generated key so the published figure and the corrected
    one cannot drift apart, and so a reader can see what changed.
    """
    from common.numbers import write_numbers
    theta, Km, ke, lm, cv = draw()
    seeds = 100 + np.arange(theta.size)
    A = np.array([sd_for(t, k, e, l, s)
                  for t, k, e, l, s in zip(theta, Km, ke, lm, seeds)])
    C = np.array([sd_for(t, k, e, l, s, c)
                  for t, k, e, l, s, c in zip(theta, Km, ke, lm, seeds, cv)])

    def frac(v):
        return 100.0 * np.mean(np.where(np.isfinite(v), v, np.inf) <= TARGET)

    N = {
        "cohort_n": theta.size,
        "cov_exact_pct": frac(A),
        "cov_drawn_pct": frac(C),
        "cohort_med_exact": float(np.median(A[np.isfinite(A)])),
        "cohort_med_drawn": float(np.median(C[np.isfinite(C)])),
        "cohort_q1_drawn": float(np.percentile(C[np.isfinite(C)], 25)),
        "cohort_q3_drawn": float(np.percentile(C[np.isfinite(C)], 75)),
        "cv_floor_pct": 100.0 * CV_FLOOR,
        "cv_hi_pct": 100.0 * CV_HI,
        "above_floor_pct": 100.0 * float(np.mean(cv > CV_FLOOR)),
        "unreachable_pct": 100.0 * float(np.mean(
            (cv > CV_FLOOR) & np.isfinite(C) & (C > TARGET))),
    }
    fmt = {k: "{:.1f}" for k in N}
    fmt.update({"cohort_n": "{:.0f}",
                "cohort_med_exact": "{:.3f}", "cohort_med_drawn": "{:.3f}",
                "cohort_q1_drawn": "{:.3f}", "cohort_q3_drawn": "{:.3f}",
                "cv_floor_pct": "{:.1f}", "cv_hi_pct": "{:.1f}"})
    write_numbers("cohort", N, fmt)
    return N


def main():
    print(__doc__)
    bar = "=" * 92
    theta, Km, ke, lm, cv = draw()
    n = theta.size
    seeds = 100 + np.arange(n)

    print(bar)
    print("1. THE FOUR ARMS")
    print(bar)
    print(f"  {n} patients, the same four published ranges and the same "
          f"per-patient seeds\n")
    print(f"  {'arm':<44} {'median':>7}  {'quartiles':<16} "
          f"{'meet 0.20':>9} {'unident':>6}")

    A = np.array([sd_for(t, k, e, l, s)
                  for t, k, e, l, s in zip(theta, Km, ke, lm, seeds)])
    fA, mA = summarise("A  increment exact (Section 5.5 as published)", A)

    B = np.array([sd_for(t, k, e, l, s, CV_MEAN)
                  for t, k, e, l, s in zip(theta, Km, ke, lm, seeds)])
    fB, mB = summarise(f"B  increment at the published mean, {100*CV_MEAN:.2f}%", B)

    C = np.array([sd_for(t, k, e, l, s, c)
                  for t, k, e, l, s, c in zip(theta, Km, ke, lm, seeds, cv)])
    fC, mC = summarise("C  increment CV drawn, 0.39%-50.84% log-uniform", C)

    D = np.array([sd_for(P.theta, k, e, l, s, c)
                  for k, e, l, s, c in zip(Km, ke, lm, seeds, cv)])
    fD, mD = summarise("D  as C but theta held at the illustrative 3000", D)

    print("\n" + bar)
    print("2. WHAT THE INCREMENT COSTS, ATTRIBUTED")
    print(bar)
    print(f"  coverage, published footing                  {fA:5.1f}%")
    print(f"  at the published MEAN repeatability          {fB:5.1f}%   "
          f"({fB - fA:+.1f} points)")
    print(f"  with the repeatability DRAWN per patient     {fC:5.1f}%   "
          f"({fC - fA:+.1f} points)")
    print(f"  the median precision moves {mA:.4f} -> {mC:.4f}, "
          f"a factor of {mC / mA:.2f}")
    print()
    frac_above = 100.0 * np.mean(cv > CV_FLOOR)
    print(f"  {frac_above:.1f}% of the cohort draws an increment CV above the "
          f"{100*CV_FLOOR:.1f}% floor,")
    print(f"  beyond which no number of monthly visits reaches the gate at "
          f"all.")
    both = 100.0 * np.mean((cv > CV_FLOOR) & np.isfinite(C) & (C > TARGET))
    print(f"  {both:.1f}% of the cohort both draws above the floor and misses "
          f"the target.")

    print("\n" + bar)
    print("3. COVERAGE AGAINST THE REPEATABILITY, WHICH IS THE PROTOCOL'S "
          "CHOICE")
    print(bar)
    print("  If Tier 1 reports a repeatability and the protocol enrols on it,")
    print("  this is the coverage that follows.\n")
    print(f"  {'increment CV':>14} {'median sd':>11} {'meet 0.20':>11}")
    for c in (None, 0.01, 0.05, 0.10, CV_MEAN, 0.187, 0.25, 0.51):
        v = np.array([sd_for(t, k, e, l, s, c)
                      for t, k, e, l, s in zip(theta, Km, ke, lm, seeds)])
        f = np.isfinite(v)
        nm = "exact" if c is None else f"{100*c:.2f}%"
        print(f"  {nm:>14} {np.median(v[f]):11.4f} "
              f"{100*np.mean(v[f] <= TARGET):10.1f}%")

    print("\n" + bar)
    print("R1. RE-CHECK: arm A reproduces Section 5.5")
    print(bar)
    ok = np.isfinite(A)
    print(f"  median {np.median(A[ok]):.4f} against the published 0.192")
    print(f"  quartiles {np.percentile(A[ok], 25):.4f} and "
          f"{np.percentile(A[ok], 75):.4f} against the published 0.132 and "
          f"0.329")
    print(f"  fraction meeting the target {fA:.1f}%, which the paper calls "
          f"'just over half'")
    print("\nR2. RE-CHECK: the fifth column at a vanishing CV returns arm A")
    #  1e-3, not 1e-6: the kappa prior enters the information matrix as
    #  1/cv^2, so a probe at 1e-6 puts 1e12 on the diagonal and the rank test
    #  below, which is relative to the largest eigenvalue, then rejects
    #  perfectly good matrices. 1e-3 is still a thousand times tighter than
    #  any repeatability the literature reports.
    v = np.array([sd_for(t, k, e, l, s, 1e-3)
                  for t, k, e, l, s in zip(theta[:200], Km[:200], ke[:200],
                                           lm[:200], seeds[:200])])
    #  a patient the published footing already cannot identify returns inf in
    #  both arms, and inf/inf is not a disagreement. Compare where both are
    #  finite, and report the infinities as a count rather than folding them
    #  into a maximum that would then always be inf.
    a = A[:200]
    both = np.isfinite(v) & np.isfinite(a)
    d = float(np.max(np.abs(v[both] - a[both]) / a[both])) if both.any() else 0.0
    print(f"  worst relative difference over the {both.sum()} patients finite "
          f"in both: {d:.2e}")
    print(f"  the two arms agree on which patients are identifiable at all "
          f"in {int(np.sum(np.isfinite(v) == np.isfinite(a)))} of 200")


if __name__ == "__main__":
    main()
    cohort_numbers()
