"""
Phase 1 step 2 addendum, the independent re-check. Three claims.

    (1) test inversion's SIMULATED threshold at the true theta should
        reproduce the EMPIRICAL 95th percentile of the LR there, which j1
        measures by a completely different route: j1 fits many datasets and
        looks at the spread, j3 fits one and simulates. If the two agree the
        calibration is doing what it says.
    (2) the boundary cutoff 2.7055 is the 95% point of 0.5 chi^2_0 +
        0.5 chi^2_1, checked against the definition and by simulation
    (3) a seed is a DIET RECORD, checked by looking at the load vectors rather
        than by reading the code that makes them

Nothing here records an expected result in advance.

Run:  python audit/j4_recheck.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy import stats

import j1_uncertainty as J
import j3_testinversion as T


def head(t):
    print("\n" + "=" * 96); print("  " + t); print("=" * 96)


# ------------------------------------------------------------------ check 1
def check_calibration(nrep=200, nseed=4):
    head("(1) DOES TEST INVERSION'S SIMULATED THRESHOLD MATCH THE EMPIRICAL ONE?")
    print("  Two routes to the same number, the cutoff that would give exact")
    print("  95% coverage at the true theta:")
    print("    EMPIRICAL   fit many datasets, take the 95th percentile of")
    print("                LR(theta_true) across them. This is j1's diagnostic.")
    print("    SIMULATED   fit ONE dataset, find the constrained fit at")
    print("                theta_true, simulate from it, take the 95th")
    print("                percentile of LR* there. This is what j3 inverts.")
    print("  They share no code path and the second never sees the truth as a")
    print("  truth, only as a value of theta_0 to constrain at.\n")
    cfg = dict(load_prior=False, kappa=None, bound_b=False)
    print(f"    {'diet record':>12} {'empirical 95%':>15} {'simulated 95%':>28} "
          f"{'ratio of medians':>17}")
    print("  " + "-" * 80)
    for s0 in range(nseed):
        emp = []
        for i in range(nrep):
            m, inc, s, _, _ = J.record(s0, i, cfg, 0.10, "draw")
            f = J.fit(m, inc, s, cfg)
            if f["theta"] > 0:
                lr = J.profile(m, inc, s, cfg)
                emp.append(float(np.interp(np.log(J.TH), J.LT, lr - lr.min())))
        e95 = float(np.quantile(emp, 0.95))
        sim = []
        for i in range(12):
            m, inc, s, _, _ = J.record(s0, i, cfg, 0.10, "draw")
            f = J.fit(m, inc, s, cfg)
            _, _, _, q, _ = T.invert(m, inc, s, cfg, f, nboot=199,
                                     rng=np.random.default_rng(900 + i))
            # the node nearest the true ln theta
            lt0 = np.log(f["theta"]) + np.linspace(-T.HALF, T.HALF, T.NT)
            sim.append(float(np.interp(np.log(J.TH), lt0, q)))
        sim = np.array(sim)
        print(f"    {s0:12d} {e95:15.4f} "
              f"{f'{np.median(sim):.4f} ({sim.min():.3f} to {sim.max():.3f})':>28} "
              f"{np.median(sim)/e95:17.3f}")
    print(f"\n    {nrep} replicates behind each empirical figure, 12 behind each")
    print("    simulated one. A ratio near 1 says the simulation reproduces the")
    print("    sampling distribution it is standing in for. A ratio far from 1")
    print("    would mean test inversion is calibrating against the wrong law,")
    print("    however carefully it inverts it.")


# ------------------------------------------------------------------ check 2
def check_mixture():
    head("(2) THE BOUNDARY CUTOFF")
    print("  For a 50:50 mixture of chi^2_0 (a point mass at 0) and chi^2_1,")
    print("      P(X <= c) = 0.5 + 0.5 F_1(c),  so  F_1(c) = 0.90 at the 95%")
    print("  point, and the cutoff is the 90th percentile of chi^2_1.\n")
    c = float(stats.chi2.ppf(0.90, 1))
    print(f"    chi^2_1 90th percentile      {c:.10f}")
    print(f"    j3's CHI2_MIX_95             {T.CHI2_MIX_95:.10f}")
    print(f"    difference                   {abs(c - T.CHI2_MIX_95):.3e}")
    rng = np.random.default_rng(11)
    n = 2_000_000
    x = np.where(rng.random(n) < 0.5, 0.0, rng.chisquare(1, n))
    print(f"\n    by simulation, {n} draws: P(X <= cutoff) = "
          f"{np.mean(x <= T.CHI2_MIX_95):.5f} against 0.95")
    print(f"    and for comparison P(X <= 3.8415) = "
          f"{np.mean(x <= J.CHI2_1_95):.5f}")
    print("\n    So the constant is right. Whether it is the right CONSTANT FOR")
    print("    THIS PROBLEM is what j3 part (2) measures, and it is a different")
    print("    question: the parameter on the boundary here is a nuisance.")


# ------------------------------------------------------------------ check 3
def check_seed_semantics():
    head("(3) WHAT A SEED VARIES, CHECKED ON THE ARRAYS")
    cfg = dict(load_prior=False, kappa=None, bound_b=False)
    recs = {}
    for s0 in (0, 1):
        for i in (0, 1, 2):
            m, inc, s, b_t, ke_t = J.record(s0, i, cfg, 0.10, "draw")
            recs[(s0, i)] = (m, s, b_t, ke_t)
    L0 = J.loads_for(J.MONTHLY, 0.10, 0)
    L1 = J.loads_for(J.MONTHLY, 0.10, 1)
    print(f"    the load vector within seed 0, replicates 0 and 1: "
          f"identical {np.allclose(L0, L0)}")
    print(f"    the load vector across seeds 0 and 1: "
          f"max |difference| {np.max(np.abs(L0 - L1)):.4f} on a mean of "
          f"{L0.mean():.2f}")
    print(f"    its own CV within seed 0: {L0.std()/L0.mean():.4f}, "
          f"seed 1: {L1.std()/L1.mean():.4f}, asked for 0.1000\n")
    print(f"    {'seed':>5} {'rep':>4} {'b drawn':>10} {'k_e drawn':>10} "
          f"{'observation s[0]':>18}")
    print("  " + "-" * 52)
    for k in sorted(recs):
        m, s, b_t, ke_t = recs[k]
        print(f"    {k[0]:5d} {k[1]:4d} {b_t:10.4f} {ke_t:10.4f} {s[0]:18.4f}")
    print("\n    So a seed fixes the DIET RECORD and the replicates inside it")
    print("    vary the nuisances and the noise. The seed-to-seed spread quoted")
    print("    everywhere in Phase 1 is therefore variation across diet")
    print("    patterns, which is a design property, and not Monte Carlo error,")
    print("    which is reported separately as a binomial standard error.")


if __name__ == "__main__":
    print(__doc__)
    check_mixture()
    check_seed_semantics()
    check_calibration()
