"""
10. THE LOAD-ACCURACY x CAPACITY GRID, A BOX NOBODY HAS EVALUATED.

Scoping only. Nothing here reaches the paper.

Every load-accuracy sweep in this repository runs at nominal capacity
(`table_load_precision`, Section 5.3). Every capacity sweep runs at the
specification's own 5% load accuracy (`n4_theta_scope`, Section 6.2). The cell
where the two cross has never been computed, and the flux identity says it
should not be flat.

THE PROPAGATION, WHICH IS WHY THE QUESTION IS WORTH ASKING. The identity is

    d s / d ln theta  =  theta g(m)  =  L - k_e m,

a difference of two large numbers at low capacity. At nominal the difference
is 322.7 - 50 = 272.7, so a 5% error on L, which is 16.1, is about 6% of the
quantity being estimated. At 5% of nominal it is 322.7 - 270 = 52.7, and the
same 16.1 is 31%. Six times more expensive -- or, read the other way,
tightening the load record should buy about six times more in the disease
range than at nominal.

That is a single-parameter propagation and it holds everything else fixed. The
full information matrix need not agree: theta, K_m and k_e are correlated, and
the correlations can eat an apparent gain. So this measures it rather than
arguing it.

WHAT "LOAD ACCURACY" MEANS HERE, and it is two things, not one. Measurement C2
established that the paper's own sweep moves both at once: the systematic
offset, which enters as the prior on b and runs along the fibre, and the
per-visit scatter, which enters only through sigma. At nominal capacity 98% of
the degradation is the scatter and 2% is the offset. Part 1 moves them
together, which is the paper's own convention and the comparable number; Part 2
splits them at every capacity, because whether that 98/2 split survives into
the disease range is exactly what decides which study would buy anything.

Convention A throughout, the paper's own: the recorded load L is held at its
theta = 3000 value while the marker rises, so a column of this grid is one
diet across capacities.

Nothing here records an expected result in advance.

Run:  python audit/t8_load_by_capacity.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import model_A.run as LOC
from common.models import ParamsA
from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, SEEDS, delta, load, median_spread)

KM, KE, TH0 = 100.0, 5.0, 3000.0
LB = load(0.0)
KE_PRIOR_LOG = 0.8
GATE = 0.20
MONTHLY = np.arange(0.0, 730.0, 30.4)
INC = np.asarray(delta(MONTHLY), float)
CAPS = (1.00, 0.50, 0.20, 0.10, 0.05, 0.02)
ACCS = (0.01, 0.02, 0.05, 0.10, 0.20)


def level(theta, inc, L):
    B = KE * KM + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * L * KM)) / (2.0 * KE)


def sd_of(theta, seed, bias, scatter):
    """
    sd(ln theta) with the two halves of the load record set independently.

    `bias` is the prior width on the constant offset, as a fraction of the
    load; `scatter` is the per-visit relative error, which reaches the fit
    only through sigma. bias = scatter is the paper's own sweep.
    """
    L = LB * (1.0 + LOAD_CV *
              np.random.default_rng(seed).standard_normal(MONTHLY.size))
    m = level(theta, INC, L)
    sig = LOC.sigma_total(scatter, ASSAY_CV,
                          p=ParamsA(theta=theta, Km=KM, ke=KE),
                          m_bar=level(theta, 0.0, LB))
    pri_b = np.inf if bias is None else (bias * LB) ** 2
    F = np.diag([0.0, 0.0, 0.0 if np.isinf(pri_b) else 1.0 / pri_b,
                 1.0 / KE_PRIOR_LOG ** 2])
    for mi, di in zip(m, INC):
        gm = mi / (KM + mi)
        q = np.array([theta * gm, -(theta + di) * KM * mi / (KM + mi) ** 2,
                      1.0, KE * mi])
        F = F + np.outer(q, q) / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over(theta, bias, scatter):
    return median_spread([sd_of(theta, s, bias, scatter) for s in SEEDS])[0]


def enz_share(theta):
    m = level(theta, 0.0, LB)
    return (LB - KE * m) / LB, LB - KE * m


def loadcap_numbers():
    """
    The load-accuracy x capacity results Sections 8 and 11 quote.

    WHY THIS FUNCTION NOW EXISTS. See REVIEW_LOG.md R81: this script's
    docstring said nothing of it reached the paper, and it owns route 5's
    factors and Section 8's capacity split.
    """
    from common.numbers import write_numbers
    r = {}
    for cap in CAPS:
        th = 3000.0 * cap
        r[cap] = over(th, 0.20, 0.20) / over(th, 0.01, 0.01)
    dis = [v for c, v in r.items() if c <= 0.20]
    N = {"loadcap_ratio_nominal": r[1.00],
         "loadcap_ratio_dis_lo": min(dis), "loadcap_ratio_dis_hi": max(dis)}
    for cap, tag in ((1.00, "nominal"), (0.50, "half"), (0.20, "fifth")):
        th = 3000.0 * cap
        b0 = over(th, LOAD_FRAC, LOAD_FRAC)
        bo = over(th, 0.20, LOAD_FRAC)
        sc = over(th, LOAD_FRAC, 0.20)
        N["loadcap_share_" + tag] = 100.0 * (bo - b0) / ((bo - b0) + (sc - b0))
    N["loadcap_perfect_5pct"] = over(150.0, 1e-6, 1e-6)
    #  each format is the precision the PAPER prints, so a macro renders the
    #  literal it replaces exactly
    fmt = {k: "{:.2f}" for k in N}
    fmt.update({"loadcap_share_nominal": "{:.1f}",
                "loadcap_share_half": "{:.1f}", "loadcap_share_fifth": "{:.1f}",
                "loadcap_perfect_5pct": "{:.3f}"})
    write_numbers("loadcap", N, fmt)
    return N


def main():
    print(__doc__)
    bar = "=" * 84

    print(bar)
    print("0. THE PROPAGATION THIS IS TESTING, RE-DERIVED")
    print(bar)
    print(f"  {'theta':>10} {'% nominal':>10} {'m':>9} {'L - k_e m':>11} "
          f"{'5% of L as a % of it':>22}")
    for f in CAPS:
        th = TH0 * f
        sh, enz = enz_share(th)
        print(f"  {th:10.0f} {100*f:9.0f}% {level(th, 0.0, LB):9.3f} "
              f"{enz:11.2f} {100*LOAD_FRAC*LB/enz:21.1f}%")

    print("\n" + bar)
    print("1. THE GRID, BOTH CHANNELS MOVED TOGETHER (THE PAPER'S CONVENTION)")
    print(bar)
    print(f"  sd(ln theta), medians over {len(SEEDS)} intake records, monthly "
          f"design\n")
    print(f"  {'theta':>12} " + "".join(f"{f'{100*a:.0f}%':>12}" for a in ACCS))
    grid = {}
    for f in CAPS:
        th = TH0 * f
        row = [f"  {f'{100*f:.0f}% of nom':>12}"]
        for a in ACCS:
            v = over(th, a, a)
            grid[(f, a)] = v
            mark = "*" if v < GATE else " "
            row.append(f"{v:11.4f}{mark}")
        print("".join(row))
    print("\n  * meets the 0.20 gate")

    print("\n" + bar)
    print("2. WHAT TIGHTENING THE RECORD BUYS, BY CAPACITY")
    print(bar)
    print("  the ratio sd(at 20%) / sd(at 1%): how much the whole sweep is")
    print("  worth at each capacity. If the propagation holds this rises as")
    print("  capacity falls.\n")
    print(f"  {'theta':>12} {'20% -> 1%':>12} {'20% -> 5%':>12} "
          f"{'5% -> 1%':>12}")
    for f in CAPS:
        g = grid
        print(f"  {f'{100*f:.0f}% of nom':>12} "
              f"{g[(f,0.20)]/g[(f,0.01)]:12.3f} "
              f"{g[(f,0.20)]/g[(f,0.05)]:12.3f} "
              f"{g[(f,0.05)]/g[(f,0.01)]:12.3f}")
    base = grid[(1.00, 0.20)] / grid[(1.00, 0.01)]
    print(f"\n  against the nominal row, the whole sweep is worth")
    for f in CAPS:
        r = (grid[(f, 0.20)] / grid[(f, 0.01)]) / base
        print(f"    {100*f:5.0f}% of nominal: {r:6.2f} times as much")

    print("\n" + bar)
    print("3. THE SPLIT AT EACH CAPACITY: WHICH CHANNEL CARRIES IT")
    print(bar)
    print("  C2 found 98% scatter / 2% offset at nominal. Whether that")
    print("  survives into the disease range decides which study to run.\n")
    print(f"  {'theta':>12} {'5/5 baseline':>13} {'offset at 20%':>14} "
          f"{'scatter at 20%':>15} {'offset share':>13}")
    for f in CAPS:
        th = TH0 * f
        b0 = over(th, LOAD_FRAC, LOAD_FRAC)
        bo = over(th, 0.20, LOAD_FRAC)
        sc = over(th, LOAD_FRAC, 0.20)
        lb, ls = bo - b0, sc - b0
        share = 100.0 * lb / (lb + ls) if (lb + ls) > 0 else float("nan")
        print(f"  {f'{100*f:.0f}% of nom':>12} {b0:13.4f} {bo:14.4f} "
              f"{sc:15.4f} {share:12.1f}%")

    print("\n" + bar)
    print("4. DOES ANY LOAD ACCURACY OPEN ANY PART OF THE DISEASE RANGE?")
    print(bar)
    for f in (0.10, 0.05, 0.02):
        th = TH0 * f
        v = [over(th, a, a) for a in (0.001, 0.005, 0.01)]
        print(f"  theta = {th:6.0f} ({100*f:.0f}% of nominal): "
              f"0.1% {v[0]:.4f}, 0.5% {v[1]:.4f}, 1% {v[2]:.4f}")
    print("\n  and with the load record PERFECT in both channels:")
    for f in CAPS:
        th = TH0 * f
        v = over(th, 1e-6, 1e-6)
        print(f"    {100*f:5.0f}% of nominal: {v:9.4f}"
              + ("   * meets the gate" if v < GATE else ""))

    print("\n" + bar)
    print("R1. RE-CHECK: the nominal row reproduces tab_D_loadprec")
    print(bar)
    for a, want in ((0.01, 0.1590), (0.05, 0.1699), (0.10, 0.1991),
                    (0.20, 0.2876)):
        got = grid.get((1.00, a)) or over(TH0, a, a)
        print(f"  {100*a:4.0f}%: {got:.4f} against the committed {want:.4f} "
              f"({abs(got-want)/want:.2%} apart)")
    print("\nR2. RE-CHECK: the 5% column reproduces the scope table")
    for f, want in ((0.10, 1.5830), (0.05, 3.0673)):
        got = grid[(f, 0.05)]
        print(f"  {100*f:3.0f}% of nominal: {got:.4f} against the committed "
              f"{want:.4f} ({abs(got-want)/want:.2%} apart)")


    #  the numbers the paper cites, written where the build can see them
    loadcap_numbers()

if __name__ == "__main__":
    main()
