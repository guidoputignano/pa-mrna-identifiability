"""
11. MEASURING THE LOAD RATHER THAN RECORDING IT.

Scoping only. Nothing here reaches the paper.

The paper assumes throughout that the substrate load is RECORDED: a diary,
with a systematic offset b carried as a free parameter and a per-visit scatter
setting sigma. It never considers that the load can be MEASURED. Whole-body
propionate turnover is measurable by isotope dilution, and the paper already
cites published figures for it (250-400 umol/kg/day, `propionate_sources`), so
this is an existing technique rather than a proposal.

THREE THINGS CHANGE AT ONCE, and they are not independent.

  (1) b disappears. It is not tightened, it is REMOVED: there is no unrecorded
      component to carry, so the parameter vector drops from four to three and
      the degeneracy of eq:fibrecl does not exist. That fibre is defined by
      trading b against theta and K_m; with no b there is nothing to trade.
  (2) sigma stops being a diet record's accuracy and becomes a tracer
      measurement's precision, which is a different quantity with a different
      achievable value.
  (3) The steady record closes on itself: with L known, A = K_m L gives K_m,
      B gives theta and C gives k_e, in the three-term form of Section 3.

WHAT IS HELD FIXED so the comparison is about the load and nothing else: the
same monthly schedule, the same dose train, the same k_e log-prior of 0.8, the
same assay CV. The k_e prior is KEPT even here, because nothing in a load
measurement supplies k_e and dropping it would credit the tracer with
information it does not carry.

AND THE TRADE IS PRICED, because a tracer study is not something a clinic runs
twenty-five times. Part 3 asks how few measured-load visits match the monthly
diary design, at nominal and in the disease range, and whether the ranking
reverses.

Nothing here records an expected result in advance.

Run:  python audit/t9_measured_load.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import model_A.run as LOC
from common.models import ParamsA
from design.run import (ASSAY_CV, LOAD_CV, LOAD_FRAC, SEEDS, delta, load,
                        median_spread)

KM, KE, TH0 = 100.0, 5.0, 3000.0
LB = load(0.0)
KE_PRIOR_LOG = 0.8
GATE = 0.20
MONTH = 30.4
CAPS = (1.00, 0.50, 0.20, 0.10, 0.05, 0.02)


def level(theta, inc, L):
    B = KE * KM + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * L * KM)) / (2.0 * KE)


def sd_of(theta, seed, n_visits=25, measured=False, prec=None,
          scatter=LOAD_FRAC, ke_prior=KE_PRIOR_LOG):
    """
    sd(ln theta) on n monthly visits.

    `measured=False` is the paper's own footing: the load is recorded, b is a
    free parameter with a prior of `scatter` * L, and sigma carries the same
    fraction as the ordinate noise.

    `measured=True` removes b from the parameter vector entirely -- not a
    tighter prior on it, no column at all -- and sets sigma from `prec`, the
    tracer's own precision. That is the whole of what measuring the load buys,
    and it is why the two cannot be compared by tightening a prior.
    """
    t = np.arange(0.0, n_visits * MONTH, MONTH)[:int(n_visits)]
    inc = np.asarray(delta(t), float)
    L = LB * (1.0 + LOAD_CV *
              np.random.default_rng(seed).standard_normal(t.size))
    m = level(theta, inc, L)
    frac = prec if measured else scatter
    sig = LOC.sigma_total(frac, ASSAY_CV,
                          p=ParamsA(theta=theta, Km=KM, ke=KE),
                          m_bar=level(theta, 0.0, LB))
    if measured:
        pri = [np.inf, np.inf, ke_prior ** 2]
    else:
        pri = [np.inf, np.inf, (scatter * LB) ** 2, ke_prior ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    for mi, di in zip(m, inc):
        gm = mi / (KM + mi)
        q = [theta * gm, -(theta + di) * KM * mi / (KM + mi) ** 2]
        if not measured:
            q.append(1.0)
        q.append(KE * mi)
        q = np.asarray(q)
        F = F + np.outer(q, q) / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over(theta, **kw):
    return median_spread([sd_of(theta, s, **kw) for s in SEEDS])


def visits_to_gate(theta, nmax=600, **kw):
    for n in range(3, nmax + 1):
        if over(theta, n_visits=n, **kw)[0] < GATE:
            return n
    return None


def measured_numbers():
    """
    Route 6's figures, generated rather than typed.

    WHY THIS FUNCTION NOW EXISTS. This script's docstring said "Scoping only.
    Nothing here reaches the paper." It owns the whole of Section 11's route 6
    and a sentence in Section 6, stated as prose literals that nothing in the
    build compared with anything (REVIEW_LOG.md R81, R82).
    """
    from common.numbers import write_numbers
    TH = 3000.0
    tgt = over(TH)[0]
    match = None
    for n in range(3, 61):
        if over(TH, n_visits=n, measured=True, prec=0.02)[0] <= tgt:
            match = n
            break
    th2 = 3000.0 * 0.02
    tgt2 = over(th2)[0]
    match2 = None
    for n in range(3, 61):
        if over(th2, n_visits=n, measured=True, prec=0.02)[0] <= tgt2:
            match2 = n
            break
    pinned = over(TH)[0]
    gone = over(TH, measured=True, prec=LOAD_FRAC)[0]
    N = {"meas_match_nominal": float(match),
         "meas_gate_nominal": float(visits_to_gate(TH, measured=True,
                                                   prec=0.02, nmax=120)),
         "meas_match_2pct": float(match2),
         "meas_exact_5pct": over(150.0, measured=True, prec=1e-6)[0],
         "meas_b_only_gain_pct": 100.0 * (1.0 - gone / pinned)}
    fmt = {"meas_match_nominal": "{:.0f}", "meas_gate_nominal": "{:.0f}",
           "meas_match_2pct": "{:.0f}", "meas_exact_5pct": "{:.3f}",
           "meas_b_only_gain_pct": "{:.2f}"}
    write_numbers("measured", N, fmt)
    return N


def main():
    print(__doc__)
    bar = "=" * 86

    print(bar)
    print("1. THE RECOMMENDED RECORD, RECORDED AGAINST MEASURED")
    print(bar)
    print(f"  twenty-five monthly visits, medians over {len(SEEDS)} records\n")
    print(f"  {'theta':>13} {'recorded, 5%':>14} " +
          "".join(f"{f'measured {100*p:.0f}%':>16}" for p in (0.02, 0.01)))
    rec = {}
    for f in CAPS:
        th = TH0 * f
        r = over(th)[0]
        rec[f] = r
        row = [f"  {f'{100*f:.0f}% of nom':>13} {r:13.4f}"
               + ("*" if r < GATE else " ")]
        for p in (0.02, 0.01):
            v = over(th, measured=True, prec=p)[0]
            row.append(f"{v:15.4f}" + ("*" if v < GATE else " "))
        print("".join(row))
    print("\n  * meets the 0.20 gate")

    print("\n" + bar)
    print("2. HOW FAR DOWN DOES IT REACH, AND ON WHICH OF THE TWO CHANGES")
    print(bar)
    print("  b removed and sigma tightened are separable, and separating them")
    print("  says which one is doing the work.\n")
    print(f"  {'theta':>13} {'recorded 5%':>12} {'b free, sig 1%':>15} "
          f"{'b gone, sig 5%':>15} {'b gone, sig 1%':>15}")
    for f in CAPS:
        th = TH0 * f
        a = rec[f]
        b = over(th, scatter=0.01)[0]
        c = over(th, measured=True, prec=LOAD_FRAC)[0]
        d = over(th, measured=True, prec=0.01)[0]
        print(f"  {f'{100*f:.0f}% of nom':>13} {a:12.4f} {b:15.4f} "
              f"{c:15.4f} {d:15.4f}")

    print("\n" + bar)
    print("3. THE TRADE: HOW FEW MEASURED VISITS MATCH THE MONTHLY DIARY")
    print(bar)
    print("  A tracer study is not run twenty-five times. This asks how many")
    print("  measured-load visits reach the same precision the recorded")
    print("  twenty-five deliver, and whether the answer moves with capacity.\n")
    print(f"  {'theta':>13} {'recorded n=25':>14} "
          f"{'measured n to match':>21} {'measured n to gate':>20}")
    for f in CAPS:
        th = TH0 * f
        target = rec[f]
        match = None
        for n in range(3, 61):
            if over(th, n_visits=n, measured=True, prec=0.02)[0] <= target:
                match = n
                break
        g = visits_to_gate(th, measured=True, prec=0.02, nmax=120)
        print(f"  {f'{100*f:.0f}% of nom':>13} {target:14.4f} "
              f"{(str(match) if match else '> 60'):>21} "
              f"{(str(g) if g else '> 120'):>20}")

    print("\n" + bar)
    print("4. WHAT IT WOULD TAKE, IF ANYTHING, IN THE DISEASE RANGE")
    print(bar)
    for f in (0.10, 0.05, 0.02):
        th = TH0 * f
        print(f"  theta = {th:6.0f} ({100*f:.0f}% of nominal)")
        for p, nm in ((0.02, "tracer to 2%"), (0.01, "to 1%"),
                      (0.001, "to 0.1%"), (1e-6, "exact")):
            v = over(th, measured=True, prec=p)
            print(f"     {nm:>14}, n=25: {v[0]:9.4f} "
                  f"[{v[1]:.4f}, {v[2]:.4f}]"
                  + ("   * GATE MET" if v[0] < GATE else ""))
        n = visits_to_gate(th, measured=True, prec=0.01, nmax=600)
        print(f"     at 1% precision the gate needs "
              + (f"{n} monthly visits, {n/12:.1f} years" if n
                 else "more than 600 monthly visits"))
        print()

    print(bar)
    print("5. AND WITHOUT THE k_e PRIOR, WHICH NO LOAD MEASUREMENT SUPPLIES")
    print(bar)
    print("  Measuring the load says nothing about the alternative channel, so")
    print("  the k_e prior is kept above. This is what it is worth here.\n")
    for f in (1.00, 0.05):
        th = TH0 * f
        a = over(th, measured=True, prec=0.02)[0]
        b = over(th, measured=True, prec=0.02, ke_prior=1e9)[0]
        print(f"  {100*f:5.0f}% of nominal: {a:9.4f} with the prior, "
              f"{b:11.4f} without, a factor of {b/a:.1f}")

    print("\n" + bar)
    print("R1. RE-CHECK: the recorded column reproduces the committed numbers")
    print(bar)
    for f, want in ((1.00, 0.1699), (0.10, 1.5830), (0.05, 3.0673)):
        print(f"  {100*f:5.0f}% of nominal: {rec[f]:.4f} against the "
              f"committed {want:.4f} ({abs(rec[f]-want)/want:.2%} apart)")
    print("\nR2. RE-CHECK: removing b is NOT the same as pinning it tight")
    for f in (1.00, 0.05):
        th = TH0 * f
        gone = over(th, measured=True, prec=LOAD_FRAC)[0]
        tiny = over(th, scatter=LOAD_FRAC)[0]
        pinned = median_spread([sd_of(th, s, scatter=LOAD_FRAC) for s in SEEDS])[0]
        print(f"  {100*f:5.0f}% of nominal: b removed {gone:.4f}, "
              f"b free at the same sigma {pinned:.4f}")


    #  the numbers the paper cites, written where the build can see them
    measured_numbers()

if __name__ == "__main__":
    main()
