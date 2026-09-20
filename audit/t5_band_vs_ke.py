"""
3. THE WORKING BAND AGAINST THE TRUE k_e, NOT AGAINST HOW WELL IT IS KNOWN.

Scoping only. Nothing here reaches the paper.

Section 6.2 reports a working band, theta in [2472, 18814], and Section 6.2's
own mechanism says why: the information the record carries about capacity is
the enzymatic flux,

    d s / d ln theta  =  theta g(m)  =  L - k_e m,

so the enzymatic share is 1 - k_e m / L. That share, and therefore the band's
LOCATION, depends on the VALUE of k_e and not only on how well it is known.
Every existing number prices the second question: the factor of 2.67 the gate
needs, and the 0.388 that an exactly known k_e still leaves at 5% of nominal.
The first question has not been asked, and it is the one on which the paper's
most important result rests, because k_e is the one constant the paper states
that nothing in the literature supplies.

The footing is Section 6.2's own, convention A: the recorded load L is held at
its theta = 3000, k_e = 5 value while the marker rises, so the comparison is
across capacities at a fixed diet. The load is NOT re-solved at each k_e
either, for the same reason: re-solving it would move the operating point and
the enzymatic share together and confound the two.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

import model_A.run as LOC
from common.models import ParamsA
from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, SEEDS, delta, load, median_spread)

KM, TH0, KE0 = 100.0, 3000.0, 5.0
LB = load(0.0)
KE_PRIOR_LOG = 0.8
GATE = 0.20
MONTHLY = np.arange(0.0, 730.0, 30.4)
INC = np.asarray(delta(MONTHLY), float)


def level(theta, ke, inc, L):
    B = ke * KM + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * L * KM)) / (2.0 * ke)


def share(theta, ke, L=LB):
    """The enzymatic share of the load at the untreated balance, 1 - k_e m/L."""
    m = level(theta, ke, 0.0, L)
    return 1.0 - ke * m / L


def sd_of(theta, ke, seed):
    L = LB * (1.0 + LOAD_CV *
              np.random.default_rng(seed).standard_normal(MONTHLY.size))
    m = level(theta, ke, INC, L)
    sig = LOC.sigma_total(LOAD_FRAC, ASSAY_CV,
                          p=ParamsA(theta=theta, Km=KM, ke=ke),
                          m_bar=level(theta, ke, 0.0, LB))
    F = np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * LB) ** 2,
                 1.0 / KE_PRIOR_LOG ** 2])
    for mi, di in zip(m, INC):
        gm = mi / (KM + mi)
        q = np.array([theta * gm, -(theta + di) * KM * mi / (KM + mi) ** 2,
                      1.0, ke * mi])
        F = F + np.outer(q, q) / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over(theta, ke):
    return median_spread([sd_of(theta, ke, s) for s in SEEDS])[0]


def band(ke, lo=100.0, hi=1e6):
    """
    The interval of theta on which the median meets the gate.

    sd is not monotone in theta, so the gate is crossed twice and each edge is
    bracketed separately around the interior minimum.
    """
    grid = np.geomspace(lo, hi, 90)
    vals = np.array([over(t, ke) for t in grid])
    if not np.any(vals < GATE):
        return None, None, float(np.min(vals)), float(grid[np.argmin(vals)])
    i = int(np.argmin(vals))
    f = lambda t: over(t, ke) - GATE
    j = int(np.where(vals[:i] >= GATE)[0][-1]) if np.any(vals[:i] >= GATE) else None
    k = (i + int(np.where(vals[i:] >= GATE)[0][0])
         if np.any(vals[i:] >= GATE) else None)
    e_lo = brentq(f, grid[j], grid[j + 1], xtol=1.0) if j is not None else None
    e_hi = brentq(f, grid[k - 1], grid[k], xtol=1.0) if k is not None else None
    return e_lo, e_hi, float(vals[i]), float(grid[i])


def main():
    print(__doc__)
    bar = "=" * 78
    KES = (1.25, 2.5, 5.0, 10.0, 20.0)

    print(bar)
    print("1. THE ENZYMATIC SHARE, WHICH IS THE MECHANISM, AGAINST k_e")
    print(bar)
    print(f"  the share is 1 - k_e m/L at the untreated balance, L held at "
          f"{LB:.1f}\n")
    print(f"  {'k_e':>8} {'x nominal':>10} " +
          "".join(f"{lab:>12}" for lab in
                  ("100%", "20%", "10%", "5%", "2%")))
    for ke in KES:
        row = [f"  {ke:8.2f} {ke / KE0:10.2f}"]
        for frac in (1.0, 0.2, 0.1, 0.05, 0.02):
            row.append(f"{100 * share(TH0 * frac, ke):11.1f}%")
        print("".join(row))

    print("\n" + bar)
    print("2. THE WORKING BAND, RE-SOLVED AT EACH k_e")
    print(bar)
    print(f"  medians over {len(SEEDS)} intake records on the recommended "
          f"monthly design\n")
    print(f"  {'k_e':>8} {'x nominal':>10} {'band lo':>10} {'band hi':>10} "
          f"{'as % of 3000':>22} {'best sd':>9} {'fold':>7}")
    out = {}
    for ke in KES:
        lo, hi, best, at = band(ke)
        out[ke] = (lo, hi, best, at)
        if lo is None and hi is None:
            print(f"  {ke:8.2f} {ke / KE0:10.2f} {'---':>10} {'---':>10} "
                  f"{'the gate is never met':>22} {best:9.4f} {'---':>7}")
        else:
            pct = (f"{100 * lo / TH0:.0f}%-{100 * hi / TH0:.0f}%"
                   if hi else f"{100 * lo / TH0:.0f}%-open")
            fold = f"{hi / lo:.1f}" if hi else "open"
            print(f"  {ke:8.2f} {ke / KE0:10.2f} {lo:10.0f} "
                  f"{(hi if hi else float('nan')):10.0f} {pct:>22} "
                  f"{best:9.4f} {fold:>7}")

    print("\n" + bar)
    print("3. WHAT IT WOULD TAKE FOR THE DISEASE TO SIT INSIDE THE BAND")
    print(bar)
    print("  the disease is theta <= 150, which is 5% of nominal\n")
    for ke in KES:
        s = over(150.0, ke)
        lo = out[ke][0]
        gap = (150.0 / lo) if lo else float("nan")
        print(f"  k_e = {ke:6.2f}: sd at theta = 150 is {s:9.4f}, and the "
              f"lower edge is a factor of {1 / gap:6.1f} above it"
              if lo else
              f"  k_e = {ke:6.2f}: sd at theta = 150 is {s:9.4f}, and there "
              f"is no band at all")
    print()
    print("  and the k_e at which theta = 150 would itself meet the gate:")
    try:
        k_star = brentq(lambda k: over(150.0, k) - GATE, 1e-4, 5.0, xtol=1e-4)
        print(f"    k_e = {k_star:.4f}, which is {KE0 / k_star:.0f} times "
              f"smaller than the value the paper carries")
    except ValueError:
        v_lo, v_hi = over(150.0, 1e-4), over(150.0, 5.0)
        print(f"    none in [1e-4, 5]: the sd runs {v_lo:.4f} to {v_hi:.4f} "
              f"over that range, so no value of k_e opens the disease range")

    print("\n" + bar)
    print("R1. RE-CHECK: at k_e = 5 this reproduces the published band")
    print(bar)
    lo, hi, best, at = out[5.0]
    print(f"  band [{lo:.0f}, {hi:.0f}] against the published [2472, 18814]")
    print(f"  best sd {best:.4f} at theta = {at:.0f} against the published "
          f"0.123 at about twice nominal")
    print(f"  share at 5% of nominal {100 * share(150.0, 5.0):.1f}% against "
          f"the published 16.3%")


if __name__ == "__main__":
    main()
