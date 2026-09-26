"""
Section 8's load-record decomposition, on the paper's OWN convention.

WHY THIS EXISTS. The provenance sweep reported these figures as reproducing
"within a few per cent, consistent with a seed convention". By this project's
standards that is not a verdict, it is the shape of the disagreements that have
already cost three rounds: sample against population standard deviation, the
t = 0 convention, and the two nuisance conventions. So the convention is
identified here rather than described.

IT IS THE t = 0 CONVENTION, AGAIN. `n4_theta_scope.MONTHLY` is
`arange(0.0, 730.0, 30.4)`: 25 draws STARTING AT t = 0, one of them before the
first infusion. The sweep's check used `arange(30.4, 30.4 * 26, 30.4)`, which
drops the t = 0 draw and runs 30 days past the two-year window. On the paper's
grid `n4.over_seeds(theta=3000)` returns 0.1699, which is the paper's 0.170.

This file mirrors `n4.sd_ln_theta` exactly, differing in one respect only: the
load record's accuracy enters `sigma` and the offset prior through SEPARATE
arguments, because Section 8's whole point is to separate them and n4's own
function ties them to one `load_frac`.

Run:  python audit/q18_section8.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import n4_theta_scope as N4

KM, KE, LOC = N4.KM0, N4.KE, N4.LOC
TH_NOM = 3000.0
SEEDS = N4.SEEDS


def sd(theta, seed, sig_frac, b_frac, times=None, load_cv=0.10,
       ke_prior=0.8, D=None):
    """`n4.sd_ln_theta`, with sigma's load term and the offset prior split."""
    times = N4.MONTHLY if times is None else times
    D = N4.D0 if D is None else D
    L, m_bar = N4.anchor(theta, KM, KE, "A")
    p = N4.ParamsA(theta=theta, Km=KM, ke=KE)
    sigma = float(LOC.sigma_total(sig_frac, N4.ASSAY_CV, p=p, m_bar=m_bar))
    S_B = b_frac * L
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        loads = np.full(times.size, L)
    inc = N4.delta_at(times, D)
    B = KE * KM + theta + inc - loads
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * loads * KM)) / (2.0 * KE)
    g = m / (KM + m)
    Q = np.stack([theta * g, -(theta + inc) * KM * m / (KM + m) ** 2,
                  np.ones_like(m), KE * m], axis=1)
    F = (np.diag([0.0, 0.0, 1.0 / S_B ** 2, 1.0 / ke_prior ** 2])
         + Q.T @ Q / sigma ** 2)
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def cell(theta, sig, b, **kw):
    v = np.array([sd(theta, s, sig, b, **kw) for s in SEEDS])
    v = v[np.isfinite(v)]
    return float(np.median(v))


def banner(t):
    print("\n" + "=" * 96)
    print("  " + t)
    print("=" * 96)


def main():
    print("=" * 96)
    print("  SECTION 8's LOAD DECOMPOSITION, ON THE PAPER'S OWN CONVENTION")
    print("=" * 96)
    banner("0. THE CONVENTION, PINNED BEFORE ANYTHING IS COMPARED")
    ref = N4.over_seeds(theta=TH_NOM)[0]
    mine = cell(TH_NOM, N4.LOAD_FRAC, N4.LOAD_FRAC)
    print(f"    n4.MONTHLY: {N4.MONTHLY.size} draws, first {N4.MONTHLY[0]:.1f} d,"
          f" last {N4.MONTHLY[-1]:.1f} d, so a t = 0 draw IS included")
    print(f"    n4.over_seeds(theta=3000)  {ref:.4f}   the paper's 0.170")
    print(f"    this file, same settings   {mine:.4f}   "
          f"{'IDENTICAL' if abs(mine - ref) < 5e-5 else 'DIFFERS'}")
    print("    LOAD_FRAC = %.2f, ASSAY_CV = %.2f, load_cv = 0.10, seeds = 24"
          % (N4.LOAD_FRAC, N4.ASSAY_CV))

    banner("1. THE 20% ROW, SEPARATED (the paper's own numbers beside)")
    base = cell(TH_NOM, 0.05, 0.05)
    both = cell(TH_NOM, 0.20, 0.20)
    sc = cell(TH_NOM, 0.20, 0.05)
    off = cell(TH_NOM, 0.05, 0.20)
    tot = both - base
    rows = [("baseline, both channels at 5%", base, 0.170),
            ("both channels at 20%", both, None),
            ("per-visit scatter alone", sc, 0.283),
            ("systematic offset alone", off, 0.172)]
    print(f"    {'quantity':34s} {'paper':>8s} {'here':>8s} {'diff':>8s}")
    for lab, v, paper in rows:
        d = "" if paper is None else f"{100 * (v / paper - 1):+7.1f}%"
        p = "" if paper is None else f"{paper:8.3f}"
        print(f"    {lab:34s} {p:>8s} {v:8.3f} {d:>8s}")
    print(f"\n    scatter share {100 * (sc - base) / tot:.1f}%, "
          f"offset share {100 * (off - base) / tot:.1f}%   "
          f"paper says 98% and 2%")

    banner("2. THE OFFSET's SHARE BY CAPACITY (paper: 1.8 / 9.9 / 26.2)")
    print(f"    {'capacity':>12s} {'paper':>8s} {'here':>8s}")
    for cap, lab, paper in ((1.0, "nominal", 1.8), (0.5, "half", 9.9),
                            (0.2, "a fifth", 26.2), (0.1, "a tenth", None)):
        th = TH_NOM * cap
        b0, bt = cell(th, 0.05, 0.05), cell(th, 0.20, 0.20)
        of = cell(th, 0.05, 0.20)
        v = 100 * (of - b0) / (bt - b0)
        print(f"    {lab:>12s} {('%.1f' % paper) if paper else '':>8s} "
              f"{v:8.1f}")

    banner("3. THE THRESHOLD THE PAPER STATES AS A REQUIREMENT: 10.4%")
    print("""    "holding the offset at 5%, the scatter must be tighter than 10.4%".
    Bisected on the same convention.""")
    lo, hi = 0.01, 0.60
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        if cell(TH_NOM, mid, 0.05) <= 0.20:
            lo = mid
        else:
            hi = mid
    print(f"\n    paper 10.4%   here {100 * lo:.2f}%   "
          f"difference {100 * (lo / 0.104 - 1):+.1f}%")
    print("""
    And the companion claim, that with the scatter held at 5% the target is met
    at every offset from 1% to 500%:""")
    bad = [f"{100 * b:g}%" for b in (0.01, 0.05, 0.20, 1.0, 2.0, 5.0)
           if cell(TH_NOM, 0.05, b) > 0.20]
    print(f"    offsets failing the 0.20 target: "
          f"{', '.join(bad) if bad else 'none, the claim holds'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
