"""
Phase 9, closing 3(a): solve the sigma budget by FITTING, not by the Fisher bound.

3a.9 confirmed there is no needle: at every capacity the fitted robust spread at
sigma* stays within about 50% of the Fisher reading. That settles the question it
was asked. It also shows the budget is OPTIMISTIC: a reader takes sigma* as the
noise at which the design delivers 0.20, and in the worst cell it delivers 0.30.

So this script does not report a Fisher solution with a caveat attached. It
solves the thing the reader wants: the largest sigma at which the FITTED robust
spread reaches 0.20. That is the house rule taken one step on, from confirming a
Fisher figure that decides a verdict to solving for one that IS the verdict.

THE FISHER VALUE IS KEPT BESIDE IT, because the ratio between them is the
quantity a reader of the earlier tables needs in order to correct them.

METHOD. The robust spread is not guaranteed monotone in sigma and a bisection on
an assumed bracket can return the bracket's own endpoint, which an earlier
version of this solver did on three cells before the check caught it. So the
crossing is found by SCANNING a log grid and bracketing the last crossing, and a
cell whose gate is never met, or met everywhere on the scan, is reported as such
rather than given a number.

Nothing here records an expected result in advance.

Run:  python audit/q12_fitted_budget.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import q6_sigma_budget as S6

GATE = S6.GATE
NUISANCE = "centre"          #: the per-patient measurement, the conservative one
KE_CASES = (0.8, 0.054)
SCAN = np.geomspace(0.02, 2.0, 25)


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


def robust_at(theta, times, sigma, s_b, ke):
    return S6.mc_at_sigma(theta, times, sigma, s_b, ke_sd=ke,
                          nuisance=NUISANCE)["robust"]


def fitted_star(theta, times, s_b, ke, fisher):
    """
    The largest sigma whose fitted robust spread meets the gate.

    Returns (value, note, grid, values). `value` is None when the scan cannot
    bracket a crossing, and the note says which way it failed.
    """
    grid = fisher * SCAN
    vals = np.array([robust_at(theta, times, g, s_b, ke) for g in grid])
    ok = vals <= GATE
    if not ok.any():
        return None, "gate never met on the scan", grid, vals
    if ok.all():
        return None, "gate met even at 2x Fisher", grid, vals
    j = int(np.nonzero(ok)[0].max())
    if j == len(grid) - 1:
        return float(grid[j]), "at the scan's top", grid, vals
    lo, hi = float(grid[j]), float(grid[j + 1])
    for _ in range(7):
        mid = float(np.sqrt(lo * hi))
        if robust_at(theta, times, mid, s_b, ke) <= GATE:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi)), "", grid, vals


def main():
    print("=" * 104)
    print("  THE sigma BUDGET SOLVED BY FITTING")
    print("=" * 104)
    print(f"\n  gate {GATE}, nuisance convention `{NUISANCE}`, 96 records per "
          f"evaluation, scan {SCAN[0]:.2f}x to {SCAN[-1]:.2f}x the Fisher value")
    banner("THE FITTED BUDGET AGAINST THE FISHER ONE")
    print(f"    {'capacity':18s} {'schedule':16s} {'k_e':>6s} {'Fisher':>9s} "
          f"{'fitted':>9s} {'fitted/Fisher':>14s}  note")
    out = {}
    for cname, theta in S6.CAPS.items():
        L = S6.anchor(theta)[0]
        s_b = 0.10 * L
        for name, times in S6.SCHED.items():
            for ke in KE_CASES:
                f = S6.sigma_star(theta, times, s_b, ke_prior=ke)
                v, note, grid, vals = fitted_star(theta, times, s_b, ke, f)
                out[(cname, name, ke)] = (f, v, note)
                if v is None:
                    print(f"    {cname:18s} {name:16s} {ke:6.3f} {f:9.2f} "
                          f"{'none':>9s} {'n/a':>14s}  {note}", flush=True)
                    k = max(len(grid) // 6, 1)
                    print(f"        robust at sigma/F "
                          f"{np.round(grid[::k] / f, 3)} = "
                          f"{np.round(vals[::k], 3)}", flush=True)
                else:
                    print(f"    {cname:18s} {name:16s} {ke:6.3f} {f:9.2f} "
                          f"{v:9.2f} {v / f:14.3f}  {note}", flush=True)
    banner("WHAT THE CORRECTION IS, PER k_e PRIOR")
    for ke in KE_CASES:
        r = [v / f for (f, v, _) in
             (out[k] for k in out if k[2] == ke) if v is not None]
        n_none = sum(1 for k in out if k[2] == ke and out[k][1] is None)
        if r:
            print(f"    k_e {ke:5.3f}: fitted/Fisher runs {min(r):.3f} to "
                  f"{max(r):.3f} over {len(r)} cells, {n_none} with no "
                  f"solution")
        else:
            print(f"    k_e {ke:5.3f}: no cell has a solution")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
