"""
Phase 7 follow-up, X22: what assuming K_m = 100 costs.

`src/common/params.py:19-20` states the rule this repository is built on:

    An `assumed` row is not a defect on its own. An `assumed` row with no bound
    on what it costs is.

`audit/PHASE7_TABLE2_PROVENANCE.md` X22 records that six rows carry no cost
bound and that `K_m` is the expensive one, because the paper's central negative
result is stated at `m_bar/K_m = 0.1` and the excitation measure depends on
`m/K_m`, not on `m/m_bar`. The caption's protection, that every design quantity
is reported relative to the patient's own operating level, does not reach it.
This script supplies the missing bound.

WHAT IS HELD FIXED, and why. The patient stays managed at `m_bar = 10`. So the
recorded load is re-solved at every `K_m`,

    L(K_m) = theta m_bar/(K_m + m_bar) + k_e m_bar,

rather than held at its `K_m = 100` value. Holding `L` instead would move the
operating point as well as the curvature and would confound the two. The load
prior `S_B = LOAD_FRAC * L` therefore moves with it, and so does `sigma`, whose
assay term carries `b = fit_ab(p, m_bar)[1]` and is `K_m`-dependent. All three
are recomputed per `K_m`; holding any of them would understate the cost.

Everything else is the paper's own footing: `theta = 3000`, `k_e = 5`,
`D = 6000`, `T = 14 d`, `lambda = ln2/3.5`, the monthly record of
`np.arange(0, 730, 30.4)`, the `ln k_e` log-prior at `0.8`, the load-bias prior
at `LOAD_FRAC`, and the information matrix in `(ln theta, ln K_m, b, ln k_e)`.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from common.models import ParamsA
from design.run import delta, DOSE_INTERVAL, LOAD_FRAC, ASSAY_CV
import model_A.run as LOC

TH, KE, M_BAR = 3000.0, 5.0, 10.0
KM0 = 100.0
GATE = 0.20
MONTHLY = np.arange(0.0, 730.0, 30.4)


def footing(Km):
    """Load, ordinate sd and load-bias prior at this K_m, patient held at m_bar."""
    p = ParamsA(theta=TH, Km=Km, ke=KE)
    L = TH * M_BAR / (Km + M_BAR) + KE * M_BAR
    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV, p=p, m_bar=M_BAR)
    return L, float(sigma), LOAD_FRAC * L


def levels(times, loads, Km):
    """The quasi-steady marker, positive root of the balance."""
    inc = np.asarray(delta(times), float)
    B = KE * Km + TH + inc - loads
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * loads * Km)) / (2.0 * KE)


def columns(m, inc, Km):
    """d(ordinate)/d(ln theta, ln K_m, b, ln k_e), g_design_common's form."""
    g = m / (Km + m)
    return np.stack([TH * g,
                     -(TH + inc) * Km * m / (Km + m) ** 2,
                     np.ones_like(m),
                     KE * m], axis=1)


def sd_ln_theta(Km, load_cv=0.10, seed=0, ke_prior=0.8, cols=columns):
    L, sigma, S_B = footing(Km)
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), MONTHLY.size))
    else:
        loads = np.full(MONTHLY.size, L)
    m = levels(MONTHLY, loads, Km)
    inc = np.asarray(delta(MONTHLY), float)
    Q = cols(m, inc, Km)
    F = (np.diag([0.0, 0.0, 1.0 / S_B ** 2, 1.0 / ke_prior ** 2])
         + Q.T @ Q / sigma ** 2)
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over_seeds(Km, load_cv=0.10, seeds=range(24), **kw):
    v = np.array([sd_ln_theta(Km, load_cv, s, **kw) for s in seeds], float)
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


# ---------------------------------------------------------------- re-check
def columns_numeric(m, inc, Km, h=1e-4):
    """
    The same Jacobian by CENTRAL DIFFERENCES of the ordinate itself.

    The ordinate is L + b = (theta + Delta) m/(K_m+m) + k_e m evaluated at the
    OBSERVED levels, so its derivatives in the log coordinates are what the
    filter needs. This shares no algebra with `columns`, which was differentiated
    by hand, so a hand-differentiation slip cannot pass both.
    """
    def ordinate(lt, lk, b, lke):
        t, K, k = np.exp(lt), np.exp(lk), np.exp(lke)
        return (t + inc) * m / (K + m) + k * m + b

    base = np.array([np.log(TH), np.log(Km), 0.0, np.log(KE)])
    out = []
    for j in range(4):
        h0 = h * max(abs(base[j]), 1.0)
        d = []
        for hh in (h0, h0 / 2.0):
            pp, pm = base.copy(), base.copy()
            pp[j] += hh
            pm[j] -= hh
            d.append((ordinate(*pp) - ordinate(*pm)) / (2 * hh))
        # Richardson: (4 D(h/2) - D(h))/3 is O(h^4). A plain central difference
        # bottoms out near 1e-7 relative here, because the ordinate is of order
        # 300 and round-off enters as eps*|f|/h. That is the difference
        # quotient's floor, not a disagreement, and the first version of this
        # check mistook one for the other. Extrapolating raises the reference
        # rather than loosening the threshold.
        out.append((4.0 * d[1] - d[0]) / 3.0)
    return np.stack(out, axis=1)


def km_numbers():
    """
    The K_m sweep as generated numbers.

    Section 2.5 and Section 5.3 both carry this sweep, each with four
    hand-typed literals, and the registry's K_m row carried a fifth set. One
    computation, one set of keys.
    """
    from common.numbers import write_numbers
    from scipy.optimize import brentq
    GATE = 0.20
    lo_km, hi_km = 5.0, 1000.0
    sd_lo = over_seeds(lo_km)[0]
    sd_hi = over_seeds(hi_km)[0]
    cross = brentq(lambda k: over_seeds(k)[0] - GATE, 100.0, 400.0, xtol=1e-3)
    N = {"km_sweep_lo": sd_lo, "km_sweep_hi": sd_hi,
         "km_sweep_fold": sd_hi / sd_lo,
         "km_crossing": cross, "km_headroom": cross / 100.0,
         "km_invitro_sd": over_seeds(290.0)[0]}
    fmt = {"km_sweep_lo": "{:.4f}", "km_sweep_hi": "{:.4f}",
           "km_sweep_fold": "{:.1f}", "km_crossing": "{:.0f}",
           "km_headroom": "{:.2f}", "km_invitro_sd": "{:.4f}"}
    write_numbers("kmcost", N, fmt)
    return N


def recheck():
    print("\n  RE-CHECK: the analytic columns against central differences of the")
    print("  ordinate, at every K_m in the sweep. Same filter, independent Jacobian.")
    worst_j, worst_sd, floor = 0.0, 0.0, 0.0
    h = 1e-4
    for Km in KM_GRID:
        L, sigma, S_B = footing(Km)
        m = levels(MONTHLY, np.full(MONTHLY.size, L), Km)
        inc = np.asarray(delta(MONTHLY), float)
        A, N = columns(m, inc, Km), columns_numeric(m, inc, Km)
        sc = np.maximum(np.abs(A).max(axis=0), 1e-300)
        worst_j = max(worst_j, float(np.max(np.abs(A - N) / sc)))
        # the difference quotient's own round-off floor, eps |f| / h, relative
        # to the same column scale. The threshold is DERIVED from this rather
        # than chosen, so the check cannot be passed by loosening it.
        ordin = np.abs((TH + inc) * m / (Km + m) + KE * m)
        floor = max(floor, float(np.max(np.finfo(float).eps * ordin.max()
                                        / h / sc)))
        a = sd_ln_theta(Km, 0.0, cols=columns)
        n = sd_ln_theta(Km, 0.0, cols=columns_numeric)
        worst_sd = max(worst_sd, abs(a / n - 1.0))
    bound = 20.0 * floor
    ok = worst_j < bound and worst_sd < bound
    print(f"      round-off floor of the quotient  {floor:.3e}")
    print(f"      threshold, 20x that floor        {bound:.3e}")
    print(f"      worst column disagreement        {worst_j:.3e}")
    print(f"      worst sd(ln theta) disagreement  {worst_sd:.3e}")
    print(f"      [{'PASS' if ok else 'FAIL'}]")
    return ok


KM_GRID = np.array([5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 290.0, 500.0, 1000.0])


def main():
    km_numbers()
    print("=" * 78)
    print("  X22. WHAT ASSUMING K_m = 100 COSTS")
    print("=" * 78)
    print(__doc__.split("Nothing here")[0].strip()[-600:])
    print()
    print(f"  Monthly record, n = {MONTHLY.size}, load CV 10%, median over 24 seeds")
    print(f"  with the [min, max] beside it. Target is sd(ln theta) <= {GATE}.\n")
    print(f"      {'K_m':>7} {'m_bar/K_m':>10} {'L':>9} {'sigma':>8} "
          f"{'sd(ln theta)':>13} {'[min, max]':>20} {'vs K_m=100':>11}")
    print("  " + "-" * 84)
    base = over_seeds(KM0)[0]
    rows = []
    for Km in KM_GRID:
        L, sigma, _ = footing(Km)
        md, lo, hi = over_seeds(Km)
        rows.append((Km, md))
        mark = " <-- assumed" if Km == KM0 else ""
        print(f"      {Km:7.0f} {M_BAR/Km:10.3f} {L:9.2f} {sigma:8.3f} "
              f"{md:13.4f} [{lo:8.4f},{hi:8.4f}] {md/base:10.2f}x{mark}")

    fin = [(k, v) for k, v in rows if np.isfinite(v)]
    worst = max(fin, key=lambda r: r[1])
    best = min(fin, key=lambda r: r[1])
    pass_gate = [k for k, v in fin if v <= GATE]
    print(f"\n  Over a {KM_GRID.max()/KM_GRID.min():.0f}-fold range of K_m, "
          f"sd(ln theta) runs {best[1]:.4f} to {worst[1]:.4f}, a factor of "
          f"{worst[1]/best[1]:.2f}.")
    if pass_gate:
        print(f"  The {GATE} target is met for K_m from {min(pass_gate):.0f} to "
              f"{max(pass_gate):.0f}, that is m_bar/K_m from "
              f"{M_BAR/max(pass_gate):.3f} to {M_BAR/min(pass_gate):.2f}.")
    fails = [k for k, v in fin if v > GATE]
    if fails:
        print(f"  It FAILS at K_m = {', '.join(f'{k:.0f}' for k in fails)}.")
    print()
    recheck()


if __name__ == "__main__":
    main()
