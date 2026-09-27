"""
The (n, k) trade curve: how many monthly draws a record needs at each number
of breath measurements.

WHY THIS REPLACES audit/u3_joint_fisher.py's TABLE. That table fixes the
record at 25 monthly draws and asks how many breath tests reach the gate. It
therefore reads as "how many breath tests INSTEAD of a design", and the
decision-relevant object is the other one: a protocol chooses both, so what
it needs is the curve of (draws, breath tests) pairs that reach the gate.
Two designs on that curve are exchangeable; a point off it is not a design.

It also puts the k = 1 result where it belongs. One breath reading
contributes exactly nothing at every capacity, because the unknown recovery
factor absorbs it, and on a trade curve that shows up as the k = 1 column
being IDENTICAL to the k = 0 column rather than merely close to it. That is
the cleanest statement of why the recovery factor forces pairs: the first
reading of a pair buys nothing at all, so breath tests are bought two at a
time.

FOOTING. As audit/u3_joint_fisher.py, whose sd_joint this imports so the two
cannot drift: theta with the recorded load held fixed at L = 322.727 umol/L/d
and the operating level, sigma and the load-bias prior re-solved at each
capacity, K_m = 100, k_e = 5/d with a 0.8 log-prior, D = 6000 on a 14 d train
with a 3.5 d half-life, load-bias prior 0.05 L, gate 0.20, constant recorded
load for the headline grid. Breath readings alternate the dose train's trough
and peak and are attached to visits the record already contains, so they cost
an assay and not a draw.

Nothing here records an expected result in advance.

Run:  python audit/u11_trade_curve.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import u3_joint_fisher as U3

LB, KM, KE = U3.LB, U3.KM, U3.KE
GATE = U3.GATE
MONTH = 30.4
YEAR = 365.25
NMAX = 4000


def rows_with_kappa(theta, m, inc, breath=False):
    """
    Design rows in (ln theta, ln K_m, b, ln k_e, ln c, ln kappa).

    THE BREATH BLOCK'S kappa ENTRY IS NOT THE MARKER BLOCK'S, and writing it
    as if it were is the error this function exists to avoid. The marker's
    flux responds to the dose scale DIRECTLY, d(flux)/d ln kappa = Delta g(m),
    a quantity of order 200 in flux units. The breath reading responds only
    THROUGH THE LEVEL, so

        d ln F / d ln kappa = (d ln F / dm)(dm / d ln kappa)
                            = (-k_e/S) * (-Delta g(m)/J) = k_e Delta g /(S J),

    which is of order 1e-3. Appending the marker column to the breath rows
    instead gave two breath tests the leverage of a hundred draws and made
    the Tier 1 record collapse from 185 draws to 11.
    """
    m = np.asarray(m, dtype=float)
    inc = np.asarray(inc, dtype=float)
    g = m / (KM + m)
    if not breath:
        base = U3.marker_rows(theta, m, inc)          # 5 columns
        return np.column_stack([base, inc * g])
    base = U3.breath_rows(theta, m, inc, modelled=True)
    S = theta + KE * (KM + m)
    J = (theta + inc) * KM / (KM + m) ** 2 + KE
    d_lkap = KE * inc * g / (S * J)
    return np.column_stack([base, d_lkap])


def sd_joint_k(theta, m, inc, sigma, mb, incb, breath_cv, kappa_cv=None):
    """u3's joint filter, optionally carrying the increment's scale."""
    if kappa_cv is None:
        return U3.sd_joint(theta, m, inc, sigma, mb, incb, breath_cv)
    Q = rows_with_kappa(theta, m, inc)
    pri = [np.inf, np.inf, (U3.LOAD_FRAC * LB) ** 2, U3.KE_PRIOR ** 2,
           np.inf, kappa_cv ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    if len(mb):
        R = rows_with_kappa(theta, mb, incb, breath=True)
        F = F + R.T @ R / breath_cv ** 2
    else:
        keep = [0, 1, 2, 3, 5]                 # no breath data, no ln c
        F = F[np.ix_(keep, keep)]
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def sd_at(theta, n, k, breath_cv, sigma, kappa_cv=None, load_cv=0.0, seed=0):
    """sd(ln theta) for n monthly draws plus k breath readings."""
    m, inc = U3.record(theta, MONTH * np.arange(n, dtype=float), load_cv, seed)
    incs, _ = U3.breath_design(k)
    mb = U3.level(theta, incs, LB) if k else np.array([])
    return sd_joint_k(theta, m, inc, sigma, mb, incs, breath_cv, kappa_cv)


#: audit/t2_tier1_margin.py's own record convention, reused verbatim so the
#: k = 0 column of the Tier 1 arm reproduces its 118 rather than approximating
#: it. Two things differ from u3's: the dose train is taken at STEADY STATE
#: (an offset of fifty intervals) rather than from the first dose, and the
#: intake is Gaussian rather than lognormal. Together they move the count from
#: 118 to 141, which is a footing difference and not a disagreement.
T_DOSE_T2 = 50.0 * U3.DOSE_INTERVAL
PEAK_T2 = U3.D_AMP / (1.0 - np.exp(-U3.LAMBDA * U3.DOSE_INTERVAL))
LOAD_CV_T2 = 0.10


def levels_t2(theta, n, seed):
    t = MONTH * np.arange(n, dtype=float)
    tau = np.mod(t - T_DOSE_T2, U3.DOSE_INTERVAL)
    inc = PEAK_T2 * np.exp(-U3.LAMBDA * tau)
    Ls = LB * (1.0 + LOAD_CV_T2 *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    return m, inc


def sd_med(theta, n, k, breath_cv, sigma, kappa_cv=None, seeds=range(24)):
    """Median over intake records on t2_tier1_margin's own convention."""
    incs, _ = U3.breath_design(k)
    mb = U3.level(theta, incs, LB) if k else np.array([])
    v = []
    for sd in seeds:
        m, inc = levels_t2(theta, n, sd)
        v.append(sd_joint_k(theta, m, inc, sigma, mb, incs, breath_cv,
                            kappa_cv))
    v = np.array(v, float)
    f = v[np.isfinite(v)]
    return float(np.median(f)) if f.size else np.inf


def smallest_n(theta, k, breath_cv, sigma, kappa_cv=None, nmax=NMAX,
               fn=None):
    """Fewest monthly draws reaching the gate at this k, by doubling+bisect."""
    fn = fn or sd_at

    def ok(n):
        return fn(theta, n, k, breath_cv, sigma, kappa_cv) <= GATE
    if ok(4):
        lo = 4
        while lo > 4 and ok(lo - 1):
            lo -= 1
        return 4
    hi = 8
    while hi <= nmax and not ok(hi):
        hi *= 2
    if hi > nmax:
        return None
    lo = hi // 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ======================================================================
def part_a():
    print("=" * 88)
    print("PART A.  the trade curve: monthly draws needed at each k")
    print("=" * 88)
    print("\n  Each cell is the fewest monthly draws reaching sd(ln theta) =")
    print("  0.20 with that many breath readings. A dash means the gate is")
    print("  not reached at 4000 draws. Years are the record's span.\n")
    KS = (0, 1, 2, 4, 8, 16, 32, 64, 128)
    for cv in (0.05, 0.10, 0.14):
        print(f"  breath CV {cv:.0%} per reading, increment known exactly")
        print(f"    {'theta':>7} " + "".join(f"{f'k={k}':>11}" for k in KS))
        print("    " + "-" * (8 + 11 * len(KS)))
        for theta in (3000.0, 600.0, 300.0):
            sig = U3.sigma_at(theta, float(U3.level(theta, 0.0, LB)))
            cells = []
            for k in KS:
                n = smallest_n(theta, k, cv, sig)
                cells.append("-" if n is None
                             else f"{n} ({n * MONTH / YEAR:.1f}y)")
            print(f"    {theta:7.0f} " + "".join(f"{c:>11}" for c in cells))
        print()


def part_b():
    print("=" * 88)
    print("PART B.  the same on the Tier 1 footing, increment at 14.02%")
    print("=" * 88)
    print("\n  This is the arm that needs 118 draws with no breath tests, so")
    print("  it is the arm where the trade is worth making. Medians over")
    print("  eight intake records at the specification's 10% load CV, which")
    print("  is t2_tier1_margin's footing, so the k = 0 column should")
    print("  reproduce its 118.\n")
    KS = (0, 1, 2, 4, 8, 16, 32, 64, 128)
    for cv in (0.05, 0.10, 0.14):
        print(f"  breath CV {cv:.0%} per reading")
        print(f"    {'theta':>7} " + "".join(f"{f'k={k}':>11}" for k in KS))
        print("    " + "-" * (8 + 11 * len(KS)))
        for theta in (3000.0, 600.0, 300.0):
            sig = U3.sigma_at(theta, float(U3.level(theta, 0.0, LB)))
            cells = []
            for k in KS:
                n = smallest_n(theta, k, cv, sig, kappa_cv=0.1402,
                               fn=sd_med)
                cells.append("-" if n is None
                             else f"{n} ({n * MONTH / YEAR:.1f}y)")
            print(f"    {theta:7.0f} " + "".join(f"{c:>11}" for c in cells))
        print()


def part_c():
    print("=" * 88)
    print("PART C.  k = 1 buys nothing, exactly, and that is why pairs")
    print("=" * 88)
    print("\n  If the recovery factor absorbs a single reading exactly, the")
    print("  k = 0 and k = 1 columns are not merely close: they are the same")
    print("  number to machine precision, at every capacity, every breath CV")
    print("  and every record length. Tested rather than asserted.\n")
    worst = 0.0
    print(f"    {'theta':>7} {'n':>5} {'CV':>5} {'k=0':>12} {'k=1':>12} "
          f"{'difference':>12} {'k=2':>12}")
    print("    " + "-" * 70)
    for theta in (3000.0, 600.0, 300.0):
        sig = U3.sigma_at(theta, float(U3.level(theta, 0.0, LB)))
        for n in (13, 25, 118):
            for cv in (0.05, 0.14):
                a = sd_at(theta, n, 0, cv, sig)
                b = sd_at(theta, n, 1, cv, sig)
                c = sd_at(theta, n, 2, cv, sig)
                worst = max(worst, abs(a - b))
                print(f"    {theta:7.0f} {n:5d} {cv:5.0%} {a:12.6f} "
                      f"{b:12.6f} {abs(a - b):12.2e} {c:12.6f}")
    print(f"\n  largest difference between k = 0 and k = 1: {worst:.3e}")
    print()
    print("  So breath tests are bought TWO AT A TIME. The first reading of")
    print("  a pair establishes nothing because the subject's recovery")
    print("  factor is free, and only the contrast between two readings")
    print("  carries information. A protocol that books an odd number books")
    print("  one it cannot use.")
    return worst


def part_d():
    print()
    print("=" * 88)
    print("PART D.  where the curve is worth walking")
    print("=" * 88)
    print("\n  Tier 1 footing, 5% breath CV. Reading the curve as a protocol")
    print("  choice: what each pair of breath tests saves in draws.\n")
    for theta in (3000.0, 600.0, 300.0):
        sig = U3.sigma_at(theta, float(U3.level(theta, 0.0, LB)))
        prev = None
        print(f"  theta = {theta:.0f}")
        print(f"    {'k':>5} {'draws':>7} {'years':>7} {'draws saved':>13} "
              f"{'per breath test':>17}")
        print("    " + "-" * 54)
        for k in (0, 2, 4, 8, 16, 32, 64):
            n = smallest_n(theta, k, 0.05, sig, kappa_cv=0.1402, fn=sd_med)
            if n is None:
                print(f"    {k:5d} {'-':>7}")
                continue
            if prev is None:
                print(f"    {k:5d} {n:7d} {n * MONTH / YEAR:7.1f}")
            else:
                dk = k - prev[0]
                dn = prev[1] - n
                print(f"    {k:5d} {n:7d} {n * MONTH / YEAR:7.1f} "
                      f"{dn:13d} {dn / dk:17.1f}")
            prev = (k, n)
        print()


if __name__ == "__main__":
    part_a()
    part_b()
    w = part_c()
    part_d()
    print("=" * 88)
    print(f"k=0 vs k=1 verdict: largest difference {w:.3e}")
    print("=" * 88)
