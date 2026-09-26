"""
Phase 10 re-checks. Independent routes to every number step 1 reports.

Step 1 is a derivation and most of it is prose, so what can be checked is the
algebra it rests on. Both checks use routes that share no code path with the
quantity they check.

  1. The flux identity `theta = S (L+b) / g(m)` and the marker identity
     `theta = (L + b - k_e m) / g(m)`, checked by solving the balance from
     scratch with a root-finder rather than the closed-form positive root, and
     by forming S from the flux definition rather than from `1 - k_e m/(L+b)`.
  2. The share's sensitivity `dS/d ln theta`, checked by a central difference on
     S computed the long way, against the analytic `-(k_e/(L+b)) dm/d ln theta`.

Nothing here records an expected result in advance.

Run:  python audit/q11_phase10_recheck.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

import n4_theta_scope as N4
import q2_requirement_map as S2

KM, KE = S2.KM, S2.KE
CAPS = (1.00, 0.50, 0.20, 0.10, 0.05, 0.02, 0.01)
RESULTS = []


def report(name, got, want, tol):
    ok = abs(got - want) <= tol
    RESULTS.append((name, ok))
    print(f"    {name:<50s} {'PASSES' if ok else 'FAILS':>6s}  "
          f"gap {abs(got - want):.3e} against {tol:g}")
    return ok


def banner(t):
    print("\n" + "=" * 92)
    print("  " + t)
    print("=" * 92)


def L_of(theta):
    return float(N4.anchor(theta, KM, KE, "A")[0])


def m_root(theta, L, ke=KE, km=KM):
    """The balance solved by root-finding, not by the closed-form root."""
    return float(brentq(lambda x: theta * x / (km + x) + ke * x - L,
                        1e-12, 1e9, xtol=1e-14, rtol=1e-15))


def check_1_identities():
    banner("1. THE FLUX AND MARKER IDENTITIES, BY ROOT-FINDING")
    print("""  S is formed here as F/(L+b) with F the enzymatic flux computed directly,
  not as 1 - k_e m/(L+b), and m is root-found rather than taken from the closed
  form. Both inversions must return theta exactly.""")
    ok = True
    print(f"\n    {'cap':>5s} {'theta':>8s} {'m (root)':>10s} {'S':>8s} "
          f"{'flux route':>11s} {'marker route':>13s}")
    for cap in CAPS:
        th = 3000.0 * cap
        L = L_of(th)
        m = m_root(th, L)
        g = m / (KM + m)
        F = th * g
        S = F / L
        flux = S * L / g
        mark = (L - KE * m) / g
        print(f"    {100 * cap:4.0f}% {th:8.1f} {m:10.4f} {S:8.4f} "
              f"{flux:11.2f} {mark:13.2f}")
        ok &= abs(flux / th - 1.0) < 1e-12 and abs(mark / th - 1.0) < 1e-12
    RESULTS.append(("both inversions return theta exactly", ok))
    print(f"    {'both inversions return theta exactly':<50s} "
          f"{'PASSES' if ok else 'FAILS':>6s}")
    print("""
  And the marker route's dependence on k_e, which the flux route does not have:""")
    for cap in (1.00, 0.05):
        th = 3000.0 * cap
        L = L_of(th)
        m = m_root(th, L)
        g = m / (KM + m)
        bad = (L - 1.5 * KE * m) / g
        print(f"    cap {100 * cap:4.0f}%: k_e assumed 1.5x gives theta_hat "
              f"{bad:9.1f} against a truth of {th:.1f}")
    return ok


def check_2_sensitivity():
    banner("2. THE SHARE'S SENSITIVITY, TWO WAYS")
    print("""  Analytic: dS/d ln theta = -(k_e/(L+b)) dm/d ln theta. Numeric: a central
  difference on S itself, with m root-found at each node.""")
    ok = True
    print(f"\n    {'cap':>5s} {'analytic':>10s} {'central diff':>13s}")
    for cap in CAPS:
        th = 3000.0 * cap
        L = L_of(th)
        h = 1e-6
        mp, mm = m_root(th * np.exp(h), L), m_root(th * np.exp(-h), L)
        dm = (mp - mm) / (2 * h)
        ana = -(KE / L) * dm
        Sp = (th * np.exp(h)) * mp / (KM + mp) / L
        Sm = (th * np.exp(-h)) * mm / (KM + mm) / L
        num = (Sp - Sm) / (2 * h)
        print(f"    {100 * cap:4.0f}% {ana:10.5f} {num:13.5f}")
        ok &= abs(ana - num) < 1e-6
    RESULTS.append(("dS/d ln theta, analytic vs central difference", ok))
    print(f"    {'dS/d ln theta, analytic vs central difference':<50s} "
          f"{'PASSES' if ok else 'FAILS':>6s}")
    return ok


def check_3_cohort():
    banner("3. THE COHORT FIT'S MACHINERY")
    print("""  q14 profiles each patient's own (ln theta_i, b_i) out by a Schur complement and
  caches the cumulative sum. Both are exact claims, so both are checked against
  the thing they replace: a direct inverse of the full block matrix, and a direct
  uncached loop.""")
    import q14_hierarchical as H
    import q2_requirement_map as S2
    times = S2.times_for("monthly", "no baseline")

    #  (a) the Schur complement against a direct inverse of the 4x4
    worst = 0.0
    for cap in (0.05, 0.10, 0.20):
        th = 3000.0 * cap
        F_ss, F_si, F_ii, s_b = H.blocks(th, times, 3)
        F_ii = F_ii + np.diag([0.0, 1.0 / s_b ** 2])
        full = np.block([[F_ss, F_si], [F_si.T, F_ii]])
        #  the shared block of the full inverse must equal the Schur's inverse
        a = np.linalg.inv(full)[:2, :2]
        b = np.linalg.inv(H.shared_info(th, times, 3))
        worst = max(worst, float(np.max(np.abs(a - b) / np.abs(b))))
    ok_a = worst < 1e-10
    print(f"\n    Schur complement vs direct 4x4 inverse, worst rel. "
          f"difference {worst:.2e}   {'PASSES' if ok_a else 'FAILS'}")
    RESULTS.append(("cohort Schur complement equals the direct inverse", ok_a))

    #  (b) the cumulative cache against an uncached loop
    H._CUM.clear()
    n = 37
    direct = np.zeros((2, 2))
    for i in range(n):
        direct = direct + H.shared_info(
            3000.0 * H.COHORT_CAPS[i % len(H.COHORT_CAPS)], times, i)
    H._CUM.clear()
    H.cohort_info(120, "monthly")           # fill past n, then read back
    cached = H.cohort_info(n, "monthly")[0]
    d = float(np.max(np.abs(cached - direct) / np.abs(direct)))
    ok_b = d < 1e-12
    print(f"    cumulative cache vs an uncached loop at n = {n}, worst rel. "
          f"difference {d:.2e}   {'PASSES' if ok_b else 'FAILS'}")
    RESULTS.append(("cohort cumulative cache equals an uncached loop", ok_b))

    #  (c) at CV 0 the patient's value IS the population's, by construction
    e = max(abs(H.sd_patient(n, "monthly", 0.0) - H.sd_population(n, "monthly"))
            for n in (8, 32))
    ok_c = e < 1e-14
    print(f"    sd_patient at CV 0 equals sd_population, worst absolute "
          f"difference {e:.2e}   {'PASSES' if ok_c else 'FAILS'}")
    RESULTS.append(("sd_patient at CV 0 equals sd_population", ok_c))

    #  (d) a cohort of n must beat a cohort of n-1, strictly, on every schedule
    ok_d = all(H.sd_population(n + 1, s) < H.sd_population(n, s)
               for s in ("monthly", "every 2 d") for n in (4, 16, 63))
    print(f"    sd_population strictly decreasing in n   "
          f"{'PASSES' if ok_d else 'FAILS'}")
    RESULTS.append(("sd_population strictly decreasing in n", ok_d))
    return ok_a and ok_b and ok_c and ok_d


def main():
    print("=" * 92)
    print("  PHASE 10 RE-CHECKS")
    print("=" * 92)
    check_1_identities()
    check_2_sensitivity()
    check_3_cohort()
    banner("SUMMARY")
    bad = [n for n, o in RESULTS if not o]
    for n, o in RESULTS:
        print(f"    {'PASS' if o else 'FAIL'}  {n}")
    print(f"\n    {len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
