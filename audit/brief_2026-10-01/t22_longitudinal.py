"""
Round-3 add-on. Two questions about realistic use.

(1) SAMPLING AT REALISTIC k_e: sd(ln theta) at 5% vs draws/year, holding
    sd(ln k_e) = 0.06 (what the safe transient delivers with L0 to 1%, R-B).
(2) LONGITUDINAL PRECISION: sd(ln theta_2 - ln theta_1) within one patient,
    where the shared k_e and load-calibration errors cancel. This is the right
    quantity if the use is MONITORING rather than classification.
"""
import numpy as np

import core as C

HORIZON = 730.0  # two-year programme


def sd_rate(n, ske):
    def one(s):
        return C.sd_th(C.fim(C.record(C.TH_5, seed=s, times=np.linspace(0, HORIZON, n)),
                             noise="perdraw", ke_prior=ske))
    return C.over_seeds(one)[0]


def data_fim(r):
    """record FIM with NO shared priors (data only; K_m already unprioured)."""
    return C.fim(r, noise="perdraw", ke_prior=0.0, b_prior=np.inf, km_prior=None)


def longitudinal(n1, n2, s1, s2, ke_prior=0.8, bfrac=0.05,
                 drift_ke=0.0, drift_b=0.0):
    """
    Two monthly series of ONE patient. Parameters
        z = (ln th1, ln th2, ln K_m, b1, b2, ln ke1, ln ke2),
    with K_m shared exactly; b and ln ke shared but allowed to DRIFT between the
    two visits by independent increments of sd drift_b (flux) and drift_ke (log).
    drift = 0 pins them equal (perfect cancellation). Returns
    (sd(ln th2 - ln th1), cross-sectional sd(ln th1) in this joint design).
    """
    r1 = C.record(C.TH_5, seed=s1, times=np.linspace(0, HORIZON, n1))
    r2 = C.record(C.TH_5, seed=s2, times=np.linspace(0, HORIZON, n2))
    L = r1["L"]
    F1, F2 = data_fim(r1), data_fim(r2)
    Z = np.zeros((7, 7))
    #            lnth1 lnKm  b1   lnke1          lnth2 lnKm  b2   lnke2
    idx1 = [0, 2, 3, 5]
    idx2 = [1, 2, 4, 6]
    for a in range(4):
        for b_ in range(4):
            Z[idx1[a], idx1[b_]] += F1[a, b_]
            Z[idx2[a], idx2[b_]] += F2[a, b_]

    def couple(i, j, sd):          # (x_i - x_j) ~ N(0, sd^2)
        if sd <= 0:
            g = 1e12               # pin equal
        else:
            g = 1.0 / sd ** 2
        Z[i, i] += g; Z[j, j] += g; Z[i, j] -= g; Z[j, i] -= g

    couple(3, 4, drift_b * L)      # b1,b2 share, drift in flux units
    couple(5, 6, drift_ke)         # ln ke1, ln ke2 share, drift in log units
    # population priors on the shared quantities (one each, split across the pair
    # via a tight common-mean prior is unnecessary: put the pop prior on each)
    Z[5, 5] += ke_prior ** -2
    Z[6, 6] += ke_prior ** -2
    Z[3, 3] += (bfrac * L) ** -2
    Z[4, 4] += (bfrac * L) ** -2
    P = np.linalg.inv(Z)
    e = np.array([-1, 1, 0, 0, 0, 0, 0.])
    return float(np.sqrt(e @ P @ e)), float(np.sqrt(P[0, 0]))


def long_med(n, **kw):
    return float(np.median([longitudinal(n, n, 2 * s, 2 * s + 1, **kw)[0]
                            for s in range(8)]))


def main():
    print("=" * 92)
    print("  ADD-ON 1&2. SAMPLING AT REALISTIC k_e, AND LONGITUDINAL PRECISION")
    print("=" * 92)

    print("\n  (1) sd(ln theta) AT 5% vs SAMPLING, with sd(ln k_e) = 0.06 (feasible),")
    print("      over a two-year programme (median over 24 seeds):")
    print(f"      {'protocol':>22}  {'draws/2yr':>9}  {'k_e=0.06':>9}  {'k_e exact (ref)':>15}")
    for n, lab in ((25, "monthly (~25)"), (52, "fortnightly (~52)"),
                   (104, "weekly (~104)"), (365, "~daily/2d (~365)"), (730, "daily (~730)")):
        print(f"      {lab:>22}  {n:>9}  {sd_rate(n, 0.06):>9.3f}  {sd_rate(n, 1e-4):>15.3f}")
    print("      The k_e = 0.06 floor (|beta|*0.06 ~ 0.30) is approached by daily; monthly")
    print("      sits above it at ~0.47, so sampling DOES buy ~1/3 (0.47 -> 0.31). None")
    print("      reaches 0.20: a feasible monthly protocol lands at ~0.47, weekly ~0.38.")

    print("\n  (2) LONGITUDINAL CONTRAST sd(ln theta_2 - ln theta_1), ONE patient.")
    print("      The shared k_e, load-bias b and K_m cancel in the difference, so the")
    print("      contrast is limited only by the INDEPENDENT per-visit noise -- far below")
    print("      the cross-sectional bar. (median over seed pairs; k_e prior 0.8 as usual)")
    print(f"      {'draws per visit':>16}  {'sd(d ln theta)':>14}  {'cross-sectional sd(ln theta)':>28}")
    for n in (25, 52, 104, 365):
        d = long_med(n)
        x = float(np.median([longitudinal(n, n, 2 * s, 2 * s + 1)[1] for s in range(8)]))
        print(f"      {n:>16}  {d:>14.3f}  {x:>28.3f}")
    print("      So monitoring a CHANGE is 0.26 (monthly) / 0.13 (weekly) -- the weekly")
    print("      contrast is BELOW 0.20, which no cross-sectional design reaches. This is")
    print("      the one regime where the marker is clinically usable at target precision.")

    print("\n  (2b) BUT ONLY IF k_e AND THE LOAD CALIBRATION ARE STABLE BETWEEN VISITS.")
    print("      The cancellation is exact only if the shared nuisances do not drift. The")
    print("      'interval between estimates' enters ONLY through drift (with perfect")
    print("      sharing the contrast is interval-independent). Sensitivity, n=25/visit:")
    print(f"      {'k_e drift (log units)':>22}  {'sd(d ln theta)':>14}")
    for dk in (0.0, 0.02, 0.05, 0.10, 0.20):
        print(f"      {dk:>22.2f}  {long_med(25, drift_ke=dk):>14.3f}")
    print(f"      {'load-bias drift (frac L)':>22}  {'sd(d ln theta)':>14}")
    for db in (0.0, 0.02, 0.05, 0.10):
        print(f"      {db:>22.2f}  {long_med(25, drift_b=db):>14.3f}")
    print("      The advantage survives drift up to ~0.05 log units in k_e (contrast")
    print("      ~0.34); beyond that it erodes toward the cross-sectional value. Whether")
    print("      k_e is that stable over the monitoring interval is a clinical assumption,")
    print("      not something the repository can settle.")


if __name__ == "__main__":
    main()
