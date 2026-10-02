"""
Round-3 close: stress-test the longitudinal contrast (the project's main result).

(1) The contrast under IMPERFECT sharing: k_e drift, load-calibration drift, and
    an antibiotic course between the two windows. When does it still deliver?
(2) What real capacity change it can DETECT at 80% power, and the minimum
    programme length (the contrast is set by draws-per-window, not calendar span).
"""
import numpy as np

import core as C

# 80%-power two-sided 5% test: detectable effect = (z_{.975}+z_{.8}) * sd
ZPOW = 1.959963985 + 0.841621234


def data_fim(r):
    return C.fim(r, noise="perdraw", ke_prior=0.0, b_prior=np.inf, km_prior=None)


def contrast(n, win, s1, s2, ke_prior=0.8, bfrac=0.05, drift_ke=0.0,
             drift_b=0.0, delta2=0.0):
    """sd(ln th2 - ln th1) for two windows of n draws over `win` days; and, if
    delta2>0, the linear-response fake contrast from an antibiotic load drop of
    delta2 on window 2 that the model does not know about."""
    r1 = C.record(C.TH_5, seed=s1, times=np.linspace(0, win, n))
    r2 = C.record(C.TH_5, seed=s2, times=np.linspace(0, win, n))
    L = r1["L"]
    F1, F2 = data_fim(r1), data_fim(r2)
    Z = np.zeros((7, 7))
    i1, i2 = [0, 2, 3, 5], [1, 2, 4, 6]   # lnth1|lnth2, lnKm(shared), b1|b2, lnke1|lnke2
    for a in range(4):
        for b in range(4):
            Z[i1[a], i1[b]] += F1[a, b]
            Z[i2[a], i2[b]] += F2[a, b]

    def cp(i, j, sd):
        g = 1e12 if sd <= 0 else 1.0 / sd ** 2
        Z[i, i] += g; Z[j, j] += g; Z[i, j] -= g; Z[j, i] -= g
    cp(3, 4, drift_b * L)
    cp(5, 6, drift_ke)
    Z[5, 5] += ke_prior ** -2; Z[6, 6] += ke_prior ** -2
    Z[3, 3] += (bfrac * L) ** -2; Z[4, 4] += (bfrac * L) ** -2
    P = np.linalg.inv(Z)
    e = np.array([-1, 1, 0, 0, 0, 0, 0.])
    sd = float(np.sqrt(e @ P @ e))
    fake = None
    if delta2:
        m2 = C.m_ss(C.TH_5, C.KM0, C.KE, r2["loads"] * (1 - delta2), r2["inc"])
        pred = (C.TH_5 + r2["inc"]) * m2 / (C.KM0 + m2) + C.KE * m2
        res = r2["loads"] - pred
        sg2 = C.sigmas(r2, "perdraw")
        Q2 = C.flux_rows(dict(theta=C.TH_5, Km=C.KM0, ke=C.KE, m=m2, inc=r2["inc"], L=L))
        score = np.zeros(7)
        g2 = (Q2.T / sg2 ** 2) @ res
        for a in range(4):
            score[i2[a]] += g2[a]
        fake = float((P @ score)[1] - (P @ score)[0])
    return sd, fake


def med_sd(n, win, **kw):
    return float(np.median([contrast(n, win, 2 * s, 2 * s + 1, **kw)[0] for s in range(8)]))


def med_fake(n, win, delta2):
    return float(np.median([contrast(n, win, 2 * s, 2 * s + 1, delta2=delta2)[1] for s in range(8)]))


def pct(lnd):
    return f"{np.exp(lnd) - 1:+.0%}"


def main():
    print("=" * 92)
    print("  CLOSE. STRESS-TESTING THE LONGITUDINAL CONTRAST")
    print("=" * 92)

    print("\n  (1a) THE CONTRAST DEPENDS ON DRAWS-PER-WINDOW, NOT CALENDAR SPAN.")
    print(f"      {'window':>8}  {'monthly':>18}  {'weekly':>16}  {'daily':>16}")
    for win, wl in ((180, "6 months"), (365, "1 year"), (730, "2 years")):
        cells = []
        for cad in (30.4, 7, 1):
            n = max(int(round(win / cad)), 3)
            cells.append(f"n={n:<4d} {med_sd(n, win):.3f}")
        print(f"      {wl:>8}  {cells[0]:>18}  {cells[1]:>16}  {cells[2]:>16}")
    print("      So two SHORT dense windows work: 6-month daily (n=180) gives 0.10, and")
    print("      ~45 draws/window reaches 0.20. You do not need years of weekly draws;")
    print("      a ~6-week daily burst per window suffices. Calendar span barely matters")
    print("      at fixed n -- but a longer GAP between windows allows more drift (1b).")

    print("\n  (1b) IMPERFECT SHARING: where the contrast degrades past 0.20.")
    print("      k_e drift between windows (weekly n=104 base 0.13; monthly n=25 base 0.26):")
    print(f"      {'drift (log units)':>18}  {'weekly':>8}  {'monthly':>8}")
    for dk in (0.0, 0.02, 0.05, 0.10, 0.20):
        print(f"      {dk:>18.2f}  {med_sd(104, 730, drift_ke=dk):>8.3f}  {med_sd(25, 730, drift_ke=dk):>8.3f}")
    print("      load-calibration (b) drift between windows:")
    print(f"      {'drift (frac L)':>18}  {'weekly':>8}  {'monthly':>8}")
    for db in (0.0, 0.01, 0.02, 0.05):
        print(f"      {db:>18.2f}  {med_sd(104, 730, drift_b=db):>8.3f}  {med_sd(25, 730, drift_b=db):>8.3f}")
    print("      The weekly contrast holds below 0.20 up to ~0.06 log units of k_e drift")
    print("      and ~0.015 L of calibration drift; beyond that it degrades toward the")
    print("      cross-sectional bar. PHYSIOLOGY: k_e is the marker's clearance (renal +")
    print("      non-enzymatic). Over a week it is effectively fixed (<0.01); over a year")
    print("      in a growing child renal clearance scales with size and 0.05-0.1 is")
    print("      plausible. The load calibration b drifts with any diet/formula/adherence")
    print("      change. So the contrast is reliable for SHORT intervals (weeks-months)")
    print("      on a stable regimen, and degrades over a year or across a regimen change.")

    print("\n  (1c) AN ANTIBIOTIC COURSE STARTING BETWEEN WINDOWS IS CATASTROPHIC.")
    print("      A load drop of delta on window 2 that the record does not capture, as a")
    print("      fake capacity change (linear response; low markers at the nominal")
    print("      recorded load read as HIGHER capacity):")
    print(f"      {'antibiotic depth':>18}  {'fake d ln theta':>16}  {'apparent capacity change':>24}")
    for d in (0.10, 0.20, 0.41):
        f = med_fake(104, 730, d)
        print(f"      {d:>17.0%}  {f:>+16.2f}  {pct(f):>24}")
    print("      Even a 10% unrecorded load drop fakes ~+1.1 ln theta (a ~3x apparent")
    print("      rise), ~8x the 0.13 weekly precision; at the realistic 41% metronidazole")
    print("      depth it is ~+2.9 (a ~20x apparent rise). It would MASK any real decline")
    print("      and MIMIC a large improvement. Antibiotic status (and any load-regimen")
    print("      change) MUST be recorded and held constant across the two windows.")

    print("\n  (2) WHAT REAL CHANGE IS DETECTABLE AT 80% POWER (two-sided 5%).")
    print("      Minimum detectable |d ln theta| = 2.80 * sd(contrast); as a rise/fall:")
    print(f"      {'design':>28}  {'contrast sd':>11}  {'min detectable rise':>19}  {'fall':>6}")
    for lab, n, win in (("monthly, 2-yr windows", 25, 730),
                        ("weekly, 1-yr windows", 52, 365),
                        ("weekly, 2-yr windows", 104, 730),
                        ("daily, 6-month windows", 180, 180),
                        ("daily, 2-yr windows", 730, 730)):
        sd = med_sd(n, win)
        mdc = ZPOW * sd
        print(f"      {lab:>28}  {sd:>11.3f}  {pct(mdc):>19}  {pct(-mdc):>6}")
    print("      So weekly windows detect a ~+45% rise (or ~-30% fall); a 6-month daily")
    print("      programme detects ~+30%/-24%; monthly windows only a doubling (~+110%).")
    print("      A 50% capacity rise (ln 1.5 = 0.405) needs contrast sd <= 0.145 -- weekly")
    print("      windows or denser. None of this is absolute capacity: it is CHANGE, and")
    print("      only valid if k_e, the load calibration and antibiotic status are held")
    print("      fixed between the two windows.")


if __name__ == "__main__":
    main()
