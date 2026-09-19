"""
m3 -- can Extended Data Table 3 tell us WHEN it was sampled?

OPEN_QUESTIONS.md section 2 calls the sampling time of Extended Data Table 3
the single largest swing in the design: a trough reading puts monthly sampling
at 10-12 visits to sd(ln theta) = 0.20, a peak reading at 42 or never.

The table's title is "Biomarker percentage change from baseline after the third
dose" and the Methods define baseline as predose. The phase of the post-dose
reading is not stated. But the table reports cohort 1 (0.30 mg/kg Q3W) apart
from cohorts 2-5 (Q2W), and the two candidate conventions differ in what that
should do:

  TROUGH convention  (the sample is the next predose): Q3W reads 21 days after
      its dose and Q2W 14 days after theirs. On a 3.5-day half-life those are
      2^-6 and 2^-4 of one dose, so the Q3W arm should show a SMALLER fall.
  FIXED-CLOCK convention (a set number of hours or days after dose 3): both arms
      sit at the same phase and should show the SAME fall.

THE POWER OF THIS TEST IS STATED AND PRINTED BEFORE THE DATA ARE TOUCHED.
That ordering is deliberate: a discriminator whose power is not established
first will always seem to have found something.

Data entered from the Nature 2024 PDF (PMC11156579 publishes these tables as
images). Mean % change from baseline after dose 3, SD in brackets.
"""
import sys
import numpy as np
sys.path.insert(0, 'src')
import matplotlib
matplotlib.use('Agg')
import design.run as D

P = D.P
L0 = D.load(0.0)
LAM = D.LAMBDA

COHORT = {                       # name: (mg/kg, interval_days, n)
    1: (0.30, 21.0, 4),
    2: (0.30, 14.0, 3),
    3: (0.45, 14.0, 3),
    4: (0.60, 14.0, 3),
    5: (0.90, 14.0, 3),
}
TAB = {   # marker: {cohort: (mean %, sd %)}
    '3-HP': {1: (-29, 25), 2: (-3, 26), 3: (-20, 24), 4: (+5, 57), 5: (-20, 41)},
    '2-MC': {1: (-1, 31), 2: (-2, 25), 3: (-22, 43), 4: (-14, 39), 5: (-17, 37)},
    'C3':   {1: (-26, 22), 2: (-1, 23), 3: (-17, 23), 4: (+7, 42), 5: (-11, 23)},
    'n-PG': {1: (-38, 27), 2: (-36, 18), 3: (-35, 54), 4: (+49, 170), 5: (+19, 108)},
}


def marker(inc):
    B = P.ke * P.Km + P.theta + inc - L0
    return (-B + np.sqrt(B ** 2 + 4.0 * P.ke * L0 * P.Km)) / (2.0 * P.ke)


def pct_change(inc):
    m0, m = marker(0.0), marker(inc)
    return 100.0 * (m - m0) / m0


def accumulated(amp, interval, t):
    """Increment at time t from doses at 0, interval, 2*interval."""
    return sum(amp * np.exp(-LAM * (t - k * interval)) for k in range(3)
               if t >= k * interval)


def welch(m1, s1, n1, m2, s2, n2):
    se = np.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2)
    df = se ** 4 / ((s1 ** 2 / n1) ** 2 / (n1 - 1) + (s2 ** 2 / n2) ** 2 / (n2 - 1))
    return m1 - m2, se, df


def t_crit(df, p=0.975):
    from scipy.stats import t
    return float(t.ppf(p, df))


def main():
    print(__doc__)
    print("=" * 92)
    print("PART A. THE POWER OF THE TEST, ESTABLISHED BEFORE LOOKING AT THE MEANS")
    print("=" * 92)
    print()
    print("A1. The only contrast free of the dose confound is cohort 1 against cohort 2:")
    print("    both are 0.30 mg/kg, and they differ only in interval (21 d against 14 d).")
    print("    Cohorts 3-5 differ from cohort 1 in dose AND interval AND, per Table 1, age.")
    print()
    for mk in TAB:
        _, s1 = TAB[mk][1]; _, s2 = TAB[mk][2]
        _, se, df = welch(0, s1, COHORT[1][2], 0, s2, COHORT[2][2])
        mdd = (t_crit(df) + t_crit(df, 0.80)) * se       # 80% power, two-sided 0.05
        print(f"    {mk:5} SE of the difference {se:6.1f} pp, df {df:4.1f}, "
              f"smallest difference detectable at 80% power {mdd:6.1f} pp")
    print()
    print("A2. What each convention predicts for that same contrast, swept over D_AMP")
    print("    because D_AMP = 6000 is an assumed constant this repository does not move.")
    print()
    print(f"    {'D_AMP':>8} {'Q3W fall %':>12} {'Q2W fall %':>12} {'predicted Q3W - Q2W':>22}")
    print("    trough convention (post-dose sample is the next predose):")
    for amp in (1500.0, 3000.0, 6000.0, 12000.0, 24000.0):
        q3 = pct_change(accumulated(amp, 21.0, 63.0))
        q2 = pct_change(accumulated(amp, 14.0, 42.0))
        print(f"    {amp:8.0f} {q3:12.1f} {q2:12.1f} {q3 - q2:22.1f}")
    print("    fixed-clock convention (one day after dose 3):")
    for amp in (6000.0,):
        q3 = pct_change(accumulated(amp, 21.0, 43.0))
        q2 = pct_change(accumulated(amp, 14.0, 29.0))
        print(f"    {amp:8.0f} {q3:12.1f} {q2:12.1f} {q3 - q2:22.1f}")
    print()
    print("=" * 92)
    print("PART B. THE DATA")
    print("=" * 92)
    print()
    print(f"    {'marker':6} {'c1 Q3W':>10} {'c2 Q2W':>10} {'difference':>12} {'t':>7} {'p':>8}")
    from scipy.stats import t as tdist
    for mk in TAB:
        m1, s1 = TAB[mk][1]; m2, s2 = TAB[mk][2]
        d, se, df = welch(m1, s1, COHORT[1][2], m2, s2, COHORT[2][2])
        tv = d / se
        p = 2 * (1 - tdist.cdf(abs(tv), df))
        print(f"    {mk:6} {m1:+6.0f}({s1:3.0f}) {m2:+6.0f}({s2:3.0f}) "
              f"{d:+12.1f} {tv:7.2f} {p:8.3f}")
    print()
    print("PART C. Is there a dose trend among the four Q2W cohorts at all?")
    print("    (if the falls do not order by dose, cohort-to-cohort differences are")
    print("     not reading dose and cannot be read as reading interval either)")
    print()
    doses = np.array([COHORT[c][0] for c in (2, 3, 4, 5)])
    for mk in TAB:
        means = np.array([TAB[mk][c][0] for c in (2, 3, 4, 5)], dtype=float)
        sds = np.array([TAB[mk][c][1] for c in (2, 3, 4, 5)], dtype=float)
        w = 3.0 / sds ** 2
        b = (np.sum(w * (doses - np.average(doses, weights=w))
                    * (means - np.average(means, weights=w)))
             / np.sum(w * (doses - np.average(doses, weights=w)) ** 2))
        r = np.corrcoef(doses, means)[0, 1]
        print(f"    {mk:6} weighted slope {b:+8.1f} pp per mg/kg, "
              f"unweighted r = {r:+.2f} over four dose levels")


def part_d():
    """The comparison that IS informative: the near-peak reading is far from
    every observed cohort mean, so it can be tested even though the two
    conventions above cannot be separated from each other."""
    from scipy.stats import t as tdist
    print()
    print("PART D. Can a near-peak reading be excluded?")
    print("    A sample one day after dose 3 predicts the fall below, at the assumed")
    print("    D_AMP = 6000. Each cohort mean is tested against it (one sample, n per")
    print("    cohort). This is conditional on D_AMP and is not a free test.")
    print()
    pred = pct_change(accumulated(6000.0, 14.0, 29.0))
    print(f"    predicted fall one day after dose 3, Q2W: {pred:.1f}%")
    print()
    print(f"    {'marker':6} {'cohort':>7} {'mean(sd)':>12} {'n':>3} {'t vs peak':>11} {'p':>8}")
    for mk in TAB:
        for c in (1, 2, 3, 4, 5):
            m, sd = TAB[mk][c]
            n = COHORT[c][2]
            tv = (m - pred) / (sd / np.sqrt(n))
            pv = 2 * (1 - tdist.cdf(abs(tv), n - 1))
            print(f"    {mk:6} {c:7d} {str(m) + '(' + str(sd) + ')':>12} {n:3d} "
                  f"{tv:11.2f} {pv:8.3f}")


if __name__ == "__main__":
    main()
    part_d()
