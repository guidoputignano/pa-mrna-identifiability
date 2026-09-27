"""
The fast is withdrawn, and the substitute is priced: can suppressing
gut-derived propionate serve as the trough perturbation for a relaxation burst?

WHY THE FAST IS WITHDRAWN. audit/u14_dynamic_burst.py priced a 50% "fast" as a
load REDUCTION. In propionic acidemia it is not one. 01_biology.tex:151-156
already says so in the paper's own words: the dangerous events are triggered by
"infection, FASTING, surgery, trauma", and "the patient stops eating, the body
enters catabolism, and skeletal muscle protein is degraded for energy... the
propionate load rises sharply precisely when the patient is eating nothing."
That is why emergency management gives intravenous glucose to suppress
catabolism rather than restricting protein (01_biology.tex:160).

AND THE PAPER HAS ALREADY PRICED IT. 07b_decision.tex:55-62 parametrises the
recorded load as L_base + rho (s_true - L_base), so rho <= 0 is exactly "crisis
protocol, where prescribed intake is cut while catabolic muscle breakdown
supplies the load". At rho = -0.25, protein stopped, the sign error is 100% and
the median theta_hat/theta is 0.009: the estimate collapses by two orders of
magnitude. u14's fast arm assumed a genuine reduction that the record tracks,
which is rho = 1, and PA's response to a prescribed cut is rho < 0.

So u14's fast numbers are not merely unsupported, they contradict an existing
paper result about the same intervention. They are withdrawn.

WHAT THE SUBSTITUTE HAS TO BE. A load reduction that does NOT provoke
catabolism, so that the true load genuinely falls and the record tracks it.
Suppressing gut-derived propionate removes substrate at source with no
catabolic stress, and it is routine management.

THIS IS A DIFFERENT USE FROM THE ONE ALREADY CLOSED, and the earlier negative
does not transfer. audit/u2_breath_functional.py part E found a load step
WORSE than the dose cycle for a breath RATIO, because a ratio wants the widest
marker-level contrast and the dose train supplies a larger one. Here the step
is not widening a contrast: it is the input whose RELAXATION carries k_e, at
the trough where k_e owns the relaxation rate (R113). Different quantity,
different design, different answer required.

THE AMPLITUDE IS ALREADY REGISTERED, AND TWO DIFFERENT NUMBERS ARE IN PLAY.
  src/common/params.py, the s_base row: the partition of whole-body propionate
  production is measured, gut flora "at least 22.2 +/- 8.5% by the fall under
  metronidazole". That is a share of PRODUCTION, and it bounds what suppressing
  the gut can remove.
  audit/m2_antibiotic_cycle.py, from Thompson 1990 (PMID 2226555, 4 PA and 5
  MMA children, verified in that file): on metronidazole urinary propionate
  metabolites fell a mean of 41% (range 12 to 76) and plasma propionate from
  45.0 to 25.1 umol/L. That is a fall in EXCRETED METABOLITES and in the
  marker, not a share of production.
The two are different objects and the amplitude is swept across both rather
than one being chosen.

THE OBJECTION THAT DECIDES IT, AND IT IS ABOUT TIME AND NOT AMPLITUDE. A
relaxation is resolvable only if the input moves faster than the relaxation it
drives. The marker's relaxation time at the trough is 2.65 h at a tenth of
nominal (R113). Metronidazole is dosed over one to two weeks
(Baumgartner 2014, in m2's docstring) and flora suppression develops over
days. PART C sweeps the input's onset time constant for that reason.

THAT OBJECTION IS THE ONE THIS FILE EXPECTED TO BE DECISIVE, AND THE ANSWER
TOOK TWO CORRECTIONS TO REACH. Both are left visible rather than absorbed.

  First reading (this docstring, as written above): the onset decides it.
  Second reading: refuted, the amplitude moves sd(ln k_e) about 2.4x more
    over each input's registered range, and the onset's ordering is not even
    monotone. PART E computed that and it was arithmetically correct.
  Third reading (PART E as it now stands): BOTH were wrong as stated,
    because each holds in one regime only. At a 5% assay and one admission
    the record dominates the posterior and the amplitude moves it more; at a
    1% assay and sixteen admissions the burst dominates and the ONSET binds.
    A "which input dominates" claim is not well posed without naming the
    regime, and a span ratio is evidence only where it was computed.

The second reading was also computed on a DEFECTIVE load-bias sensitivity,
fixed in burst_rows_ramp below and recorded at R119, which made the burst
look 2.6x less informative than it is and so hid the regime where the onset
binds. PART F is the finite-difference guard that caught it and is now
permanent.

FOOTING as audit/u14_dynamic_burst.py, whose sensitivity system this imports.

Nothing here records an expected result in advance.

Run:  python audit/u15_gut_suppression.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.integrate import solve_ivp

import u14_dynamic_burst as U14

KM, KE, LB = U14.KM, U14.KE, U14.LB
HOUR, TROUGH = U14.HOUR, U14.TROUGH
LAMBDA = U14.LAMBDA
TIER1_N, SEEDS = U14.TIER1_N, U14.SEEDS
FRONTIER = U14.FRONTIER
PUBLISHED_KAPPA_CV = U14.PUBLISHED_KAPPA_CV

#: src/common/params.py, s_base row: gut flora's share of production
GUT_SHARE, GUT_SD = 0.222, 0.085
#: audit/m2_antibiotic_cycle.py, Thompson 1990: the fall in metabolites
THOMPSON_FALL, THOMPSON_RANGE = 0.41, (0.12, 0.76)


def burst_rows_ramp(theta, times_d, amp, tau_in_d, assay_cv=U14.ASSAY_CV,
                    Km=KM, ke=KE, Ltot=None):
    """
    Rows of d(ln m)/d(params) for a load suppression with a FINITE onset.

    L(t) = L (1 - amp (1 - exp(-t/tau_in))), so tau_in = 0 is the ideal step
    u14 assumed and a large tau_in is an input slower than the relaxation it
    drives. The increment is held at the trough throughout, which is where
    k_e owns the relaxation rate.
    """
    Ltot = LB if Ltot is None else Ltot
    m0 = float(U14.qss(theta, TROUGH, Ltot, Km, ke))
    J0 = (theta + TROUGH) * Km / (Km + m0) ** 2 + ke
    y0 = [m0,
          -theta * (m0 / (Km + m0)) / J0,
          (theta + TROUGH) * Km * m0 / (Km + m0) ** 2 / J0,
          1.0 / J0,
          -ke * m0 / J0]

    def frac(t):
        return 1.0 if tau_in_d <= 0 else (1.0 - np.exp(-t / tau_in_d))

    def Lt(t):
        return Ltot * (1.0 - amp * frac(t))

    def rhs(t, y):
        m, qt, qK, qb, qk = y
        m = max(m, 1e-14)
        th = theta + TROUGH * np.exp(-LAMBDA * t)
        J = th * Km / (Km + m) ** 2 + ke
        #  THE LOAD-BIAS FORCING CARRIES THE INTERVENTION FACTOR, and this
        #  is where an earlier version of this file was wrong (R119). u14's
        #  burst_rows uses a forcing of exactly 1.0, which is correct for
        #  what u14 models: a prescribed change of KNOWN ABSOLUTE size, with
        #  the unknown bias added after it, so dL/db = 1. A gut-derived
        #  share is not that. It removes a PROPORTION of the load actually
        #  being produced, so if true production is higher than recorded the
        #  gut share of it is proportionally higher too:
        #      L(t) = (Ltot + b)(1 - amp f(t))   =>   dL/db = 1 - amp f(t).
        #  Inheriting 1.0 here made the burst look 2.6x LESS informative
        #  than it is, and the finite-difference check in part_f is what
        #  catches it: every other column agreed to 3e-9 and this one to
        #  2.2e-1.
        return [Lt(t) - th * m / (Km + m) - ke * m,
                -theta * (m / (Km + m)) - J * qt,
                th * Km * m / (Km + m) ** 2 - J * qK,
                (1.0 - amp * frac(t)) - J * qb,
                -ke * m - J * qk]

    T = float(max(times_d))
    sol = solve_ivp(rhs, (0.0, T + 1e-9), y0, dense_output=True,
                    rtol=1e-11, atol=1e-14, method="LSODA")
    Y = sol.sol(np.asarray(times_d, dtype=float))
    m = Y[0]
    return (np.stack([Y[1] / m, Y[2] / m, Y[3] / m, Y[4] / m], axis=1)
            / assay_cv), m


def sd_ke(theta, times, amp, tau_in, burst_cv=U14.ASSAY_CV, n_bursts=1,
          n_rec=TIER1_N):
    """Median sd(ln k_e) from the Tier 1 record plus the ramped burst."""
    R, _ = burst_rows_ramp(theta, times, amp, tau_in, burst_cv)
    R = np.column_stack([R, np.zeros(R.shape[0])])      # no kappa in a burst
    out = []
    for s in SEEDS:
        F = U14.fisher(theta, n_rec, s, kappa_cv=PUBLISHED_KAPPA_CV)
        out.append(U14.sd_of(F + n_bursts * (R.T @ R), 3))
    v = np.array(out, float)
    f = v[np.isfinite(v)]
    return float(np.median(f)) if f.size else np.inf


def qss_floor(theta):
    """sd(ln k_e) from the record alone: what a slow input must converge to."""
    v = np.array([U14.sd_of(U14.fisher(theta, TIER1_N, s,
                                       kappa_cv=PUBLISHED_KAPPA_CV), 3)
                  for s in SEEDS], float)
    return float(np.median(v[np.isfinite(v)]))


# ======================================================================
def part_a():
    print("=" * 88)
    print("PART A.  the fast is withdrawn, and the paper had already priced it")
    print("=" * 88)
    print("""
  01_biology.tex:151-156 states the mechanism in the paper's own words: the
  dangerous events are triggered by "infection, FASTING, surgery, trauma",
  and "the patient stops eating, the body enters catabolism, and skeletal
  muscle protein is degraded for energy... the propionate load rises sharply
  precisely when the patient is eating nothing". 01_biology.tex:160 gives the
  consequence for management: intravenous glucose to suppress catabolism.

  07b_decision.tex:55-62 then PRICES a prescribed cut. Writing the recorded
  load as L_base + rho (s_true - L_base), rho <= 0 is "crisis protocol, where
  prescribed intake is cut while catabolic muscle breakdown supplies the
  load", and at rho = -0.25, protein stopped, the sign error is 100% with a
  median theta_hat/theta of 0.009.

  audit/u14's fast arm assumed a genuine reduction the record tracks, which
  is rho = 1. So its fast numbers do not merely lack support: they contradict
  an existing result about the same intervention. WITHDRAWN, and the 0.1121
  at a tenth of nominal is not a safe number and is not carried forward.

  WHAT THE PAPER DOES NOT HAVE is the coupling in the IDENTIFICATION model.
  eq:local carries L as an input dialable in either direction. The catabolic
  response lives in the prose and, as rho, in a RECORDING asymmetry; it is
  not a physiological coupling from a prescribed cut to a rise in the true
  load. For everything the paper currently claims that is sufficient, because
  its load variation is spontaneous day-to-day intake at 10% and its
  decompensations are modelled as load spikes. For any PRESCRIBED reduction
  it is not sufficient, and rho is the hook a coupling would attach to.
""")


def part_b():
    print("=" * 88)
    print("PART B.  the substitute at an IDEAL step, swept over the amplitude")
    print("=" * 88)
    print("\n  24 hourly draws at the trough, Tier 1 record, one admission.")
    print("  The amplitude's two registered anchors are marked.\n")
    tt = HOUR * np.arange(0, 25, dtype=float)
    print(f"    {'amplitude':>10} {'source':>34}" +
          "".join(f"{f'theta={t:.0f}':>12}" for t in (600, 300)))
    print("    " + "-" * 70)
    marks = {0.222: "params.py gut share of production",
             0.137: "that share minus one sd",
             0.307: "that share plus one sd",
             0.41: "Thompson 1990 metabolite fall",
             0.12: "Thompson range, low",
             0.76: "Thompson range, high"}
    for amp in (0.10, 0.12, 0.137, 0.222, 0.307, 0.41, 0.50, 0.76):
        cells = []
        for theta in (600.0, 300.0):
            v = sd_ke(theta, tt, amp, 0.0)
            fr = FRONTIER[theta]
            cells.append(f"{v:.4f}" + ("*" if fr and v <= fr else ""))
        print(f"    {amp:10.1%} {marks.get(amp, ''):>34}"
              + "".join(f"{c:>12}" for c in cells))
    print("\n  * reaches the frontier. Record-only floors: "
          f"{qss_floor(600.0):.4f} at 600, {qss_floor(300.0):.4f} at 300.")


def part_c():
    print()
    print("=" * 88)
    print("PART C.  the onset time constant, which is what decides the route")
    print("=" * 88)
    print("\n  A relaxation is resolvable only if the input moves faster than")
    print("  the relaxation it drives. The marker's relaxation time at the")
    print("  trough is 2.65 h at a tenth of nominal. Amplitude held at the")
    print("  registered gut share of 22.2%. The burst spans four onset time")
    print("  constants or 24 h, whichever is longer, at 24 draws.\n")
    print(f"    {'onset tau':>12} {'burst span':>11}" +
          "".join(f"{f'theta={t:.0f}':>13}" for t in (600, 300))
          + f"{'  what it is':>14}")
    print("    " + "-" * 74)
    labels = {0.0: "ideal step", 1 / 24: "1 h", 3 / 24: "3 h", 6 / 24: "6 h",
              0.5: "12 h", 1.0: "1 day", 3.0: "3 days", 7.0: "1 week"}
    for tau in (0.0, 1 / 24, 3 / 24, 6 / 24, 0.5, 1.0, 3.0, 7.0):
        span = max(1.0, 4.0 * tau)
        tt = np.linspace(0.0, span, 25)
        cells = []
        for theta in (600.0, 300.0):
            v = sd_ke(theta, tt, GUT_SHARE, tau)
            fr = FRONTIER[theta]
            cells.append(f"{v:.4f}" + ("*" if fr and v <= fr else ""))
        note = ("metronidazole" if tau >= 1.0 else "")
        print(f"    {labels[tau]:>12} {span:10.1f}d"
              + "".join(f"{c:>13}" for c in cells) + f"  {note}")
    print()
    print(f"  Record-only floors, which a slow input must converge to: "
          f"{qss_floor(600.0):.4f} and {qss_floor(300.0):.4f}.")
    print()
    print("  The direction is as expected: a slow input converges to the")
    print("  floor, because it leaves the marker on the quasi-steady locus")
    print("  where u14 PART A found no k_e information at any density. The")
    print("  MAGNITUDE is not. PART E measures it against the amplitude's")
    print("  instead of asserting which one decides, because this docstring")
    print("  asserted the onset and the numbers above do not support it.")


def part_d():
    print()
    print("=" * 88)
    print("PART D.  what would have to be true for the route to work")
    print("=" * 88)
    print("\n  Best case: a 1% burst assay and sixteen admissions, at each")
    print("  amplitude and onset. WHICH AMPLITUDE IS DEFENSIBLE MATTERS HERE.")
    print("  params.py's 22.2% is a share of PRODUCTION, so it bounds what")
    print("  suppressing the gut can remove from L. Thompson's 41% is a fall")
    print("  in urinary METABOLITES, which is not a production share, so it")
    print("  is shown for contrast and is not a justified step amplitude.\n")
    print(f"    {'amplitude':>10} {'onset':>11}"
          + "".join(f"{f'theta={t:.0f}':>13}" for t in (600, 300))
          + "   what it rests on")
    print("    " + "-" * 79)
    rest = {0.222: "production share, defensible",
            0.307: "that share plus one sd",
            0.41: "metabolite fall, NOT a share"}
    for amp in (0.222, 0.307, 0.41):
        for tau, lab in ((6 / 24, "6 h"), (1.0, "1 day")):
            ttx = np.linspace(0.0, max(1.0, 4.0 * tau), 25)
            cells = []
            for theta in (600.0, 300.0):
                v = sd_ke(theta, ttx, amp, tau, burst_cv=0.01, n_bursts=16)
                fr = FRONTIER[theta]
                cells.append(f"{v:.4f}" + ("*" if fr and v <= fr else ""))
            print(f"    {amp:10.1%} {lab:>11}" + "".join(f"{c:>13}" for c in cells)
                  + f"   {rest[amp] if tau == 6 / 24 else ''}")
    print(f"\n  * reaches it. Frontiers {FRONTIER[600.0]:.4f} at 600, "
          f"{FRONTIER[300.0]:.4f} at 300.")
    print()
    print("  THE SLOWEST ONSET THAT STILL REACHES THE FRONTIER, which is the")
    print("  actionable number, since the amplitude is capped by biology and")
    print("  the onset is a property of the agent chosen:")
    print()
    print(f"    {'amplitude':>10} {'theta':>7} {'frontier':>9} "
          f"{'slowest onset that reaches it':>31}")
    print("    " + "-" * 61)
    for amp in (0.222, 0.307, 0.41):
        for theta in (600.0, 300.0):
            fr = FRONTIER[theta]
            lo, hi = 0.0, 7.0
            if sd_ke(theta, np.linspace(0.0, 1.0, 25), amp, 0.0,
                     burst_cv=0.01, n_bursts=16) > fr:
                print(f"    {amp:10.1%} {theta:7.0f} {fr:9.4f} "
                      f"{'none, not even an ideal step':>31}")
                continue
            for _ in range(11):
                mid = 0.5 * (lo + hi)
                v = sd_ke(theta, np.linspace(0.0, max(1.0, 4.0 * mid), 25),
                          amp, mid, burst_cv=0.01, n_bursts=16)
                lo, hi = (mid, hi) if v <= fr else (lo, mid)
            tau = 0.5 * (lo + hi)
            print(f"    {amp:10.1%} {theta:7.0f} {fr:9.4f} "
                  f"{f'{24 * tau:.1f} h':>31}")
    print()
    print("  READING. At the DEFENSIBLE amplitude the route DOES reach both")
    print("  frontiers in this best case, which the pre-correction version of")
    print("  this file denied (R119). What it needs is an onset measured in")
    print("  hours. Metronidazole's is measured in days, and at a one-day")
    print("  onset the same best case misses the tighter frontier by 2.0x.")
    print("  So the binding requirement on the INTERVENTION is its time")
    print("  constant, and the binding requirement on the STUDY is the 1%")
    print("  assay and sixteen admissions that put it in this regime at all.")


def part_e():
    print()
    print("=" * 88)
    print("PART E.  which input actually decides it, measured rather than asserted")
    print("=" * 88)
    print("\n  This script's docstring said the onset decides the route.")
    print("  Rather than assert either way, the two inputs are swept over")
    print("  each one's registered range IN BOTH REGIMES, at a tenth of")
    print("  nominal where the requirement is tightest. The relative-span")
    print("  measure is in log units so the two are comparable.\n")
    taus = (0.0, 1 / 24, 3 / 24, 6 / 24, 0.5, 1.0, 3.0, 7.0)
    amps = (0.12, 0.222, 0.307, 0.41)
    print(f"    {'regime':>26} {'amplitude span':>15} {'onset span':>12} "
          f"{'which binds':>22}")
    print("    " + "-" * 78)
    for cv, nb, lab in ((U14.ASSAY_CV, 1, "5% assay, 1 admission"),
                        (0.01, 16, "1% assay, 16 admissions")):
        av = [sd_ke(300.0, HOUR * np.arange(0, 25, dtype=float), a, 0.0,
                    burst_cv=cv, n_bursts=nb) for a in amps]
        tv = [sd_ke(300.0, np.linspace(0.0, max(1.0, 4.0 * t), 25),
                    GUT_SHARE, t, burst_cv=cv, n_bursts=nb) for t in taus]
        #  ratio of largest to smallest, which compares fairly across
        #  regimes whose absolute levels differ by an order of magnitude
        ar, tr = max(av) / min(av), max(tv) / min(tv)
        who = "amplitude" if ar > tr else "ONSET"
        print(f"    {lab:>26} {f'x{ar:.2f}':>15} {f'x{tr:.2f}':>12} "
              f"{f'{who}, by {max(ar, tr) / min(ar, tr):.1f}x':>22}")
    print()
    print("  READING, AND IT IS NOT THE ONE EITHER EARLIER VERSION GAVE.")
    print("  WHICH INPUT BINDS DEPENDS ON THE REGIME. With a clinic assay")
    print("  and one admission the record dominates the posterior, both")
    print("  inputs move it little, and the amplitude moves it more. With a")
    print("  1% assay and sixteen admissions the burst dominates and the")
    print("  ONSET binds: at the registered amplitude a 6 h onset reaches")
    print("  the frontier and a one-day onset misses it by 2.0x.")
    print()
    print("  So the docstring's expectation was right in the regime that")
    print("  matters for the verdict and wrong in the one that is cheap to")
    print("  run, and the intermediate 'corrected reading' this file carried")
    print("  before R119 was wrong in the regime that matters. Neither")
    print("  single claim was available: the answer is regime-dependent.")


def part_f():
    """
    Every sensitivity column against a central finite difference of the
    level trajectory alone, which shares no line with the sensitivity ODE.

    THIS EXISTS BECAUSE IT CAUGHT A REAL ERROR. The load-bias column was
    inherited from u14, whose forcing of 1.0 is correct for a prescribed
    change of known absolute size and wrong for a proportional suppression.
    Three columns agreed to 3e-9 and that one to 2.2e-1, and the consequence
    on sd(ln k_e) was a factor of 2.6. This is the third time in this audit
    that a sensitivity disagreed with its own forward model (R101 twice, and
    now R119), so the check is permanent rather than ad hoc.
    """
    print()
    print("=" * 88)
    print("PART F.  every sensitivity column against a finite difference")
    print("=" * 88)
    names = ["dln m/dln theta", "dln m/dln Km", "dln m/db", "dln m/dln ke"]
    print()
    print(f"    {'case':>22} " + " ".join(f"{n:>16}" for n in names))
    print("    " + "-" * 90)
    ok = True
    for amp, tau, lab in ((0.222, 0.0, "22.2%, ideal step"),
                          (0.222, 6 / 24, "22.2%, 6 h onset"),
                          (0.410, 1.0, "41%, 1 day onset")):
        tt = np.linspace(0.0, max(1.0, 4.0 * tau), 25)
        R, _ = burst_rows_ramp(300.0, tt, amp, tau, assay_cv=1.0)

        def lvl(th, Km, ke, Lt):
            m0 = float(U14.qss(th, TROUGH, Lt, Km, ke))
            J0 = (th + TROUGH) * Km / (Km + m0) ** 2 + ke

            def rhs(t, y):
                m = max(float(y[0]), 1e-14)
                thx = th + TROUGH * np.exp(-LAMBDA * t)
                f = 1.0 if tau <= 0 else (1.0 - np.exp(-t / tau))
                return [Lt * (1.0 - amp * f)
                        - thx * m / (Km + m) - ke * m]

            s = solve_ivp(rhs, (0.0, float(max(tt)) + 1e-9), [m0],
                          dense_output=True, rtol=1e-12, atol=1e-15,
                          method="LSODA")
            del J0
            return np.log(s.sol(tt)[0])

        #  h = 1e-4 rather than 1e-6: at rtol 1e-12 the integrator's own
        #  noise divided by 2h is the floor, so too small a step is as bad
        #  as too large. 1e-4 puts truncation near 1e-8 and roundoff near
        #  5e-9. The load step is relative for the same reason.
        h, hL = 1e-4, 1e-4 * LB
        fd = np.stack([
            (lvl(300.0 * np.exp(h), KM, KE, LB)
             - lvl(300.0 * np.exp(-h), KM, KE, LB)) / (2 * h),
            (lvl(300.0, KM * np.exp(h), KE, LB)
             - lvl(300.0, KM * np.exp(-h), KE, LB)) / (2 * h),
            (lvl(300.0, KM, KE, LB + hL)
             - lvl(300.0, KM, KE, LB - hL)) / (2 * hL),
            (lvl(300.0, KM, KE * np.exp(h), LB)
             - lvl(300.0, KM, KE * np.exp(-h), LB)) / (2 * h)], axis=1)
        cells = []
        for j in range(4):
            sc = max(np.max(np.abs(R[:, j])), np.max(np.abs(fd[:, j])), 1e-30)
            rel = float(np.max(np.abs(R[:, j] - fd[:, j])) / sc)
            ok &= rel < 1e-6
            cells.append(f"{rel:.2e}")
        print(f"    {lab:>22} " + " ".join(f"{c:>16}" for c in cells))
    print()
    print(f"  max relative disagreement under 1e-6 in every column: "
          f"{'YES' if ok else 'NO'}")
    print("  The finite difference perturbs the PRE-SUPPRESSION total load,")
    print("  so the b column agreeing is exactly the statement that the")
    print("  bias is suppressed along with the load it is a bias on.")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
