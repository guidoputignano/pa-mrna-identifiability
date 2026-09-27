"""
The elective glucose-protected withdrawal, and the load-record term that may
not apply during it.

THE MANOEUVRE. 01_biology.tex:160 gives the standard emergency regimen as
intravenous glucose to suppress catabolism, together with protein withdrawal.
Run ELECTIVELY rather than in response to illness, that is a load reduction
with three properties no earlier candidate had all of:

  fast onset          glucose acts on a circulatory and metabolic timescale of
                      hours, not on the bacterial-turnover timescale that makes
                      metronidazole slow (R118).
  no catabolic penalty  the glucose is there precisely to suppress catabolism,
                      which is the objection that withdrew the fast (R116).
  a MEASURED load     the load is delivered intravenously and recorded from the
                      infusion, not reconstructed from a diet diary.

THE THIRD PROPERTY IS THE ONE WORTH THE MOST, and it is the point of this file.
LOAD_FRAC, the 5% load-record error, is half of sigma at the nominal footing.
It enters the burst calculation of u15 in exactly one place: the prior on the
load-bias parameter b, which R119 showed is the near-collinear direction that
absorbs the k_e signal. If the burst's load is measured rather than recalled,
that nuisance is not the diet record's and is not 5% wide. That changes sigma,
which is the currency the whole requirement frontier is stated in, rather than
merely changing the amplitude.

WHY b MUST BECOME TWO PARAMETERS. u15 gives the burst and the record a SINGLE
load-bias b, which is correct while both loads come from the same diet record.
Under an infused load they are different quantities: the diet record's bias
still applies to the 118 quasi-steady visits and says nothing about what a pump
delivered. So the design carries b_rec, with its 5% prior, and b_burst, whose
prior width is swept from the same 5% down to exactly known.

THE AMPLITUDE, from the register, by an UNBOUNDED search of src/common/params.py
for the partition (5 matches, all read). The s_base row:

    "amino-acid catabolism a mean of 51.7% (range 24.5 to 66.4%), gut flora at
     least 22.2 +/- 8.5% by the fall under metronidazole, and roughly 30%
     unaccounted for and attributed to odd-chain fatty acids; A PROTEIN-
     RESTRICTED DIET LOWERS THE FIRST TWO DELIBERATELY."

So the register itself says the two shares are separately reducible, which is
the additivity this file prices. Glucose with protein withdrawal targets the
first, 51.7%; an osmotic bowel preparation removes the fermenting biomass and
its substrate on a transit timescale and targets the second, 22.2%. Together
73.9%. Each is an UPPER BOUND on what that arm can remove, since obligate
protein turnover and residual flora are not zero, so every amplitude is swept
below its anchor rather than assumed reached.

FOOTING as audit/u15_gut_suppression.py, whose corrected burst sensitivities
this imports. PART F is the anchor: with b shared, this file must reproduce u15
exactly.

Nothing here records an expected result in advance.

Run:  python audit/u17_glucose_burst.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import u14_dynamic_burst as U14
import u15_gut_suppression as U15

KM, KE, LB = U14.KM, U14.KE, U14.LB
HOUR, TROUGH = U14.HOUR, U14.TROUGH
TIER1_N, SEEDS = U14.TIER1_N, U14.SEEDS
FRONTIER = U14.FRONTIER
KE_PRIOR = U14.KE_PRIOR
PUBLISHED_KAPPA_CV = U14.PUBLISHED_KAPPA_CV
SB = U14.LOAD_FRAC * LB          # the 5% load-record error, in flux units

#: src/common/params.py, s_base row, read in full
CATAB_SHARE, CATAB_RANGE = 0.517, (0.245, 0.664)
GUT_SHARE, GUT_SD = 0.222, 0.085
BOTH = CATAB_SHARE + GUT_SHARE


def fisher_split(theta, seed, times, amp, tau_in, burst_cv, n_bursts,
                 f_burst, n_rec=TIER1_N, share_b=False):
    """
    Joint information with the record's load bias and the BURST's load bias as
    separate parameters.

    Order: (ln theta, ln Km, b_rec, ln ke, ln kappa, b_burst).

    f_burst is the burst load-bias prior as a fraction of the diet record's
    5%. f_burst = 0 means an exactly known infused load, and the column is
    dropped rather than given an infinite precision.

    share_b = True collapses b_burst back onto b_rec, which is u15's model and
    is the anchor PART F checks.
    """
    Q, m, inc = U14.record_rows(theta, n_rec, seed)
    g = m / (KM + m)
    sig = U14.sigma_flux(theta)
    #  record: (lnth, lnKm, b_rec, lnke, lnkappa), nothing about the infusion
    Qr = np.column_stack([Q, inc * g])
    R, _ = U15.burst_rows_ramp(theta, times, amp, tau_in, burst_cv)
    #  burst: (lnth, lnKm, b_BURST, lnke), and no kappa sensitivity worth
    #  carrying, as R110 established for a burst row
    drop = share_b or f_burst == 0.0
    npar = 5 if drop else 6
    Qr = np.pad(Qr, ((0, 0), (0, npar - 5)))
    Rb = np.zeros((R.shape[0], npar))
    Rb[:, 0], Rb[:, 1], Rb[:, 3] = R[:, 0], R[:, 1], R[:, 3]
    if share_b:
        Rb[:, 2] = R[:, 2]              # u15: one shared load bias
    elif f_burst > 0.0:
        Rb[:, 5] = R[:, 2]              # its own, with its own prior
    #  f_burst == 0 and not share_b: the infused load is known, no column
    pri = [0.0, 0.0, 1.0 / SB ** 2, 1.0 / KE_PRIOR ** 2,
           1.0 / PUBLISHED_KAPPA_CV ** 2]
    if npar == 6:
        pri.append(1.0 / (f_burst * SB) ** 2)
    F = np.diag(pri) + Qr.T @ Qr / sig ** 2 + n_bursts * (Rb.T @ Rb)
    return F


def sd_ke(theta, times, amp, tau_in, burst_cv=U14.ASSAY_CV, n_bursts=1,
          f_burst=1.0, share_b=False):
    """Median over record seeds of sd(ln k_e)."""
    v = []
    for s in SEEDS:
        F = fisher_split(theta, s, times, amp, tau_in, burst_cv, n_bursts,
                         f_burst, share_b=share_b)
        v.append(U14.sd_of(F, 3))
    a = np.array(v, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else np.inf


#: An admission is an admission. u15 let the burst window grow as four onset
#: time constants, which silently hands a SLOW input a huge free window: at a
#: one-week onset that is 28 days and 24 draws spanning two dose cycles, which
#: is a record extension at a better assay and not a relaxation measurement at
#: all. R118 already noticed the symptom, that the slow end sat under the
#: record-only floor instead of converging to it. Fixing the window at one
#: inpatient day is what makes the onset sweep measure the onset.
WINDOW_D = 1.0


def win(tau, span=None):
    """Sample times over a FIXED admission window, 24 draws."""
    return np.linspace(0.0, WINDOW_D if span is None else span, 25)


def slowest_onset(theta, amp, burst_cv, n_bursts, f_burst, hi=14.0,
                  fixed_window=True):
    """Slowest onset time constant that still reaches the frontier, in hours."""
    fr = FRONTIER[theta]
    if fr is None:
        return None
    span = (lambda t: WINDOW_D) if fixed_window else \
        (lambda t: max(WINDOW_D, 4.0 * t))
    if sd_ke(theta, win(0.0, span(0.0)), amp, 0.0, burst_cv,
             n_bursts, f_burst) > fr:
        return float("nan")
    lo = 0.0
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        v = sd_ke(theta, win(mid, span(mid)), amp, mid,
                  burst_cv, n_bursts, f_burst)
        lo, hi = (mid, hi) if v <= fr else (lo, mid)
    out = 24.0 * 0.5 * (lo + hi)
    return float("inf") if out > 24.0 * 13.0 else out


# ======================================================================
def part_a():
    print("=" * 88)
    print("PART A.  the amplitude anchors, read from the register in full")
    print("=" * 88)
    print(f"""
  src/common/params.py, s_base row, quoted rather than paraphrased:

    "amino-acid catabolism a mean of 51.7% (range 24.5 to 66.4%), gut flora
     at least 22.2 +/- 8.5% by the fall under metronidazole, and roughly 30%
     unaccounted for ...; a protein-restricted diet lowers the first two
     deliberately."

  So the register states the additivity rather than this file assuming it.

    arm                                        anchor        what it removes
    ----------------------------------------------------------------------
    glucose + protein withdrawal              {CATAB_SHARE:6.1%}  amino-acid catabolism
      its range                        {CATAB_RANGE[0]:5.1%} to {CATAB_RANGE[1]:.1%}
    osmotic bowel preparation                 {GUT_SHARE:6.1%}  gut flora
      plus or minus one sd                     {GUT_SD:6.1%}
    both, the register's "first two"          {BOTH:6.1%}

  Each anchor is an UPPER BOUND on its arm: obligate protein turnover and
  residual flora are not zero. Every table below sweeps under its anchor.
""")


def part_b():
    print("=" * 88)
    print("PART B.  the prize: the burst's load-record term, swept to zero")
    print("=" * 88)
    print("\n  24 hourly draws, a 6 h onset, Tier 1 record, one admission at a")
    print("  5% assay. f is the burst load-bias prior as a fraction of the diet")
    print("  record's 5%; f = 0 is an infused load known exactly.\n")
    tt = np.linspace(0.0, 1.0, 25)
    fs = (1.0, 0.5, 0.2, 0.1, 0.0)
    print(f"    {'arm':>30} {'theta':>6}" +
          "".join(f"{f'f={x:g}':>9}" for x in fs) + f"{'  f=0 vs f=1':>13}")
    print("    " + "-" * 92)
    for lab, amp in (("bowel prep alone", GUT_SHARE),
                     ("glucose + withdrawal", CATAB_SHARE),
                     ("both", BOTH)):
        for theta in (600.0, 300.0):
            vals = [sd_ke(theta, tt, amp, 6 / 24, f_burst=x) for x in fs]
            cells = "".join(f"{v:9.4f}" for v in vals)
            print(f"    {lab:>30} {theta:6.0f}{cells}"
                  f"{vals[0] / vals[-1]:12.2f}x")
    print()
    print(f"  Record-only floors: {U15.qss_floor(600.0):.4f} at 600, "
          f"{U15.qss_floor(300.0):.4f} at 300.")
    print("  Frontiers: 0.1567 at 600, 0.0633 at 300.")


def part_c():
    print()
    print("=" * 88)
    print("PART C.  the best case, and whether the margin now clears the spread")
    print("=" * 88)
    print("\n  A 1% burst assay and sixteen admissions, 6 h onset. The ratio to")
    print("  the frontier is what matters: R124 put the method spread between")
    print("  u15's Fisher and u16's Monte Carlo at 9% to 33%, so a margin of")
    print("  less than about 1.33x is not a result.\n")
    tt = np.linspace(0.0, 1.0, 25)
    print(f"    {'arm':>30} {'theta':>6} {'f=1':>9} {'f=0':>9} "
          f"{'frontier':>9} {'margin at f=0':>14}")
    print("    " + "-" * 82)
    for lab, amp in (("bowel prep alone", GUT_SHARE),
                     ("glucose + withdrawal", CATAB_SHARE),
                     ("both", BOTH)):
        for theta in (600.0, 300.0):
            a = sd_ke(theta, tt, amp, 6 / 24, 0.01, 16, 1.0)
            b = sd_ke(theta, tt, amp, 6 / 24, 0.01, 16, 0.0)
            fr = FRONTIER[theta]
            mark = "" if fr / b < 1.33 else " *"
            print(f"    {lab:>30} {theta:6.0f} {a:9.4f} {b:9.4f} "
                  f"{fr:9.4f} {f'{fr / b:.2f}x{mark}':>14}")
    print()
    print("  * margin exceeds the method spread, so the cell is a result and")
    print("    not merely a number that lands on the right side of a line.")


def part_d():
    print()
    print("=" * 88)
    print("PART D.  what the measured load does to the required onset")
    print("=" * 88)
    print("\n  The slowest onset time constant that still reaches the frontier,")
    print("  at a 1% assay and sixteen admissions. This is the number that says")
    print("  which manoeuvre is fast enough.\n")
    print("  The window is a FIXED one-day admission, which is the correction")
    print("  this file makes to u15: letting the window grow with the onset")
    print("  hands a slow input a free record extension (see WINDOW_D).\n")
    fmt = (lambda x: "never" if x is None or (isinstance(x, float) and np.isnan(x))
           else ("beyond 13 d" if np.isinf(x) else f"{x:.1f} h"))
    print(f"    {'arm':>30} {'theta':>6} {'fixed 1 d window':>18} "
          f"{'u15 growing window':>20}")
    print("    " + "-" * 78)
    for lab, amp in (("bowel prep alone", GUT_SHARE),
                     ("glucose + withdrawal", CATAB_SHARE),
                     ("both", BOTH)):
        for theta in (600.0, 300.0):
            a = slowest_onset(theta, amp, 0.01, 16, 0.0, fixed_window=True)
            b = slowest_onset(theta, amp, 0.01, 16, 0.0, fixed_window=False)
            print(f"    {lab:>30} {theta:6.0f} {fmt(a):>18} {fmt(b):>20}")
    print()
    print("  The two columns differ by more than an order of magnitude where")
    print("  the growing window runs to weeks, which is the artefact, not a")
    print("  result. The left column is the requirement on the intervention.")
    print()
    print("  For scale: glucose acts over hours, an osmotic preparation over a")
    print("  gut transit, and metronidazole over days (R118), which is a")
    print("  property of bacterial turnover and not of the route of")
    print("  administration, so no formulation change reaches the top rows.")


def part_g():
    print()
    print("=" * 88)
    print("PART G.  where the load-record prize actually is, since PART B says")
    print("         it is not in the burst")
    print("=" * 88)
    print("""
  PART B found removing the burst's load-record term worth 1.00x to 1.02x.
  The diagnosis is structural rather than numerical: LOAD_FRAC enters the
  burst in ONE place, the prior on its load bias, because a burst observation
  is a MARKER LEVEL with assay noise and carries no load measurement at all.
  It is the RECORD whose sigma is hypot(load term, assay term), and the
  record is still the ordinary outpatient diet record whatever the admission
  does. So the premise that the 5% error is half of sigma is right, and it is
  half of the RECORD's sigma, which an infused burst does not touch.

  The adjacent question is then worth asking: what WOULD the load-record term
  be worth if it were removed where it actually lives?
""")
    tt = win(6 / 24)
    print(f"    {'design':>44} {'theta=600':>10} {'theta=300':>10}")
    print("    " + "-" * 68)
    #  "measured exactly" means b is KNOWN, so its prior precision goes to
    #  INFINITY. Setting the precision to zero instead makes b a free
    #  parameter with no prior, which is less information and not more, and
    #  it showed up immediately as a frac = 0 row worse than frac = 0.2.
    #  A tiny width is used rather than dropping the column, so that the
    #  ln k_e index stays put, and two of them are shown to prove it is
    #  converged rather than merely small.
    for lab, frac in (("as published, 5% load-record error", 1.0),
                      ("halved, 2.5%", 0.5),
                      ("a fifth, 1%", 0.2),
                      ("a hundredth, 0.05%", 0.01),
                      ("known to 1e-4 of it", 1e-4),
                      ("known to 1e-6 of it (convergence check)", 1e-6)):
        cells = []
        for theta in (600.0, 300.0):
            v = []
            for s in SEEDS:
                Q, m, inc = U14.record_rows(theta, TIER1_N, s)
                g = m / (KM + m)
                #  sigma_flux is hypot(load term, assay term). Scale only the
                #  load term, rebuilding the assay term from its definition
                #  rather than subtracting under a square root.
                mbar = float(U14.qss(theta, 0.0, LB))
                slope = theta * KM / (KM + mbar) ** 2 + KE
                sig = float(np.hypot(frac * SB,
                                     slope * U14.ASSAY_CV * mbar))
                Qr = np.column_stack([Q, inc * g])
                R, _ = U15.burst_rows_ramp(theta, tt, CATAB_SHARE, 6 / 24, 0.01)
                Rb = np.zeros((R.shape[0], 5))
                Rb[:, 0], Rb[:, 1], Rb[:, 3] = R[:, 0], R[:, 1], R[:, 3]
                pri = [0.0, 0.0, 1.0 / (frac * SB) ** 2,
                       1.0 / KE_PRIOR ** 2, 1.0 / PUBLISHED_KAPPA_CV ** 2]
                F = np.diag(pri) + Qr.T @ Qr / sig ** 2 + 16 * (Rb.T @ Rb)
                v.append(U14.sd_of(F, 3))
            a = np.array(v, float)
            a = a[np.isfinite(a)]
            cells.append(f"{np.median(a):10.4f}" if a.size else f"{'inf':>10}")
        print(f"    {lab:>44} " + " ".join(cells))
    print()
    print("  This is the glucose arm's burst held fixed while the RECORD's")
    print("  load term is swept, so the column shows what the 5% is worth")
    print("  where it lives: both in the record's per-visit sigma and in the")
    print("  width of the load-bias prior, which are the only two places it")
    print("  appears anywhere in the design.")


def part_h():
    print()
    print("=" * 88)
    print("PART H.  the amplitude has an INTERIOR optimum, so the two arms do")
    print("         not add")
    print("=" * 88)
    print("\n  PART E's column is not monotone. Located here on a fine sweep,")
    print("  1% assay, sixteen admissions, 6 h onset, infused load.\n")
    tt = win(6 / 24)
    grid = np.arange(0.20, 0.86, 0.05)
    print(f"    {'amplitude':>10} {'theta=600':>11} {'theta=300':>11}")
    print("    " + "-" * 36)
    best = {}
    for theta in (600.0, 300.0):
        vals = [sd_ke(theta, tt, a, 6 / 24, 0.01, 16, 0.0) for a in grid]
        best[theta] = (grid[int(np.argmin(vals))], min(vals))
    for i, a in enumerate(grid):
        row = []
        for theta in (600.0, 300.0):
            v = sd_ke(theta, tt, a, 6 / 24, 0.01, 16, 0.0)
            row.append(f"{v:11.4f}" + ("*" if abs(a - best[theta][0]) < 1e-9
                                       else " "))
        print(f"    {a:10.1%} " + "".join(row))
    print()
    for theta in (600.0, 300.0):
        print(f"    optimum at theta = {theta:.0f}: "
              f"{best[theta][0]:.0%}, giving {best[theta][1]:.4f}")
    print(f"""
  READING. Suppressing more load is not better. Below the optimum the marker
  falls far enough off enzyme saturation that the relaxation rate J rises and
  k_e's share of it falls, which is the same quantity R113 used to reject the
  post-dose placement. So the register's additivity, "a protein-restricted
  diet lowers the first two deliberately", is real as BIOLOGY and is the
  wrong thing to want as DESIGN: glucose plus withdrawal at {CATAB_SHARE:.1%} beats
  glucose plus withdrawal plus bowel preparation at {BOTH:.1%}. Adding the second
  arm makes the measurement worse.""")


def part_e():
    print()
    print("=" * 88)
    print("PART E.  the amplitude swept under each anchor, at the infused load")
    print("=" * 88)
    print("\n  Neither arm reaches its anchor in practice, so the question is")
    print("  how far under it the verdict survives. 1% assay, sixteen")
    print("  admissions, 6 h onset, f = 0.\n")
    tt = np.linspace(0.0, 1.0, 25)
    print(f"    {'amplitude':>10} {'theta=600':>10} {'margin':>8} "
          f"{'theta=300':>10} {'margin':>8}  what it is")
    print("    " + "-" * 76)
    marks = {CATAB_RANGE[0]: "catabolism range, low",
             GUT_SHARE: "gut share alone",
             CATAB_SHARE: "catabolism mean",
             CATAB_RANGE[1]: "catabolism range, high",
             BOTH: "both arms, the register's 'first two'"}
    for amp in sorted({0.10, 0.15, CATAB_RANGE[0], GUT_SHARE, 0.35,
                       CATAB_SHARE, CATAB_RANGE[1], BOTH}):
        row = []
        for theta in (600.0, 300.0):
            v = sd_ke(theta, tt, amp, 6 / 24, 0.01, 16, 0.0)
            row += [f"{v:10.4f}", f"{FRONTIER[theta] / v:7.2f}x"]
        print(f"    {amp:10.1%} " + " ".join(row) + f"  {marks.get(amp, '')}")
    print()
    print("  A margin above 1.33x clears R124's method spread.")


def part_f():
    print()
    print("=" * 88)
    print("PART F.  the anchor: with b shared, this file must reproduce u15")
    print("=" * 88)
    print("\n  If the split-bias construction is right, collapsing b_burst onto")
    print("  b_rec must return u15's own numbers exactly, not approximately.\n")
    worst = 0.0
    print(f"    {'case':>34} {'u15':>10} {'u17, b shared':>15} {'diff':>11}")
    print("    " + "-" * 74)
    for lab, amp, tau, cv, nb in (
            ("22.2%, 6 h, 5% assay, 1 adm", GUT_SHARE, 6 / 24, U14.ASSAY_CV, 1),
            ("22.2%, 6 h, 1% assay, 16 adm", GUT_SHARE, 6 / 24, 0.01, 16),
            ("51.7%, ideal step, 5%, 1 adm", CATAB_SHARE, 0.0, U14.ASSAY_CV, 1),
            ("73.9%, 1 day, 1% assay, 4 adm", BOTH, 1.0, 0.01, 4)):
        tt = np.linspace(0.0, max(1.0, 4.0 * tau), 25)
        a = U15.sd_ke(300.0, tt, amp, tau, burst_cv=cv, n_bursts=nb)
        b = sd_ke(300.0, tt, amp, tau, cv, nb, 1.0, share_b=True)
        worst = max(worst, abs(a - b))
        print(f"    {lab:>34} {a:10.6f} {b:15.6f} {abs(a - b):11.2e}")
    print()
    print(f"  largest disagreement {worst:.2e}, "
          f"{'PASS' if worst < 1e-9 else 'FAIL'}")
    print("  The remaining difference between the two files is therefore the")
    print("  one modelling change under test, and nothing else.")


if __name__ == "__main__":
    part_a()
    part_f()
    part_b()
    part_c()
    part_d()
    part_e()
    part_g()
    part_h()
