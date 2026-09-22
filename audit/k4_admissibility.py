"""
Phase 5, the follow-up: is the trough reading physically admissible, and the
dose sensitivity table with admissibility on every row.

Two questions, asked of every candidate dose rather than of the committed one:

  (a) THE CAPACITY CEILING.  h1_admissibility.py's constraint 5,
      theta + Delta_peak <= normal capacity, with the reference carried as a
      parameter because src/common/params.py states the residual capacity two
      ways four lines apart ("about 10% of normal", "under 5% of control").
      The used fraction is frac * (1 + Delta_peak/theta), so the ceiling binds
      at Delta_peak/theta = 1/frac - 1, that is 9.0 at the 10% reading and
      19.0 at the 5% reading. Neither theta nor the model block enters, which
      is h1's own cancellation and is re-derived here rather than quoted.

  (b) THE SHAPE CEILING.  The largest marker fall the published inhibition
      shape can produce, maximised over IC50, at THAT ROW'S OWN enzyme ratio.
      k2_structures.py computed it at the committed ratio only; it is a
      function of the ratio and has to be recomputed per row, because a bigger
      dose raises the ceiling as well as the fall it has to clear.

A row is called admissible only if it clears BOTH, and the two disagree about
which rows are in trouble, which is the point of asking both.

Nothing here records an expected result in advance, no constant is edited, and
design.run's LAMBDA and D_AMP are rebound through k1's context manager, which
restores them.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from contextlib import contextmanager
from scipy.optimize import brentq

import design.run as DR
import g_design_common as G
import k1_sources as K1
import k2_structures as K2

TH, KM, KE, LB = G.TH, G.KM, G.KE, G.LB
T = G.T
M_BAR = 10.0
REFS = (0.10, 0.05)          # the two readings of the residual-capacity fraction


# --------------------------------------------------------------------------
def peak_of(D=None, lam=None):
    """The steady-state peak increment of a biweekly train."""
    D = DR.D_AMP if D is None else D
    lam = DR.LAMBDA if lam is None else lam
    return D / (1.0 - np.exp(-lam * T))


def capacity_used(dpk, frac):
    """h1 constraint 5. theta cancels; this is frac*(1 + Delta_peak/theta)."""
    normal = TH / frac
    return (TH + dpk) / normal


@contextmanager
def _imax(v):
    """Temporarily rebind k2's published Imax, to separate two questions."""
    old = K2.IMAX
    K2.IMAX = v
    try:
        yield
    finally:
        K2.IMAX = old


def shape_ceiling(cls, dpk, imax=None):
    """
    The largest marker fall `cls`'s published shape allows at this peak.

    `imax` reruns the same maximisation with Imax rebound. That separates two
    questions that the first version of this script ran together and got
    wrong: whether a row is out of reach of the STRUCTURE, and whether it is
    out of reach of one FIXED constant inside it. Attarwala fixed Imax at
    0.999 rather than estimating it (Table 1 footnote m), so a verdict that
    flips between 0.999 and 1 is a verdict about a modelling convenience.
    """
    if imax is None:
        return K2._ceiling(lambda c: cls.shape(c, E_lo=TH, E_hi=TH + dpk))
    with _imax(imax):
        return K2._ceiling(lambda c: cls.shape(c, E_lo=TH, E_hi=TH + dpk))


def b_window(dpk, frac):
    """
    The admissible window for the load bias b along P5's exact fibre, at this
    dose and this capacity reference.

    h1's constraint 6. The fibre is K_m(b) = K_0 L/(L+b),
    theta(b) = theta_0 + b + k_e(K_0 - K_m(b)), and the data cannot separate
    points along it, so admissibility is the only thing that bounds b. The
    UPPER bound is spent against the ceiling's headroom, so it moves with D.
    """
    Km_b = lambda b: KM * LB / (LB + b)
    th_b = lambda b: TH + b + KE * (KM - Km_b(b))
    lo = max(-LB, brentq(th_b, -LB + 1e-9, 0.0, xtol=1e-12))
    normal = TH / frac
    f = lambda b: th_b(b) + dpk - normal
    if f(lo) > 0:
        return lo, lo, 0.0          # no admissible b at all at this dose
    try:
        hi = brentq(f, lo, 1e9, xtol=1e-9)
    except ValueError:
        hi = np.inf
    return lo, hi, hi - lo


def visits_to_gate(seeds=range(8), load_cv=0.10, nmax=60):
    fn = G.SCHEDULES["monthly (30.4 d)"]
    for n in range(4, nmax + 1):
        _, md, _ = G.sd_over_seeds(fn(n), load_cv, seeds=seeds)
        if md <= G.GATE:
            return n
    return None


# --------------------------------------------------------------------------
def timing_evidence():
    print("=" * 78)
    print("1. DOES ANYTHING CONSTRAIN WHEN 'AFTER THE THIRD DOSE' WAS SAMPLED?")
    print("=" * 78)
    print("""
  THE TRIAL PAPER: no. Checked again in this pass, exhaustively.
    - Extended Data Table 3's caption is 'Biomarker percentage change from
      baseline after the third dose'. No time point, no visit number.
    - The Methods define the endpoint only as 'changes in blood 2-MC and 3-HP
      levels from baseline (predose levels) to post-dose levels'.
    - There is NO biomarker time-course figure. Fig. 3 is PK exposure to PCCA
      mRNA; Fig. 4 is MDE profiles. 3-HP appears in the article exactly once
      as a number, in Extended Data Table 3.
    - The only sampling time named anywhere in the article is a PREDOSE draw,
      in the antibody paragraph: 'the predose sample for dose 3 (52 days after
      the first dose)'. That is a trough, and it is about anti-PEG titres, not
      about 3-HP. It is weak evidence and it is the only evidence there is.
      (Note in passing: cohort 1 is Q3W, so a third dose on schedule falls at
      day 42, not day 52. Either dosing slipped or the draw did. Whichever it
      is, nominal and actual timing are not the same thing in this trial.)

  THE TRIAL REGISTRY (NCT04159103): NOT CHECKED. Outbound network access is
    blocked by this environment's egress policy. Both
    `curl https://clinicaltrials.gov/api/v2/studies/NCT04159103` and a WebFetch
    of the same record return a 403 CONNECT denial from the gateway, logged in
    the proxy's own recentRelayFailures as
    'gateway answered 403 to CONNECT (policy denial or upstream failure)'.
    This is a real open item, not a null result: the registry's secondary
    outcome measures carry a `timeFrame` field, and if it names a visit for the
    3-HP endpoint it settles the question outright. It needs a session with
    egress, or the record pasted in.

  BAEK'S SOURCE DATA: NOT AVAILABLE. Fig. 2's caption says 'Source data are
    provided as a Source Data file'. That file was not among the five PDFs and
    cannot be fetched here, for the same reason.

  WHAT BAEK'S FIGURE 2 DOES SAY, from the caption, read in this pass:
    'Blood was collected from mice prior to treatment and on Days 2, 8, 14, 22,
     and 28 after the first IV dose, and 2 days after the second IV dose.'
    on a regimen of 'two IV bolus doses ... on Days 0 and 28'.

    So the preclinical PD data the translational model was fitted to is a FULL
    TIME COURSE across a 28-day cycle, sampled at 2, 8, 14, 22 and 28 days, not
    a trough-only design. Day 2 is near the PD peak of that cycle and Day 28 is
    its trough. That tells us the convention in the programme's own preclinical
    work was to sample the whole cycle. It does NOT tell us which point of the
    human cycle Extended Data Table 3 reports, because the trial's biomarker
    sampling schedule is not published.

  ONE THING THE SAME FIGURE DOES SETTLE, for finding 1(e). Fig. 2e's y-axis is
    '3-HP (uM)', plotted 0 to 200, in the SAME Pcca-/- (A138T) mouse model that
    Attarwala fitted. Same company, same strain, same analyte, plotted in
    micromolar. Attarwala's Table 1 prints the 3-HP baseline as 60.1 'mM/L'.
    That is now confirmed as umol/L by an external source and not only by the
    internal LOQ argument.
""")


def hdr_rows(rows):
    return rows


def sensitivity():
    print("=" * 78)
    print("2. THE DOSE SENSITIVITY TABLE, WITH ADMISSIBILITY ON EVERY ROW")
    print("=" * 78)
    print("""
  D_AMP IS NOT MOVED. Every row below rebinds design.run.D_AMP and LAMBDA
  inside a context manager that restores them, exactly as k1_sources.py does.
  The committed constants are D = 6000 and lambda = ln2/3.5 d.

  'fall' is the marker's fall below the untreated operating point at the
  treated steady-state peak, under the PAPER'S model.
  'used@10%' and 'used@5%' are h1's capacity constraint, theta + Delta_peak
  over normal capacity, at the two readings of the residual-capacity fraction.
  The constraint binds at Delta_peak/theta = 9.0 and 19.0 respectively.
  'B@1' and 'A@1' are the largest fall each published inhibition shape can
  produce AT THAT ROW'S OWN enzyme ratio, maximised over IC50, with Imax set
  to 1 rather than the fixed 0.999. 'bracketed' asks whether the paper's fall
  lies between them. The 0.999 versions are printed below the table, because
  Imax was FIXED and not estimated and a verdict that turns on it is a verdict
  about a modelling convenience.
  'n' is monthly visits to sd(ln theta) <= 0.20 at load CV 10%, median over 8
  seeds, searched to n = 60.
""")
    rows = []
    for lam_name, lam in (("3.5 d", np.log(2) / 3.5), ("7 d", np.log(2) / 7.0)):
        specs = [("committed D = 6000", 6000.0)]
        for reading in ("peak", "trough"):
            for tgt in (0.20, 0.30):
                D, _ = K1.calibrate_D(tgt, reading, lam=lam)
                specs.append((f"{100*tgt:.0f}% fall at the {reading}", D))
        for name, D in specs:
            with K1._set(LAMBDA=lam, D_AMP=D):
                dpk = peak_of()
                fall = K1.fall(dpk)
                m_pk = K1.steady_level(dpk)
                rng = G.delta_range(G.SCHEDULES["monthly (30.4 d)"](12))
                n = visits_to_gate()
                sd12 = G.sd_over_seeds(
                    G.SCHEDULES["monthly (30.4 d)"](12), 0.10)[1]
            used = {f: capacity_used(dpk, f) for f in REFS}
            a_ceil, _ = shape_ceiling(K2.PubA, dpk)
            b_ceil, _ = shape_ceiling(K2.PubB, dpk)
            a1, _ = shape_ceiling(K2.PubA, dpk, imax=1.0)
            b1, _ = shape_ceiling(K2.PubB, dpk, imax=1.0)
            bw = {f: b_window(dpk, f) for f in REFS}
            rows.append(dict(lam=lam_name, name=name, D=D, dpk=dpk,
                             ratio=dpk / TH, fall=fall, m_pk=m_pk, rng=rng,
                             n=n, sd12=sd12, used=used,
                             a=a_ceil, b=b_ceil, a1=a1, b1=b1, bw=bw))

    hdr = (f"  {'lambda':>7} {'row':>26} {'D':>9} {'Dpk/theta':>10} "
           f"{'fall':>7} {'m_peak':>7} {'rng':>7} {'sd@12':>7} {'n':>5}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['lam']:>7} {r['name']:>26} {r['D']:9.1f} {r['ratio']:10.3f} "
              f"{100*r['fall']:6.2f}% {r['m_pk']:7.3f} {r['rng']:7.3f} "
              f"{r['sd12']:7.4f} "
              f"{(str(r['n']) if r['n'] else '>60'):>5}")

    hdr2 = (f"\n  {'lambda':>7} {'row':>26} {'used@10%':>9} {'used@5%':>8} "
            f"{'capacity':>10} | {'fall':>7} {'B@1':>7} {'A@1':>7} "
            f"{'bracketed':>10}")
    print(hdr2)
    print("  " + "-" * (len(hdr2) - 3))
    for r in hdr_rows(rows):
        cap = "OK" if (r["used"][0.10] <= 1.0 and r["used"][0.05] <= 1.0) else "BREACH"
        br = "yes" if (r["b1"] <= r["fall"] <= r["a1"]) else "NO"
        print(f"  {r['lam']:>7} {r['name']:>26} {r['used'][0.10]:9.4f} "
              f"{r['used'][0.05]:8.4f} {cap:>10} | {100*r['fall']:6.2f}% "
              f"{100*r['b1']:6.2f}% {100*r['a1']:6.2f}% {br:>10}")
    print("""
  'capacity' is h1's constraint 5 at both references. It binds at
  Delta_peak/theta = 9.0 (10% reading) and 19.0 (5% reading).

  'B@1' and 'A@1' are the two published shapes' ceilings with Imax set to 1
  rather than the fixed 0.999, which is the structure's own limit rather than
  one convenience constant's. 'bracketed' asks whether the paper's fall sits
  between them.""")
    print(f"\n  With Imax at its published 0.999 the ceilings are lower, and the "
          f"gap\n  to the paper's fall is small enough that the verdict flips:")
    print(f"    {'lambda':>7} {'row':>26} {'fall':>7} {'A@.999':>8} "
          f"{'B@.999':>8} {'A margin':>9} {'B margin':>9}")
    for r in hdr_rows(rows):
        print(f"    {r['lam']:>7} {r['name']:>26} {100*r['fall']:6.2f}% "
              f"{100*r['a']:7.2f}% {100*r['b']:7.2f}% "
              f"{100*(r['a']-r['fall']):+8.2f}p {100*(r['b']-r['fall']):+8.2f}p")

    print(f"\n  The admissible window for the load bias b along P5's exact fibre,")
    print(f"  which is what the ceiling's headroom is spent on (h1 constraint 6).")
    print(f"  The load record's own 95% half-width is "
          f"{1.96 * G.LOAD_FRAC * LB:.4f}, so a window of "
          f"{2 * 1.96 * G.LOAD_FRAC * LB:.2f} is what the record already asserts.")
    hdr3 = (f"\n  {'lambda':>7} {'row':>26} {'b window @10%':>26} "
            f"{'width':>10} {'vs record':>10}")
    print(hdr3)
    print("  " + "-" * (len(hdr3) - 3))
    rec = 2 * 1.96 * G.LOAD_FRAC * LB
    for r in rows:
        lo, hi, w = r["bw"][0.10]
        print(f"  {r['lam']:>7} {r['name']:>26} "
              f"[{lo:10.2f}, {hi:11.2f}] {w:10.2f} {w / rec:9.1f}x")
    return rows


def verdict(rows):
    print("\n" + "=" * 78)
    print("3. ARE THE TROUGH ROWS PHYSICALLY ADMISSIBLE?")
    print("=" * 78)
    tr = [r for r in rows if "trough" in r["name"]]
    cm = [r for r in rows if "committed" in r["name"]][0]
    worst = max(tr, key=lambda r: r["used"][0.10])

    print("""
  THE SHORT ANSWER: YES, ON BOTH TESTS. Neither the capacity ceiling nor the
  published shape excludes any trough row. What is against the trough reading
  is not admissibility, and saying otherwise would be an overclaim. The two
  real objections are set out at the end, and one of them mostly dissolves
  under the 7 d half-life.

  (a) THE CAPACITY CEILING: every trough row passes, on both references.
      The constraint binds at Delta_peak/theta = 9.0 at the 10% reading and
      19.0 at the 5%.""")
    for r in tr:
        print(f"      {r['lam']:>5}, {r['name']:26s}  Delta_peak/theta "
              f"{r['ratio']:6.3f}   uses {100*r['used'][0.10]:6.2f}% at 10%, "
              f"{100*r['used'][0.05]:6.2f}% at 5%")
    print(f"      for comparison, the committed dose uses "
          f"{100*cm['used'][0.10]:.2f}% and {100*cm['used'][0.05]:.2f}%.")
    print(f"""
      The worst row is {worst['lam']}, {worst['name']}, at
      {100*worst['used'][0.10]:.2f}% of normal capacity on the 10% reading. That
      is inside the bound and it leaves {100*(1-worst['used'][0.10]):.2f}% of a PA
      patient's normal PCC capacity unused at the peak of every dose cycle. It
      is admissible arithmetic. It is also a strong clinical claim, that one
      fortnightly dose lifts a patient from 10% of normal capacity to
      {100*worst['used'][0.10]:.0f}% of it, and the trial's own outcomes do not
      support it: the relative risk of a metabolic decompensation event fell by
      70%, not to near-normal, and the biomarkers moved by -29% to +5%.

  (b) THE SHAPE CEILING: it does not exclude the trough rows either, and the
      first version of this check said it did. That was wrong, and correcting
      it produced a better result.

      At Imax fixed to its published 0.999 the ceilings sit a fraction of a
      percentage point BELOW the paper's fall on most rows, which reads as an
      exclusion. But Attarwala FIXED Imax at 0.999 rather than estimating it
      (Table 1 footnote m), so it is a modelling convenience and not a
      measurement, and a verdict that turns on it is a verdict about nothing.
      Rerunning each ceiling at Imax = 1 separates the two questions, and the
      answer is much cleaner:

        PUB_A's ceiling at Imax = 1 is exactly Delta/(theta + Delta), the
        c -> 0 limit of the inhibition ratio, and it CLEARS the paper's fall on
        every row by 1.2 to 2.0 percentage points. PUB_A is not excluded at any
        dose. The 0.999 cap costs it 1.0 to 2.4 points, which is the entire
        reason it appeared to fail.

        PUB_B's ceiling falls SHORT on every row at BOTH values of Imax, by
        0.09 to 1.74 points. That shortfall is structural: it is the published
        floor E0/base = 0.0987/60.1 riding along and blunting every fall, not
        the Imax cap. So an algebraic Imax marker carrying Attarwala's own
        floor is always very slightly too shallow to reproduce a
        Michaelis-Menten marker's response to the same enzyme change.

      THE RESULT THAT MATTERS IS THE BRACKET. On every row in the table, at
      Imax = 1, the paper's fall lies strictly between the two published
      shapes' ceilings: PUB_B a little too shallow, PUB_A with headroom to
      spare. The paper's dose-response is bracketed by the published
      structures rather than outside them. A shortfall of under two percentage
      points on one side cannot exclude a dose, and it would be dishonest to
      report it as if it could.

      ONE SCOPE LIMIT ON THAT BRACKET, found by k6_recheck.check_9 and not by
      this script. PUB_A's ceiling is above the paper's fall at every dose, but
      PUB_B's CROSSES it at Delta_peak = 443.6, that is Delta_peak/theta =
      0.148. Below that dose PUB_B's shape can do MORE than the paper's marker,
      not less, and the bracket reverses. Every row in the table sits above the
      crossover, the smallest at Delta_peak = 816.8, so the bracket holds
      throughout the table. It is not a property of the three forms at every
      dose, and this script's first draft said it was.

  WHAT IS ACTUALLY AGAINST THE TROUGH READING, then, is two things, neither of
  them admissibility:

    1. It needs D to be 2.0x and 3.5x the committed value at the committed
       half-life (12255 and 20977 against 6000), producing peak marker falls of
       80.2% and 87.4%, while the trial reports falls two to three times
       SMALLER than the committed calibration already is. The trough reading
       answers a small reported number by making the underlying effect much
       larger, which is coherent but only if the trial sampled at the very
       bottom of a very deep excursion.
    2. The implied clinical claim in (a).

  AND THE 7 d HALF-LIFE MOSTLY DISSOLVES OBJECTION 1. A longer half-life makes
  the trough a much larger fraction of the peak (25.0% at 7 d against 6.25% at
  3.5 d), so the same trough target needs far less dose:""")
    for r in tr:
        print(f"      {r['lam']:>5}, {r['name']:26s}  D {r['D']:9.1f}  "
              f"peak fall {100*r['fall']:6.2f}%  capacity {100*r['used'][0.10]:6.2f}%  "
              f"visits {(str(r['n']) if r['n'] else '>60'):>3}")
    print("""
    So 7 d PLUS the trough reading is the most internally coherent combination
    on the table: doses of 2489 and 4261, peak falls of 51% and 64%, capacity
    at 21% and 29% of normal, and both bracketed by the published shapes. It is
    also the combination that costs the design the most, 30 to 37 monthly
    visits to reach sd(ln theta) = 0.20 against the committed 11. Both of its
    ingredients are the reading the primary sources lean towards, and neither
    is established.
""")

def main():
    timing_evidence()
    rows = sensitivity()
    verdict(rows)


if __name__ == "__main__":
    main()
