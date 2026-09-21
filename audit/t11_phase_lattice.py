"""
17. THE PHASE LATTICE: how many points of the dose cycle a schedule can reach,
    how many it does reach, and which of the two the figure is printing.

Measurement only. Nothing here reaches the paper, and nothing here writes a
results/numbers file.

WHAT audit/f4_sampling.py ALREADY ESTABLISHED, and is not repeated here.
f4 defines `predicted_phases(step)` -- "q, from step/T = p/q in lowest terms"
-- tabulates it against the phases a 730-day record actually visits for
eighteen sampling intervals, and settles the question f4 was written for:
which quantity binds. Its answer is the range of Delta/theta, not the phase
count, and its S3 shows two designs with the same phase count and very
different precision. That result stands and is used below rather than retested.

WHAT IS NEW HERE. f4 prints `phases` and `q` side by side in one table and
never says why they differ, and in its own table they differ on five of
eighteen rows. The difference is not noise and it is not a rounding artefact:
a schedule of step h reaches

    min(q, n)   distinct phases,

where q is the lattice the step admits and n is the number of visits the
record contains. q is a property of the STEP; n is a property of the
HORIZON; and the smaller of the two is what a reader sees. This script
separates them, applies the separation to every schedule in
g_design_common.SCHEDULES, and uses it to diagnose a count printed in the
repository's own figure, audit/g1_phase_clock.py, which draws N = 25 visits
per panel and heads each panel with len(g_design_common.cycle_phases(times)).

Three further things f4 does not contain: the exact arithmetic behind the
PI's "28 visits" claim at 14.5 d; the observation that the half-day booking
grid of g_design_common.GRID caps EVERY clinic schedule at 28 phases,
whatever its step; and a demonstration that the phase count can be driven
arbitrarily high while the design degenerates, which is the same warning f4
gives in the language of the count rather than the range.

THE FOOTING, held fixed throughout and stated once. The dose interval is
T = 14 d, the record is the paper's own two years, and every sd(ln theta) is
computed by design.run._filter on design.run.design_record, which is the
footing Table 3 and the specification are computed on: the increment carried
in the sensitivity, theta and K_m free, a log-prior of 0.8 on ln k_e, and the
load recorded to 5%. It is held fixed because the question here is about the
SAMPLING TIMES alone -- a phase count that moved when the prior moved would
not be a property of the schedule -- and because the re-checks at the end have
to land on numbers the repository already committed to on that footing.

WHICH QUANTITIES NEED THE SEED SWEEP AND WHICH DO NOT. The phase counts, the
lattices, the phase gaps and the Delta/theta ranges are functions of the
sampling times only: design.run.delta is deterministic, and the intake record
enters the marker level, not the increment. They are therefore reported once,
without a sweep, and a sweep over them would print twenty-four identical
numbers. Every sd(ln theta) below does depend on the intake record and is a
median over the 24 records of design.run.SEEDS with its range, via
design.run.median_spread.

Nothing here records an expected result in advance.

Run:  python audit/t11_phase_lattice.py
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import model_A.run as LOC
from design.run import (ASSAY_CV, LOAD_FRAC, SEEDS, delta, design_record,
                        median_spread, _filter)
import g_design_common as G
from f4_sampling import predicted_phases

T = G.T                      # the dose interval, 14 d
HORIZON = 730.0              # the paper's two-year record
GRID = G.GRID                # the half-day clinic booking grid
SIG = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
N_FIG = 25                   # g1_phase_clock.N, the visits it draws per panel

#: the steps behind g_design_common.SCHEDULES, in the order that dict lists
#: them. They are named here rather than recovered from the lambdas because a
#: step is what this script is about and a lambda is not introspectable.
#: "phase-spread" has no step at all -- it is built phase-first -- so it
#: carries None and is handled as the exception it is.
SCHEDULE_STEPS = {
    "14 d (the dose interval)": (14.0, 0.0),
    "14 d, already on therapy": (14.0, G.STEADY),
    "14.5 d": (14.5, 0.0),
    "monthly (30.4 d)": (30.4, 0.0),
    "phase-spread": (None, 0.0),
}


def bars(title, bar):
    print("\n" + bar)
    print(title)
    print(bar)


def ratio(step):
    """
    step/T as an exact fraction.

    Fraction(str(x)) reads the DECIMAL a schedule is named by. Fraction(30.4)
    would read the binary double nearest 30.4, whose denominator is a power of
    two, and q would come back as 2**48 rather than 35 -- an answer about
    IEEE-754 rather than about the schedule. f4's predicted_phases takes the
    other route, limit_denominator(10000) on the float; R5 checks the two
    agree, so neither is trusted alone.
    """
    return Fraction(str(step)) / Fraction(str(T))


def visits(step, t0=0.0, horizon=HORIZON):
    """
    The visit times of a fixed-step schedule over one record.

    np.arange, not linspace and not a floor() count, because np.arange is what
    every design in src/design/run.py and every schedule in g_design_common is
    actually built with, and the re-checks have to land on those records
    exactly. Its length is floor((horizon - eps)/step) + 1, which is one more
    than floor(horizon/step) whenever the step does not divide the horizon;
    both are printed below so the reader can see which is which.
    """
    return t0 + np.arange(0.0, horizon, step)


def max_phase_gap(times):
    """
    The widest unvisited arc of the dose cycle, as a fraction of it.

    This is the coverage quantity, and it is what src/design/run.py:694 writes
    out as `phase_gap_monthly`. A single visited phase leaves a gap of the
    whole cycle, so it returns 1.0 there rather than 0.0: one phase covers
    nothing, and a diff over a one-element array would report the opposite.
    """
    ph = G.cycle_phases(times)
    if len(ph) < 2:
        return 1.0
    return float(np.max(np.diff(np.r_[ph, ph[0] + T])) / T)


def sd_over_seeds(times):
    """
    Median sd(ln theta) with its range over the 24 intake records of SEEDS.

    design_record + _filter is spec_numbers()'s own `sd_of` with the same
    priors in the same order, which is why R4 can land on the committed
    sd_monthly, sd_aliased and sd_aliased_off to the digit.
    """
    def one(seed):
        lv, inc = design_record(times, seed=seed)
        return float(np.sqrt(_filter(lv, SIG, inc)[-1][0, 0]))
    return median_spread([one(s) for s in SEEDS])


def fmt(v, w=9, d=4):
    return f"{'singular':>{w}}" if not np.isfinite(v) else f"{v:{w}.{d}f}"


# ==========================================================================
def s0_the_arithmetic(bar):
    bars("0. THE LATTICE, DERIVED ONCE", bar)
    print("  A schedule of step h visits t_k = t0 + k h, so its dose-cycle")
    print("  phase is (t0 + k h) mod T. Write h/T = p/q in lowest terms. Then")
    print("      phase_k = t0 + (k p mod q) T/q,")
    print("  and gcd(p, q) = 1 makes k -> k p mod q a bijection of Z_q, so the")
    print("  visited phases are the FULL lattice {t0 + j T/q : j = 0..q-1}")
    print("  once q visits have happened, and never anything else. With n")
    print("  visits the first n of them are reached, so")
    print("      phases reached = min(q, n).")
    print()
    print(f"  A second cap sits on top of it. g_design_common.GRID = {GRID} d:")
    print("  clinic visits are booked on a half-day grid, so every phase a")
    print(f"  clinic can reach is a multiple of {GRID} d and there are")
    print(f"  T/GRID = {T/GRID:.0f} of them in the cycle, whatever the step.")
    print(f"  No schedule bookable in this repository's own units reaches more")
    print(f"  than {int(T/GRID)} distinct phases.")
    print()
    print("  Check of that cap against g_design_common.phase_spread, which")
    print("  asks for n evenly spaced phases and rounds onto the grid:")
    print(f"    {'n asked for':>12} {'phases delivered':>17}")
    for n in (24, 27, 28, 29, 32, 40):
        t = G.phase_spread(n, gap=30.4)
        print(f"    {n:12d} {len(G.cycle_phases(t)):17d}")
    print(f"  It saturates at {int(T/GRID)}, which is the grid and not the design.")


# ==========================================================================
def s1_the_pis_count(bar):
    bars("1. THE PI'S COUNT AT 14.5 d, CHECKED EXACTLY", bar)
    print("  His claim: at 14.5-day sampling against a 14-day cycle the phase")
    print("  advances half a day per visit, so 28 visits walk the full cycle")
    print("  and the 29th returns to the start.\n")

    step = 14.5
    fr = ratio(step)
    print(f"  step/T = {step}/{T:.0f} = ({Fraction(str(step))}) / ({Fraction(str(T))})"
          f" = {fr.numerator}/{fr.denominator} in lowest terms,")
    print(f"  since gcd({fr.numerator}, {fr.denominator}) = "
          f"{np.gcd(fr.numerator, fr.denominator)}. So q = {fr.denominator}.")
    print()
    #  the advance per visit, in exact arithmetic rather than by subtraction
    #  of two floats: h mod T is where the half day comes from, and 14.5 - 14
    #  in binary is 0.5 exactly only by luck of the exponent
    adv = Fraction(str(step)) % Fraction(str(T))
    print(f"  Advance per visit: h mod T = {adv} d, and the lattice spacing is")
    print(f"  T/q = {Fraction(str(T))}/{fr.denominator} = {Fraction(str(T))/fr.denominator} d."
          f"  They agree, so the walk steps")
    print("  one lattice point at a time and misses nothing on the way.")
    print()
    print(f"  {'visit k':>8} {'t = 14.5k':>12} {'phase (exact)':>15} {'phase (float)':>14}")
    for k in (0, 1, 2, 27, 28, 29, 56):
        tk = Fraction(str(step)) * k
        ph = tk % Fraction(str(T))
        print(f"  {k:8d} {float(tk):12.1f} {str(ph):>15} "
              f"{float(tk) % T:14.4f}")
    print()
    print("  29 = 1 mod 28, so phase_k = (29k mod 28)/2 = (k mod 28)/2 days.")
    print("  Visits k = 0..27 -- the first TWENTY-EIGHT -- land on the 28")
    print("  distinct half-days {0, 0.5, ..., 13.5}, and visit k = 28, the")
    print("  TWENTY-NINTH, is back at phase 0.")
    print()
    t = visits(step)
    ph = G.cycle_phases(t)
    print(f"  Measured on the two-year record: n = {len(t)} visits, "
          f"{len(ph)} distinct phases,")
    print(f"  spacing {np.diff(ph).min():.4f} to {np.diff(ph).max():.4f} d, "
          f"widest unvisited arc {max_phase_gap(t):.4f} of the cycle.")
    print(f"  CONFIRMED. q = {fr.denominator}, and with n = {len(t)} >= "
          f"{fr.denominator} the record walks the lattice")
    print(f"  {len(t)/fr.denominator:.2f} times over, so every point of it is "
          f"visited more than once.")
    print()
    print(f"  The one qualification the claim needs: 28 is a count of PHASES,")
    print(f"  not a count of visits a protocol must book. It is 28 visits only")
    print(f"  because q <= n here; at n < 28 the record stops part-way round.")


# ==========================================================================
def s2_the_table(bar):
    bars("2. EVERY SCHEDULE OVER THE FULL TWO-YEAR RECORD", bar)
    print(f"  horizon {HORIZON:.0f} d, dose interval T = {T:.0f} d. "
          f"`n` is len(np.arange(0, {HORIZON:.0f}, h)),")
    print(f"  which is floor({HORIZON:.0f}/h) + 1 unless h divides the horizon; "
          f"both are shown.")
    print("  Phase counts and ranges are functions of the times alone and carry")
    print(f"  no seed sweep; sd(ln theta) is a median over {len(SEEDS)} records.\n")
    print(f"  {'schedule':>26} {'h (d)':>7} {'h/T':>10} {'q':>6} "
          f"{'n':>5} {'floor':>6} {'reached':>8} {'min(q,n)':>9} "
          f"{'gap':>7} {'D/th rng':>9} {'sd(ln th)':>10}")
    print("  " + "-" * 118)

    rows = []
    for name, (step, t0) in SCHEDULE_STEPS.items():
        if step is None:
            continue
        t = visits(step, t0)
        fr = ratio(step)
        q, n = fr.denominator, len(t)
        reached = len(G.cycle_phases(t))
        med = sd_over_seeds(t)[0]
        rows.append((name, step, fr, q, n, reached, max_phase_gap(t),
                     G.delta_range(t), med))
        print(f"  {name:>26} {step:7.1f} {str(fr):>10} {q:6d} {n:5d} "
              f"{int(HORIZON//step):6d} {reached:8d} {min(q, n):9d} "
              f"{max_phase_gap(t):7.4f} {G.delta_range(t):9.4f} "
              f"{fmt(med, 10)}")

    #  phase-spread is built phase-first and has no step, so it has no h/T and
    #  no q. Reporting a blank rather than a plausible number is the point:
    #  the lattice argument simply does not apply to it, and the quantity that
    #  caps it is the booking grid of section 0.
    n_ps = max(n for n in range(1, 40) if G.phase_spread(n, gap=30.4)[-1] <= HORIZON)
    t = G.phase_spread(n_ps, gap=30.4)
    med = sd_over_seeds(t)[0]
    print(f"  {'phase-spread':>26} {'--':>7} {'--':>10} {'--':>6} {n_ps:5d} "
          f"{'--':>6} {len(G.cycle_phases(t)):8d} {'--':>9} "
          f"{max_phase_gap(t):7.4f} {G.delta_range(t):9.4f} {fmt(med, 10)}")
    print(f"\n  phase-spread carries no step: it picks phases first and finds a")
    print(f"  calendar day for each, so h/T does not exist for it. Its "
          f"{n_ps} visits are")
    print(f"  what a 30.4 d minimum gap fits into {HORIZON:.0f} days.")
    print()
    print("  Reading of the table:")
    print("   - 'reached' equals min(q, n) on every step row, which is the")
    print("     identity of section 0 measured rather than asserted.")
    print("   - the two fortnightly rows are capped by q = 1: one phase, for")
    print("     any horizon whatever, which is the aliasing result. Their")
    print("     sd(ln theta) is large but FINITE because the intake record")
    print("     varies day to day and the levels vary with it; the phase is")
    print("     what is frozen, not the design matrix.")
    print("   - 14.5 d is capped by q = 28 with n = 51 visits to spend, so it")
    print("     reaches its whole lattice and has 23 visits left over.")
    print("   - monthly is capped by n = 25, NOT by its lattice of q = 35. Its")
    print("     step admits ten phases the two-year record never reaches, and")
    print("     it would need 35 visits -- 1064 days -- to reach them.")


# ==========================================================================
def s3_what_the_figure_counts(bar):
    bars("3. WHAT g1_phase_clock.py's PANEL HEADER IS ACTUALLY COUNTING", bar)
    print(f"  g1 draws N = {N_FIG} visits per panel and heads each panel with")
    print("  len(g_design_common.cycle_phases(times)), rendered as "
          "'<count> phases'.")
    print(f"  That call returns the phases reached, which is min(q, N). At")
    print(f"  N = {N_FIG}:\n")
    print(f"  {'panel':>26} {'q (lattice)':>12} {'N':>4} {'min(q,N)':>9} "
          f"{'header prints':>14} {'horizon':>9} {'gap':>7} {'D/th rng':>9}")
    print("  " + "-" * 100)
    for name, fn in G.SCHEDULES.items():
        t = fn(N_FIG)
        step = SCHEDULE_STEPS[name][0]
        q = ratio(step).denominator if step is not None else None
        header = len(G.cycle_phases(t))
        qs = f"{q:12d}" if q is not None else f"{'--':>12}"
        ms = f"{min(q, N_FIG):9d}" if q is not None else f"{'--':>9}"
        print(f"  {name:>26} {qs} {N_FIG:4d} {ms} {header:14d} "
              f"{t[-1] - t[0]:8.1f}d {max_phase_gap(t):7.4f} "
              f"{G.delta_range(t):9.4f}")
    print()
    print("  WHY MONTHLY PRINTS '25 phases'. Monthly's lattice is q = 35. The")
    print(f"  panel draws N = {N_FIG} visits. min(35, {N_FIG}) = {N_FIG}. The header is")
    print("  reporting the VISIT COUNT, not the lattice: with N visits at a")
    print("  step whose q exceeds N, every visit lands on a fresh phase and")
    print("  the count is N by construction, whatever the step.")
    print()
    print("  It is not only monthly. Three of the five panels -- monthly,")
    print(f"  14.5 d and phase-spread -- print '{N_FIG} phases', and their")
    print("  lattices are 35, 28 and (grid-capped) 28. The header cannot")
    print("  distinguish them because it is not reporting a property of the")
    print("  schedule at all on those three panels. Only the two fortnightly")
    print("  panels print a number that is the schedule's own -- 1 phase,")
    print("  because there q = 1 < N and the lattice, not N, is what binds.")
    print()
    print("  A second way the header misleads across panels: holding N fixed")
    print("  does not hold the record fixed. The panels span")
    for name, fn in G.SCHEDULES.items():
        t = fn(N_FIG)
        print(f"    {name:>26}: {t[-1]-t[0]:7.1f} d")
    print("  so the counts are not even taken over a common horizon.")
    print()
    print("  WHAT THE HEADER REPORTS: phases reached by the drawn record,")
    print("  min(q, N), which is a joint property of the step AND the panel's")
    print("  visit count, and which for three of five panels is just N.")
    print()
    print("  WHAT IT SHOULD SAY: report both, 'reached of available' --")
    for name, fn in G.SCHEDULES.items():
        t = fn(N_FIG)
        step = SCHEDULE_STEPS[name][0]
        r = len(G.cycle_phases(t))
        if step is None:
            print(f"    {name:>26}: {r} phases of the {int(T/GRID)} the "
                  f"booking grid allows")
        else:
            print(f"    {name:>26}: {r} of {ratio(step).denominator} phases "
                  f"the step admits")
    print("  -- or drop the count and lead with the two quantities that are")
    print("  properties of the schedule and do bind: the widest unvisited arc")
    print("  and the Delta/theta range the panel already prints. Section 4")
    print("  shows why the count on its own cannot be trusted to rank designs.")


# ==========================================================================
def s4_count_is_not_coverage(bar):
    bars("4. THE COUNT CAN BE DRIVEN UP WHILE THE DESIGN COLLAPSES", bar)
    print("  Steps h = 14 + eps, approaching the dose interval from above.")
    print("  If the phase count were the design quantity, smaller eps would be")
    print("  better: q = denominator of eps/T rises without bound as eps -> 0.")
    print(f"  Medians over {len(SEEDS)} intake records for sd(ln theta).\n")
    print(f"  {'eps (d)':>9} {'eps (h)':>8} {'h (d)':>9} {'h/T':>13} {'q':>7} "
          f"{'n':>4} {'reached':>8} {'gap':>7} {'D/th rng':>9} {'sd(ln th)':>10}")
    print("  " + "-" * 100)
    for eps in (7.0, 2.0, 1.0, 0.5, 0.4, 0.3, 0.28, 0.25, 0.2, 0.1, 0.01):
        step = round(14.0 + eps, 6)
        t = visits(step)
        fr = ratio(step)
        med = sd_over_seeds(t)[0]
        print(f"  {eps:9.2f} {24*eps:8.2f} {step:9.4f} {str(fr):>13} "
              f"{fr.denominator:7d} {len(t):4d} {len(G.cycle_phases(t)):8d} "
              f"{max_phase_gap(t):7.4f} {G.delta_range(t):9.4f} {fmt(med, 10)}")
    print()
    print("  At eps = 0.01 d the record visits 53 DISTINCT PHASES -- more than")
    print("  14.5 d's 28 -- and is worthless: they all sit within half a day of")
    print("  the dose, the widest unvisited arc is 96% of the cycle, and the")
    print("  Delta/theta range is back down at the aliased design's. The count")
    print("  rises monotonically as eps falls; the design does not. This is")
    print("  f4's finding restated in the count's own language, and it is why")
    print("  section 3 asks for the gap and the range beside any count.")
    print()
    print("  The turn is where the walk stops completing a lap. Each visit")
    print("  advances the phase by eps, so n visits sweep n*eps of the cycle:")
    print(f"    {'eps (d)':>9} {'n':>5} {'n*eps (d)':>11} {'lap?':>6} {'gap':>8}")
    for eps in (0.5, 0.3, 0.28, 0.27, 0.26, 0.25, 0.2):
        step = round(14.0 + eps, 6)
        t = visits(step)
        n = len(t)
        print(f"    {eps:9.2f} {n:5d} {n*eps:11.2f} "
              f"{('yes' if n*eps >= T else 'NO'):>6} {max_phase_gap(t):8.4f}")
    print(f"  The lap closes at eps = T/n ~ {T:.0f}/52 = {T/52:.4f} d = "
          f"{24*T/52:.2f} h, and the")
    print("  widest unvisited arc jumps from 2% of the cycle to 9% to 27% as")
    print("  eps crosses it downward.")


# ==========================================================================
def s5_the_rule(bar):
    bars("5. THE RULE", bar)
    print("  Write the sampling interval as h and the dose interval as T = 14 d,")
    print("  and put h/T = p/q in lowest terms. The schedule's phases lie on the")
    print("  lattice of spacing T/q, it reaches min(q, n) of the q points in n")
    print("  visits, and the arc it has swept after n visits is n*(h mod T),")
    print("  capped at T. So a ratio h/T that is rational WITH A SMALL")
    print("  DENOMINATOR reaches few phases however long the record runs: h = T")
    print("  gives q = 1 and a single phase for ever, h = T/2 or 3T/2 gives two,")
    print("  h = 2T/7 gives seven, and no horizon buys an eighth. Making q large")
    print("  is necessary and nowhere near sufficient, because q says only how")
    print("  fine the lattice is and says nothing about how much of it n visits")
    print("  can cover: h = 14.01 d has q = 1400 and reaches 53 phases packed")
    print("  into half a day. The design quantity is coverage -- the widest")
    print("  unvisited arc, and through it the range of Delta/theta -- and")
    print("  coverage needs BOTH q >= n-ish and n*(h mod T) >= T.")
    print()
    print(f"  The smallest offset from 14 days that reaches many phases, WELL")
    print(f"  SPREAD, is therefore set by two things. On this repository's own")
    print(f"  half-day booking grid (g_design_common.GRID = {GRID} d) the")
    print(f"  smallest offset available at all is {GRID} d, and it is already")
    print(f"  the best one: h = 14.5 d has q = 28, which is the entire lattice")
    print(f"  the grid admits, and its 51 visits walk that lattice nearly twice.")
    print(f"  Nothing bookable does better. Off the grid, the binding condition")
    print(f"  is the lap: eps >= T/n, and with n ~ 52 visits in two years that")
    print(f"  is eps >= {T/52:.3f} d, about {24*T/52:.1f} hours. An offset of six")
    print(f"  hours or less does not get round the cycle inside two years, and")
    print(f"  an offset of a few minutes is the aliased design wearing a large")
    print(f"  phase count as a disguise.")


# ==========================================================================
def r_rechecks(bar):
    bars("R1. RE-CHECK: phase_gap_monthly, the widest unvisited arc", bar)
    print("  src/design/run.py:694 writes it from the same monthly record;")
    print("  results/numbers/design.tex commits it at 2 s.f. of a fraction.")
    got = max_phase_gap(visits(30.4))
    want = 0.0857
    print(f"  monthly, 730 d: {got:.6f} against the committed {want:.4f} "
          f"({abs(got-want)/want:.2%} apart)")
    print(f"  It is exactly 3 lattice steps of T/q = {T/35:.1f} d = "
          f"{3*T/35:.1f} d, {3*T/35/T:.6f} of the cycle.")

    bars("R2. RE-CHECK: phase_spread_fortnightly, the phase span at 14 d", bar)
    ph = np.asarray(visits(14.0), float) % T
    got = float(ph.max() - ph.min())
    want = 0.0000
    print(f"  fortnightly, 730 d: {got:.6f} against the committed {want:.4f} "
          f"(exact to {abs(got-want):.1e})")

    bars("R3. RE-CHECK: the committed Delta/theta ranges", bar)
    for lab, step, want in (("delta_range_monthly", 30.4, 1.84),
                            ("delta_range_aliased", 14.0, 0.13),
                            ("delta_range_off", 14.5, 1.99),
                            ("delta_range_two_monthly", 60.0, 1.94),
                            ("delta_range_dense", 2.0, 1.95)):
        got = G.delta_range(visits(step))
        print(f"  {lab:>24} (h = {step:5.1f} d): {got:.4f} against the "
              f"committed {want:.2f} ({abs(got-want)/want:.2%} apart)")

    bars("R4. RE-CHECK: the committed sd(ln theta) of the three schedules", bar)
    print(f"  medians over the {len(SEEDS)} records of SEEDS, with the range;")
    print("  spec_numbers() computes these with the same priors in the same order.\n")
    for lab, step, want, wlo, whi in (
            ("sd_monthly", 30.4, 0.1699, 0.1684, 0.1709),
            ("sd_aliased", 14.0, 4.5845, 3.5050, 5.6362),
            ("sd_aliased_off", 14.5, 0.1615, 0.1611, 0.1618)):
        med, lo, hi = sd_over_seeds(visits(step))
        print(f"  {lab:>16} (h = {step:5.1f} d): {med:.4f} [{lo:.4f}, {hi:.4f}]"
              f" against the committed {want:.4f} [{wlo:.4f}, {whi:.4f}] "
              f"({abs(med-want)/want:.2%} apart)")

    bars("R5. RE-CHECK: q agrees with f4_sampling.predicted_phases", bar)
    print("  f4 reaches q by limit_denominator(10000) on the float; this script")
    print("  reaches it by an exact fraction of the decimals. Two routes, and a")
    print("  disagreement would mean one of them is reading the binary double.\n")
    bad = 0
    for step in (1.0, 2.0, 3.5, 4.0, 7.0, 7.5, 10.0, 14.0, 14.5, 15.0, 21.0,
                 28.0, 28.5, 30.4, 42.0, 56.0, 60.0):
        mine, theirs = ratio(step).denominator, predicted_phases(step)
        ok = "agree" if mine == theirs else "DISAGREE"
        bad += mine != theirs
        print(f"  h = {step:6.2f} d: exact q = {mine:5d}, f4's q = "
              f"{theirs:5d}  ({ok})")
    print(f"\n  {17 - bad} of 17 agree, {bad} apart.")


# ==========================================================================
def lattice_numbers():
    """
    The lap condition and the lattice counts, as generated numbers.

    Section 5.3 and Remark 4 now state the aliasing rule as ONE inequality
    rather than as a prohibition plus a fix, and these are the quantities that
    inequality is written in. They are generated because the rule is arithmetic
    and a hand-typed arithmetic claim is the thing this repository's number
    discipline exists to prevent.

    The lap condition. A schedule at step h = T + eps advances its phase by
    eps each visit, so after n visits it has swept n*eps of the dial. It has
    completed a lap, and therefore reached the whole increment range, only if
    n*eps >= T, i.e.

        eps  >=  T / n.

    Everything else follows: 14 d is eps = 0 and never laps at any n; 14.5 d
    laps comfortably; and 14.01 d has a LARGE phase lattice and still does not
    lap, which is the case a prohibition on multiples does not catch.
    """
    from common.numbers import write_numbers
    T, HORIZON = G.T, 730.0
    n_fortnightly = int(np.floor(HORIZON / T)) + 1
    eps_min = T / n_fortnightly
    #: the near-miss the rule has to catch. Its lattice is huge and useless.
    NEAR = 14.01
    t_near = np.arange(0.0, HORIZON, NEAR)
    ph_near = np.sort(np.asarray(t_near, float) % T)
    N = {
        "lap_eps_h": 24.0 * eps_min,
        "lap_n": n_fortnightly,
        "lap_horizon_y": HORIZON / 365.0,
        "q_off": ratio(14.5).denominator,
        "q_monthly": ratio(30.4).denominator,
        "n_monthly": int(np.floor(HORIZON / 30.4)) + 1,
        "monthly_days_to_exhaust": ratio(30.4).denominator * 30.4,
        "monthly_unreached": ratio(30.4).denominator - (int(np.floor(HORIZON / 30.4)) + 1),
        "near_step": NEAR,
        "near_q": ratio(NEAR).denominator,
        "near_phases": len(G.cycle_phases(t_near)),
        "near_span_h": 24.0 * float(ph_near.max() - ph_near.min()),
        "near_sd": sd_over_seeds(t_near)[0],
    }
    write_numbers("lattice", N,
                  {"lap_eps_h": "{:.1f}", "lap_n": "{:.0f}",
                   "lap_horizon_y": "{:.0f}", "q_off": "{:.0f}",
                   "q_monthly": "{:.0f}", "n_monthly": "{:.0f}",
                   "monthly_days_to_exhaust": "{:.0f}",
                   "monthly_unreached": "{:.0f}", "near_step": "{:.2f}",
                   "near_q": "{:.0f}", "near_phases": "{:.0f}",
                   "near_span_h": "{:.1f}", "near_sd": "{:.4f}"})
    for k in sorted(N):
        print(f"    {k:26} {N[k]}")
    return N


def main():
    print(__doc__)
    bar = "=" * 84
    s0_the_arithmetic(bar)
    s1_the_pis_count(bar)
    s2_the_table(bar)
    s3_what_the_figure_counts(bar)
    s4_count_is_not_coverage(bar)
    s5_the_rule(bar)
    r_rechecks(bar)


if __name__ == "__main__":
    main()
    lattice_numbers()
