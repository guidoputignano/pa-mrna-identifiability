"""
Phase 5, point 4: what Baek's Source Data could and could not discriminate.

WRITTEN BEFORE THE FILE ARRIVES, DELIBERATELY. The instruction is to report
what the data can and cannot settle before any analysis of it, so everything
here is computed from the DESIGN (its schedule, its group size, the absence of
a diet record) and from the simulators already in this repository. No mouse
datum enters. When refs/data/ holds the file, this file is the yardstick the
fit is read against, and it cannot be retuned afterwards.

THE DESIGN, from Baek et al. 2024 Fig. 2's caption, read in Phase 5 step 0:

    "Pcca-/- (A138T) mice received 2 IV bolus doses of Tris-sucrose (control)
     or 0.5 mg/kg mRNA-3927 on Days 0 and 28 (n = 6 mice/group) ... Blood was
     collected from mice prior to treatment and on Days 2, 8, 14, 22, and 28
     after the first IV dose, and 2 days after the second IV dose."

So: one treated arm, n = 6, a pre-dose sample and five post-dose samples across
one 28-day cycle, plus one sample two days into the second cycle. Several cells
are n = 2 or n = 5 rather than 6, per the same caption. No diet record exists,
and a laboratory colony is fed a fixed chow, so the substrate load is constant
by construction.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import design.run as DR
import g_design_common as G
import k1_sources as K1
import k2_structures as K2

TH, LB = G.TH, G.LB
#: the mouse sampling grid, days after the first dose. The Day 0 sample is
#: "prior to treatment", so it is evaluated PRE-dose; a first version put it at
#: t = 0+ and every shape statistic below came out at 100%, because an
#: algebraic marker is already at its minimum the instant the enzyme peaks
#: while a turnover marker has not moved yet. That is an artefact of placing a
#: pre-dose sample after the dose, not a discrimination.
#:
#: Day 28 carries the same ambiguity as the human Extended Data Table 3, in
#: miniature: doses are on Days 0 and 28, and the caption's "on ... Day 28" does
#: not say which side of the second dose it falls. Taken here as the trough,
#: which is the natural reading of a Day-28 sample in the FIRST cycle, and
#: flagged because it is the same unstated convention twice.
MOUSE_PRE = -1e-6
MOUSE_POST = np.array([2.0, 8.0, 14.0, 22.0, 28.0 - 1e-6, 30.0])
MOUSE_DAYS = np.concatenate([[MOUSE_PRE], MOUSE_POST])
MOUSE_CYCLE = 28.0
N_MICE = 6
IIV_BASELINE = 0.511          # Attarwala Table 1, 3-HP baseline, between MICE


def head(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


def _mouse_profiles(c_A, c_B, lam=None, dose=None, npts=60001):
    """
    PUB_A and PUB_B on the mouse regimen: untreated at t = 0-, doses at Days 0
    and 28, read over [0, 30].

    The initial condition is the UNTREATED steady state. A first version of
    this function started PUB_A at the steady state of the POST-dose enzyme,
    which put its nadir at t = 0 and made every statistic below meaningless.
    """
    lam = DR.LAMBDA if lam is None else lam
    dose = DR.D_AMP if dose is None else dose
    t = np.linspace(MOUSE_PRE, 30.0, npts)
    E = np.full(t.size, float(TH))
    for k in (0.0, MOUSE_CYCLE):
        on = t >= k
        E[on] += dose * np.exp(-lam * (t[on] - k))

    gA, gB = K2.PubA(c_A), K2.PubB(c_B)
    fA = K2._inhib(E, c_A)
    mA = np.empty_like(t)
    mA[0] = LB * K2._inhib(float(TH), c_A) / gA.kout     # UNTREATED
    dt = t[1] - t[0]
    ex = np.exp(-gA.kout * dt)
    for i in range(1, t.size):                            # exponential integrator
        mA[i] = mA[i - 1] * ex + LB * fA[i] * (1.0 - ex) / gA.kout
    mB = np.asarray(gB.steady(E), float)
    return t, E, mA, mB, gA, gB


def what_it_can():
    head("WHAT THE SOURCE DATA CAN DISCRIMINATE")
    print("""
  ONE THING, AND IT IS THE RIGHT ONE: whether the marker has a STATE.

  Step 3 established that the whole difference between the two published
  structures is a marker time constant. PUB_A's marker is a turnover state and
  its transient supplies the fourth number that makes theta identifiable;
  PUB_B's marker is algebraic in the enzyme and has no transient at all. A
  TIME COURSE is exactly the instrument that separates those two, and a
  cross-section is not. Baek's figure is a time course.

  The observable signature of a marker state is LAG: PUB_B's marker reaches
  its minimum at the moment the enzyme peaks, which is the dose itself, while
  PUB_A's lags by an amount set by its own half-life. Simulated on the mouse
  schedule at OUR operating point, which is a statement about the SCHEDULE and
  not about the mouse:
""")
    print(f"      {'calibration':12s} {'branch':6s} {'IC50':>10} "
          f"{'marker t1/2':>12} {'lag to nadir':>13} {'shape gap':>10} "
          f"{'at sampled days':>16}")
    print("  " + "-" * 84)
    rows = []
    for tgt in (0.30, 0.60):
        cB = float(K2._solve_shape(K2.PubB.shape, tgt)[0])
        for c, tag in zip(K2._solve_shape(K2.PubA.shape, tgt),
                          ("slow", "fast")):
            t, E, mA, mB, gA, gB = _mouse_profiles(c, cB)
            win = (t > 0.0) & (t <= MOUSE_CYCLE)
            lag = float(t[win][int(np.argmin(mA[win]))])
            zA = (mA - mA.min()) / np.ptp(mA)
            zB = (mB - mB.min()) / np.ptp(mB)
            gap_all = float(np.max(np.abs(zA - zB)))
            sA = np.interp(MOUSE_POST, t, zA)
            sB = np.interp(MOUSE_POST, t, zB)
            gap_s = float(np.max(np.abs(sA - sB)))
            rows.append((tgt, tag, c, np.log(2) / gA.kout, lag, gap_all, gap_s))
            print(f"      {100*tgt:9.0f}%   {tag:6s} {c:10.4g} "
                  f"{np.log(2)/gA.kout:10.4g} d {lag:11.3g} d "
                  f"{100*gap_all:9.1f}% {100*gap_s:15.1f}%")
    slowgaps = " and ".join(f"{100*r[6]:.1f}%"
                            for r in rows if r[1] == "slow")
    fastgaps = " and ".join(f"{100*r[6]:.1f}%"
                            for r in rows if r[1] == "fast")
    se_best = 100 * IIV_BASELINE / np.sqrt(N_MICE)
    se_thin = 100 * IIV_BASELINE / np.sqrt(2)
    print(f"""
  'shape gap' is the largest difference between the two NORMALISED profiles,
  so it is free of any scale or baseline difference between the structures; the
  last column restricts it to the seven days actually sampled, which is what a
  fit would see.

  Read against the noise: Attarwala's own between-mouse variability on the
  3-HP baseline is {100*IIV_BASELINE:.1f}%, and n = {N_MICE} per group, with n = 2 or 5 in
  several cells per the caption. The standard error on a group mean is
  therefore about {se_best:.1f}% at best and {se_thin:.1f}% in the thinnest cells, BEFORE assay
  error. A shape gap has to clear that to be a discrimination rather than a
  preference, and only the SLOW branch does: {slowgaps} against the fast
  branch's {fastgaps}.
""")
    return rows


def what_it_cannot(rows):
    head("WHAT IT CANNOT DISCRIMINATE")
    fastest = min((r for r in rows if r[1] == "fast"), key=lambda r: r[3])
    fast_hl, fast_lag, fast_gap = fastest[3], fastest[4], 100 * fastest[6]
    fast_min = 24 * 60 * fast_hl
    se = 100 * IIV_BASELINE / np.sqrt(N_MICE)
    print(f"""
  1. IT CANNOT IDENTIFY theta UNDER PUB_B, AND THIS IS STRUCTURAL, NOT A
     POWER PROBLEM. A laboratory colony on fixed chow has a CONSTANT substrate
     load. Step 3b's rank table at a constant load gives PUB_B rank 3 of 5, a
     deficit of TWO. No number of mice and no number of timepoints changes
     that: it is a property of the map, not of the sample size. The same table
     gives the paper's own model rank 3 of 4 at a constant load, so this
     dataset cannot close the paper's fibre either.

  2. IT CANNOT SETTLE THE HUMAN SAMPLING-TIME QUESTION, which is the one open
     item that moves the design numbers. It is a different species, a
     different regimen (0.5 mg/kg Q4W against 0.30 to 0.90 mg/kg Q2W or Q3W),
     and it says nothing about when the TRIAL drew its post-dose samples.

  3. IT CANNOT TRANSFER AN OPERATING POINT. The Pcca-/- (A138T) mouse carries
     about 2% of normal liver PCC activity. The repository's theta is set at a
     different fraction, stated two ways in params.py (10% and 5%), and the
     ratio m_bar/K_m = 0.1 that every excitation statement is made against is
     ours. Nothing in a mouse time course fixes either.

  4. IT PROBABLY CANNOT SEPARATE INDIVIDUAL MICE. The figure reports mean and
     SD per timepoint. If the Source Data holds only those summaries rather
     than per-animal series, then with {100*IIV_BASELINE:.0f}% between-mouse variability the
     within-animal shape, which is what a turnover model is identified from,
     is not recoverable at all. THIS IS THE FIRST THING TO CHECK WHEN THE FILE
     ARRIVES, and it decides whether point 1 of the previous section is even
     available.

  5. ON THE FAST BRANCH IT CANNOT SEE THE STATE. The fast calibration branch
     of PUB_A has a marker half-life of {fast_hl:.4g} d, about {fast_min:.0f} minutes, and
     lags the enzyme by only {fast_lag:.3g} d. Sampled on a grid whose finest spacing
     is 2 days, such a marker is indistinguishable from an instantaneous one:
     its shape gap against PUB_B at the sampled days is {fast_gap:.1f}%, against a
     standard error on a group mean of {se:.1f}% at best. So a null result on this
     dataset would rule out the SLOW branch of PUB_A and say nothing about the
     fast one, and a null result is therefore NOT evidence for PUB_B.
""")


def the_standing_rule():
    head("THE STANDING RULE FOR WHEN THE FILE ARRIVES")
    print("""
  DO NOT FIT IT. The instruction is to report what it can and cannot
  discriminate BEFORE any analysis, and this file is that report. When
  refs/data/ holds the Source Data the first pass is descriptive only:

    a. Is it per-animal or per-timepoint summaries? (Section above, item 4.
       If summaries, stop and say so: most of the value is gone.)
    b. How many animals and how many timepoints actually carry a 3-HP value,
       against the caption's n = 6 with cells at n = 2 and n = 5?
    c. Is the pre-treatment sample present for every animal, since a turnover
       model's baseline is what its transient is measured against?
    d. What are the units and the assay LOQ, against the 25 umol/L that
       Attarwala's text quotes and the 0 to 200 uM axis of Fig. 2e?

  Then report, and wait. Nothing is fitted until you say so.
""")


def main():
    print(__doc__.strip())
    rows = what_it_can()
    what_it_cannot(rows)
    the_standing_rule()


if __name__ == "__main__":
    main()
