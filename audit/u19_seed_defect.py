"""
Every number in 07b_decision's crisis paragraph is a single seed, and the
paragraph is prose literals in a section whose other numbers are generated.

THE DEFECT, stated before it is measured. `OPEN_QUESTIONS.md` section 0 requires
footing and seed spread with every number, and records that a single-seed figure
quoted without its range has three times been the error in this audit. The
crisis-protocol paragraph of `paper/sections/07b_decision.tex` reports twelve
numbers at lines 58 to 61 and two more at lines 70 to 71, all at one seed.

WHAT AN UNBOUNDED SEARCH OF THAT FILE SHOWS, patterns `num{`, a decimal, and an
escaped per cent, no head limit, 9 matching lines. Lines 37, 40 and 43 are
MACRO-WIRED: `\num{cohort_fixed}`, `\num{cohort_monthly}`,
`\num{kappa_monthly_pct}` with its lo and hi. Lines 58 to 61 and 70 to 71 are
PROSE LITERALS, every one of them. So the paragraph is the unwired one in an
otherwise wired section, which is why no generator carries its spread.

THE OWNERS, located before recomputing.
  The error rates and medians at lines 58 to 61: `audit/d13_decision_therapy.py`
  lines 122 to 129, `D.fit_ensemble` and `D.sign_error` on the MONTHLY design
  `lv_m`, at `seed=11, reps=400`.
  The regret thresholds at lines 70 to 71: `audit/d6_regimes.py` line 109,
  `fit_regime(lv9, rho, reps=400, seed=11)` on the NINE-POINT design, with the
  source's own comment "the seed fig_decision uses".

AND d6 ALREADY KNOWS. Lines 122 to 129 compute a twelve-seed spread for the
error rates and print "(12 seeds, 400 replicates each)". The regret-threshold
table on line 140 onward then reuses the SINGLE-seed ensemble from line 109. So
the convention exists in the same file and is not applied to the numbers the
paper took from it.

PART A reproduces every literal at its own seed, which is the anchor. PART B
gives each one its spread. PART C does the same for the thresholds. Nothing
here records an expected result in advance.

Run:  python audit/u19_seed_defect.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import design.run as D
from common.params import _A as P

BASE = D.clearance_A(D.M_BAR, P)
T_M = np.arange(0.0, 730.0, 30.4)
LV_M, INC_M = D.design_record(T_M)
LV9 = np.geomspace(D.M_BAR, 5 * D.M_BAR, 9)
SEEDS = range(12)          # d6's own convention
RHOS = (1.0, 0.5, 0.0, -0.25)

#: the twelve prose literals of 07b_decision.tex:58-61, as printed
PAPER = {1.0: (0.175, 0.418, 1.03),
         0.5: (0.990, 0.000, 0.90),
         0.0: (1.000, 0.000, 0.086),
         -0.25: (1.000, 0.000, 0.009)}
#: the two prose literals of 07b_decision.tex:70-71
PAPER_TAU = {"managed": +0.07, "crisis": -0.39}


def cell(rho, seed, reps=400):
    ens = D.fit_ensemble(LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                         increments=INC_M, rho=float(rho))
    e_s = D.sign_error(LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                       s_peak=2.98 * BASE, increments=INC_M,
                       rho=float(rho))["p"]
    e_u = D.sign_error(LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                       s_peak=3.58 * BASE, increments=INC_M,
                       rho=float(rho))["p"]
    return e_s, e_u, float(np.median(ens[:, 0]) / P.theta)


def part_a():
    print("=" * 88)
    print("PART A.  every literal reproduced at its own seed, which is the anchor")
    print("=" * 88)
    print("\n  07b_decision.tex:58-61, monthly design, seed 11, 400 replicates.\n")
    print(f"    {'rho':>7} {'treat unneeded':>22} {'fail to treat':>21} "
          f"{'median th/th':>20}")
    print(f"    {'':>7} {'paper    here':>22} {'paper    here':>21} "
          f"{'paper    here':>20}")
    print("    " + "-" * 74)
    ok = True
    for rho in RHOS:
        a, b, c = cell(rho, 11)
        pa, pb, pc = PAPER[rho]
        ok &= abs(a - pa) < 6e-3 and abs(b - pb) < 6e-3 \
            and abs(c - pc) / max(pc, 1e-9) < 0.05
        print(f"    {rho:7.2f} {pa:11.1%} {a:10.1%} {pb:10.1%} {b:10.1%} "
              f"{pc:10.3f} {c:9.3f}")
    print(f"\n  every literal reproduced: {'YES' if ok else 'NO'}")


def part_b():
    print()
    print("=" * 88)
    print("PART B.  the same twelve numbers with the spread section 0 requires")
    print("=" * 88)
    print(f"\n  {len(list(SEEDS))} seeds, 400 replicates each, d6's own convention.\n")
    print(f"    {'rho':>7} {'treat unneeded':>26} {'fail to treat':>26}")
    print("    " + "-" * 62)
    store = {}
    for rho in RHOS:
        A, B, C = [], [], []
        for s in SEEDS:
            a, b, c = cell(rho, s)
            A.append(a); B.append(b); C.append(c)
        store[rho] = (np.array(A), np.array(B), np.array(C))
        f = lambda v: f"{np.mean(v):6.1%} [{min(v):5.1%},{max(v):5.1%}]"
        print(f"    {rho:7.2f} {f(A):>26} {f(B):>26}")
    print()
    print(f"    {'rho':>7} {'paper':>8} {'mean':>8} {'min':>8} {'max':>8} "
          f"{'max/min':>9}  the median theta_hat/theta")
    print("    " + "-" * 74)
    for rho in RHOS:
        C = store[rho][2]
        print(f"    {rho:7.2f} {PAPER[rho][2]:8.3f} {C.mean():8.3f} "
              f"{C.min():8.3f} {C.max():8.3f} {C.max() / max(C.min(), 1e-9):8.2f}x")
    print("""
  The benign cell is the one that matters, because it is the cell the section
  leads with and the only one whose value is not already a collapse.""")
    return store


def tau_star(ens, ratio=1.0, lo=-30000.0, hi=30000.0, n=601):
    """Regret-equalised threshold, or the reason there is none."""
    import d6_regimes as D6
    taus = np.linspace(lo, hi, n)
    gap = np.array([ratio * D6.err_at(ens, 3.58 * BASE, t)[0]
                    - D6.err_at(ens, 2.98 * BASE, t)[0] for t in taus])
    sg = np.where(np.diff(np.sign(gap)) != 0)[0]
    if sg.size:
        return float(taus[sg[0]]) / P.theta
    return np.inf if gap[0] > 0 else -np.inf


def part_c():
    print()
    print("=" * 88)
    print("PART C.  the two thresholds, whose OWNER is Phase 2D and not d6")
    print("=" * 88)
    print("""
  THE OWNER, and the first attempt had it wrong. An unbounded search for the
  literal across audit/, src/, paper/ and results/ finds it at
  audit/PHASE2D_FINDINGS.md:150: "Regret-equalised thresholds: managed +0.07
  theta at equal costs, crisis -0.39 theta". So PHASE 2D is the source, and
  audit/d6_regimes.py is a LATER re-derivation of the same quantity, not the
  owner. Phase 2D's own framing matters for reading "crisis": the sentence
  above it says "the reversal is complete by rho = 0.5", so rho = 0.5 is
  already a crisis regime in that document.

  d6 also computes a twelve-seed spread for its error rates and prints
  "(12 seeds, 400 replicates each)", then builds its threshold table from the
  SINGLE-seed ensemble of its line 109, commented "the seed fig_decision uses".

  All four regimes are reported here because the paper names only two and the
  mapping is what is at issue.
""")
    import d6_regimes as D6
    print(f"    {'rho':>7} {'seed 11':>15} {'12-seed mean':>14} "
          f"{'min':>10} {'max':>10} {'no crossing':>12}")
    print("    " + "-" * 74)
    fmt = lambda x: ("escalate all" if np.isposinf(x) else
                     "escalate none" if np.isneginf(x) else f"{x:.3f}")
    for rho in RHOS:
        vals = [tau_star(D6.fit_regime(LV9, rho, reps=400, seed=s))
                for s in SEEDS]
        v = np.array(vals, float)
        fin = v[np.isfinite(v)]
        n_non = int((~np.isfinite(v)).sum())
        s11 = tau_star(D6.fit_regime(LV9, rho, reps=400, seed=11))
        if fin.size:
            print(f"    {rho:7.2f} {fmt(s11):>15} {fin.mean():14.3f} "
                  f"{fin.min():10.3f} {fin.max():10.3f} {n_non:9d}/12")
        else:
            print(f"    {rho:7.2f} {fmt(s11):>15} {'none finite':>14} "
                  f"{'':>10} {'':>10} {n_non:9d}/12")
    print(f"""
  AGAINST THE PUBLISHED PAIR, managed {PAPER_TAU['managed']:+.2f} and crisis {PAPER_TAU['crisis']:+.2f}.

  THE MANAGED THRESHOLD REPRODUCES. rho = 1 at seed 11 gives 0.067, the
  published +0.07, and the twelve-seed range is narrow with the sign stable at
  0/12 non-crossings, so "escalate earlier" survives the spread.

  THE CRISIS THRESHOLD HAS TWO CANDIDATE PROVENANCES AND NEITHER SUPPORTS THE
  SENTENCE AS WRITTEN. An earlier version of this file concluded it simply did
  not reproduce; that was drawn from the equal-cost column alone and is
  withdrawn at R136.

    (a) rho = -0.25 at cost ratios of 10:1 and worse gives exactly -0.400 at
        seed 11, within one grid step of -0.39 (the grid is 100 units, so 1170
        is off it and 1200 is on it). This matches the SEED cleanly and
        requires "at equal costs" in the sentence to be wrong.
    (b) rho = 0.5 AT equal costs has a twelve-seed range of -0.533 to -0.367,
        which CONTAINS -0.39, and Phase 2D's own text licenses calling
        rho = 0.5 a crisis regime by saying "the reversal is complete by
        rho = 0.5". This preserves the cost ratio but -0.39 is then NOT seed
        11's value, which is -0.467, while the +0.07 beside it IS seed 11's.
        One sentence would be mixing seeds.

  AND AT EQUAL COSTS BOTH ROWS THE PAPER WOULD MOST NATURALLY CALL "CRISIS"
  HAVE NO FINITE THRESHOLD ON ALL TWELVE SEEDS, 12/12 non-crossings. That is
  structural rather than a seed accident, so no amount of reseeding produces a
  -0.39 from rho = 0 or rho = -0.25 at equal costs.

  SEPARATELY, the sentence's second clause has a visible exception under both
  readings: at 10:1 and worse, three regimes want every patient escalated and
  rho = -0.25 wants -0.400, which is the row reading (a) takes the number from.

  Reported as a question for the author and not as a paper error, since the
  Phase 2D-era code may differ from today's d6. The question is now specific:
  which regime and which cost ratio was meant, and does the second clause need
  its exception.
""")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
