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


def part_c():
    print()
    print("=" * 88)
    print("PART C.  the regret thresholds, whose owner knows about seeds and does")
    print("         not use them here")
    print("=" * 88)
    print("""
  audit/d6_regimes.py computes a twelve-seed spread for its error rates and
  prints "(12 seeds, 400 replicates each)", then builds the regret-threshold
  table from the SINGLE-seed ensemble of its line 109. These are the two
  numbers 07b_decision.tex:70-71 took from that table.
""")
    import d6_regimes as D6
    print(f"    {'regime':>10} {'paper':>8} {'seed 11':>9} "
          f"{'mean':>8} {'min':>8} {'max':>8}")
    print("    " + "-" * 56)
    for lab, rho in (("managed", 1.0), ("crisis", -0.25)):
        vals = []
        for s in SEEDS:
            ens = D6.fit_regime(LV9, rho, reps=400, seed=s)
            taus = np.linspace(-30000.0, 30000.0, 601)
            #  d6's published amplitudes, its own line 98, and the same
            #  baseline. They are locals in d6 so they are restated here
            #  rather than imported, and 2.98 / 3.58 are the values
            #  07b_decision itself quotes.
            gap = np.array([1.0 * D6.err_at(ens, 3.58 * BASE, t)[0]
                            - D6.err_at(ens, 2.98 * BASE, t)[0]
                            for t in taus])
            sg = np.where(np.diff(np.sign(gap)) != 0)[0]
            vals.append(float(taus[sg[0]]) / P.theta if sg.size else np.nan)
        v = np.array(vals, float)
        fin = v[np.isfinite(v)]
        s11 = vals[11] if len(vals) > 11 else np.nan
        print(f"    {lab:>10} {PAPER_TAU[lab]:8.2f} {s11:9.2f} "
              f"{fin.mean():8.2f} {fin.min():8.2f} {fin.max():8.2f}")
    print("""
  If these spreads straddle zero the qualitative claim changes, because the
  sign of tau* is what "escalate earlier" versus "escalate later" means.""")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
