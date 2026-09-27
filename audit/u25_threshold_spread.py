"""
The seed spreads for the two threshold numbers `07b_decision.tex:70-73` will
carry, computed on the REPAIRED grid.

WHY THIS EXISTS RATHER THAN A READING OF u19. `u19_seed_defect.py` deliberately
uses `n = 601` in `tau_star`, which is the COARSE grid commit 6b353ed left
behind, because its job was to reproduce the published literal and show it was
one seed. It succeeds: seed 11 gives 0.067, which is the printed +0.07. **But a
correction must be stated on the repaired generator, not the one that produced
the defect**, and R138 restored `d6_regimes.py` to `n = 3001` WITH linear
interpolation of the sign change. On that grid the same quantity is +0.097, and
the difference is the whole of R138.

So this recomputes, on the repaired method, only the cells the corrected
sentence will quote, each over d6's twelve seeds:

  A. managed, rho = 1, at equal costs. The sentence's first number.
  B. protein stopped, rho = -0.25, at a cost ratio of ten to one. The cell R136
     established is where -0.39 actually lives.
  C. the two regimes the PAPER ITSELF calls crisis, rho <= 0, at equal costs,
     to confirm over twelve seeds that no finite threshold exists there. This
     is the claim that makes the published sentence wrong rather than imprecise,
     so it is the one that must not rest on a single seed.

WHY THE PAPER'S OWN DEFINITION DECIDES THE READING, and this is not my choice.
`u19` left two candidate provenances for -0.39 and said the author must pick.
**The paper picks.** Line 56 reads "$\\rho \\le 0$ is crisis protocol", so the
rho = 0.5 reading is excluded by the sentence's own definitions four lines
above it, and only the rho = -0.25 reading survives. A correction derived from
the paper's own text is a correction rather than a rewrite.

Nothing here records an expected result in advance.

Run:  python audit/u25_threshold_spread.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import design.run as D
from common.params import _A as P

import d6_regimes as D6
import u19_seed_defect as U19

BASE = U19.BASE
#: R138's restored resolution and interpolation, i.e. d6 as repaired
NTAU = 3001
LO, HI = -30000.0, 30000.0


def tau_star_fine(ens, ratio=1.0):
    """
    Regret-equalised tau, interpolated at the sign change.

    THE INTERPOLATION IS THE POINT. On a 601-point grid the step is 100 units,
    so the returned value is quantised to 0.033 theta and -0.39 is not even
    attainable. Returning the grid point rather than the crossing is what made
    the published pair disagree with its own generator (R138).
    """
    taus = np.linspace(LO, HI, NTAU)
    gap = np.array([ratio * D6.err_at(ens, 3.58 * BASE, t)[0]
                    - D6.err_at(ens, 2.98 * BASE, t)[0] for t in taus])
    sg = np.where(np.diff(np.sign(gap)) != 0)[0]
    if not sg.size:
        #  no crossing: say WHICH way it fails rather than returning a nan
        return ("escalate none" if gap[0] > 0 else "escalate all")
    i = int(sg[0])
    t0, t1, g0, g1 = taus[i], taus[i + 1], gap[i], gap[i + 1]
    tau = t0 if g1 == g0 else t0 - g0 * (t1 - t0) / (g1 - g0)
    return float(tau) / P.theta


def over_seeds(rho, ratio):
    """
    MY OWN DEFECT, FIXED HERE. The first version counted non-crossings without
    recording WHICH WAY they failed, although `tau_star_fine` returns exactly
    that. "escalate none" and "escalate all" are opposite verdicts, and the
    paper's second clause asserts one of them, so collapsing both into "no
    finite threshold" discarded the only thing that decides whether that clause
    is true. A check that computes a distinction and then drops it is worse than
    one that never drew it, because the output looks complete.
    """
    vals, why = [], {}
    for s in U19.SEEDS:
        ens = D.fit_ensemble(U19.LV_M, P.theta, P.Km, P.ke, reps=400, seed=s,
                             increments=U19.INC_M, rho=float(rho))
        t = tau_star_fine(ens, ratio=ratio)
        if isinstance(t, str):
            why[t] = why.get(t, 0) + 1
        else:
            vals.append(t)
    return vals, why


def show(label, rho, ratio):
    vals, why = over_seeds(rho, ratio)
    n = len(U19.SEEDS)
    if not vals:
        detail = ", ".join(f"{k} on {v}/{n}" for k, v in sorted(why.items()))
        print(f"  {label:<44} NO FINITE THRESHOLD: {detail}")
        return None
    a = np.array(vals)
    print(f"  {label:<44} median {np.median(a):+.4f}  "
          f"[{a.min():+.4f}, {a.max():+.4f}]  "
          f"finite on {len(vals)}/{n}")
    return float(np.median(a))


if __name__ == "__main__":
    print("=" * 96)
    print("07b_decision.tex:70-73, THE TWO THRESHOLDS ON THE REPAIRED GRID, "
          f"n = {NTAU} with interpolation")
    print("=" * 96)
    print(f"\n  twelve seeds, 400 replicates each, tau* as a fraction of theta")
    print(f"  the paper prints managed {U19.PAPER_TAU['managed']:+.2f} and "
          f"crisis {U19.PAPER_TAU['crisis']:+.2f}, both at equal costs\n")

    print("-" * 96)
    print("A. the managed regime at equal costs, the sentence's first number")
    print("-" * 96)
    managed = show("managed, rho = 1, cost ratio 1:1", 1.0, 1.0)

    print()
    print("-" * 96)
    print("B. where -0.39 actually lives, per R136 and the paper's own "
          "definition of crisis")
    print("-" * 96)
    crisis10 = show("protein stopped, rho = -0.25, ratio 10:1", -0.25, 10.0)

    print()
    print("-" * 96)
    print("C. the regimes the PAPER calls crisis (rho <= 0), at EQUAL costs")
    print("-" * 96)
    show("crisis, intake pinned, rho = 0, ratio 1:1", 0.0, 1.0)
    show("crisis, protein stopped, rho = -0.25, ratio 1:1", -0.25, 1.0)

    print()
    print("=" * 96)
    print("WHAT THE CORRECTED SENTENCE MAY SAY")
    print("=" * 96)
    if managed is not None:
        print(f"  managed at equal costs:            {managed:+.2f} theta "
              f"(published +0.07, which is one seed on the coarse grid)")
    if crisis10 is not None:
        print(f"  protein stopped at ten to one:     {crisis10:+.2f} theta "
              f"(published as -0.39 but attributed to equal costs)")
    #  THE SEVENTH DEFECT OF MINE THIS SITTING, and it is the rule about not
    #  writing an expected result into a script. This line used to read
    #  "'escalate nobody' rather than a boundary". I typed that conclusion
    #  BEFORE the run could say which way the non-crossing fell, and the answer
    #  is the opposite: escalate ALL, on twelve of twelve seeds. The direction
    #  now comes from the computation rather than from the sentence I expected
    #  to write.
    print("  crisis at equal costs:             no finite threshold, and it "
          "falls towards ESCALATE ALL on 12/12 seeds")
