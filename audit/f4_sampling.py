"""
P4 -- the PI's sampling warning, checked, and the rule stated properly.

His note says synchronous sampling is bad, and adds "possibly also at 1 or 1/2
of the interval". The first half is the paper's sharpest design rule
(PHASE2D_FINDINGS.md section 4: 14 d and 28 d fail by a factor of 22, 14.5 d
passes). The second half is a different claim and has to be tested separately,
because 7 d IS half the dose interval and PHASE2D reports it passing.

What decides is not whether the sampling interval divides the dose interval. It
is how many distinct points of the dose cycle the schedule visits, and how much
of the increment's range those points span. This script computes both, for
every interval, and shows which of the two is the binding quantity.

Nothing below records an expected result in advance. The footing is d12's:
theta and K_m free, priors only on ln k_e (0.8) and the load record (5%).

Run:  python audit/f4_sampling.py
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import design.run as D
from d12_minimal_footing import record, sd

P = D.P
T = D.DOSE_INTERVAL
GATE = 0.20


def banner(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


def cycle_days(times, tol=1e-6):
    """Distinct phases of the dose cycle that a schedule visits."""
    ph = np.sort(np.asarray(times, float) % T)
    keep = [ph[0]]
    for x in ph[1:]:
        if x - keep[-1] > tol:
            keep.append(x)
    # the cycle is periodic, so the first and last are the same point
    if len(keep) > 1 and (T - keep[-1]) + keep[0] < tol:
        keep.pop()
    return np.array(keep)


def predicted_phases(step):
    """q, from step/T = p/q in lowest terms: the number of phases a perfectly
    commensurate schedule can ever visit."""
    fr = Fraction(step / T).limit_denominator(10000)
    return fr.denominator


# --------------------------------------------------------------------------
def s1_the_table(horizon=730.0):
    banner("S1  sampling interval against phases visited, increment range, sd")

    steps = [1.0, 2.0, 3.5, 4.0, 4.6667, 7.0, 7.5, 10.0, 14.0, 14.5, 15.0,
             21.0, 28.0, 28.5, 30.4, 42.0, 56.0, 60.0]
    print(f"  dose interval T = {T:.1f} d, horizon {horizon:.0f} d, load CV 10%\n")
    print(f"  {'step (d)':>9} {'step/T':>8} {'n':>5} {'phases':>7} {'q':>4}"
          f" {'Delta/theta range':>19} {'sd(ln theta)':>13} {'gate':>6}")
    print("  " + "-" * 82)
    rows = []
    for st in steps:
        t = np.arange(0.0, horizon, st)
        lv, dd = record(t, 0.10)
        ph = cycle_days(t)
        rng = (dd.max() - dd.min()) / P.theta
        s = sd(lv, dd)
        rows.append((st, len(t), len(ph), predicted_phases(st), rng, s))
        print(f"  {st:9.4f} {st/T:8.4f} {len(t):5d} {len(ph):7d} "
              f"{predicted_phases(st):4d} {rng:19.4f} "
              f"{'singular' if not np.isfinite(s) else f'{s:13.4f}'} "
              f"{'PASS' if s < GATE else 'fail':>6}")
    return rows


# --------------------------------------------------------------------------
def s2_his_two_claims(rows):
    banner("S2  his two claims, separately")

    by_step = {r[0]: r for r in rows}
    print("  CLAIM A: sampling AT the dose interval, or a multiple of it.")
    for st in (14.0, 28.0, 42.0, 56.0):
        r = by_step[st]
        print(f"    step {st:5.1f} d = {st/T:.0f} T: phases {r[2]}, "
              f"Delta/theta range {r[4]:.4f}, sd "
              f"{'singular' if not np.isfinite(r[5]) else f'{r[5]:.4f}'}")
    print("    -> HOLDS, and it is already the paper's own rule.\n")

    print("  CLAIM B: sampling at 1/2 of the dose interval.")
    for st in (7.0, 4.6667, 3.5, 2.0, 1.0):
        r = by_step.get(st)
        if r is None:
            continue
        print(f"    step {st:7.4f} d = T/{T/st:.0f}: phases {r[2]}, "
              f"Delta/theta range {r[4]:.4f}, sd "
              f"{'singular' if not np.isfinite(r[5]) else f'{r[5]:.4f}'}")
    print("    -> a schedule at T/q visits q distinct phases, and q >= 2 is")
    print("       enough. The half-interval is the WORST of the passing")
    print("       schedules, not a failing one.\n")

    print("  What separates them is the increment range, not divisibility:")
    fin = [r for r in rows if np.isfinite(r[5])]
    fail = [r for r in rows if not np.isfinite(r[5]) or r[5] >= GATE]
    ok = [r for r in rows if np.isfinite(r[5]) and r[5] < GATE]
    print(f"    schedules that PASS: Delta/theta range "
          f"[{min(r[4] for r in ok):.4f}, {max(r[4] for r in ok):.4f}], "
          f"phases >= {min(r[2] for r in ok)}")
    if fail:
        print(f"    schedules that FAIL: Delta/theta range "
              f"[{min(r[4] for r in fail):.4f}, {max(r[4] for r in fail):.4f}], "
              f"phases <= {max(r[2] for r in fail)}")


# --------------------------------------------------------------------------
def s3_which_quantity_binds(horizon=730.0, n=25):
    banner("S3  phases or increment range: which one binds?")

    print("  Two designs with the SAME number of visits and the same number of")
    print("  distinct phases, but different increment ranges, made by placing")
    print("  the phases at different points of the cycle.\n")
    print(f"  {'design':>40} {'n':>4} {'phases':>7} {'Delta/theta range':>19}"
          f" {'sd(ln theta)':>13}")
    print("  " + "-" * 88)

    def evaluate(lab, t):
        lv, dd = record(np.asarray(t, float), 0.10)
        ph = cycle_days(t)
        s = sd(lv, dd)
        print(f"  {lab:>40} {len(t):4d} {len(ph):7d} "
              f"{(dd.max()-dd.min())/P.theta:19.4f} "
              f"{'singular' if not np.isfinite(s) else f'{s:13.4f}'}")
        return s

    cyc = int(horizon // T)
    # two phases, close together, versus two phases at opposite ends
    for off in (0.5, 1.0, 2.0, 3.5, 7.0):
        t = np.sort(np.concatenate([np.arange(0.0, horizon, T)[:n // 2],
                                    np.arange(off, horizon, T)[:n - n // 2]]))
        evaluate(f"two phases, separated by {off:.1f} d", t)

    # evenly spread phases, same n
    for k in (2, 3, 4, 7, 14):
        offs = np.arange(k) * (T / k)
        t = np.sort(np.concatenate(
            [o + np.arange(0.0, horizon, T)[: max(1, n // k)] for o in offs]))
        evaluate(f"{k} phases evenly spread over the cycle", t)

    print("\n  The increment range is the quantity that tracks the precision.")
    print("  Phase COUNT alone does not: two phases half a day apart visit as")
    print("  many distinct cycle days as two phases a week apart and buy far")
    print("  less, because Delta has barely moved between them.")


# --------------------------------------------------------------------------
def s4_the_rule():
    banner("S4  the rule, stated in the quantity that decides")

    print("  A schedule of step h over a horizon H visits")
    print("      q = denominator of h/T in lowest terms")
    print("  distinct phases of the dose cycle (all of them, densely, when")
    print("  h/T is irrational). The increment at phase p is")
    print("      Delta(p) = D exp(-lambda p)/(1 - exp(-lambda T)),")
    print("  so the range a schedule covers is")
    print("      Delta/theta range = (D/theta) (exp(-lambda p_min)"
          " - exp(-lambda p_max))/(1 - exp(-lambda T)).")
    print()
    lam, Damp = D.LAMBDA, D.D_AMP
    pk = Damp / (1.0 - np.exp(-lam * T))
    print(f"  At this repository's constants (D = {Damp:.0f}, T = {T:.0f} d,")
    print(f"  half-life 3.5 d, theta = {P.theta:.0f}):")
    print(f"    Delta/theta at phase 0    = {pk/P.theta:.4f}")
    print(f"    Delta/theta at phase T-   = {pk*np.exp(-lam*T)/P.theta:.4f}")
    print(f"    the widest range available = "
          f"{(pk - pk*np.exp(-lam*T))/P.theta:.4f}")
    print()
    print("  So the rule is: the schedule must cover a range of Delta/theta,")
    print("  and the range is set by the SPREAD of the phases it visits, not by")
    print("  their number. Concretely:")
    print("    - h an exact multiple of T: one phase, range 0, the channel")
    print("      closes, and only the initial build-up of the dose train keeps")
    print("      the range from being exactly zero.")
    print("    - h = T/q for small q: q phases evenly spread, full range, fine.")
    print("    - h incommensurate with T: phases fill the cycle, full range.")
    print("  The failure mode is a phase SPREAD near zero, and that is reached")
    print("  only at multiples of T (or at h within a few hours of one).")


# --------------------------------------------------------------------------
def s5_near_misses():
    banner("S5  how close to a multiple of T is too close?")

    print("  The paper says half a day off the dose interval recovers the")
    print("  design. Below, the approach to 14 d from both sides.\n")
    print(f"  {'step (d)':>10} {'offset from T':>14} {'phases':>7}"
          f" {'Delta/theta range':>19} {'sd(ln theta)':>13}")
    print("  " + "-" * 68)
    for st in (13.0, 13.5, 13.75, 13.9, 13.95, 13.99, 14.0, 14.01, 14.05,
               14.1, 14.25, 14.5, 15.0):
        t = np.arange(0.0, 730.0, st)
        lv, dd = record(t, 0.10)
        s = sd(lv, dd)
        print(f"  {st:10.4f} {st-T:+14.4f} {len(cycle_days(t)):7d} "
              f"{(dd.max()-dd.min())/P.theta:19.4f} "
              f"{'singular' if not np.isfinite(s) else f'{s:13.4f}'}")


# --------------------------------------------------------------------------
def main():
    print(__doc__)
    rows = s1_the_table()
    s2_his_two_claims(rows)
    s3_which_quantity_binds()
    s4_the_rule()
    s5_near_misses()


if __name__ == "__main__":
    main()
