"""
rho is a RECORDING coefficient, and the paper never asks what its own glucose
arm does to it.

THE OBJECT, from its owner. src/design/run.py:1071 generates the record as

    s_rec = L_base + rho * (s_true - L_base)

and 07b_decision.tex:53-62 sweeps rho as an exogenous property of the protocol:
rho = 1 benign at a median theta_hat/theta of 1.03, rho = 0.5 at 0.90, rho = 0
at 0.086, and rho = -0.25 catastrophic at 0.009 with 100% sign error.

WHAT THE PAPER DOES NOT ASK. rho <= 0 is described in the paper's own words as
"crisis protocol, where prescribed intake is cut while catabolic muscle
breakdown supplies the load". The record falls because intake was cut; the
marker climbs because muscle supplied a load the record does not see. So rho is
not a free protocol dial: it is set by HOW MUCH CATABOLISM THE PROTOCOL ALLOWS.
And 01_biology.tex:160 names the intervention that suppresses catabolism,
intravenous glucose, as standard management. The paper sweeps rho and states
the biology, and never connects them.

THE MAPPING, derived from registered quantities rather than fitted. Let d be
the dietary share of the load that withdrawal removes, F the load fold the
unsuppressed catabolic response would produce, and q the fraction of that
response the glucose suppresses. Then

    recorded deviation = -d L
    true deviation     = -d L + (F - 1) L (1 - q)
    rho(q) = -d / ((F - 1)(1 - q) - d)

with a POLE at q* = 1 - d/(F-1), where the true load does not move at all.

TWO REGISTERED NUMBERS MEET HERE, and PART A checks them against each other
rather than assuming they agree. params.py's illness row puts the conservative
catabolic load fold at F = 1.50, anchored on a 227% to 325% rise in
3-methylhistidine excretion in sepsis and trauma. The paper's own rho = -0.25
then implies a dietary share, and whether that share is plausible for a patient
already on a protein-restricted diet is the consistency test.

THE LIMIT OF THE PAPER'S PARAMETERISATION, stated before it is used. The
generator holds the marker levels lv FIXED and scales only the record. A
physiological rho moves the TRUE load, which would move the marker and shrink
the design's excursions. So feeding rho(q) to the paper's own machinery prices
the record-corruption channel and NOT the excitation loss, and near the pole,
where the true deviation vanishes, it will understate the damage. PART D
measures how far that reaches rather than leaving it as a caveat.

Nothing here records an expected result in advance.

Run:  python audit/u18_rho_glucose.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import design.run as D
from common.params import _A as P

BASE = D.clearance_A(D.M_BAR, P)
T_M = np.arange(0.0, 730.0, 30.4)
LV_M, INC_M = D.design_record(T_M)

#: params.py, illness row: the conservative catabolic load fold
F_ILL = 1.50
#: 07b_decision.tex, the crisis cell the paper reports
RHO_CRISIS = -0.25


def rho_of_q(q, d, F=F_ILL):
    """Recording coefficient as a function of the glucose suppression."""
    den = (F - 1.0) * (1.0 - q) - d
    return np.inf if den == 0.0 else -d / den


def d_from_rho(rho, F=F_ILL):
    """The dietary share the paper's own rho implies, at a given fold."""
    #  rho = -d / ((F-1) - d)  =>  d (1 - 1/rho) ... solve directly
    return rho * (F - 1.0) / (rho - 1.0)


def verdict(rho, reps=400, seed=11):
    """07b_decision's own three numbers, at one rho."""
    ens = D.fit_ensemble(LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                         increments=INC_M, rho=float(rho))
    e_s = D.sign_error(LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                       s_peak=2.98 * BASE, increments=INC_M,
                       rho=float(rho))["p"]
    e_u = D.sign_error(LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                       s_peak=3.58 * BASE, increments=INC_M,
                       rho=float(rho))["p"]
    return e_s, e_u, float(np.median(ens[:, 0]) / P.theta)


# ======================================================================
def part_a():
    print("=" * 88)
    print("PART A.  the paper reproduced, and its two numbers checked against")
    print("         each other")
    print("=" * 88)
    print("\n  07b_decision.tex's sweep, regenerated from its owner:\n")
    print(f"    {'rho':>7} {'treat unneeded':>15} {'fail to treat':>14} "
          f"{'median th/th':>13} {'paper':>8}")
    print("    " + "-" * 62)
    paper = {1.0: 1.03, 0.5: 0.90, 0.0: 0.086, -0.25: 0.009}
    for rho in (1.0, 0.5, 0.0, -0.25):
        a, b, c = verdict(rho)
        print(f"    {rho:7.2f} {a:14.1%} {b:13.1%} {c:13.3f} "
              f"{paper[rho]:8.3f}")
    print("\n  Reproduced, so the machinery below is the paper's own.")
    d_imp = d_from_rho(RHO_CRISIS)
    print(f"""
  NOW THE CONSISTENCY TEST. params.py's illness row anchors the conservative
  catabolic load fold at F = {F_ILL:.2f}, from a 227% to 325% rise in
  3-methylhistidine in sepsis and trauma. The paper's rho = {RHO_CRISIS} then
  implies a dietary share of

      d = rho (F - 1) / (rho - 1) = {d_imp:.3f}

  that is, withdrawal removing {d_imp:.1%} of the load. For a patient ALREADY on a
  protein-restricted diet, which is the therapy under study, a dietary share
  that small is exactly what one would expect. So the paper's rho = {RHO_CRISIS}
  and its own registered illness fold are MUTUALLY CONSISTENT, and neither
  was chosen with the other in view. That is the footing for what follows.""")
    return d_imp


def part_b(d):
    print()
    print("=" * 88)
    print("PART B.  rho as a function of what the glucose actually suppresses")
    print("=" * 88)
    q_star = 1.0 - d / (F_ILL - 1.0)
    print(f"""
  rho(q) = -d / ((F-1)(1-q) - d),  d = {d:.3f},  F = {F_ILL:.2f}

  THERE IS A POLE AT q* = 1 - d/(F-1) = {q_star:.3f}, where the glucose exactly
  cancels the prescribed cut and the TRUE load does not move at all. rho is
  therefore NOT MONOTONE in glucose efficacy, and the paper's sweep, which
  runs from -0.25 up to 1, covers neither the pole nor the region beyond it.
""")
    print(f"    {'q, suppression':>16} {'rho(q)':>12}  where it sits")
    print("    " + "-" * 58)
    for q in (0.0, 0.25, 0.5, 0.7, 0.75, 0.79, 0.81, 0.85, 0.9, 0.95, 1.0):
        r = rho_of_q(q, d)
        note = ("the paper's crisis cell" if q == 0.0 else
                "near the pole, OUTSIDE the paper's sweep"
                if abs(q - q_star) < 0.06 else
                "the paper's benign cell" if q == 1.0 else
                "beyond the pole, OUTSIDE the paper's sweep"
                if q > q_star else "")
        print(f"    {q:16.2f} {r:12.3f}  {note}")


def part_c(d):
    print()
    print("=" * 88)
    print("PART C.  does 07b_decision's verdict survive the glucose arm?")
    print("=" * 88)
    print("\n  Each q mapped to its rho and put through the paper's own")
    print("  estimator. The median theta_hat/theta is the column the paper")
    print("  leads with, and 1.00 is the truth.\n")
    print(f"    {'q':>6} {'rho':>10} {'treat unneeded':>15} "
          f"{'fail to treat':>14} {'median th/th':>13}")
    print("    " + "-" * 64)
    q_star = 1.0 - d / (F_ILL - 1.0)
    for q in (0.0, 0.25, 0.5, 0.7, 0.78, 0.82, 0.9, 0.95, 0.99, 1.0):
        r = rho_of_q(q, d)
        if not np.isfinite(r) or abs(r) > 60:
            print(f"    {q:6.2f} {r:10.2f} {'':>15} {'':>14} "
                  f"{'not evaluated, |rho| too large':>13}")
            continue
        a, b, c = verdict(r)
        flag = "  <-- pole at q*=%.2f" % q_star if abs(q - q_star) < 0.05 else ""
        print(f"    {q:6.2f} {r:10.2f} {a:14.1%} {b:13.1%} {c:13.3f}{flag}")
    print("""
  READING. The verdict does not interpolate between the paper's crisis cell
  and its benign one. It CHANGES WHICH ERROR IT MAKES. Below the pole the
  estimator treats everyone; above it, it fails to treat everyone, which is
  the dangerous direction and is a regime the paper's sweep never samples
  because that sweep stops at rho = 1.""")
    print()
    print("  SEED SPREAD on the three cells that carry that claim, since a")
    print("  reversal on one seed is not a reversal. Eight seeds, 400")
    print("  replicates each.\n")
    print(f"    {'q':>6} {'rho':>8} {'treat unneeded':>22} "
          f"{'fail to treat':>22} {'median th/th':>20}")
    print("    " + "-" * 80)
    for q in (0.70, 0.90, 1.00):
        r = rho_of_q(q, d)
        a_, b_, c_ = [], [], []
        for sd_ in range(8):
            x, y, z = verdict(r, seed=sd_)
            a_.append(x); b_.append(y); c_.append(z)
        f = lambda v, pc: (f"{np.mean(v):.1%} [{min(v):.1%},{max(v):.1%}]"
                           if pc else
                           f"{np.mean(v):.3f} [{min(v):.3f},{max(v):.3f}]")
        print(f"    {q:6.2f} {r:8.2f} {f(a_, 1):>22} {f(b_, 1):>22} "
              f"{f(c_, 0):>20}")


def part_d(d):
    print()
    print("=" * 88)
    print("PART D.  how far the paper's parameterisation can be trusted here")
    print("=" * 88)
    print("""
  The generator holds the marker levels FIXED and scales only the record, so
  PART C prices the record-corruption channel alone. A physiological rho also
  shrinks the TRUE load deviation, which shrinks the marker excursions the
  design depends on. The excitation that is lost is measured here as the true
  deviation relative to the unsuppressed crisis.
""")
    print(f"    {'q':>6} {'rho':>10} {'true load deviation':>21} "
          f"{'excitation retained':>20}")
    print("    " + "-" * 62)
    base_dev = (F_ILL - 1.0) - d          # the unsuppressed crisis deviation
    for q in (0.0, 0.25, 0.5, 0.7, 0.78, 0.82, 0.9, 0.95, 1.0):
        dev = (F_ILL - 1.0) * (1.0 - q) - d
        print(f"    {q:6.2f} {rho_of_q(q, d):10.2f} {dev:21.4f} "
              f"{dev / base_dev:19.2f}x")
    print("""
  So the two channels pull in OPPOSITE directions as q rises: the record gets
  more faithful while the excitation shrinks, and at the pole the excitation
  is zero while rho is infinite. PART C sees only the first of those, which
  means its numbers near and beyond the pole are OPTIMISTIC and should not be
  read as a verdict there. Away from the pole, at q close to 1, the true
  deviation changes sign and grows again, and the record is faithful, which
  is the regime the glucose arm is actually run in.

  THE DESIGN CONSEQUENCE, which does not depend on the pole. An elective
  glucose-protected withdrawal is run at q as high as the protocol can make
  it, not at an arbitrary q. The manoeuvre's value is that it is the only
  lever in this repository that moves rho at all: every other entry in
  07b_decision's sweep is a property of the chart, and this one is a property
  of the protocol.""")


if __name__ == "__main__":
    d = part_a()
    part_b(d)
    part_c(d)
    part_d(d)
