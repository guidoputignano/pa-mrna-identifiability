"""
Round 19: can the on-therapy record settle WHICH SIDE of the 21-day crossing a
patient is on?

R-O and R-P said the record "reads lambda at 7-18 sigma" and therefore settles the
two sides of the 3.5452 d crossing "a fortiori". The 7-18 sigma is the separation
between a 3.6 d and a 7 d half-life, ln(7/3.6) = 0.665. The distance from the
measured 3.6 d to the crossing is ln(3.6/3.5452) = 0.015, about 43 times smaller.
So placing a value near 3.6 d on one side of the crossing is much HARDER than
telling 3.6 d from 7 d, not easier. This script prices how much harder.

Input: the per-patient sd(ln lambda) table that t46_timing.py prints (FIM, by
record length and cadence). It is parsed from out/t46_timing.txt, not retyped, so
the numbers here move if t46 is rerun.

Outputs, per patient:
  * z = ln(t_true / 3.5452) / sd(ln lambda) at the measured 3.6 d and at the
    round-17 central estimate 4.3 d
  * the 2-sigma "dead zone": true half-lives a single patient's record cannot
    place on either side of the crossing
  * the first record length at which a 4.3 d truth is placed at 2 sigma
  * how many patients a POOLED estimate would need to place a 3.6 d truth at
    2 sigma, if (and only if) lambda were the same in every patient
"""
import re
from pathlib import Path

import numpy as np

CROSS = 3.5452          # 21-day crossing (R-N, t48), Model C intrinsic theta_crit
MEASURED = 3.6          # Chandler & Ballard 1985, human fibroblast PCC
CENTRAL = 4.3           # round-17 compartment-anchored liver estimate (R-P, t50)
MOUSE = 7.0             # Attarwala 2023 effective, the comparator of R-N's 7-18 sigma
Z = 2.0

HERE = Path(__file__).resolve().parent


def t46_table():
    """{cadence: {months: sd(ln lambda)}} from t46's ITEM 2 table."""
    txt = (HERE / "out" / "t46_timing.txt").read_text()
    block = txt.split("ITEM 2.")[1].split("spec reached by")[0]
    rows = re.findall(r"^\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$", block, re.M)
    if not rows:
        raise RuntimeError("could not parse the sd(ln lambda) table from out/t46_timing.txt")
    out = {"monthly": {}, "fortnightly": {}, "weekly": {}}
    for mo, a, b, c in rows:
        out["monthly"][int(mo)] = float(a)
        out["fortnightly"][int(mo)] = float(b)
        out["weekly"][int(mo)] = float(c)
    return out


def main():
    tab = t46_table()
    g_meas = np.log(MEASURED / CROSS)
    g_cent = np.log(CENTRAL / CROSS)
    g_mouse = np.log(MOUSE / MEASURED)
    print("=" * 92)
    print("  ROUND 19. CAN THE RECORD SETTLE WHICH SIDE OF THE 3.5452 d CROSSING? (per patient)")
    print("=" * 92)
    print(f"""
  log-distances:  3.6 d vs 7 d (R-N's 7-18 sigma)      ln = {g_mouse:.3f}
                  3.6 d vs the 3.5452 d crossing       ln = {g_meas:.4f}   ({g_mouse / g_meas:.0f}x smaller)
                  4.3 d vs the 3.5452 d crossing       ln = {g_cent:.3f}
  sd(ln lambda) per patient: parsed from out/t46_timing.txt (FIM).
""")
    print(f"  {'cadence':>12} {'months':>6} {'sd(ln lam)':>10} {'z @3.6d':>8} {'z @4.3d':>8}"
          f" {'z 3.6 vs 7':>10}   2-sigma dead zone (true half-life, d)")
    for cad in ("monthly", "fortnightly", "weekly"):
        for mo, s in sorted(tab[cad].items()):
            lo, hi = CROSS * np.exp(-Z * s), CROSS * np.exp(Z * s)
            print(f"  {cad:>12} {mo:>6} {s:>10.3f} {g_meas / s:>8.2f} {g_cent / s:>8.2f}"
                  f" {g_mouse / s:>10.1f}   {lo:.2f} - {hi:.2f}")
    best = min(min(v.values()) for v in tab.values())
    print(f"""
  Best per-patient precision in the table: sd(ln lambda) = {best:.3f}  ->  z at 3.6 d = {g_meas / best:.2f}.
  A single patient's record never places a 3.6 d truth on a side of the crossing: even
  after 5 years of weekly draws it sits at {g_meas / best:.1f} sigma, inside the dead zone.
""")
    print("  First record length at which a 4.3 d truth is placed at >= 2 sigma, per patient:")
    for cad in ("monthly", "fortnightly", "weekly"):
        hit = [mo for mo, s in sorted(tab[cad].items()) if g_cent / s >= Z]
        print(f"    {cad:>12}: {str(hit[0]) + ' months' if hit else 'not within 60 months'}"
              f"   (z = {g_cent / tab[cad][hit[0]]:.2f})" if hit else "")
    print("\n  Pooling (ONLY if lambda is common to every patient, which nothing measured shows):")
    for cad, mo in (("monthly", 12), ("weekly", 12), ("weekly", 24)):
        s = tab[cad][mo]
        n = int(np.ceil((Z * s / g_meas) ** 2))
        print(f"    {cad:>8} draws, {mo:>2} months: {n:>4} patients to place a 3.6 d truth at 2 sigma")
    print(f"""
  READING. R-N's 7-18 sigma is right for its question (3.6 d or 7 d?). It does not carry
  to the crossing: the crossing is {g_mouse / g_meas:.0f}x closer, so "a fortiori" runs the wrong way.
  Per patient the record settles the side only when the true liver half-life is clear of
  the crossing: about 6 months of draws if it is near the 4.3 d central estimate, never
  if it is near the measured 3.6 d. A pooled trial estimate could do it with a few dozen
  patients, but only if lambda does not vary between patients, which is unmeasured. So
  the two emails still matter, and the on-therapy record is still the best route, but
  it settles the knife-edge only if the answer is not on the knife-edge.
""")


if __name__ == "__main__":
    main()
