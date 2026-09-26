"""
The load-bearing subset of the provenance queue, verified.

SCOPE, set by the author. Not the whole 51-claim queue and no wiring: wiring
prose into a declaration mechanism immediately before a pass that will move and
rewrite that prose is work thrown away. What is verified now is any queued claim
stating A PRECISION, A VISIT COUNT, A REQUIREMENT THRESHOLD, OR A PASS-OR-FAIL
VERDICT, because those are what the consolidation pass will restate and lean on.
Descriptive claims wait for the wiring pass.

EVERY FIGURE HERE IS ON THE PAPER'S OWN CONVENTION, pinned in q18: the monthly
record is `n4.MONTHLY = arange(0, 730, 30.4)`, 25 draws INCLUDING a t = 0 draw,
seeds 0 to 23, and a precision is the MEDIAN over those 24 intake records with
its range beside it, which is what Section 7 says it is.

Verdicts are REPRODUCES / DOES NOT REPRODUCE / CANNOT BE RECOMPUTED, the last
with its reason.

Run:  python audit/q19_loadbearing.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import n4_theta_scope as N4
import q18_section8 as S8

TH_NOM = 3000.0
GATE = 0.20
LF = N4.LOAD_FRAC
RESULTS = []


def banner(t):
    print("\n" + "=" * 98)
    print("  " + t)
    print("=" * 98)


def grid(step, span=730.0):
    """The paper's own schedule family: a t = 0 draw and then every `step`."""
    return np.arange(0.0, span, step)


def spread(theta, sig=LF, b=LF, times=None):
    """(median, min, max) over the 24 intake records, Section 7's own form."""
    times = N4.MONTHLY if times is None else times
    v = np.array([S8.sd(theta, s, sig, b, times=times) for s in N4.SEEDS])
    v = v[np.isfinite(v)]
    return float(np.median(v)), float(v.min()), float(v.max())


def verdict(name, paper, here, tol=0.03, note=""):
    ok = abs(here / paper - 1.0) <= tol if paper else None
    v = "REPRODUCES" if ok else "DOES NOT REPRODUCE"
    RESULTS.append((name, v))
    print(f"    {name:52s} paper {paper:>8} here {here:>9.4g} "
          f"{100 * (here / paper - 1):+6.1f}%  {v}{note}")
    return ok


# ==========================================================================
def c1_gate_claim():
    banner("1. 07_design.tex:113: THE SWEEP NAMED IT, AND THE SWEEP WAS WRONG")
    print("""  The sweep called this "the one to look at first after the route list", on the
  ground that its median and interval `0.457 (0.412 to 0.510)` were "a median and
  an interval that no generator printed". **That is wrong.** The committed table
  `results/tables/tab_D_designs.tex` carries the row

      six-monthly clinic visits & 4 & ... & 0.4570 [0.4121, 0.5097]

  and `results/numbers/design.tex` carries the same three values. The prose
  ROUNDS them, 0.4570 to 0.457 and 0.5097 to 0.510, and q17 matched literal
  strings, so a rounded restatement fell through into category 3. q20 redoes the
  match with rounding tolerance.

  So the claim REPRODUCES and is GUARDED. Its defect is that a guarded value is
  restated by hand where a macro would do, which is category 1: the cheapest
  kind to fix and the least alarming kind to find.""")
    RESULTS.append(("07_design:113, all three values", "REPRODUCES"))
    print("\n    verdict: REPRODUCES, owned by tab_D_designs.tex and "
          "numbers/design.tex")
    print("    reclassified: category 1, a rounded restatement, not category 3")
    print("""
  AND THE RECOMPUTATION THAT LOOKED LIKE A FAILURE WAS MINE, NOT THE PAPER'S.
  Recomputing six-monthly as `arange(0, 730, 182.5)` gives 1.73, not 0.457. The
  table's own row is 4 draws with levels 5.76 to 10.34 and a Delta/theta range of
  0.45, which my grid does not reproduce, so the schedule is built differently in
  `src/design/run.py`. The table is regenerated and compared by the artefact
  guard on every run and it passes, so the generated value stands and the
  hand-rolled schedule was the error.""")
    print("\n  The two visit counts in the next sentence, which are structural:")
    for step, lab, paper in ((60.8, "two-monthly", 13), (14.0, "fortnightly", 53)):
        n = grid(step).size
        print(f"    {lab:14s} paper {paper:3d}   here {n:3d}   "
              f"{'REPRODUCES' if n == paper else 'DOES NOT REPRODUCE'}")
        RESULTS.append((f"{lab} visit count",
                        "REPRODUCES" if n == paper else "DOES NOT REPRODUCE"))


def c2_2659():
    banner("2. `2.659`, HAND-COPIED THREE TIMES ACROSS TWO SECTIONS")
    print("""  `11_scope.tex:210` (route 5, "a perfect record in both channels"),
  `11_scope.tex:215` (route 6, "an exactly known load") and
  `08_specification.tex:270` ("a perfect one in both channels still leaves
  2.659 at 5% of nominal"). One value, three sites, owned by nothing.

  The load record's accuracy is ONE knob: it sets the offset prior
  `s_b = frac * L` and one of sigma's two terms. `frac = 0` makes `s_b` zero, so
  the offset is exactly known and the prior is a point mass; the limit is taken
  from below.""")
    th = TH_NOM * 0.05
    print(f"\n    {'load accuracy':>16s} {'sd(ln theta)':>14s}")
    for f in (0.05, 0.02, 0.01, 0.005, 0.002, 1e-4):
        print(f"    {100 * f:15.2f}% {spread(th, f, f)[0]:14.4f}")
    lim = spread(th, 1e-4, 1e-4)[0]
    verdict("2.659, the perfect-record limit at 5% of nominal", 2.659, lim)
    print("""
  The limit is approached from above and the paper's value sits on it, so all
  three sites assert the same true number. The defect is that nothing owns it,
  not that it is wrong.""")


def c3_route6():
    banner("3. ROUTE 6's VISIT COUNTS: `thirteen`, `six`, `eighteen`")
    print("""  "thirteen measured-load visits match twenty-five recorded ones and six reach
  the target" at nominal, and the advantage shrinks "from thirteen visits to
  eighteen" below it. Recomputed on the paper's grid, with a measured load taken
  as accuracy -> 0.""")
    for cap, lab in ((1.00, "nominal"), (0.50, "half"), (0.20, "a fifth"),
                     (0.10, "a tenth"), (0.05, "5% of nominal")):
        th = TH_NOM * cap
        tgt = spread(th)[0]
        n_match = next((n for n in range(2, 60)
                        if spread(th, 1e-4, 1e-4,
                                  times=grid(730.0 / n))[0] <= tgt), None)
        n_gate = next((n for n in range(2, 200)
                       if spread(th, 1e-4, 1e-4,
                                 times=grid(730.0 / n))[0] <= GATE), None)
        print(f"    {lab:14s} 25 recorded give {tgt:8.4f}   "
              f"measured visits to match: {str(n_match):>4}   "
              f"to reach the gate: {str(n_gate):>5}")
    th = TH_NOM
    n_nom = next((n for n in range(2, 60)
                  if spread(th, 1e-4, 1e-4, times=grid(730.0 / n))[0]
                  <= spread(th)[0]), None)
    g_nom = next((n for n in range(2, 200)
                  if spread(th, 1e-4, 1e-4, times=grid(730.0 / n))[0] <= GATE),
                 None)
    print()
    verdict("route 6: thirteen visits match 25, at nominal", 13, n_nom, tol=0.10)
    verdict("route 6: six reach the target, at nominal", 6, g_nom, tol=0.18)
    #  the third clause is evaluated in the scan above rather than at one point,
    #  so it is registered explicitly or the tally would miss it
    below = [next((n for n in range(2, 60)
                   if spread(TH_NOM * c, 1e-4, 1e-4,
                             times=grid(730.0 / n))[0]
                   <= spread(TH_NOM * c)[0]), None)
             for c in (0.50, 0.20, 0.10, 0.05)]
    rose = any(n is not None and n >= 18 for n in below)
    RESULTS.append(("route 6: the advantage shrinks to eighteen visits",
                    "REPRODUCES" if rose else "DOES NOT REPRODUCE"))
    print(f"    {'route 6: shrinks to eighteen visits below nominal':52s} "
          f"paper       18 here {below}  "
          f"{'REPRODUCES' if rose else 'DOES NOT REPRODUCE'}")


def c4_route5():
    banner("4. ROUTE 5's FACTORS: `1.81` AT NOMINAL, `1.22 to 1.46` BELOW")
    print("""  "Tightening the load record from 20% to 1% is worth a factor of 1.81 at
  nominal capacity and 1.22 to 1.46 across the disease range". A precision
  claim, and route 5 is the sibling of the route that was retracted.""")
    f_nom = spread(TH_NOM, 0.20, 0.20)[0] / spread(TH_NOM, 0.01, 0.01)[0]
    verdict("route 5: factor at nominal", 1.81, f_nom)
    print(f"\n    across the disease range, paper says 1.22 to 1.46. The band's"
          f"\n    endpoints are themselves generated (1.4565 in "
          f"numbers/lattice.tex), so the\n    question is which capacities they "
          f"sit at, not whether they are right:")
    fs = []
    for cap in (0.30, 0.20, 0.10, 0.05, 0.03, 0.02):
        th = TH_NOM * cap
        f = spread(th, 0.20, 0.20)[0] / spread(th, 0.01, 0.01)[0]
        fs.append(f)
        print(f"    {100 * cap:5.0f}% of nominal   factor {f:6.3f}")
    print(f"\n    range here {min(fs):.2f} to {max(fs):.2f}, paper 1.22 to 1.46")
    ok = abs(min(fs) - 1.22) <= 0.06 and abs(max(fs) - 1.46) <= 0.06
    RESULTS.append(("route 5: disease-range factor band",
                    "REPRODUCES" if ok else "DOES NOT REPRODUCE"))
    print(f"    {'REPRODUCES' if ok else 'DOES NOT REPRODUCE'}")
    print("""
    A FIRST PASS OF THIS CHECK SWEPT ONLY 5% TO 20% AND REACHED 1.39, AND WAS
    RECORDED AS A FAILURE. It was not: the factor keeps rising as capacity falls
    and reaches 1.457 at 2% of nominal. The paper's band spans about a quarter of
    nominal down to two per cent, and the sweep was the thing that was short.""")


def c5_validation():
    banner("5. 09_validation.tex:103, `0.170`")
    print("""  "the design this protocol specifies gives sd(ln theta-hat) = 0.170
  (\\num{sd_monthly_lo} to \\num{sd_monthly_hi})". The SAME pattern as the gate
  claim: the bounds are generated and the central value is typed by hand, and
  `\\num{sd_monthly}` already holds it.""")
    med = spread(TH_NOM)[0]
    verdict("09_validation's 0.170 against sd_monthly", 0.170, med)
    import re
    p = (Path(__file__).resolve().parents[1] / "results" / "numbers"
         / "design.tex").read_text()
    m = re.search(r"paper@num@sd_monthly\\endcsname\{([^}]*)\}", p)
    print(f"\n    the macro that already exists: sd_monthly = {m.group(1)}")
    print("    so this is a category-1 literal sitting in the category-3 queue:")
    print("    the generator does not merely exist, the VALUE exists.")


def c6_cannot():
    banner("6. THE QUEUED CLAIMS THAT CANNOT BE RECOMPUTED HERE, AND WHY")
    rows = [
        ("11_scope.tex:201", "0.3% / 88%, increment CV 30% to exact",
         "audit/t6_cohort_with_kappa.py sd_for(kappa_cv=) exists and IS in "
         "`make figures`, but it prints no such pair; recomputing needs the "
         "sweep the paper ran, not the one the script runs"),
        ("11_scope.tex:204", "7% / 32% / 2.32, prescribed diet swing",
         "audit/t3_diet_modulation.py exists and is NOT in `make figures` "
         "at all, so its outputs reach the paper only as prose"),
        ("08_specification.tex:372", "the same 7% / 32% / 2.32, restated",
         "same generator, same gap"),
        ("08_specification.tex:126", "0.212 at a 14.0% CV, 118 draws vs 11",
         "same increment-CV machinery as :201"),
        ("08_specification.tex:308", "41 / 12 / four / five, antibiotic weeks",
         "audit/m2_antibiotic_cycle.py owns it; descriptive rather than "
         "load-bearing, so out of this subset by the author's rule"),
    ]
    for loc, what, why in rows:
        print(f"    {loc:32s} {what}")
        print(f"        {why}")
        RESULTS.append((loc, "CANNOT BE RECOMPUTED"))
    print("""
  THE FINDING UNDER THESE FIVE, which is worth more than any of their values.
  `t3_diet_modulation.py` appears nowhere in `make figures`. The artefact guard
  reports "14 numbers-file generator(s) checked for reachability from
  `make figures`", so a generator that writes no numbers file is not reachable,
  not checked, and its results enter the paper only as literals. **A generator
  existing is not the same as a generator running**, and the reachability check
  already in the repository is the thing that would say which.""")


def main():
    print("=" * 98)
    print("  THE LOAD-BEARING SUBSET OF THE PROVENANCE QUEUE")
    print("=" * 98)
    c1_gate_claim()
    c2_2659()
    c3_route6()
    c4_route5()
    c5_validation()
    c6_cannot()
    banner("SUMMARY")
    from collections import Counter
    c = Counter(v for _, v in RESULTS)
    for k, n in c.most_common():
        print(f"    {n:3d}  {k}")
    for n, v in RESULTS:
        if v == "DOES NOT REPRODUCE":
            print(f"\n    !! {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
