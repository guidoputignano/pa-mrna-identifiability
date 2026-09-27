"""
The register's own `k_e` requirement against the frontier, at the capacities
where the clinical question is live.

THE CLAIM UNDER AUDIT. `src/common/params.py:206`, the `k_e` row of the
register that GENERATES the paper's parameter table, says of the parameter it
calls "the priority gap":

    "Section 7 needs it within a factor of about three, and no independent
     measurement in humans supplies it."

`audit/u9_ke_frontier.py` computes the same quantity at other capacities and
gets 0.1567 and 0.0633 in log units. A factor of three is a log-sd of
ln 3 = 1.0986. So the question is whether the register understates its own
requirement, and by how much.

WHAT IS AND IS NOT AT ISSUE. The number is not wrong. `07_design.tex:60` is
already WIRED to the generated macro and prints 2.67, which this script
reproduces from the generator; the register's "about three" is that value
rounded, well inside the repository's own rounding convention. What is at
issue is that the requirement is stated with NO CAPACITY QUALIFIER, in the
row a reader consults to learn what the priority gap needs, while the
generator that owns it evaluates it at nominal capacity alone.

THE OWNER, located before recomputing this time. `src/design/run.py`
`spec_numbers()` solves

    brentq(lambda x: sd_of(monthly, ke_prior=x, seed=s) - 0.20, 1e-3, 3.0)

over seeds and reports the median. `monthly` is the 25-draw record and
`sd_of` runs at the module's default `P.theta`, which is 3000. So the
published figure is the nominal-capacity requirement on the recommended
record, and that is exactly `audit/u9_ke_frontier.py`'s own "recommended
record, theta = 3000" cell. PART A checks that the two agree, which is the
anchor for everything after.

Nothing here records an expected result in advance.

Run:  python audit/u13_ke_requirement_audit.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import u9_ke_frontier as U9

LN3 = float(np.log(3.0))


def part_a():
    print("=" * 84)
    print("PART A.  the published requirement, from its own generator")
    print("=" * 84)
    import design.run as D
    N = D.spec_numbers()
    pub_log = float(N["ke_needed"])
    pub_fac = float(np.exp(pub_log))
    print(f"\n  src/design/run.py spec_numbers(), theta = {D.P.theta:.0f}:")
    print(f"    ke_needed, as a log sd   {pub_log:.4f}  "
          f"({N['ke_needed_lo']:.4f} to {N['ke_needed_hi']:.4f})")
    print(f"    as a factor              {pub_fac:.3f}")
    print(f"    the register's prose     'about three', ln 3 = {LN3:.4f}")
    print()
    mine, note = U9.frontier(3000.0, 25, None)
    print(f"  audit/u9_ke_frontier.py, same capacity, same record length:")
    print(f"    frontier, as a log sd    "
          f"{(mine if mine is not None else float('nan')):.4f}")
    if mine:
        print(f"    difference from the generator  "
              f"{abs(mine - pub_log):.4f}  ({100 * abs(mine - pub_log) / pub_log:.2f}%)")
    print("\n  Two code paths, one quantity. That is the anchor.")
    return pub_log, pub_fac


def part_b(pub_log, pub_fac):
    print()
    print("=" * 84)
    print("PART B.  the same requirement at the capacities the question is about")
    print("=" * 84)
    print("\n  The generator evaluates at nominal. The frontier evaluates")
    print("  everywhere. Tier 1 record, which is the arm that needs a decade.\n")
    print(f"    {'theta':>7} {'% nominal':>10} {'frontier log sd':>16} "
          f"{'as a factor':>12} {'vs ln 3':>9} {'vs the generator':>17}")
    print("    " + "-" * 76)
    rows = []
    for theta in (3000.0, 600.0, 300.0, 150.0):
        s, note = U9.frontier(theta, U9.TIER1_N, U9.PUBLISHED_KAPPA_CV)
        if s is None:
            print(f"    {theta:7.0f} {100 * theta / 3000:9.1f}% "
                  f"{'0.8 suffices':>16}")
            continue
        if isinstance(s, float) and np.isnan(s):
            print(f"    {theta:7.0f} {100 * theta / 3000:9.1f}% "
                  f"{'unreachable':>16}")
            rows.append((theta, None))
            continue
        print(f"    {theta:7.0f} {100 * theta / 3000:9.1f}% {s:16.4f} "
              f"{np.exp(s):12.3f} {LN3 / s:8.1f}x {pub_log / s:16.1f}x")
        rows.append((theta, s))
    print()
    print("  The last two columns are how many times TIGHTER the requirement")
    print("  is than the register states it, anchored two ways: on the prose")
    print(f"  'about three' (ln 3 = {LN3:.4f}) and on the generated value")
    print(f"  ({pub_fac:.3f}, log sd {pub_log:.4f}).")
    return rows


def part_c(rows, pub_log):
    print()
    print("=" * 84)
    print("PART C.  what the correction is, stated precisely")
    print("=" * 84)
    print("""
  IT IS NOT THAT THE NUMBER IS WRONG. The generator is correct, the macro is
  wired, and 07_design.tex:60 prints it with its seed range. The register's
  "about three" is that value rounded.

  IT IS THAT THE REQUIREMENT HAS A CAPACITY AXIS AND IS STATED WITHOUT ONE,
  in the row a reader consults to learn what the parameter the paper calls
  its priority gap actually needs. The same row already says the design does
  not reach the disease range, and the two statements do not connect. The
  frontier connects them: the disease range IS reachable, and what it needs
  is the number in PART B.

  THE READER THIS MISLEADS is the one costing the measurement. "Within a
  factor of about three" reads as a modest ask that an analogue or a
  literature bound might satisfy. The requirement where the clinical
  question lives does not read that way at all.
""")
    live = [(t, s) for t, s in rows if s is not None]
    if live:
        lo = min(LN3 / s for _, s in live)
        hi = max(LN3 / s for _, s in live)
        print(f"  Across the live capacities the understatement runs "
              f"{lo:.1f}x to {hi:.1f}x on the prose anchor,")
        lo2 = min(pub_log / s for _, s in live)
        hi2 = max(pub_log / s for _, s in live)
        print(f"  and {lo2:.1f}x to {hi2:.1f}x on the generated one.")
    print()
    print("  AND ONE ROW BELOW THE LIVE RANGE THE AXIS CHANGES CHARACTER.")
    print("  At a twentieth of nominal there is no frontier at all, because")
    print("  the increment's own prior floors the design above the gate")
    print("  whatever k_e is. So the requirement is not simply steeper as")
    print("  capacity falls: it is a MEASUREMENT requirement between a fifth")
    print("  and a tenth of nominal and a STRUCTURAL limit below that, and a")
    print("  single unqualified factor cannot express either.")


if __name__ == "__main__":
    a, f = part_a()
    r = part_b(a, f)
    part_c(r, a)
