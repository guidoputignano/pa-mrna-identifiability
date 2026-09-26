"""
Rounding defeats the category-1 match, and it hid a whole class.

THE BUG, found by checking the sweep's own headline claim. `07_design.tex:113`
prints `0.457 ($0.412$ to $0.510$)` and was classified category 3, prose
asserting a computed result no generator owns. It is not. The committed table
`results/tables/tab_D_designs.tex` carries

    six-monthly clinic visits & 4 & ... & 0.4570 \\scriptsize[0.4121, 0.5097]

so the value IS generated and IS guarded, and the prose merely ROUNDS it:
0.4570 -> 0.457, 0.4121 -> 0.412, 0.5097 -> 0.510. The classifier matched
literal strings, so every rounded restatement fell through into category 3.

That matters in both directions. It inflates category 3 with claims that are
actually checked, and it understates category 1, which is the cheap-to-fix pile.
It also means the sweep's own "look at this one first" was wrong: that claim
reproduces, and its defect is DUPLICATION of a guarded value, not an unowned
assertion.

This file redoes the match with rounding tolerance and reports what moves.

Run:  python audit/q20_rounding.py
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CAT3 = ROOT / "audit" / "q17_category3.txt"


def generated():
    """Every value a generator owns, as floats with their source."""
    out = []
    for p in (RESULTS / "numbers").glob("*.tex"):
        for m in re.finditer(r"\\endcsname\{([^}]*)\}", p.read_text()):
            for t in re.findall(r"\d+(?:\.\d+)?", m.group(1)):
                out.append((float(t), f"numbers/{p.name}"))
    for p in (RESULTS / "tables").glob("*.tex"):
        for t in re.findall(r"\d+(?:\.\d+)?", p.read_text()):
            out.append((float(t), f"tables/{p.name}"))
    return out


#: A rounding match is only evidence when the literal is precise enough that a
#: coincidence is unlikely. `0.4570 -> 0.457` is a restatement; `0.45 -> 0.5` is
#: two different round numbers meeting, and `0.1696 -> 0` is nothing at all. The
#: threshold is THREE significant digits, below which a match is not claimed.
MIN_SIG = 3


def sigfigs(lit):
    t = lit.lstrip("0.").replace(".", "").rstrip()
    return len(t.lstrip("0"))


def rounds_to(g, lit):
    """Does the generated value `g` round to the literal `lit` as printed?"""
    if sigfigs(lit) < MIN_SIG:
        return False
    if "." in lit:
        d = len(lit.split(".")[1])
        return f"{g:.{d}f}" == lit
    try:
        return f"{g:.0f}" == lit and abs(g - float(lit)) < 0.5
    except ValueError:
        return False


def main():
    gen = generated()
    print("=" * 96)
    print("  ROUNDED RESTATEMENTS OF GENERATED VALUES, MISSED BY AN EXACT MATCH")
    print("=" * 96)
    rows, moved = [], []
    for line in CAT3.read_text().splitlines():
        loc, val, why, ctx = line.split("\t", 3)
        v = val.replace(",", "")
        if not re.match(r"^\d+(\.\d+)?$", v):
            rows.append((loc, val, None))
            continue
        #  `0.4570` and `0.457` are EQUAL as floats but differ as strings, so
        #  they fell through the exact-string test in q17 and would have fallen
        #  through here too if float equality were excluded. It is not.
        hit = next(((g, src) for g, src in gen if rounds_to(g, v)), None)
        rows.append((loc, val, hit))
        if hit:
            moved.append((loc, val, hit))
    loose = sum(1 for line in CAT3.read_text().splitlines()
                if sigfigs(line.split("\t")[1].replace(",", "")) < MIN_SIG)
    print(f"\n  {len(rows)} category-3 literals examined, "
          f"{len(moved)} are a ROUNDED form of a value a generator already owns")
    print(f"  at {MIN_SIG} or more significant digits. Those belong in "
          f"category 1, not category 3.")
    print(f"  {loose} literals carry fewer than {MIN_SIG} significant digits, "
          f"where a rounding match\n  would be a coincidence and none is "
          f"claimed.\n")
    by = Counter(src for _, _, (g, src) in moved)
    print("  Where the owning value lives:")
    for k, n in by.most_common(12):
        print(f"    {n:5d}  results/{k}")
    print("\n  The moved literals, by file:")
    grp = defaultdict(list)
    for loc, val, (g, src) in moved:
        grp[loc.rsplit(":", 1)[0]].append((loc, val, g, src))
    for f in sorted(grp):
        print(f"\n    {f}")
        for loc, val, g, src in sorted(grp[f]):
            print(f"      {loc.rsplit(':',1)[1]:>5s}  prose {val:<9s} "
                  f"generated {g:<10g} in {src}")
    print(f"""
  THE HEADLINE THE SWEEP GOT WRONG. `07_design.tex:113` was named as the claim
  to look at first, "a median and an interval that no generator printed". All
  three of its numbers are in `tables/tab_D_designs.tex` at four decimals. It
  reproduces, it is guarded, and its defect is that a guarded value is restated
  by hand where a macro would do. That is category 1: the cheapest kind to fix
  and the least alarming kind to find.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
