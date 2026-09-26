"""
Classify every unguarded numeric literal in paper/ into the author's four
categories, and isolate category 3.

  1  GENERATED-ELSEWHERE  the value also exists in results/numbers or a
     committed table, so a generator owns the quantity and the prose retypes it
  2  STRUCTURAL           equation constant, drawing coordinate, count, year
  3  ASSERTING-COMPUTED   prose-only and stating a computed outcome
  4  HARMLESS             prose-only, asserting nothing computed

SPELLED-OUT NUMERALS ARE INCLUDED, because route 6's `thirteen` and `eighteen`
are words and no scan over digits would ever have seen them. That is the
sharpest lesson of R75 for this sweep: the hole is wider than a regex.

WHAT CATEGORY 3 REQUIRES, and the thresholds are deliberately strict because an
over-inclusive list is not usable and would bury the members that matter:

  * not inside a tikzpicture, which is drawing coordinates
  * not a hyphenated compound (`eight-fold`, `six-monthly`, `one-dimensional`),
    which is an adjective and not a count
  * a RESULT VERB in the sentence: buys, beats, reaches, leaves, moves, costs,
    is worth, achieves, gives, returns, falls to
  * and a numeric SHAPE that a computation produces: a decimal carrying two or
    more places, a percentage, an explicit factor or multiple, or an integer
    attached to a countable design quantity (visits, draws, records, patients,
    seeds, replicates)

Run:  python audit/q17_classify.py
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import q17_numeric_provenance as SCAN

ROOT, PAPER = SCAN.ROOT, SCAN.PAPER
RESULTS = ROOT / "results"

WORDS = ("one two three four five six seven eight nine ten eleven twelve "
         "thirteen fourteen fifteen sixteen seventeen eighteen nineteen "
         "twenty thirty forty fifty sixty seventy eighty ninety hundred "
         "twenty-five twenty-eight twenty-four").split()
WORD_RE = re.compile(r"(?<![-\w])(" + "|".join(sorted(WORDS, key=len,
                                                      reverse=True))
                     + r")(?![-\w])", re.I)
#: a verb that turns a number into an assertion about an outcome
RESULT = re.compile(
    r"\b(buys?|bought|beats?|reach(es|ed)?|leaves?|left|moves?|costs?|"
    r"worth|achieves?|gives?|returns?|falls?|rises?|improves?|degrades?|"
    r"yields?|delivers?|supplies|separates?|recovers?|matches?|"
    r"is\s+enough|suffices?|match|equal(s|led)?)\b", re.I)
#: the shape of a number a computation produced
DECIMAL = re.compile(r"^\d+\.\d{2,}$")
FACTOR = re.compile(r"\b(factor|times|fold|multiple|ratio|to one)\b", re.I)
PCT = re.compile(r"\\%|\bper cent\b|%", re.I)
COUNTABLE = re.compile(
    r"\b(visits?|draws?|records?|patients?|subjects?|seeds?|replicates?|"
    r"samples?|measurements?|studies|breath tests?|infusions?|"
    r"recorded ones|reach(es)? the target|reach(es)? the gate)\b", re.I)
BENIGN = re.compile(r"(Section|Table|Figure|Appendix|Proposition|Remark|"
                    r"Lemma|Theorem|Tier)~?\\?r?e?f?", re.I)
YEAR = re.compile(r"^(19|20)\d\d$")


def tikz_lines(text):
    """Line numbers inside a tikzpicture, whose numbers are coordinates."""
    bad, depth = set(), 0
    for n, line in enumerate(text.splitlines(), 1):
        if "\\begin{tikzpicture}" in line:
            depth += 1
        if depth:
            bad.add(n)
        if "\\end{tikzpicture}" in line:
            depth = max(depth - 1, 0)
    return bad


def generated_values():
    vals = set()
    for p in (RESULTS / "numbers").glob("*.tex"):
        for m in re.finditer(r"\\endcsname\{([^}]*)\}", p.read_text()):
            vals.update(re.findall(r"\d+(?:\.\d+)?", m.group(1)))
    for p in (RESULTS / "tables").glob("*.tex"):
        vals.update(re.findall(r"\d+(?:\.\d+)?", p.read_text()))
    return vals


def in_math(ctx, value):
    i = ctx.find(value)
    if i < 0:
        return False
    return (ctx[:i].count("$") - ctx[:i].count("\\$")) % 2 == 1


#: `one` and `two` are usually determiners ("one pool, two channels", "an exact
#: one"), so they enter only when a countable design quantity is beside them. A
#: genuine computed result of one or two would be typed as a digit anyway.
DETERMINER = {"one", "two"}


def computed_shape(v, ctx, spelled):
    """Does this literal have the shape of a computed result?"""
    if spelled:
        if v.lower() in DETERMINER:
            return bool(COUNTABLE.search(ctx))
        return bool(FACTOR.search(ctx) or COUNTABLE.search(ctx)
                    or PCT.search(ctx))
    bare = v.replace(",", "")
    if DECIMAL.match(bare):
        return True
    if PCT.search(ctx) or FACTOR.search(ctx):
        return True
    if COUNTABLE.search(ctx) and bare.isdigit() and int(bare) >= 3:
        return True
    return False


def classify(r, gen, tikz):
    v, ctx, spelled = r["value"], r["ctx"], r.get("spelled", False)
    if r["line"] in tikz.get(r["file"], ()):
        return 2, "tikz drawing coordinate"
    bare = v.replace(",", "")
    if not spelled:
        if YEAR.match(bare):
            return 2, "year"
        if BENIGN.search(ctx) and not RESULT.search(ctx):
            return 2, "cross-reference"
        if bare in gen and not in_math(ctx, v):
            return 1, ("value exists in results/, retyped in a claim"
                       if RESULT.search(ctx) else "value exists in results/")
        if in_math(ctx, v) and not RESULT.search(ctx):
            return 2, "equation constant"
    if RESULT.search(ctx) and computed_shape(v, ctx, spelled):
        return 3, ("spelled-out numeral asserting an outcome" if spelled
                   else "prose asserting a computed outcome")
    return 4, "prose, nothing computed"


def main():
    gen = generated_values()
    files = sorted(PAPER.glob("*.tex")) + sorted((PAPER / "sections").glob("*.tex"))
    tikz = {str(f.relative_to(ROOT)): tikz_lines(f.read_text()) for f in files}
    rows = SCAN.sweep()
    for f in files:
        for n, raw in enumerate(f.read_text().splitlines(), 1):
            for m in WORD_RE.finditer(SCAN.strip_comments(raw)):
                rows.append({"file": str(f.relative_to(ROOT)), "line": n,
                             "value": m.group(1), "ctx": raw.strip()[:220],
                             "spelled": True})
    print("=" * 100)
    print("  CLASSIFYING EVERY UNGUARDED LITERAL IN paper/")
    print("=" * 100)
    nd = sum(1 for r in rows if not r.get("spelled"))
    print(f"\n  {nd} digit literals and {len(rows) - nd} spelled-out numerals, "
          f"{len(rows)} in all.")
    buckets = defaultdict(list)
    for r in rows:
        c, why = classify(r, gen, tikz)
        r["why"] = why
        buckets[c].append(r)
    names = {1: "GENERATED-ELSEWHERE, retyped here", 2: "STRUCTURAL",
             3: "ASSERTING-COMPUTED, guard-invisible", 4: "HARMLESS"}
    print()
    for c in (1, 2, 3, 4):
        print(f"    {c}  {names[c]:42s} {len(buckets[c]):5d}")
    for c, nm in ((3, "category3"), (1, "category1")):
        out = ROOT / "audit" / f"q17_{nm}.txt"
        with out.open("w") as fh:
            for r in sorted(buckets[c], key=lambda r: (r["file"], r["line"])):
                fh.write(f"{r['file']}:{r['line']}\t{r['value']}\t{r['why']}\t"
                         f"{r['ctx']}\n")
        print(f"\n  Category {c} -> audit/q17_{nm}.txt ({len(buckets[c])} rows)")
    print("\n  Category 3 by file:")
    for f, k in Counter(r["file"] for r in buckets[3]).most_common():
        print(f"    {k:5d}  {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
