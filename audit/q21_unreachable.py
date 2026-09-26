"""
Scripts declaring nothing of theirs reaches the paper, that own a number in it.

WHY. R80 found `audit/t3_diet_modulation.py` outside `make figures`, so its
results reach the paper only as prose literals, and reported it as one instance.
It is not one instance. Most scripts in `audit/` DECLARE in their docstrings that
nothing of theirs reaches the paper, and that declaration is the only thing
separating "exploratory scoping work, correctly excluded from the build" from "a
generator no entry point calls, whose numbers are in the paper by hand".

THIS IS R23's FAILURE MODE ON A SECOND SURFACE. R23 was a generator no entry
point called, which froze and printed stale numbers into a TABLE. The same shape
appears here in PROSE: a generator the build never runs, whose numbers the paper
states as literals, with a docstring asserting the opposite.

METHOD, and why it is not digit matching. An earlier version of this file matched
every value a script printed against every literal in `paper/`, and drowned:
round constants (100, 365, 14.5, 600, 0.20) are printed by every script here and
a digit collision on them means nothing. Ownership is SEMANTIC. So each confirmed
case below names the QUANTITY, the paper's sentence, and the script's own printed
line, and the check is that the two are the same quantity to the precision the
paper prints. Candidates were found by reading each declaring script against the
sections its own docstring names, and each was then put to three independent
agents instructed to refute it.

Nothing is changed and `paper/` is not touched.

Run:  python audit/q21_unreachable.py
"""
import ast
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEN = [Path("/tmp/claude-0/-home-user-system-identification-pa/"
            "b7530985-7665-5e5d-a8f1-f1a9534f728b/gen"),
       Path("/tmp/claude-0/-home-user-system-identification-pa/"
            "b7530985-7665-5e5d-a8f1-f1a9534f728b/gen2")]

DECLARES = re.compile(
    r"(nothing here reaches the paper"
    r"|nothing (in this|of this)[^.]{0,40}reaches the paper"
    r"|scoping only|measurement only"
    r"|nothing here writes a results/numbers file)", re.I)

#: The confirmed ownerships, each stated as (paper site, what the paper says,
#: the script's own printed evidence). These are quantity matches read out of the
#: scripts and their captured stdout, not digit collisions.
CONFIRMED = [
    ("t4_bias_vs_scatter", [
        ("08_specification.tex:255", "per-visit scatter alone gives 0.283",
         "at 20% in the SCATTER alone   0.2829   (1.67x)"),
        ("08_specification.tex:255", "the systematic offset alone gives 0.172",
         "at 20% in the BIAS alone      0.1718   (1.01x)"),
        ("08_specification.tex:255", "against a 0.170 baseline",
         "bias to 5%:  0.1696"),
        ("08_specification.tex:255", "98% of the degradation is scatter",
         "the scatter channel accounts for 0.1133  (98%)"),
        ("08_specification.tex:256", "the target is met at every offset from "
                                     "1% to 500%",
         "the gate holds at EVERY bias from 1% to 500% (0.1691, 0.1733)"),
        ("06_excitation.tex:286", "sd runs 0.159 at a 1% recording error to "
                                  "0.288 at 20%",
         "1%  0.1589  ... 20%  0.2863"),
        ("06_excitation.tex:287", "98% of that is the per-visit scatter and "
                                  "2% the systematic offset",
         "the same 98% split, restated in a second section"),
        ("06_excitation.tex:290", "the target is met at every bias from 1% "
                                  "to 500%", "the same 1%-to-500% sweep"),
    ]),
    ("t8_load_by_capacity", [
        ("11_scope.tex:207", "worth a factor of 1.81 at nominal capacity",
         "100% of nom   1.808"),
        ("11_scope.tex:207", "and 1.22 to 1.46 across the disease range",
         "20% of nom 1.215 ... 2% of nom 1.464"),
        ("08_specification.tex:265", "the offset's share runs 1.8% at nominal, "
                                     "9.9% at half, 26.2% at a fifth",
         "offset share: 100% 1.8% | 50% 9.9% | 20% 26.2%"),
        ("11_scope.tex:210", "a perfect record in both channels leaves 2.659",
         "theta = 150 (5% of nominal): 0.1% 2.6595   [also printed by t9]"),
    ]),
    ("t9_measured_load", [
        ("11_scope.tex:213", "thirteen measured-load visits match twenty-five "
                             "recorded ones", "100% of nom  ...  13"),
        ("11_scope.tex:213", "and six reach the target",
         "100% of nom  ...  measured n to gate  6"),
        ("11_scope.tex:215", "an exactly known load leaves 2.659 at 5% of "
                             "nominal", "theta = 150: to 0.1%, n=25:  2.6595"),
        ("11_scope.tex:216", "the advantage shrinks ... from thirteen visits "
                             "to eighteen", "2% of nom  7.4406  18"),
        ("08_specification.tex:270", "a perfect one in both channels still "
                                     "leaves 2.659", "the same 2.6595"),
        ("06_excitation.tex:309", "measuring the load instead of recording it "
                                  "is worth under one per cent",
         "100% of nom: recorded 5% 0.1699 vs b gone, sig 5% 0.1693 = 0.35%"),
    ]),
    ("t3_diet_modulation", [
        ("11_scope.tex:204", "a prescribed +/-15% swing is worth 7%, an "
                             "implausible +/-50% 32%, leaving 2.32",
         "R80's finding, the first instance"),
        ("08_specification.tex:372", "the same three values restated",
         "same generator"),
        ("08_specification.tex:375", "prescribing it is worth between nothing "
                                     "and two per cent",
         "designed 0.1696 vs spontaneous 0.1698, ratio 1.001 at nominal; the "
         "upper endpoint sits at a capacity not located here"),
    ]),
]


def docstring(p):
    try:
        return ast.get_docstring(ast.parse(p.read_text())) or ""
    except SyntaxError:
        return ""


def in_make_figures(name):
    mk = (ROOT / "Makefile").read_text()
    body = mk.split("figures:", 1)[1].split("\nrefs:", 1)[0]
    return f"{name}.py" in body


def writes_numbers(name):
    """Does the script write a results/numbers file, which is what the guard's
    reachability check actually tracks?"""
    src = (ROOT / "audit" / f"{name}.py").read_text()
    return bool(re.search(r"write_numbers\(", src))


def ran(name):
    for d in GEN:
        if (d / f"{name}.log").exists():
            return True
    return False


def main():
    print("=" * 100)
    print("  SCRIPTS THAT SAY THEY DO NOT REACH THE PAPER, AND THE ONES THAT DO")
    print("=" * 100)
    scripts = sorted((ROOT / "audit").glob("*.py"))
    decl = [(p.stem, DECLARES.search(docstring(p))) for p in scripts]
    decl = [(n, m) for n, m in decl if m]
    owners = {n for n, _ in CONFIRMED}
    print(f"\n  {len(scripts)} scripts in audit/. {len(decl)} declare that "
          f"nothing of theirs reaches the paper.\n")
    print(f"  {'script':26s} {'in make figures':>16s} {'writes numbers':>15s} "
          f"{'owns a paper claim':>19s}")
    for n, _ in decl:
        print(f"  {n:26s} {str(in_make_figures(n)):>16s} "
              f"{str(writes_numbers(n)):>15s} "
              f"{('YES, ' + str(len(dict(CONFIRMED)[n])) if n in owners else 'no'):>19s}")
    print("\n" + "=" * 100)
    print("  THE LIST")
    print("=" * 100)
    print("""
  Each row: what the paper asserts, and the line the script prints. The paper
  states these as bare literals; none is a \\num{} macro.""")
    for name, claims in CONFIRMED:
        print(f"\n  audit/{name}.py")
        print(f"      in `make figures`: {in_make_figures(name)}    "
              f"writes a numbers file: {writes_numbers(name)}")
        d = docstring(ROOT / "audit" / f"{name}.py")
        m = DECLARES.search(d)
        print(f"      declares: \"...{' '.join(d[m.start():m.end()+58].split())}...\"")
        for site, says, prints in claims:
            print(f"        {site}")
            print(f"            paper:  {says}")
            print(f"            script: {prints}")
    print("\n" + "=" * 100)
    print("  WHY THE GUARD CANNOT SEE ANY OF IT")
    print("=" * 100)
    print("""  `tools/check_artefacts.py` reports "14 numbers-file generator(s) checked for
  reachability from `make figures`". That check tracks generators that WRITE a
  numbers file. None of the four above writes one, so none is in the 14, none is
  reachable, and none is checked. The guard is not failing; it was never pointed
  at this class.

  So there are two independent conditions and the paper needs both:
      1. a generator the build runs, which `make figures` decides
      2. a paper number that cites it, which `\\num{}` decides
  These four scripts satisfy neither, and their docstrings assert that the second
  is deliberate. For t4, t8 and t9 that assertion is false.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
