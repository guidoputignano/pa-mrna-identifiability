"""
Does each category-3 claim still reproduce?

THE METHOD, and its limit is stated first because it decides how the list reads.
A category-3 literal is by definition NOT written by a generator, so there is
nothing to diff it against. What CAN be done cheaply is to run the generators
that already serve the paragraph and search their printed output: these scripts
print far more than they write as `\\num{}` macros, so a literal that appears in
its own generator's stdout is a value the code still produces and that merely
was never wired to a macro. That is a real reproduction verdict.

    REPRODUCES        the value appears in the output of a generator serving
                      that section, to the precision the paper prints it
    NOT IN OUTPUT     it does not, which is not proof it is wrong: the
                      generator may simply never print that quantity. It IS
                      proof that nothing in the build checks it.
    NO GENERATOR RUN  no generator serves that paragraph at all

Route 6's `eighteen` is the case that motivates the whole exercise: it is NOT
IN OUTPUT, and q16 part E then showed by direct recomputation that it is wrong.
So NOT IN OUTPUT is the queue to work through, not a verdict of error.

Run:  python audit/q17_reproduce.py
"""
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEN = Path("/tmp/claude-0/-home-user-system-identification-pa/"
           "b7530985-7665-5e5d-a8f1-f1a9534f728b/gen")
CLAIMS = ROOT / "audit" / "q17_claims.tsv"
WORDS = {"one": "1", "two": "2", "three": "3", "four": "4", "five": "5",
         "six": "6", "seven": "7", "eight": "8", "nine": "9", "ten": "10",
         "eleven": "11", "twelve": "12", "thirteen": "13", "fourteen": "14",
         "fifteen": "15", "sixteen": "16", "seventeen": "17",
         "eighteen": "18", "nineteen": "19", "twenty": "20", "thirty": "30",
         "forty": "40", "fifty": "50", "sixty": "60", "eighty": "80",
         "hundred": "100", "twenty-five": "25", "twenty-eight": "28",
         "twenty-four": "24"}
#: which generator log serves which numbers file
SERVES = {"scope": ["n4_theta_scope"], "transient": ["t1_transient_worth"],
          "kmcost": ["n1_km_cost"], "profile": ["n6_profile_vs_crb"],
          "kesweep": ["t5_band_vs_ke"], "lattice": ["t11_phase_lattice"],
          "structure": ["t7_structure_transient"],
          "boundary": ["h1_admissibility"], "margin": ["h1_admissibility"],
          "design": ["src_design_run"], "local": ["src_model_A_run"],
          "fibre": ["src_model_A_run"], "reduction": ["src_model_C_run"],
          "therapy": ["audit_t13_ther_forward.py"],
          "cohort": ["audit_t6_cohort_with_kappa"],
          "phaseclock": ["audit_g2_visits_vs_sd"]}


def logs():
    return {p.stem: p.read_text() for p in GEN.glob("*.log")}


def seen(val, text):
    """Is this value printed, at the precision the paper uses?"""
    v = WORDS.get(val.lower(), val).replace(",", "")
    if "." in v:
        return re.search(r"(?<![\d.])" + re.escape(v) + r"(?![\d])",
                         text) is not None
    return re.search(r"(?<![\d.])" + re.escape(v) + r"(?![\d.])",
                     text) is not None


def main():
    L = logs()
    if not L:
        print("  no generator logs captured; run the generators first")
        return 1
    print("=" * 100)
    print("  DOES EACH CATEGORY-3 CLAIM REPRODUCE")
    print("=" * 100)
    print(f"\n  Generators captured: {', '.join(sorted(L))}")
    rows, verdicts = [], Counter()
    for line in CLAIMS.read_text().splitlines()[1:]:
        f, n, vals, owners, ctx = line.split("\t", 4)
        names = []
        for o in owners.split("+"):
            names += SERVES.get(o, [])
        names = [x for x in names if x in L]
        if not names:
            #  no generator serves this paragraph. Fall back to asking whether
            #  the value is printed ANYWHERE in the generated corpus, which is
            #  weaker (it may be a coincidence) but separates a literal the
            #  code still produces from one it produces nowhere at all.
            blob = "\n".join(L.values())
            miss = [x for x in vals.split(",") if not seen(x, blob)]
            if miss:
                v, note = "NO OWNER, UNPRINTED", f"missing {','.join(miss)}"
            else:
                v, note = "NO OWNER, printed somewhere", "corpus-wide match"
        else:
            blob = "\n".join(L[x] for x in names)
            hit = [x for x in vals.split(",") if seen(x, blob)]
            miss = [x for x in vals.split(",") if not seen(x, blob)]
            if not miss:
                v, note = "REPRODUCES", "+".join(names)
            elif hit:
                v, note = "PARTIAL", f"missing {','.join(miss)}"
            else:
                v, note = "NOT IN OUTPUT", "+".join(names)
        verdicts[v] += 1
        rows.append((f, n, vals, owners, v, note, ctx))
    print()
    for k, c in verdicts.most_common():
        print(f"    {c:5d}  {k}")
    out = ROOT / "audit" / "q17_verdicts.tsv"
    with out.open("w") as fh:
        fh.write("file\tline\tvalues\tgenerator\tverdict\tnote\tcontext\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"\n  Written to {out.relative_to(ROOT)}")
    print("\n  The claims a generator RAN for and did not print:")
    for f, n, vals, ow, v, note, ctx in rows:
        if v in ("NOT IN OUTPUT", "PARTIAL"):
            print(f"    {f}:{n}  [{vals}]  {note}")
            print(f"        {ctx[:150]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
