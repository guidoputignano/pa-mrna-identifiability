"""
Convert the verified unguarded literals in paper/ to \\num{} macros.

THE RULE, and it is enforced rather than intended: a replacement happens ONLY
if the macro renders EXACTLY the literal it replaces. No value changes in this
conversion. Where the generated value and the published literal differ, both are
left alone and the pair is reported, because adopting the computed value would
silently change a published number and adopting the published one would falsify
its generator.

Spelled-out numerals are NOT converted. `\\num{}` renders a digit, so turning
"thirteen" into a macro rewords the sentence, and the paper is frozen against
prose edits. Those sites are listed for the author instead. That they cannot be
wired without rewording is the same property that made them invisible to a
digit-based scan in the first place (REVIEW_LOG.md R77).

Run:  python audit/q22_wire_macros.py            # report only
      python audit/q22_wire_macros.py --apply    # perform the conversions
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUMDIR = ROOT / "results" / "numbers"

#: (section file, the exact LaTeX to replace, the macro key)
#: The literal inside the LaTeX must equal what the macro renders.
SITES = [
    ("08_specification.tex", r"$0.212$", "tier1_sd_at_pub"),
    ("08_specification.tex", r"$118$", "tier1_visits_at_pub"),
    ("08_specification.tex", r"$1.21$", "tier1_factor"),
    ("08_specification.tex", r"$17\%$", "tier1_cv_cut_pct"),
    ("08_specification.tex", r"$0.151$", "tier1_floor_10"),
    ("08_specification.tex", r"$0.160$", "tier1_floor_req"),
    ("08_specification.tex", r"$0.174$", "tier1_floor_14"),
    ("08_specification.tex", r"$0.207$", "tier1_floor_20"),
    ("08_specification.tex", r"$0.248$", "tier1_floor_30"),
    ("08_specification.tex", r"$0.283$", "bias_scatter20"),
    ("08_specification.tex", r"$0.172$", "bias_alone20"),
    ("08_specification.tex", r"$0.170$", "bias_base"),
    ("08_specification.tex", r"$98\%$", "bias_scatter_share"),
    ("08_specification.tex", r"$9.9\%$", "loadcap_share_half"),
    ("08_specification.tex", r"$26.2\%$", "loadcap_share_fifth"),
    ("08_specification.tex", r"$2.659$", "loadcap_perfect_5pct"),
    ("11_scope.tex", r"$1.81$", "loadcap_ratio_nominal"),
    ("11_scope.tex", r"$1.22$", "loadcap_ratio_dis_lo"),
    ("11_scope.tex", r"$1.46$", "loadcap_ratio_dis_hi"),
    ("11_scope.tex", r"$2.659$", "loadcap_perfect_5pct"),
    ("06_excitation.tex", r"$0.159$", "bias_tight"),
    ("06_excitation.tex", r"$98\%$", "bias_scatter_share"),
]

#: sites deliberately NOT converted, with the reason
HELD = [
    ("08_specification.tex", "0.288 / 06_excitation.tex:286",
     "MISMATCH: the generator gives 0.286, the paper says 0.288"),
    ("08_specification.tex:372, 11_scope.tex:204", "7% and 32%",
     "MISMATCH: the generator gives 6.3% and 24.1%"),
    ("11_scope.tex:213", "thirteen, six", "spelled numeral, needs a reword"),
    ("11_scope.tex:216", "thirteen, eighteen", "spelled numeral"),
    ("08_specification.tex:126", "nearly ten years", "spelled numeral"),
    ("08_specification.tex:127", "eleven", "spelled numeral"),
    ("06_excitation.tex:309", "under one per cent", "spelled numeral"),
    ("08_specification.tex:151", "18.7%", "already a generated macro"),
]


def rendered():
    out = {}
    for p in NUMDIR.glob("*.tex"):
        for m in re.finditer(
                r"csname paper@num@([A-Za-z0-9_]+)\\endcsname\{([^}]*)\}",
                p.read_text()):
            out[m.group(1)] = m.group(2)
    return out


def main():
    apply = "--apply" in sys.argv
    vals = rendered()
    print("=" * 96)
    print("  WIRING THE VERIFIED LITERALS TO GENERATED MACROS")
    print("=" * 96)
    ok, bad = [], []
    for fn, lit, key in SITES:
        want = re.sub(r"[^\d.]", "", lit)
        got = vals.get(key)
        if got is None:
            bad.append((fn, lit, key, "no such macro"))
        elif got != want:
            bad.append((fn, lit, key, f"macro renders {got}, literal is {want}"))
        else:
            ok.append((fn, lit, key))
    print(f"\n  {len(ok)} of {len(SITES)} sites verified: the macro renders the "
          f"literal exactly.")
    if bad:
        print(f"  {len(bad)} REFUSED:")
        for fn, lit, key, why in bad:
            print(f"    {fn} {lit} -> {key}: {why}")
    if apply:
        n = 0
        for fn, lit, key in ok:
            p = ROOT / "paper" / "sections" / fn
            s = p.read_text()
            if lit not in s:
                print(f"    ! {fn}: {lit} no longer present, skipped")
                continue
            s = s.replace(lit, "$\\num{" + key + "}$" if not lit.endswith(r"\%$")
                          else "$\\num{" + key + r"}\%$", 1)
            p.write_text(s)
            n += 1
        print(f"\n  applied {n} conversions")
    else:
        print("\n  report only; pass --apply to convert")
    print("\n" + "=" * 96)
    print("  HELD BACK, AND WHY")
    print("=" * 96)
    for a, b, why in HELD:
        print(f"    {a:44s} {b:22s} {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
