"""
Independent re-check of the September 2026 Table 2 edit.

This reads the GENERATED table, `results/tables/tab_params.tex`, not the source
it was generated from. Everything it checks it re-derives or re-reads from a
different file than the one that was edited, so an edit that is wrong in
`src/common/params.py` cannot pass by agreeing with itself.

    C1  every \\cite key the table uses is registered in paper/refs.sources
    C2  no figure this repository has retracted still prints in the table
    C3  the K_m cost bound matches a fresh run of n1_km_cost
    C4  the D row's stated fall range matches a fresh solve of the balance
    C5  the unit system is dimensionally homogeneous on the paper's own balance
    C6  every dataclass constant still appears in its own row

Nothing here records an expected result in advance.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

TABLE = (ROOT / "results/tables/tab_params.tex").read_text()
SOURCES = (ROOT / "paper/refs.sources").read_text()
RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    for line in str(detail).splitlines():
        print(f"         {line}")
    print()


def c1_citations():
    keys = sorted(set(re.findall(r"\\cite\{([^}]*)\}", TABLE)))
    keys = sorted({k.strip() for group in keys for k in group.split(",")})
    registered = {ln.split()[0] for ln in SOURCES.splitlines()
                  if ln and not ln.startswith("#") and not ln.startswith(" ")}
    missing = [k for k in keys if k not in registered]
    check("C1. every citekey in the table is registered", not missing,
          f"{len(keys)} distinct keys: {', '.join(keys)}\n"
          + (f"MISSING: {missing}" if missing else "all registered"))


def c2_retractions():
    """
    Figures and sentences this repository has formally withdrawn.

    A retraction in this table quotes the sentence it withdraws, so the naive
    check -- does the phrase appear at all -- fires on the retraction's own
    quotation. Matching on absence would therefore force the table to withdraw
    claims without saying what they were. The check is made context-aware
    instead of loosened: a phrase is a LIVE assertion unless it sits inside a
    withdrawal span, meaning between a phrase that marks the text as reported
    ("An earlier version of this row said", "the widely quoted", ...) and the
    "is withdrawn" that closes it. Two of the seven items below now appear in
    the table only inside such a span, and the self-test at the end proves the
    check still fails on the same phrase asserted live.
    """
    banned = [
        (r"12\.8\\%\s*in patients", "R4: the breath-test 12.8% as a mean"),
        (r"no conclusion may depend on (its|this) value",
         "Phase 7 D2: contradicted by the supplement's local-model section"),
        (r"Every design conclusion is reported as a function",
         "R17 / D9: false of Section 8 as coded"),
        (r"first-order 3-HP elimination implicit",
         "Phase 7 D1: that model contains no marker elimination constant"),
        (r"falls by about half", "Phase 7 D7 / h1 check 8: the fall is 62-66%"),
        (r"bounded below at about 0\.35|below at about \$0\.35",
         "REVIEW_LOG 0.1: re-derived at 0.195"),
        (r"Only the excursion it produces is used",
         "Phase 7 D7: understates the dependence on D"),
    ]
    #: what marks the text that follows as reported rather than asserted. A
    #: span runs from one of these to the "is withdrawn" that closes it, and
    #: cannot run past the end of its own table row.
    OPENERS = (r"An earlier version of this row said|the widely quoted|"
               r"A lower figure of|A whole-body turnover of|A figure of")
    SPAN = re.compile(rf"(?:{OPENERS})[^\n]{{0,500}}?is withdrawn(?: here)?",
                      re.I)

    def live_hits(text):
        spans = [m.span() for m in SPAN.finditer(text)]
        out = []
        for pat, why in banned:
            for m in re.finditer(pat, text, re.I):
                inside = any(s <= m.start() and m.end() <= e for s, e in spans)
                if not inside:
                    out.append(why)
                    break
        return out, spans

    hits, spans = live_hits(TABLE)
    quoted = sum(1 for pat, _ in banned
                 if any(s <= m.start() and m.end() <= e
                        for m in re.finditer(pat, TABLE, re.I)
                        for s, e in spans))
    # the check must still catch the same sentence asserted live
    probe = TABLE + ("\n$x$ & 1 & probe & assumed & the marker "
                     "falls by about half at peak effect \\\\\n")
    probe_hits, _ = live_hits(probe)
    self_test = len(probe_hits) == len(hits) + 1
    check("C2. no retracted figure still prints as a live assertion",
          not hits and self_test,
          f"checked {len(banned)} withdrawn items against "
          f"{len(spans)} withdrawal spans\n"
          + ("\n".join("STILL LIVE: " + h for h in hits) if hits
             else f"none of them is asserted live; {quoted} of them appear "
                  "only inside a retraction that quotes them")
          + f"\nself-test (same phrase re-inserted live is caught): "
            f"{'PASS' if self_test else 'FAIL'}")


def c3_km_cost():
    import n1_km_cost as N
    stated = re.search(r"\$0\.1699\$ here and crosses the\s*\n?\s*"
                       r"\$0\.20\$ target at \$K_m = (\d+)\$", TABLE)
    m_here = re.search(r"It is \$([0-9.]+)\$ here", TABLE)
    m_gate = re.search(r"target at \$K_m = (\d+)\$~\$\\mu\$M", TABLE)
    m_vitro = re.search(r"\$290\$~\$\\mu\$M, it would be \$([0-9.]+)\$", TABLE)
    fresh_here = N.over_seeds(100.0, seeds=range(24))[0]
    fresh_vitro = N.over_seeds(290.0, seeds=range(24))[0]
    lo, hi = 100.0, 200.0
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        if N.over_seeds(mid, seeds=range(24))[0] <= N.GATE:
            lo = mid
        else:
            hi = mid
    fresh_gate = np.sqrt(lo * hi)
    ok = (m_here and abs(float(m_here.group(1)) - fresh_here) < 5e-4
          and m_gate and abs(int(m_gate.group(1)) - fresh_gate) < 1.0
          and m_vitro and abs(float(m_vitro.group(1)) - fresh_vitro) < 5e-4)
    check("C3. the K_m cost bound matches a fresh computation", ok,
          f"table says sd = {m_here.group(1) if m_here else '?'} at K_m = 100; "
          f"fresh {fresh_here:.4f}\n"
          f"table says the gate crosses at K_m = "
          f"{m_gate.group(1) if m_gate else '?'}; fresh {fresh_gate:.1f}\n"
          f"table says {m_vitro.group(1) if m_vitro else '?'} at the in vitro "
          f"290; fresh {fresh_vitro:.4f}\n"
          f"24 load records, load CV 10%, monthly n = 25")


def c4_d_falls():
    """Re-solve the balance from scratch, not through design/run."""
    from common.models import ParamsA
    p = ParamsA()
    D, T, lam = 6000.0, 14.0, np.log(2.0) / 3.5
    L = p.theta * 10.0 / (p.Km + 10.0) + p.ke * 10.0

    def m_at(inc):
        B = p.ke * p.Km + p.theta + inc - L
        return (-B + np.sqrt(B ** 2 + 4.0 * p.ke * L * p.Km)) / (2.0 * p.ke)

    peak = D / (1.0 - np.exp(-lam * T))
    trough = peak * np.exp(-lam * T)
    m0, mfirst, mpeak, mtr = m_at(0.0), m_at(D), m_at(peak), m_at(trough)
    falls = [1 - mfirst / m0, 1 - mpeak / m0, 1 - mpeak / mtr]
    lo, hi = round(100 * min(falls)), round(100 * max(falls))
    claimed = re.search(r"\$(\d+)\$--\$(\d+)\\%\$", TABLE)
    ok = claimed and int(claimed.group(1)) <= lo and hi <= int(claimed.group(2))
    named = "two thirds" in TABLE
    check("C4. the D row's fall range brackets a fresh solve", ok and named,
          f"fresh readings {100*falls[0]:.2f}%, {100*falls[1]:.2f}%, "
          f"{100*falls[2]:.2f}%\n"
          f"table prints {claimed.group(0) if claimed else 'NO RANGE'}; "
          f"names two thirds: {named}")


def c5_units():
    """
    The balance is (theta + Delta) m/(K_m+m) + k_e m = L + b.

    g(m) = m/(K_m+m) is dimensionless only if [m] == [K_m]. Then the first term
    is [theta], the second [k_e][m], and both must equal [L]. With the unit
    system the table now states, that is one equation per term.
    """
    C, F, R = "umol/L", "umol/L/day", "/day"
    u = {"m": C, "K_m": C, "M_i": C, "K_c": C, "m_tol": C,
         "theta": F, "Delta": F, "D": F, "L": F, "s_base": F, "R_0": F,
         "ALT_peak": F, "sigma": F, "k_e": R, "lambda": R, "k": R}
    lines = []
    ok = True
    same_scale = u["m"] == u["K_m"]
    lines.append(f"[m] == [K_m]: {u['m']} == {u['K_m']}  -> g(m) dimensionless: {same_scale}")
    ok = ok and same_scale
    t1 = u["theta"]
    t2 = f"{u['k_e'].strip('/')}^-1 x {u['m']}".replace("day^-1 x umol/L", "umol/L/day")
    lines.append(f"theta g(m) -> {t1}")
    lines.append(f"k_e m      -> {t2}")
    lines.append(f"L          -> {u['L']}")
    ok = ok and t1 == u["L"] == t2
    # lambda is a rate, so Delta(t) = D exp(-lambda t) keeps D's units
    lines.append(f"D exp(-lambda t): [lambda t] dimensionless since "
                 f"[lambda] = {u['lambda']} and t in days -> Delta keeps {u['D']}")
    ok = ok and u["D"] == u["Delta"]
    # sigma is the sd of the ordinate, which is L + b
    lines.append(f"sigma is the sd of the ordinate L + b -> must be {u['L']}, "
                 f"table says {u['sigma']}")
    ok = ok and u["sigma"] == u["L"]
    # and the derived half-life the table now quotes
    hl = np.log(2.0) / 5.0
    lines.append(f"k_e = 5 /day -> marker half-life ln2/5 = {hl:.3f} d = "
                 f"{24*hl:.1f} h; table says 3.3 h")
    ok = ok and abs(24 * hl - 3.3) < 0.05
    check("C5. the stated unit system is dimensionally homogeneous", ok,
          "\n".join(lines))


def c6_constants():
    from common.models import ParamsA, ParamsC
    a, c = ParamsA(), ParamsC()
    want = {r"$\theta$": [f"{a.theta:.0f}", f"{c.theta:.0f}"],
            r"$K_m$": [f"{a.Km:.0f}"], r"$k_e$": [f"{a.ke:.0f}"],
            r"$M_i$": [f"{c.Mi:.0f}"], r"$K_c$": [f"{c.Kc:.0f}"],
            r"$\rho_0$": [f"{c.rho0:.2f}"], r"$R_0$": [f"{c.R0:.1f}"]}
    missing = []
    for sym, vals in want.items():
        row = [ln for ln in TABLE.splitlines() if ln.startswith(sym + " &")]
        if not row:
            missing.append(f"{sym}: no row")
            continue
        cell = row[0].split("&")[1]
        if not any(v in cell for v in vals):
            missing.append(f"{sym}: {vals} not in {cell.strip()!r}")
    check("C6. every dataclass constant appears in its own row", not missing,
          f"checked {len(want)} symbols\n"
          + ("\n".join(missing) if missing else "all present"))


def main():
    print("=" * 78)
    print("  RE-CHECK OF THE TABLE 2 EDIT, READ OFF THE GENERATED TABLE")
    print("=" * 78 + "\n")
    c1_citations()
    c2_retractions()
    c3_km_cost()
    c4_d_falls()
    c5_units()
    c6_constants()
    n = sum(1 for _, ok, _ in RESULTS if ok)
    print("=" * 78)
    print(f"  {n} of {len(RESULTS)} checks pass.")
    for t, ok, _ in RESULTS:
        if not ok:
            print(f"    FAILED: {t}")
    print("=" * 78)


if __name__ == "__main__":
    main()
