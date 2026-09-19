"""
Phase 2, item 4: physical admissibility.

A fit can return numbers that balance the data and describe no patient. This
checks the constraints that say so, states each one as a predicate the code can
be tested against, and reports which hold.

    1  positive rates                  theta, K_m, k_e, D, lambda > 0
    2  a non-negative increment        Delta(t) >= 0, and bounded by its own
                                       steady-state peak
    3  a positive, unique level        the balance's positive root exists, is
                                       unique, and m > 0 at every visit
    4  a positive total load           L + b > 0
    5  the capacity ceiling            theta + Delta_peak <= normal capacity,
                                       WITH THE REFERENCE AS A PARAMETER
    6  plausible bounds on b           what physics allows, against what the
                                       load record assumes
    7  the settling exclusion          no sample inside the post-dose window
    8  the D calibration               params.py says D was chosen so the
                                       marker falls by about two thirds at
                                       peak effect. The row was reworded in
                                       September 2026; this is now a drift
                                       check on the stated range.

On (5) the repository stated the residual capacity two ways four lines apart,
"about 10% of normal" and "under 5% of control". The two differ by a factor of
two and they are the denominator of this constraint. THE AUTHOR HAS SETTLED
THIS AT 5%, September 2026, and 5% is the reference throughout this script and
in the text. The 10% reading is still computed and printed, as a check that the
verdict does not turn on the choice: it cannot, because theta cancels out of
the used fraction. The supporting source is fibroblast enzymology rather than
in vivo hepatic capacity, so the choice is a judgement informed by it and not a
measurement of the quantity the model carries.

Footing. Local model parameters as src/common/params.py carries them, theta =
3000, K_m = 100, k_e = 5, L = load(0.0) = 322.7273, D = 6000, T = 14 d,
lambda = ln2/3.5. The fibre of P5 is used in its own coordinates. Anything
that depends on a random draw is reported over 8 seeds with the range.

Nothing here records an expected result in advance, and nothing here edits a
constant, a docstring or a paper file.

Run:  python audit/h1_admissibility.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.optimize import brentq

from common.params import _A, _C
from design.run import (P, delta, load, quasi_steady_level, D_AMP, LAMBDA,
                        DOSE_INTERVAL, M_BAR, LOAD_FRAC)
import g_design_common as G
import g3_calendar as g3

LB = load(0.0)
#: the residual-activity reference, settled by the author at 5% in
#: September 2026. See the module docstring.
REFERENCE_FRAC = 0.05
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
TROUGH = PEAK * np.exp(-LAMBDA * DOSE_INTERVAL)
RESULTS = []


def head(t):
    print("\n" + "=" * 94); print("  " + t); print("=" * 94)


def record(tag, ok, detail):
    RESULTS.append((tag, ok, detail))
    print(f"    [{'PASS' if ok else 'FAIL' if ok is False else 'n/a ':>4}] "
          f"{tag:<44} {detail}")


# ------------------------------------------------------------------ 1, 2, 3
def rates_and_increment():
    head("1-3. POSITIVE RATES, A NON-NEGATIVE INCREMENT, A POSITIVE LEVEL")
    for nm, v in (("theta", P.theta), ("K_m", P.Km), ("k_e", P.ke),
                  ("D", D_AMP), ("lambda", LAMBDA), ("L", LB)):
        record(f"{nm} > 0", bool(v > 0), f"{nm} = {v:.6g}")

    t = np.linspace(0.0, 60 * DOSE_INTERVAL, 200001)
    d = np.asarray(delta(t), float)
    record("Delta(t) >= 0 on 60 cycles", bool(d.min() >= 0.0),
           f"min {d.min():.6g}, max {d.max():.6g}")
    record("Delta(t) <= its steady-state peak", bool(d.max() <= PEAK + 1e-9),
           f"max {d.max():.6f} against peak D/(1-e^-lambda T) = {PEAK:.6f}")
    late = float(np.atleast_1d(delta(np.array([50 * DOSE_INTERVAL - 1e-9])))[0])
    record("Delta returns to its steady trough each cycle",
           bool(abs(late - TROUGH) / TROUGH < 1e-6),
           f"at the end of cycle 50, Delta = {late:.6f} against the trough "
           f"{TROUGH:.6f}; peak/trough {PEAK/TROUGH:.4f} = "
           f"2^(T/half-life) = {2**(DOSE_INTERVAL*LAMBDA/np.log(2)):.4f}")

    # the balance k_e m^2 + (k_e K_m + theta + Delta - L) m - L K_m = 0
    # has product of roots -L K_m / k_e < 0, so exactly one root is positive.
    prod = -LB * P.Km / P.ke
    record("the balance has exactly one positive root",
           bool(prod < 0),
           f"product of roots = -L K_m/k_e = {prod:.4f} < 0 at every Delta, "
           f"so the roots straddle zero")
    m = np.asarray(quasi_steady_level(t), float)
    res = (P.ke * m ** 2
           + (P.ke * P.Km + P.theta + d - LB) * m - LB * P.Km)
    record("the closed form solves the balance",
           bool(np.max(np.abs(res)) < 1e-6 * LB * P.Km),
           f"worst residual {np.max(np.abs(res)):.3e} against L K_m = {LB*P.Km:.1f}")
    record("m(t) > 0 across the cycle", bool(m.min() > 0),
           f"m runs {m.min():.4f} to {m.max():.4f}, m_bar = {M_BAR:g}")


# ---------------------------------------------------------------------- 4
def total_load():
    head("4. A POSITIVE TOTAL LOAD, L + b > 0")
    half = 1.96 * LOAD_FRAC * LB
    record("L + b > 0 at b = 0", bool(LB > 0), f"L = {LB:.4f}")
    record("L + b > 0 across the load prior's 95% interval",
           bool(LB - half > 0),
           f"b in [{-half:.4f}, {half:.4f}], so L + b in "
           f"[{LB-half:.4f}, {LB+half:.4f}]")
    print(f"\n    The prior's 95% half-width is 1.96 x {LOAD_FRAC:g} x L = "
          f"{half:.4f}, which is")
    print(f"    {half/LB:.1%} of L. Positivity is nowhere near binding there; "
          f"it binds at")
    print(f"    b = -L = {-LB:.4f}, which is {LB/half:.1f} prior half-widths away.")
    return half


# ---------------------------------------------------------------------- 5
def ceiling():
    head("5. THE CAPACITY CEILING, WITH THE REFERENCE AS A PARAMETER")
    print("  The constraint is theta + Delta_peak <= normal capacity. `normal`")
    print("  is not a parameter of this repository; it is implied by a stated")
    print("  FRACTION, and src/common/params.py states two, four lines apart,")
    print("  on the parent theta row:")
    print("      `under 5% of control` in PA  -> normal = theta / 0.05  <-- ADOPTED")
    print("      `about 10% of normal`        -> normal = theta / 0.10")
    print("  5% IS THE REFERENCE (author, September 2026). The 10% reading is")
    print("  kept below as a check that the verdict does not turn on it.\n")
    print("  A second ambiguity has to be named before the arithmetic. The")
    print("  fraction is stated on the PARENT row (theta = "
          f"{_C.theta:.0f}); D and Delta are")
    print(f"  defined in the LOCAL block (theta = {_A.theta:.0f}). Evaluating the ceiling")
    print("  at all means carrying a fraction across models, which is an")
    print("  assumption and is made explicit here rather than silently.\n")
    print(f"    {'model':>8} {'theta':>8} {'Delta_peak':>11} {'theta+peak':>11} "
          f"{'fraction':>9} {'normal':>10} {'used':>8} {'headroom':>9}")
    print("  " + "-" * 88)
    rows = []
    for mod, th, dpk in (("local", _A.theta, PEAK),
                         ("parent", _C.theta, PEAK * _C.theta / _A.theta)):
        for frac in (REFERENCE_FRAC, 0.10):
            nrm = th / frac
            used = (th + dpk) / nrm
            rows.append((mod, frac, used))
            print(f"    {mod:>8} {th:8.0f} {dpk:11.2f} {th+dpk:11.2f} "
                  f"{frac:9.0%} {nrm:10.0f} {used:8.4f} {1-used:9.4f}")
    ok = all(u <= 1.0 for _, _, u in rows)
    record("theta + Delta_peak <= normal, at BOTH references", ok,
           f"worst used fraction {max(u for _, _, u in rows):.4f}")
    print("\n    The two readings give the same used fraction within a model,")
    print("    0.3133 at 10% and 0.1567 at 5%, because theta cancels: the ratio")
    print("    is frac x (1 + Delta_peak/theta) and Delta_peak/theta is the same")
    print("    2.1333 in both blocks once Delta is scaled with theta. So the")
    print("    cross-model question does not change the ANSWER, and the 5%/10%")
    print("    question changes it by exactly a factor of two.")
    print("    The ceiling is not breached on either reading. What the ambiguity")
    print("    costs is the HEADROOM, 0.687 against 0.843, and headroom is what")
    print("    an admissibility bound on b is spent against; see 6.")
    return rows


# ---------------------------------------------------------------------- 6
def fibre_bounds(half):
    head("6. PLAUSIBLE BOUNDS ON b: WHAT PHYSICS ALLOWS, AGAINST WHAT IS ASSUMED")
    K0, TH0, KE = P.Km, P.theta, P.ke
    Km_b = lambda b: K0 * LB / (LB + b)
    th_b = lambda b: TH0 + b + KE * (K0 - Km_b(b))
    print("  Along P5's exact fibre the data cannot separate these, so the only")
    print("  thing that bounds b is admissibility. The fibre is")
    print("      K_m(b) = K_0 L/(L+b),  theta(b) = theta_0 + b + k_e(K_0 - K_m(b)),")
    print("      k_e unchanged, and m is bit-identical along it.\n")
    lo_load = -LB
    lo_theta = brentq(th_b, -LB + 1e-9, 0.0, xtol=1e-12)
    lo = max(lo_load, lo_theta)
    print(f"    K_m(b) > 0  and  L + b > 0   bind at b = {lo_load:.4f}")
    print(f"    theta(b) > 0                 binds at b = {lo_theta:.4f}")
    print(f"    so the lower bound is        b > {lo:.4f}\n")
    his = {}
    for frac in (0.10, 0.05):
        nrm = TH0 / frac
        f = lambda b: th_b(b) + PEAK - nrm
        try:
            his[frac] = brentq(f, 0.0, 1e7, xtol=1e-9)
        except ValueError:
            his[frac] = np.inf
        print(f"    theta(b) + Delta_peak <= normal at {frac:.0%} binds at "
              f"b = {his[frac]:.4f}")
    print()
    print(f"    {'reference':>12} {'admissible b':>30} {'width':>12} "
          f"{'vs prior 95% width':>20}")
    print("  " + "-" * 82)
    for frac in (0.10, 0.05):
        w = his[frac] - lo
        print(f"    {frac:12.0%} [{lo:12.4f}, {his[frac]:12.4f}] {w:12.4f} "
              f"{w/(2*half):19.1f}x")
    print(f"\n    The load record's own 95% interval is [{-half:.4f}, {half:.4f}], "
          f"width {2*half:.4f}.")
    record("physics bounds b on both sides", bool(np.isfinite(his[0.05])),
           f"b in [{lo:.1f}, {his[0.10]:.1f}] at 10%, "
           f"[{lo:.1f}, {his[0.05]:.1f}] at 5%")
    ok = (his[0.10] - lo) > 2 * half
    record("the load prior is tighter than admissibility", bool(ok),
           f"prior width {2*half:.1f} against admissible "
           f"{his[0.10]-lo:.1f} (10%) and {his[0.05]-lo:.1f} (5%)")
    print("\n    This is the result that matters for Phase 3 section 4(a).")
    print("    Admissibility DOES close the fibre, so the constant-load case is")
    print("    not unbounded in the way `exactly singular` suggests. But the")
    print(f"    interval it leaves is {(his[0.10]-lo)/(2*half):.0f}x the load record's at the 10%")
    print(f"    reference and {(his[0.05]-lo)/(2*half):.0f}x at 5%, so on these numbers the record, and")
    print("    not physics, is carrying the identification. Replacing the prior")
    print("    with the admissibility bound would not rescue a constant-load")
    print("    design; it would only make its width finite.")
    print(f"\n    And this is the concrete cost of not choosing a reference: the")
    print(f"    pass/fail answers do not move, but the admissible range of the")
    print(f"    unidentified direction does, by {(his[0.05]-lo)/(his[0.10]-lo):.1f}x. A bound that is used")
    print(f"    quantitatively, rather than as a yes or no, needs the reference")
    print(f"    settled first.")
    return lo, his


# ---------------------------------------------------------------------- 7
def settling_exclusion():
    head("7. THE SETTLING EXCLUSION")
    settle, worst = G.settling_days(tol=0.01)
    settle5, worst5 = G.settling_days(tol=0.05)
    print(f"    window {settle*24:.3f} h at a 1% tolerance, {settle5*24:.3f} h "
          f"at the repository's 5%;")
    print(f"    g6 reconciles these against the paper's 4.2 h and P1f's 3.71 h.\n")
    r = g3.reachable(g3.AFTERNOON, settle=settle)
    booked = np.asarray(G.phase_spread(12, gap=30.4,
                                       phases=g3.spread_over(r, 12)), float)

    def inside(t, shift=0.0):
        ph = (np.asarray(t, float) + shift) % DOSE_INTERVAL
        return int(np.sum(ph < settle - 1e-9))

    print("  Every schedule in g1 and g2 starts ON the dose, at t = 0, so its")
    print("  first draw is at phase 0 by construction. That is a property of")
    print("  how the comparators were written, not of the interval. The `+0.5 d`")
    print("  column is the same schedule started half a day later, and `range`")
    print("  is what that offset does to the Delta/theta coverage it buys.\n")
    print(f"    {'schedule':>28} {'n':>3} {'inside':>7} {'+0.5 d':>7} "
          f"{'range':>8} {'range +0.5 d':>13}")
    print("  " + "-" * 82)
    cases = [("the booked 12-visit design", booked)] + \
            [(k, G.SCHEDULES[k](12)) for k in G.SCHEDULES]
    fixable = []
    for name, t in cases:
        a, b = inside(t), inside(t, 0.5)
        r0, r1 = G.delta_range(t), G.delta_range(np.asarray(t, float) + 0.5)
        fixable.append(a == 0 or b == 0)
        print(f"    {name:>28} {len(t):3d} {a:7d} {b:7d} {r0:8.4f} {r1:13.4f}")
    record("the booked design is outside the window as constructed",
           inside(booked) == 0, f"{inside(booked)} of {len(booked)} inside")
    record("every collision is removed by a half-day offset", all(fixable),
           f"{sum(fixable)} of {len(cases)} schedules clean at t0 = 0 or "
           f"t0 = 0.5 d")
    full = [g3.inc_at_phase(x) for x in np.arange(0.0, DOSE_INTERVAL, G.GRID)]
    use = [g3.inc_at_phase(x) for x in np.arange(0.0, DOSE_INTERVAL, G.GRID)
           if x >= settle - 1e-9]
    fr, ur = max(full) - min(full), max(use) - min(use)
    print("\n    Two things here came out against expectation.")
    print("\n    First, the settling collision is NEVER structural, not even for")
    print("    the fortnightly schedules. Shifting the start half a day moves")
    print("    every visit from phase 0 to phase 0.5, which is 12 h against a")
    print("    1.30 h window, and it clears all six.")
    print("\n    Second, the offset is NOT free, and this corrects Phase 3. The")
    print("    last two columns lose about 9% of the Delta/theta range, and the")
    print("    reason is that the phase being given up is phase 0, the PEAK,")
    print("    which is one endpoint of the very quantity f4 found binds.")
    print(f"    Costed properly on the half-day dial: excluding the window")
    print(f"    removes {1-ur/fr:.1%} of the reachable Delta/theta range ({fr:.4f} to")
    print(f"    {ur:.4f}) while removing {100*settle/DOSE_INTERVAL:.2f}% of the dial.")
    print("    g1's figure says the window `costs a schedule almost nothing`.")
    print("    That is true of the ANGLE and false of the VALUE, and the figure")
    print("    and the Phase 3 report are corrected to say so.")
    print("\n    The booked design already paid this: its most peak-ward phase is")
    print("    0.5, worth 1.9322, not phase 0's 2.1333. That is why its range")
    print("    goes UP under the offset, at the price of one visit landing in")
    print("    the window. It is at the constrained optimum, not below it.")
    print("\n    So the settling exclusion is a real constraint that is cheap to")
    print("    satisfy, and it is NOT an argument against any schedule. In")
    print("    particular it is not a second reason to reject `fortnightly, on")
    print("    the dose`: an offset start clears the window and leaves the rank")
    print("    deficiency exactly where it was. The rank objection is the whole")
    print("    of the case against that schedule, which is worth knowing,")
    print("    because a reader offered two objections may discount both.")


# ---------------------------------------------------------------------- 8
def d_calibration():
    head("8. THE D CALIBRATION: `chosen so the marker falls by about two thirds`")
    print("  src/common/params.py's D row states D = 6000 was chosen so the")
    print("  marker falls by about two thirds at peak effect, and prints the")
    print("  computed range. The claim needs two things named, and the row now")
    print("  names both: fall FROM what, and which peak.\n")
    # the untreated level: the balance at Delta = 0
    B0 = P.ke * P.Km + P.theta - LB
    m_untreated = (-B0 + np.sqrt(B0 ** 2 + 4.0 * P.ke * LB * P.Km)) / (2.0 * P.ke)
    readings = []
    for lab, dpk in (("first dose, D = 6000", D_AMP),
                     ("steady state, D/(1-e^-lT)", PEAK)):
        Bp = P.ke * P.Km + P.theta + dpk - LB
        mp = (-Bp + np.sqrt(Bp ** 2 + 4.0 * P.ke * LB * P.Km)) / (2.0 * P.ke)
        readings.append((lab + ", from untreated", m_untreated, mp))
    Bt = P.ke * P.Km + P.theta + TROUGH - LB
    m_trough = (-Bt + np.sqrt(Bt ** 2 + 4.0 * P.ke * LB * P.Km)) / (2.0 * P.ke)
    Bp = P.ke * P.Km + P.theta + PEAK - LB
    m_peak = (-Bp + np.sqrt(Bp ** 2 + 4.0 * P.ke * LB * P.Km)) / (2.0 * P.ke)
    readings.append(("steady state, from the pre-dose trough", m_trough, m_peak))
    print(f"    {'reading':>44} {'from':>9} {'to':>9} {'fall':>9} {'half?':>7}")
    print("  " + "-" * 84)
    for lab, a, b in readings:
        fall = 1.0 - b / a
        print(f"    {lab:>44} {a:9.4f} {b:9.4f} {fall:9.2%} "
              f"{'yes' if 0.4 <= fall <= 0.6 else 'no':>7}")
    falls = [1.0 - b / a for _, a, b in readings]
    lo, hi = min(falls), max(falls)
    # The row was reworded in September 2026, from "about half" to "about two
    # thirds" with the computed range beside it. This check is now a DRIFT
    # check on that wording rather than a verdict on the old one: every reading
    # must sit inside the range the row prints, and the row must still name the
    # fraction. Two thirds is 66.7%, and the readings run 62.25% to 66.36%.
    import re as _re
    from common.params import LOCAL as _LOCAL
    row = next(r for r in _LOCAL if r.symbol == r"$D$")
    claimed = _re.search(r"\$(\d+)\$--\$(\d+)\\%\$", row.note)
    band_ok = claimed is not None and (
        int(claimed.group(1)) <= round(100 * lo) and
        round(100 * hi) <= int(claimed.group(2)))
    named_ok = "two thirds" in row.note
    ok = band_ok and named_ok and all(0.60 <= f <= 0.70 for f in falls)
    record("the D row's stated fall matches the computed one", ok,
           f"row says {claimed.group(0) if claimed else 'NO RANGE FOUND'}, "
           f"computed {100*lo:.2f}% to {100*hi:.2f}%, "
           f"names two thirds: {named_ok}")
    print(f"\n    Every reading is larger than half and close to two thirds:")
    print(f"    {100*lo:.2f}% to {100*hi:.2f}%. The row now says `about two thirds`")
    print("    and prints the range, with the two conventions that make it vary:")
    print("    which baseline the fall is read from, and which peak.")
    print("\n    What the row no longer says is that only the excursion is used.")
    print("    That was true and read as reassurance; PHASE5_ROBUSTNESS shows a")
    print("    defensible range of D moving the headline design from 10 monthly")
    print("    visits to more than 60, so the value matters and not just the")
    print("    order. The trial's own falls are smaller than this calibration.")


def write_paper_table():
    """
    The twenty checks and their verdicts, for the paper.

    Physical admissibility is a property of the design a reader cannot see
    from a precision, and a fit that balances the data while describing no
    patient is the failure the paper is about. So the checks go in, as a
    table, rather than being asserted in a sentence.
    """
    import sys
    from pathlib import Path as _P
    sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "src"))
    from common.style import write_table
    import re as _re

    #: plain-text token -> what the paper prints for it. Order matters only
    #: in that the regex below prefers the longest alternative, so
    #: "Delta_peak" cannot be eaten by "Delta".
    _SYMBOLS = {
        "theta_crit": r"$\thc$", "Delta_peak": r"$\Delta_{\mathrm{peak}}$",
        "Delta(t)": r"$\Delta(t)$", "Delta": r"$\Delta$",
        "theta": r"$\theta$", "K_m": r"$K_m$", "k_e": r"$k_e$",
        "lambda": r"$\lambda$", "m_bar": r"$\bar m$", "m(t)": r"$m(t)$",
        "sd(ln theta)": r"$\operatorname{sd}(\ln\hat\theta)$",
        "<=": r"$\le$", ">=": r"$\ge$", ">": r"$>$", "<": r"$<$",
    }
    _TOKEN = _re.compile("|".join(_re.escape(k) for k in
                                  sorted(_SYMBOLS, key=len, reverse=True)))

    def tex(t):
        """
        The checks are written for a terminal; the paper is not a terminal.

        Symbols become symbols and LaTeX's reserved characters are escaped,
        and the two passes cannot interfere: the escaping runs only on the
        text BETWEEN symbol matches, so an underscore this function inserts
        as part of $K_m$ is never escaped back out again.
        """
        def esc(x):
            for a, b in (("\\", r"\textbackslash{}"), ("%", r"\%"),
                         ("&", r"\&"), ("_", r"\_"), ("#", r"\#"),
                         ("$", r"\$"), ("^", r"\^{}"), ("~", r"\~{}")):
                x = x.replace(a, b)
            return x
        out, last = [], 0
        for m in _TOKEN.finditer(t):
            out.append(esc(t[last:m.start()]))
            out.append(_SYMBOLS[m.group(0)])
            last = m.end()
        out.append(esc(t[last:]))
        return "".join(out)

    rows = []
    for tag, ok, detail in RESULTS:
        verdict = (r"\textbf{pass}" if ok is True
                   else (r"\textbf{FAIL}" if ok is False else "---"))
        #  the detail column is evidence, not prose: keep it short enough to
        #  read across a page and let the audit record carry the rest
        d = detail if len(detail) <= 72 else detail[:69].rstrip() + r"\,\ldots"
        rows.append([tex(tag), verdict, tex(d)])
    #  the main text wants the verdicts at a glance and the supplement wants
    #  the evidence, so both are written rather than one being cut down by
    #  hand in the .tex
    half = (len(RESULTS) + 1) // 2
    pair = []
    for k in range(half):
        left = rows[k][:2]
        right = rows[k + half][:2] if k + half < len(rows) else ["", ""]
        pair.append(left + right)
    write_table("tab_S_admissible_compact",
                ["constraint", "verdict", "constraint", "verdict"],
                pair, colspec=(r">{\raggedright\arraybackslash}p{0.38\linewidth}l "
                         r">{\raggedright\arraybackslash}p{0.38\linewidth}l"),
                note=("%d checks, %d pass, %d fail; the evidence for each is "
                      "in the supplement"
                      % (len(RESULTS),
                         sum(1 for _, ok, _ in RESULTS if ok is True),
                         sum(1 for _, ok, _ in RESULTS if ok is False))))
    print(f"  wrote tab_S_admissible_compact.tex")
    write_table("tab_S_admissible",
                ["constraint", "verdict", "evidence"],
                rows,
                colspec=(r">{\raggedright\arraybackslash}p{0.34\linewidth}l "
                         r">{\raggedright\arraybackslash}p{0.44\linewidth}"),
                note=("%d checks, %d pass, %d fail"
                      % (len(RESULTS),
                         sum(1 for _, ok, _ in RESULTS if ok is True),
                         sum(1 for _, ok, _ in RESULTS if ok is False))))
    print(f"  wrote tab_S_admissible.tex ({len(RESULTS)} checks)")


# ------------------------------------------------------------------- table
def table():
    head("PASS / FAIL")
    w = max(len(t) for t, _, _ in RESULTS)
    npass = sum(1 for _, ok, _ in RESULTS if ok is True)
    nfail = sum(1 for _, ok, _ in RESULTS if ok is False)
    for t, ok, d in RESULTS:
        print(f"    {'PASS' if ok else 'FAIL':>4}  {t:<{w}}  {d}")
    print(f"\n    {npass} pass, {nfail} fail, {len(RESULTS)} checks.")
    print("\n    Every admissibility constraint proper holds, at the adopted 5%")
    print("    reference and at the 10% reading kept as a check, and the")
    print("    recommended design satisfies all of them as constructed.")
    print("\n    Check 8 was the one failure in the September 2026 run, when the")
    print("    D row said the marker `falls by about half`. The row now says")
    print("    `about two thirds` and prints the computed range, so check 8 is a")
    print("    drift check on that range rather than a verdict on the old claim.")
    print("\n    Two things the checks corrected in this audit's own earlier")
    print("    work rather than in the repository: the settling window costs")
    print("    10.1% of the Delta/theta range and not `almost nothing` (check")
    print("    7), and the settling collision is removable by a half-day offset")
    print("    for every schedule including the fortnightly ones, so it is not a")
    print("    second objection to `fortnightly, on the dose`.")


if __name__ == "__main__":
    print(__doc__)
    rates_and_increment()
    half = total_load()
    ceiling()
    fibre_bounds(half)
    settling_exclusion()
    d_calibration()
    table()
    write_paper_table()
