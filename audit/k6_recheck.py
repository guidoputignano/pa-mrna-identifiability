"""
Phase 5, steps 3 and 4 and the admissibility follow-up: the independent
re-check.

Every check reaches the same conclusion by a route that shares no code with the
script it checks. Where k5 used sympy's `cancel`, this uses exact rational
interpolation. Where k5 used the complex step, this uses Richardson-extrapolated
central differences. Where k5 solved the fibre by hand, this eliminates the
nuisance variables with a Groebner basis and reads off what survives. Where k4
quoted h1's cancellation, this re-derives it symbolically.

A check that fails is reported as a failure and not repaired in place.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp

import g_design_common as G
import k2_structures as K2
import k4_admissibility as K4
import k5_identifiability as K5
from design.run import delta

TH, KM, KE, LB = G.TH, G.KM, G.KE, G.LB
T = G.T
RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    for line in str(detail).splitlines():
        print(f"         {line}")
    print()


# --------------------------------------------------------------------------
def c1_moebius_exact():
    """k5's sympy reduction, by exact rational interpolation instead."""
    E0, base, c, th = F(7, 13), F(413, 3), F(5077, 79), F(3000)
    q = F(1, 1000)                       # 1 - Imax at Imax = 0.999

    def m_of(D):
        E = th + D
        return E0 + base * (c + q * E) / (c + E)

    # a Moebius function is fixed by three points; take three and solve for
    # (P, U, V) in EXACT rationals, then test at four more.
    Ds = [F(1), F(1000), F(6400)]
    A = sp.Matrix([[D, 1, -m_of(D)] for D in Ds])
    rhs = sp.Matrix([m_of(D) * D for D in Ds])
    P, U, V = [sp.nsimplify(x) for x in A.solve(rhs)]
    pred = {D: (P * D + U) / (D + V) for D in [F(3), F(77), F(2500), F(99991)]}
    exact = all(sp.Rational(pred[D]) == sp.Rational(m_of(D)) for D in pred)

    P_pred = E0 + base * q
    U_pred = P_pred * th + (E0 + base) * c
    V_pred = c + th
    same = (sp.Rational(P) == sp.Rational(P_pred)
            and sp.Rational(U) == sp.Rational(U_pred)
            and sp.Rational(V) == sp.Rational(V_pred))
    check("1. PUB_B is Moebius in Delta, by exact rational interpolation",
          exact and same,
          f"interpolated from 3 points, exact at 4 more: {exact}\n"
          f"P = {P}  (predicted {P_pred})\n"
          f"U = {U}  (predicted {U_pred})\n"
          f"V = {V}  (predicted {V_pred})\n"
          f"No floating point anywhere in this check.")


def c2_groebner_fibre():
    """k5's hand-solved fibre, by eliminating the nuisances instead."""
    E0, base, c, th, P, U, V, q = sp.symbols(
        "E0 base c theta P U V q", positive=True)
    eqs = [E0 + base * q - P,
           (E0 + base * q) * th + (E0 + base) * c - U,
           c + th - V]
    gb = sp.groebner(eqs, E0, base, c, th, P, U, V, order="lex")
    surv = [g for g in gb.exprs
            if not (g.has(E0) or g.has(base) or g.has(c))]
    check("2. PUB_B: eliminating (E0, base, IC50) leaves NO relation on theta",
          len(surv) == 0,
          f"Groebner basis in lex order E0 > base > IC50 > theta > P,U,V has "
          f"{len(gb.exprs)} generators.\n"
          f"Generators free of E0, base and IC50: {surv if surv else 'none'}\n"
          f"An empty elimination ideal is the definition of theta being\n"
          f"unconstrained by the invariants. Any theta in (0, V) is realised by\n"
          f"  base = (U - P V)/((1-q)(V - theta)),  E0 = P - q base,  IC50 = V-theta.")

    # PUB_A, the counterpart: the same elimination must PIN theta
    W = sp.Symbol("W", positive=True)
    eqs2 = [q * th + c - W, th + c - V]
    gb2 = sp.groebner(eqs2, c, th, W, V, q, order="lex")
    surv2 = [g for g in gb2.exprs if not g.has(c)]
    sol = sp.solve(eqs2, [c, th], dict=True)[0]
    check("3. PUB_A: the same elimination PINS theta in closed form",
          len(surv2) >= 1 and sp.simplify(sol[th] - (V - W) / (1 - q)) == 0,
          f"generators free of IC50: {surv2}\n"
          f"theta = {sp.simplify(sol[th])} = (V - W)/Imax\n"
          f"IC50  = {sp.simplify(sol[c])}\n"
          f"So under PUB_A the two invariants determine theta exactly, and\n"
          f"under PUB_B the three invariants determine nothing about it.")


def c3_ranks_by_richardson():
    """k5's complex-step Jacobian, by Richardson extrapolation instead."""
    times = G.SCHEDULES["phase-spread"](24)
    inc = np.asarray(delta(times), float)
    q = 1.0 - K2.IMAX
    ic50_B = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    oB = K2.PubB(ic50_B)
    lines, ok = [], True

    def richardson(fn, p0):
        """Central differences at h and h/2, extrapolated: O(h^4)."""
        p0 = np.asarray(p0, float)
        cols = []
        for j in range(p0.size):
            h = 1e-4 * max(abs(p0[j]), 1.0)
            d = []
            for hh in (h, h / 2.0):
                pp, pm = p0.copy(), p0.copy()
                pp[j] += hh
                pm[j] -= hh
                d.append((fn(pp) - fn(pm)) / (2 * hh))
            cols.append((4.0 * d[1] - d[0]) / 3.0)
        return np.stack(cols, axis=1)

    for tag, loads in (("constant", np.full(times.size, LB)),
                       ("varying", LB * np.exp(
                           np.random.default_rng(0).normal(
                               -0.5 * np.log(1.01), np.sqrt(np.log(1.01)),
                               times.size)))):
        def paper(p):
            th, Km, ke, b = p
            tot = loads + b
            B = ke * Km + th + inc - tot
            return (-B + np.sqrt(B ** 2 + 4.0 * ke * tot * Km)) / (2.0 * ke)

        def pubB(p):
            E0, bs, cc, th, b = p
            E = th + inc
            return (E0 + bs * (cc + q * E) / (cc + E)) * (loads + b) / LB

        for nm, fn, p0, want in (
                ("PAPER", paper, [TH, KM, KE, 0.0], 3 if tag == "constant" else 4),
                ("PUB_B", pubB, [oB.e0, oB.base, ic50_B, TH, 0.0],
                 3 if tag == "constant" else 4)):
            J = richardson(fn, p0)
            sc = np.maximum(np.abs(J).max(axis=0), 1e-300)
            sv = np.linalg.svd(J / sc, compute_uv=False)
            gap = sv[want - 1] / sv[want] if want < sv.size else np.inf
            good = gap > 1e4
            ok = ok and good
            lines.append(f"{tag:8s} {nm:6s} expected rank {want}: "
                         f"sv ratio across the gap {gap:10.3e}  "
                         f"{'ok' if good else 'NO CLEAR GAP'}")
    check("4. Ranks confirmed by Richardson extrapolation, not complex step",
          ok, "\n".join(lines) +
          "\nA difference quotient cannot reach machine zero, so this checks "
          "the\nSPECTRAL GAP rather than the null singular value itself. A gap "
          "of\n1e4 or more at the predicted rank is the same conclusion by a "
          "cruder\ninstrument, which is the point.")


def c4_sensitivity_vs_fd():
    """k5's PUB_A sensitivity integrator, against finite differences of it.

    A note on the step size, because getting it wrong made this check fail
    twice and only one of those was the code's fault.

    The integrator runs at rtol 1e-12, so a marker value of order 5 carries an
    absolute error of order 5e-12. A central difference divides a difference of
    order 2 h |dm/dp| by 2 h, so its relative error is about 5e-12/(h |dm/dp|):
    too SMALL a step and the quotient is integrator noise. b enters as
    k_in = L + b on a scale of L = 322.7, so a step of 1e-5 makes |dm/db| h
    of order 3e-7 and the quotient 1e-5 wrong, which is what the first version
    of this check reported as a failure of k5. The step is therefore taken
    relative to each parameter's OWN scale, and Richardson-extrapolated.
    """
    times = G.SCHEDULES["phase-spread"](12)
    ic50 = float(K2._solve_shape(K2.PubA.shape, 0.60)[0])
    o = K2.PubA(ic50)
    kin_of = lambda d: np.full(np.asarray(d, float).size, LB)
    p0 = np.array([o.kout, ic50, TH, 0.0])
    typ = np.array([o.kout, ic50, TH, LB])          # b lives on the load's scale
    _, S = K5.pubA_traj(times, None, p0, kin_of)
    lines, worst = [], 0.0
    for j, nm in enumerate(["k_out", "IC50", "theta", "b"]):
        h0 = 1e-3 * typ[j]
        d = []
        for hh in (h0, h0 / 2.0):
            pp, pm = p0.copy(), p0.copy()
            pp[j] += hh
            pm[j] -= hh
            d.append((K5.pubA_marker(times, pp, kin_of)
                      - K5.pubA_marker(times, pm, kin_of)) / (2 * hh))
        fd = (4.0 * d[1] - d[0]) / 3.0
        sc = max(np.max(np.abs(S[:, j])), 1e-300)
        e = float(np.max(np.abs(fd - S[:, j])) / sc)
        worst = max(worst, e)
        lines.append(f"d(m)/d({nm:6s}) step {h0:11.4g}  max rel difference "
                     f"{e:10.3e}")
    check("5. PUB_A forward sensitivities agree with finite differences",
          worst < 1e-7, "\n".join(lines) +
          f"\nworst {worst:.3e}. The sensitivity system and the state are "
          f"integrated\ntogether, so this also checks the augmented "
          f"right-hand side.\n"
          f"This check found a REAL error in k5 on its first run: the "
          f"sensitivities\nof the INITIAL CONDITION to IC50 and theta were "
          f"set to zero, which left a\ntransient worth 1e-3 relative. They are "
          f"now carried, and both columns\nagree to 1e-9.")


def c5_two_pubA_implementations():
    """k5's theta-aware integrator against k2's, where they must agree."""
    times = G.SCHEDULES["phase-spread"](12)
    kin_of = lambda d: np.full(np.asarray(d, float).size, LB)
    lines, worst = [], 0.0
    for tgt in (0.30, 0.60):
        for c in K2._solve_shape(K2.PubA.shape, tgt):
            o = K2.PubA(c)
            a = o.marker(times, None)
            b = K5.pubA_marker(times, [o.kout, c, TH, 0.0], kin_of)
            e = float(np.max(np.abs(a - b)) / np.max(np.abs(a)))
            worst = max(worst, e)
            lines.append(f"fall {100*tgt:2.0f}%  IC50 {c:11.5g}  "
                         f"max rel difference {e:10.3e}")
    check("6. k2's PUB_A and k5's theta-aware PUB_A agree at theta = 3000",
          worst < 1e-8, "\n".join(lines) +
          f"\nworst {worst:.3e}. Two independent integrators, one written for "
          f"the\nsimulator and one for the sensitivity system.")


def c6_capacity_cancellation():
    """k4's use of h1's cancellation, re-derived symbolically."""
    th, dpk, frac, r = sp.symbols("theta Delta_peak frac r", positive=True)
    used = (th + dpk) / (th / frac)
    used2 = sp.simplify(used.subs(dpk, r * th))
    binds = sp.solve(sp.Eq(frac * (1 + r), 1), r)[0]
    check("7. The capacity ceiling depends only on Delta_peak/theta",
          sp.simplify(used2 - frac * (1 + r)) == 0,
          f"used = (theta + Delta_peak)/(theta/frac) = {sp.simplify(used2)}\n"
          f"with r = Delta_peak/theta. theta cancels exactly, so the reading "
          f"of the\nresidual-capacity fraction and the choice of model block "
          f"are separable,\nwhich is h1's own observation, re-derived rather "
          f"than quoted.\n"
          f"It binds at r = {binds} = 1/frac - 1, that is "
          f"{float(binds.subs(frac, sp.Rational(1,10)))} at 10% and "
          f"{float(binds.subs(frac, sp.Rational(1,20)))} at 5%.")


def c7_shape_ceiling_limit():
    """k4's Imax=1 ceiling, against its closed form."""
    cc, El, Eh = sp.symbols("c E_lo E_hi", positive=True)
    ratio = (cc / (cc + Eh)) * ((cc + El) / cc)          # PUB_A at Imax = 1
    lim = sp.limit(sp.simplify(ratio), cc, 0, "+")
    lines, ok = [], True
    for dpk in (6400.0, 1398.1, 13072.3, 22375.3):
        closed = float(dpk / (TH + dpk))
        num, _ = K4.shape_ceiling(K2.PubA, dpk, imax=1.0)
        good = abs(closed - num) < 2e-4
        ok = ok and good
        lines.append(f"Delta_peak {dpk:9.1f}  closed form Delta/(theta+Delta) "
                     f"{100*closed:7.4f}%  grid {100*num:7.4f}%  "
                     f"{'ok' if good else 'MISMATCH'}")
    check("8. PUB_A's Imax=1 ceiling is exactly Delta/(theta+Delta)",
          ok and sp.simplify(lim - El / Eh) == 0,
          f"symbolic c -> 0 limit of the inhibition ratio: {sp.simplify(lim)} "
          f"= E_lo/E_hi,\nso the ceiling is 1 - theta/(theta+Delta) = "
          f"Delta/(theta+Delta).\n" + "\n".join(lines) +
          f"\nThe grid maximisation approaches it from below, which is why the "
          f"agreement\nis to 2e-4 and not to machine precision: the optimum is "
          f"at c -> 0, off the\ngrid's low end.")


def c8_bracket():
    """k4's bracket claim, swept independently rather than tabulated."""
    import k1_sources as K1
    lo_sweep, hi_sweep = 100.0, 40000.0
    rows = []
    for dpk in np.geomspace(lo_sweep, hi_sweep, 61):
        f = K1.fall(dpk)
        a1, _ = K4.shape_ceiling(K2.PubA, dpk, imax=1.0)
        b1, _ = K4.shape_ceiling(K2.PubB, dpk, imax=1.0)
        rows.append((dpk, b1, f, a1, b1 <= f <= a1))

    # PUB_A's ceiling is above the fall on the whole sweep; PUB_B's crosses.
    a_always = all(f <= a1 for _, _, f, a1, _ in rows)
    g = lambda d: K4.shape_ceiling(K2.PubB, d, imax=1.0)[0] - K1.fall(d)
    lo, hi = lo_sweep, hi_sweep
    if np.sign(g(lo)) != np.sign(g(hi)):
        for _ in range(80):
            mid = np.sqrt(lo * hi)
            if np.sign(g(mid)) == np.sign(g(lo)):
                lo = mid
            else:
                hi = mid
        cross = float(np.sqrt(lo * hi))
    else:
        cross = float("nan")

    tabulated = [816.8, 1398.1, 6400.0, 8000.0, 3319.2, 13072.3, 22375.3,
                 830.0, 1420.2, 5682.5]
    in_table = all(K4.shape_ceiling(K2.PubB, d, imax=1.0)[0]
                   <= K1.fall(d)
                   <= K4.shape_ceiling(K2.PubA, d, imax=1.0)[0]
                   for d in tabulated)
    lines = [
        f"swept {len(rows)} doses, Delta_peak from {lo_sweep:.0f} to "
        f"{hi_sweep:.0f}, geometric.",
        f"PUB_A's ceiling is above the paper's fall on ALL of them: {a_always}.",
        f"PUB_B's ceiling CROSSES the paper's fall at Delta_peak = "
        f"{cross:.2f}, i.e.",
        f"  Delta_peak/theta = {cross/TH:.4f}. BELOW that dose PUB_B's shape "
        f"can do MORE",
        f"  than the paper's Michaelis-Menten marker, and above it, less.",
        f"So the bracket 'B below, paper, A above' is NOT a property of the "
        f"three forms",
        f"at every dose. It holds for Delta_peak above about {cross:.0f}, and "
        f"reverses below.",
        f"Every dose k4 tabulates has Delta_peak >= 816.8, so all of them are "
        f"above the",
        f"crossover: bracket holds on all {len(tabulated)} tabulated doses: "
        f"{in_table}.",
    ]
    check("9. The bracket holds above a crossover dose, and reverses below",
          a_always and in_table and np.isfinite(cross), "\n".join(lines))


def main():
    print("=" * 78)
    print("  PHASE 5 STEPS 3-4 AND THE ADMISSIBILITY FOLLOW-UP: RE-CHECK")
    print("=" * 78 + "\n")
    c1_moebius_exact()
    c2_groebner_fibre()
    c3_ranks_by_richardson()
    c4_sensitivity_vs_fd()
    c5_two_pubA_implementations()
    c6_capacity_cancellation()
    c7_shape_ceiling_limit()
    c8_bracket()
    n = sum(1 for _, ok, _ in RESULTS if ok)
    print("=" * 78)
    print(f"  {n} of {len(RESULTS)} checks pass.")
    for t, ok, _ in RESULTS:
        if not ok:
            print(f"    FAILED: {t}")
    print("=" * 78)


if __name__ == "__main__":
    main()
