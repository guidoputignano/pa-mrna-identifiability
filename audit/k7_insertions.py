"""
Phase 5, the follow-up to steps 3 and 4: where the load enters PUB_B, and what
the non-identifiability of theta actually rests on.

Three questions, all of them raised against the previous round's report.

  1. PUB_B AS PUBLISHED HAS NO LOAD INPUT AT ALL. Attarwala's marker equation
     is 3-HP = E0 + base(1 - Imax PCC/(IC50+PCC)), fitted in a fixed-diet
     animal colony. Every load in our version is OUR insertion. k2_structures
     chose one of several and the previous report said "a varying load does not
     help, because the load enters only as a common factor". That sentence is
     true OF THAT INSERTION. This script asks whether the conclusion survives
     the others.

         whole   m = (E0 + base f) u            <- what k2 does
         base    m = E0 + base u f              <- the load lifts only the
                                                   suppressible pool
         E0      m = E0 u + base f              <- the load lifts only the
                                                   unaffected floor
         add     m = E0 + base f + gamma u      <- the load is a separate
                                                   additive source

     with u = (L+b)/L0 and f(Delta) = (IC50 + q(theta+Delta))/(IC50+theta+Delta),
     q = 1 - Imax.

  2. THE INTUITION. "theta enters only through theta + IC50, so residual
     enzyme and potency trade off." Confirmed or corrected here, by walking
     each fibre and printing which invariants move.

  3. WHICH PARAMETER IS RESPONSIBLE. PUB_B is itself a DIRECT-RESPONSE form:
     its marker has no state, so it is an instantaneous map from enzyme to
     marker, exactly as the paper's quasi-steady balance is. The difference
     between them is therefore the SHAPE of the map and not the dynamics, and
     the question is whether the deficit comes from PUB_B's phenomenological
     parameters. Tested by removing them one at a time.

Method as in k5: the complex step, which is exact to machine precision, and
explicit construction of every fibre found. Nothing here records an expected
result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp

import g_design_common as G
import k2_structures as K2
import k5_identifiability as K5
from design.run import delta

TH, KM, KE, LB = G.TH, G.KM, G.KE, G.LB
Q = 1.0 - K2.IMAX
RESULTS = []


def head(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


def record(tag, ok, detail):
    RESULTS.append((tag, bool(ok), detail))
    print(f"\n  [{'PASS' if ok else 'FAIL'}] {tag}")
    for line in str(detail).splitlines():
        print(f"         {line}")


# --------------------------------------------------------------------------
#  the four insertions, each as a closure over (times, loads)
# --------------------------------------------------------------------------
def _f(c, th, inc):
    E = th + inc
    return (c + Q * E) / (c + E)


INSERTIONS = {
    "whole  (E0+base f)*u": (
        ["E0", "base", "IC50", "theta", "b"],
        lambda p, inc, L: (p[0] + p[1] * _f(p[2], p[3], inc)) * (L + p[4]) / LB),
    "base   E0 + base*u*f": (
        ["E0", "base", "IC50", "theta", "b"],
        lambda p, inc, L: p[0] + p[1] * (L + p[4]) / LB * _f(p[2], p[3], inc)),
    "E0     E0*u + base*f": (
        ["E0", "base", "IC50", "theta", "b"],
        lambda p, inc, L: p[0] * (L + p[4]) / LB + p[1] * _f(p[2], p[3], inc)),
    "add    E0+base f+g*u": (
        ["E0", "base", "IC50", "theta", "gamma", "b"],
        lambda p, inc, L: (p[0] + p[1] * _f(p[2], p[3], inc)
                           + p[4] * (L + p[5]) / LB)),
}

#: PUB_B exactly as published, with no load input anywhere. `b` is carried in
#: the parameter vector so that the table can show it having NO effect, which
#: is the point: a load that does not enter cannot identify anything.
AS_PUBLISHED = {
    "as published (no load)": (
        ["E0", "base", "IC50", "theta", "b"],
        lambda p, inc, L: p[0] + p[1] * _f(p[2], p[3], inc)),
}

#: the same four with the unsuppressible offset REMOVED, E0 identically zero
NO_E0 = {
    "whole, E0 = 0": (
        ["base", "IC50", "theta", "b"],
        lambda p, inc, L: p[0] * _f(p[1], p[2], inc) * (L + p[3]) / LB),
    "E0-route, E0 = 0": (
        ["base", "IC50", "theta", "b"],
        lambda p, inc, L: p[0] * _f(p[1], p[2], inc)),
}


def truth_for(names, ic50, base, e0):
    v = dict(E0=e0, base=base, IC50=ic50, theta=TH, b=0.0, gamma=1.0)
    return np.array([v[n] for n in names], float)


def ranks(times, loads, tag):
    inc = np.asarray(delta(times), float)
    ic50 = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    o = K2.PubB(ic50)
    print(f"\n  {tag}: {len(times)} visits, "
          f"{'VARYING' if np.ptp(loads) > 0 else 'CONSTANT'} recorded load.")
    print(f"      {'insertion':22s} {'params':>7} {'rank':>5} {'deficit':>8}  "
          f"{'theta in the null space?':>26}")
    print("  " + "-" * 76)
    out = {}
    for name, (names, fn) in {**AS_PUBLISHED, **INSERTIONS,
                             **NO_E0}.items():
        p0 = truth_for(names, ic50, o.base, o.e0)
        J = K5.jac_complex(lambda p: fn(p, inc, loads), p0)
        # A parameter the record cannot see AT ALL (an identically zero
        # column) is its own null direction. Scaling such a column before the
        # SVD divides by ~1e-300, so it is split out first: the rank is taken
        # on the live columns and the dead ones are added back as unit
        # vectors. Getting this wrong loses the REAL null direction, which is
        # what a first version of this block did.
        dead = [k for k in range(J.shape[1]) if np.max(np.abs(J[:, k])) == 0.0]
        live = [k for k in range(J.shape[1]) if k not in dead]
        r_live, sv, null_live = K5.rank_of(J[:, live])
        r = r_live
        vecs = []
        for k in dead:
            e = np.zeros(len(names))
            e[k] = 1.0
            vecs.append(e)
        if null_live is not None:
            for v in null_live:
                e = np.zeros(len(names))
                e[live] = v
                vecs.append(e / np.linalg.norm(e))
        null = np.stack(vecs) if vecs else None
        j = names.index("theta")
        if null is None:
            verdict = "no null space: IDENTIFIED"
        else:
            w = float(np.max(np.abs(null[:, j])))
            verdict = (f"yes, |theta| = {w:.4f}: FREE" if w > 1e-6
                       else f"no, |theta| = {w:.1e}: PINNED")
        out[name] = (r, len(names), null, names, verdict)
        print(f"      {name:22s} {len(names):7d} {r:5d} {len(names)-r:8d}  "
              f"{verdict:>26}")
    print("\n      null direction(s):")
    for name, (r, n, null, names, _) in out.items():
        if null is None:
            continue
        for v in null:
            print(f"      {name:22s} "
                  + ", ".join(f"{a}={x:+.4f}" for a, x in zip(names, v)))
    return out


def bias_by_insertion(seeds=range(8), n=12):
    """
    Fit the PAPER'S model to records generated by each insertion.

    Step 2 measured theta_hat/theta under the `whole` insertion only, so the
    headline bias figures are conditional on it. This asks whether they are.
    The ordinate, its noise and the estimator are i1_estimator's, unchanged.
    """
    import i1_estimator as I1
    head("1b. IS STEP 2's BIAS FIGURE INSERTION-ROBUST?")
    print("""
  Each generator below is one of the insertions, calibrated to the same
  operating point and the same 60% fall at the treated peak. The fit is the
  paper's model, by the same variable-projection scan and polish step 2 used.
  RMS residual is in units of the repository's own sigma; medians over 8 daily
  diet records, with the range on theta_hat/theta.
""")
    ic50 = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    o = K2.PubB(ic50)
    exp_rms = np.sqrt((n - 4) / n)
    print(f"  A well-specified fit of n = {n} points for 4 parameters has an "
          f"expected\n  RMS/sigma of sqrt(({n}-4)/{n}) = {exp_rms:.3f}. Both "
          f"residual columns are read\n  against that. 'structure' is the "
          f"noiseless fit, which is step 2's own\n  headline quantity; 'noisy' "
          f"carries the repository's ordinate noise.\n")
    print(f"      {'insertion':22s} {'schedule':26s} {'structure':>10} "
          f"{'noisy':>7} {'theta_hat/theta':>24}")
    print("  " + "-" * 92)
    out = {}
    for name, (names, fn) in {**AS_PUBLISHED, **INSERTIONS}.items():
        p0 = truth_for(names, ic50, o.base, o.e0)
        for sname in ("monthly (30.4 d)", "phase-spread"):
            times = G.SCHEDULES[sname](n)
            inc = np.asarray(delta(times), float)
            rs, r0s, ts = [], [], []
            for sd in seeds:
                path = K2.daily_loads(0.10, sd, times.min() - 50.0,
                                      times.max() + 2.0)
                L = path(times)
                m = np.asarray(fn(p0, inc, L), float)
                rng = np.random.default_rng(7_000_003 * sd + 11)
                noisy = (L * (1.0 + I1.LOAD_FRAC
                              * rng.standard_normal(times.size))
                         + I1.BLOC * I1.CV * 10.0
                         * rng.standard_normal(times.size))
                for ordn, acc in ((L.copy(), r0s), (noisy, rs)):
                    fit = I1.polish(I1.varpro(m, inc, ordn)[0],
                                    m, inc, ordn)[0]
                    r = ((fit["theta"] + inc) * m / (fit["Km"] + m)
                         + fit["ke"] * m + fit["b"] - ordn)
                    acc.append(float(np.sqrt(np.mean(r ** 2)) / G.SIG))
                    if acc is rs:
                        ts.append(fit["theta"] / TH)
            out[(name, sname)] = (float(np.median(r0s)), float(np.median(rs)),
                                  float(np.median(ts)), float(np.min(ts)),
                                  float(np.max(ts)))
            r0, md, tm, tlo, thi = out[(name, sname)]
            print(f"      {name:22s} {sname:26s} {r0:10.3f} {md:7.3f} "
                  f"{tm:9.3g} [{tlo:8.3g},{thi:9.3g}]")
    biases = [v[2] for v in out.values()]
    quiet = {k: v for k, v in out.items() if v[1] <= 1.15 * exp_rms}
    loud = {k: v for k, v in out.items() if v[1] > 1.15 * exp_rms}
    record("the bias is present under EVERY insertion; how visible it is "
           "is NOT insertion-robust",
           min(biases) > 1.05,
           "theta_hat/theta runs {:.3g} to {:.3g} across the five generators and "
           "the two\nrecommended schedules, and is above 1 in every cell. So "
           "step 2's FINDING,\nthat fitting the paper's model to a PUB_B-shaped "
           "record biases theta\nupward, survives the choice of insertion. The "
           "particular factor of 1.29\nto 2.03 was measured under `whole` and is "
           "conditional on it.\n\n"
           "WHAT IS NOT INSERTION-ROBUST IS THE ABSENCE OF A WARNING, and the "
           "previous\ndraft of this check claimed it was. Against the expected "
           "{:.3f}:\n"
           "  invisible  {}\n"
           "  partly visible  {}\n"
           "PUB_A, for comparison, gives 2.0 to 2.6 sigma. So the insertions "
           "sort into\ntwo groups: those where the load is a common factor on "
           "the whole marker or\non the suppressible pool alone are invisible, "
           "and those where it moves the\nUNSUPPRESSIBLE part are audible at "
           "about twice the expected residual.".format(
               min(biases), max(biases), exp_rms,
               ", ".join(sorted({k[0].split()[0] for k in quiet})) or "none",
               ", ".join(sorted({k[0].split()[0] for k in loud})) or "none"))
    return out


def q1_insertions(times, loads_c, loads_v):
    head("1. WHERE THE LOAD ENTERS, AND WHETHER theta SURVIVES EACH CHOICE")
    print("""
  PUB_B as published has NO load input. Every row below is our insertion, and
  the previous report's conclusion was stated for the first of them only.
""")
    a = ranks(times, loads_c, "constant load")
    b = ranks(times, loads_v, "varying load")

    free_v = [k for k, v in b.items() if "FREE" in v[4]]
    pinned_v = [k for k, v in b.items() if "FREE" not in v[4]]
    record("theta's identifiability under PUB_B DEPENDS ON THE INSERTION",
           len(free_v) > 0 and len(pinned_v) > 0,
           "At a varying recorded load:\n"
           + "\n".join(f"  {k:22s} {b[k][4]}" for k in b) +
           "\n\nSo the previous report's 'with or without load variation' was "
           "stated for\nthe `whole` insertion and does NOT hold for every "
           "plausible one.")
    return a, b


def q2_intuition(times, loads_v):
    head("2. THE INTUITION: theta and IC50 TRADE OFF, AND WHAT PAYS FOR THE REST")
    print("""
  The claim to check: "theta enters only through theta + IC50, so residual
  enzyme and potency trade off."

  The invariants of the `whole` insertion, derived in k5 step 3a, are
      P = E0 + base q ,  U = P theta + (E0+base) IC50 ,  V = IC50 + theta .
  theta appears in TWO of them, V and U, so the claim as stated is not quite
  right. Walking the fibre shows what happens to each.
""")
    q = Q
    ic50 = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    o = K2.PubB(ic50)
    P = o.e0 + o.base * q
    U = P * TH + (o.e0 + o.base) * ic50
    V = ic50 + TH
    print(f"      truth: E0 {o.e0:.6f}, base {o.base:.4f}, IC50 {ic50:.4f}, "
          f"theta {TH:.1f}")
    print(f"      P = {P:.9f}   U = {U:.4f}   V = {V:.4f}\n")
    print(f"      {'IC50':>10} {'theta':>10} {'IC50+theta':>12} "
          f"{'base':>12} {'E0':>10} {'P':>14} {'U':>12}")
    print("  " + "-" * 86)
    dV, dU, dP = 0.0, 0.0, 0.0
    for frac in (0.01, 0.1, 0.25, 0.5, ic50 / V, 0.75, 0.9, 0.99):
        R = frac * V
        th_r = V - R
        base_r = (U - P * V) / (R * K2.IMAX)
        e0_r = P - q * base_r
        if not (base_r > 0 and e0_r >= 0 and th_r > 0):
            continue
        P_r = e0_r + base_r * q
        U_r = P_r * th_r + (e0_r + base_r) * R
        dV = max(dV, abs((R + th_r) / V - 1.0))
        dU = max(dU, abs(U_r / U - 1.0))
        dP = max(dP, abs(P_r / P - 1.0))
        print(f"      {R:10.4f} {th_r:10.4f} {R + th_r:12.6f} "
              f"{base_r:12.5g} {e0_r:10.5f} {P_r:14.9f} {U_r:12.4f}")
    record("the trade-off is one-for-one in IC50 + theta, "
           "and the marker's scales pay for U and P",
           dV < 1e-12 and dU < 1e-9 and dP < 1e-9,
           f"IC50 + theta is held to {dV:.2e} relative along the whole walk.\n"
           f"U is held to {dU:.2e} and P to {dP:.2e}, but NOT by theta and "
           f"IC50:\nthey are held by `base` and `E0`, which move over three "
           f"orders of\nmagnitude to do it.\n\n"
           f"CORRECTED FORM OF THE INTUITION. theta enters two invariants, not\n"
           f"one. The reason only V = IC50 + theta binds is that the marker's\n"
           f"two free scale constants can always be re-tuned to hold U and P,\n"
           f"so the single surviving constraint is the one-for-one trade\n"
           f"between residual enzyme and potency. Under PUB_A the same V\n"
           f"exists but there is a SECOND invariant W = q theta + IC50 and no\n"
           f"free scale to absorb it, which is exactly why the trade breaks\n"
           f"and theta is pinned there.")


def q3_which_parameter(times, loads_v, ranks_v):
    head("3. PUB_B IS DIRECT-RESPONSE, AND THE DEFICIT IS ITS FREE OFFSET")
    print("""
  FIRST, THE FORM. PUB_B's marker has no state: it is
      3-HP = E0 + base (1 - Imax PCC/(IC50 + PCC)),
  an instantaneous function of the enzyme. The paper's marker is the state of
  a flux balance solved at quasi-steady state, which is also an instantaneous
  map from (Delta, L) to m. So BOTH are direct-response forms and the
  difference between them is the SHAPE of that map, not the dynamics. PUB_A is
  the only one of the three with a marker state, and that is where its fourth
  number comes from.

  SECOND, THE COUNT, side by side at a varying recorded load:

      structure  parameters                       invariants      deficit
      PAPER      theta, K_m, k_e, b          (4)  A, B, C, b  (4)      0
      PUB_B      E0, base, IC50, theta, b    (5)  P, U, V, b  (4)      1

  PUB_B has ONE more parameter and the SAME number of invariants. The extra
  parameter is a free marker scale, which the paper's model does not have at
  all: there the marker IS the state whose balance is written, so its scale is
  fixed by the load's units. PUB_B's marker is a phenomenological readout with
  a free offset E0 and a free amplitude `base`.

  THIRD, WHICH ONE. Removing E0 (setting the unsuppressible floor to zero)
  leaves four parameters and four invariants. The rank table in section 1
  reports that case directly.
""")
    for k in ("whole  (E0+base f)*u", "whole, E0 = 0"):
        r, n, null, names, verdict = ranks_v[k]
        print(f"      {k:22s} rank {r} of {n}: {verdict}")
    ok_e0 = ("FREE" in ranks_v["whole  (E0+base f)*u"][4]
             and "IDENTIFIED" in ranks_v["whole, E0 = 0"][4])

    # the other route: a KNOWN potency
    inc = np.asarray(delta(times), float)
    ic50 = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    o = K2.PubB(ic50)
    names2 = ["E0", "base", "theta", "b"]
    fn2 = lambda p, i, L: ((p[0] + p[1] * _f(ic50, p[2], i))
                           * (L + p[3]) / LB)
    p2 = np.array([o.e0, o.base, TH, 0.0])
    J2 = K5.jac_complex(lambda p: fn2(p, inc, loads_v), p2)
    r2, _, null2 = K5.rank_of(J2)
    print(f"\n      and with IC50 KNOWN instead (potency measured externally):")
    print(f"      {'whole, IC50 fixed':22s} rank {r2} of {len(names2)}: "
          f"{'no null space: IDENTIFIED' if null2 is None else 'deficit'}")
    record("the deficit is PUB_B's free unsuppressible offset E0, "
           "and either E0 = 0 or a known IC50 removes it",
           ok_e0 and null2 is None,
           "So PUB_B's non-identifiability of theta is NOT a property of\n"
           "direct-response forms in general, nor of the saturating shape. It\n"
           "is caused by carrying a free marker offset alongside a free marker\n"
           "amplitude, which is what a phenomenological readout has and a flux\n"
           "balance does not.\n\n"
           "Two independent fixes, both external to the record:\n"
           "  E0 = 0            the floor is known to be zero\n"
           "  IC50 known        the potency is measured, e.g. in vitro\n"
           "and one internal to it: put the load on `base` alone (section 1),\n"
           "which is arguably the most natural reading, since `base` is\n"
           "Attarwala's 'baseline 3-HP amenable to suppression' and E0 is\n"
           "'plasma levels not affected by mRNA-3927'. A larger substrate load\n"
           "should raise the suppressible pool and not the unaffected floor.")


def q4_class():
    head("4. THE THREE-NUMBER BOUND: THE CLASS IT IS PROVED FOR")
    print("""
  The previous report called this "the three-number bound" without naming a
  class, which reads as a general theorem. It is not one. Here is the class it
  is actually proved for, and the verification that all three structures are
  in it.

  THE CLASS. Fix the recorded load. Say a structure is FOUR-MONOMIAL if the
  set of (Delta, m) pairs it can produce is the zero set of a polynomial
  F(Delta, m) = 0 with at most four monomials, whose coefficients are functions
  of the structural parameters.

  THE BOUND, for that class. F is defined only up to an overall scale, so the
  record determines its four coefficients as a point in projective 3-space,
  that is THREE numbers, and nothing else. Any parameter vector of dimension
  n > 3 therefore has an (n-3)-dimensional fibre at a fixed load.

  This is elementary. The content is entirely in MEMBERSHIP, which is checked
  below and is claimed for these three structures ONLY.
""")
    D, m, th, Km, ke, L, b, c, E0, base, q, kout = sp.symbols(
        "Delta m theta K_m k_e L b c E0 base q k_out", positive=True)
    rows, ok = [], True

    # each relation is written as a rational expression equal to zero; the
    # numerator after cancellation is the polynomial whose monomials are
    # counted. Expanding first leaves 1/(K_m+m) factors behind, which sympy
    # rightly refuses to treat as a polynomial in m.
    for tag, expr in (
        ("PAPER, quasi-steady balance",
         (th + D) * m / (Km + m) + ke * m - (L + b)),
        ("PUB_B, algebraic marker",
         m - (E0 + base * (c + q * (th + D)) / (c + th + D))),
        ("PUB_A, at steady state",
         kout * m - (L + b) * (c + q * (th + D)) / (c + th + D)),
    ):
        num, _ = sp.fraction(sp.cancel(sp.together(expr)))
        p = sp.Poly(sp.expand(num), D, m)
        mons = sorted(p.monoms())
        rows.append((tag, mons, len(mons)))
        ok = ok and len(mons) <= 4
        names = {(0, 0): "1", (0, 1): "m", (0, 2): "m^2",
                 (1, 0): "Delta", (1, 1): "Delta*m", (1, 2): "Delta*m^2"}
        pretty = ", ".join(names.get(t, str(t)) for t in mons)
        print(f"      {tag:30s} {len(mons)} monomials: {pretty}")

    print("""
  WHERE THE CLASS STOPS, and this is the useful part. PUB_A's FULL model is
  NOT four-monomial and is not in the class at all, because its marker is a
  state: at a fixed load the pairs (Delta(t), m(t)) do not lie on any algebraic
  curve, since m depends on the history of Delta rather than on its current
  value. That is precisely where PUB_A's fourth number comes from, and it is
  why the class boundary and the mechanism are the same line.

  So the honest statement is: within the four-monomial class the bound is
  elementary, all three structures examined are in it, PUB_A leaves it as soon
  as its transient is used, and nothing is claimed about structures outside it.
""")
    record("all three structures examined are four-monomial at a fixed load",
           ok,
           "\n".join(f"{t}: {n} monomials" for t, _, n in rows) +
           "\nClaimed for these three and not beyond them.")


def main():
    times = G.SCHEDULES["phase-spread"](24)
    loads_c = np.full(times.size, LB)
    rng = np.random.default_rng(0)
    s2 = np.log(1 + 0.10 ** 2)
    loads_v = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))

    _, ranks_v = q1_insertions(times, loads_c, loads_v)
    bias_by_insertion()
    q2_intuition(times, loads_v)
    q3_which_parameter(times, loads_v, ranks_v)
    q4_class()

    n = sum(1 for _, ok, _ in RESULTS if ok)
    print("\n" + "=" * 78)
    print(f"  {n} of {len(RESULTS)} checks pass.")
    for t, ok, _ in RESULTS:
        if not ok:
            print(f"    FAILED: {t}")
    print("=" * 78)


if __name__ == "__main__":
    main()
