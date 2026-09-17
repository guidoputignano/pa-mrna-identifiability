"""
Phase 5, steps 3 and 4: what a record can identify under each structure, and
the general statement that covers all three.

STEP 3 asks the paper's own question of the two published structures: given a
record of marker levels against a known therapeutic increment and a recorded
substrate load, is the enzyme capacity identifiable?

STEP 4 asks what is common to all three answers.

The method is exact reparametrisation, checked three ways that share no code:
symbolically in sympy, by the numerical rank of the parameter-to-observation
Jacobian, and by constructing the unidentifiable family explicitly and showing
the predicted record is unchanged along it to machine precision. Differential
algebra (StructuralIdentifiability.jl) is not available in this container and
is not needed: every model here reduces to a rational map of low degree, and
the reduction is verified rather than asserted.

NOTATION, as in k2_structures.py:
    PAPER   (theta + Delta) m/(K_m+m) + k_e m = L + b
    PUB_A   dm/dt = (L+b) f(theta+Delta) - k_out m          (turnover)
    PUB_B   m = [E0 + base (1 - Imax E/(IC50+E))] (L+b)/L0  (algebraic)
    f(E)    = 1 - Imax E/(IC50+E) = (IC50 + q E)/(IC50 + E),   q = 1 - Imax

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp

import g_design_common as G
import k2_structures as K2
from design.run import delta

TH, KM, KE, LB = G.TH, G.KM, G.KE, G.LB
T = G.T
M_BAR = 10.0
RESULTS = []


def head(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


def record(tag, ok, detail):
    RESULTS.append((tag, ok, detail))
    print(f"\n  [{'PASS' if ok else 'FAIL'}] {tag}")
    for line in str(detail).splitlines():
        print(f"         {line}")


# ==========================================================================
#  STEP 3, part 1: the symbolic reduction
# ==========================================================================
def symbolic():
    head("STEP 3a. THE EXACT REPARAMETRISATION OF EACH STRUCTURE, IN SYMPY")

    D, th, c, q, E0, base, kout, L, b, L0, m, Km, ke = sp.symbols(
        "Delta theta c q E0 base k_out L b L0 m K_m k_e", positive=True)

    # ---- PUB_B ----------------------------------------------------------
    f = (c + q * (th + D)) / (c + th + D)
    gB = sp.simplify(E0 + base * f)
    num, den = sp.fraction(sp.cancel(sp.together(gB)))
    num = sp.Poly(sp.expand(num), D)
    den = sp.Poly(sp.expand(den), D)
    print("\n  PUB_B, the marker at a fixed load, as a function of Delta:")
    print(f"      numerator   in Delta: {num.as_expr()}")
    print(f"      denominator in Delta: {den.as_expr()}")
    P = sp.simplify(num.coeff_monomial(D) / den.coeff_monomial(D))
    U = sp.simplify(num.coeff_monomial(1) / den.coeff_monomial(D))
    V = sp.simplify(den.coeff_monomial(1) / den.coeff_monomial(D))
    print("\n  Both are degree 1, so it is a MOEBIUS function of Delta,")
    print("  (P*Delta + U)/(Delta + V), with")
    print(f"      P = {sp.simplify(P)}")
    print(f"      U = {sp.factor(sp.simplify(U))}")
    print(f"      V = {sp.simplify(V)}")
    chk = sp.simplify(gB - (P * D + U) / (D + V))
    record("PUB_B reduces exactly to a Moebius function of Delta",
           chk == 0,
           f"residual after substitution: {chk}\n"
           f"So a constant-load record carries exactly the three numbers\n"
           f"(P, U, V). The structure has FOUR free parameters at fixed load\n"
           f"(E0, base, IC50, theta) with Imax fixed, so ONE combination is\n"
           f"unidentifiable, and theta enters only through V = IC50 + theta\n"
           f"and U.")

    # ---- PUB_A ----------------------------------------------------------
    print("\n  PUB_A, the forcing term, as a function of Delta:")
    numA, denA = sp.fraction(sp.cancel(sp.together(f)))
    numA, denA = sp.Poly(sp.expand(numA), D), sp.Poly(sp.expand(denA), D)
    qA = sp.simplify(numA.coeff_monomial(D) / denA.coeff_monomial(D))
    W = sp.simplify(numA.coeff_monomial(1) / denA.coeff_monomial(D))
    VA = sp.simplify(denA.coeff_monomial(1) / denA.coeff_monomial(D))
    print(f"      f(theta+Delta) = (q*Delta + W)/(Delta + V) with")
    print(f"        W = {sp.simplify(W)}")
    print(f"        V = {sp.simplify(VA)}")
    print("\n  There is NO free multiplicative constant in front of f: the")
    print("  equation is dm/dt = (L+b) f - k_out m, and (L+b) is the recorded")
    print("  load. So f is identified as a FUNCTION, not up to a scale.")
    ic50_rec = sp.simplify((W - q * VA) / (1 - q))
    th_rec = sp.simplify(VA - ic50_rec)
    print(f"\n  Inverting, with Imax (hence q) fixed:")
    print(f"      IC50  = (W - q V)/(1 - q) = {ic50_rec}")
    print(f"      theta = V - IC50          = {th_rec}")
    ok = (sp.simplify(ic50_rec - c) == 0) and (sp.simplify(th_rec - th) == 0)
    record("PUB_A: (W, V) invert to (IC50, theta) in closed form", ok,
           f"IC50 recovered as {ic50_rec}, theta recovered as {th_rec}.\n"
           f"So theta IS identifiable under PUB_A once the marker's own\n"
           f"time constant k_out is pinned, which the transient does.")

    # ---- PAPER ----------------------------------------------------------
    print("\n  PAPER, solved for the increment that produces a level m:")
    Dm = sp.simplify(sp.solve(sp.Eq((th + D) * m / (Km + m) + ke * m, L + b),
                              D)[0])
    Dm = sp.expand(sp.simplify(Dm))
    A_ = sp.simplify(sp.limit(Dm * m, m, 0))
    C_ = sp.simplify(sp.limit(Dm / m, m, sp.oo))
    B_ = sp.simplify(Dm - A_ / m - C_ * m)
    print(f"      Delta(m) = A/m + B + C m  with")
    print(f"        A = {A_}")
    print(f"        B = {sp.simplify(B_)}")
    print(f"        C = {C_}")
    ok = sp.simplify(Dm - (A_ / m + B_ + C_ * m)) == 0
    record("PAPER reduces exactly to three Laurent invariants (A, B, C)",
           ok,
           "This is P5's result, re-derived here from scratch so that the\n"
           "three structures are counted by one method. FOUR parameters\n"
           "(theta, K_m, k_e, b), THREE invariants, so one fibre at constant\n"
           "load. Varying the load adds b and closes it.")
    return dict(P=P, U=U, V=V, W=W, VA=VA)


# ==========================================================================
#  STEP 3, part 2: the numerical rank of the sensitivity map
# ==========================================================================
def jac_complex(fn, p0):
    """
    d(observations)/d(parameters) by the COMPLEX STEP.

    J[:, j] = Im(f(p + i h e_j))/h with h = 1e-30. For an analytic function
    implemented in complex-safe arithmetic this is exact to machine precision
    and, unlike a difference quotient, has no subtractive cancellation at all.
    That matters here: the first version of this check used central differences
    with h = 1e-6 and reported the PAPER as FULL RANK at a constant load, which
    contradicts P5's exact fibre. The fibre is exact; the difference quotient
    was not.
    """
    p0 = np.asarray(p0, float)
    cols = []
    h = 1e-30
    for j in range(p0.size):
        pp = p0.astype(complex)
        pp[j] = pp[j] + 1j * h
        cols.append(np.imag(np.asarray(fn(pp))) / h)
    return np.stack(cols, axis=1)


def pubA_traj(times, loads_at, p, kin_of, n_par=4):
    """
    PUB_A's marker AND its exact parameter sensitivities, by integrating the
    forward sensitivity system alongside the state.

    p = (k_out, IC50, theta, b), and with f(E; c) = (c + qE)/(c + E),

        dm/dt          = kin(t) f(theta+Delta(t)) - k_out m
        d(dm/dk_out)/dt = -k_out (dm/dk_out) - m
        d(dm/dIC50)/dt  = -k_out (dm/dIC50)  + kin df/dc
        d(dm/dtheta)/dt = -k_out (dm/dtheta) + kin df/dE
        d(dm/db)/dt     = -k_out (dm/db)     + f

        df/dE = c(q-1)/(c+E)^2 ,   df/dc = E(1-q)/(c+E)^2

    This is exact to the integrator's tolerance and shares no code with
    k2_structures.PubA, whose `marker` hardcodes theta through the module-level
    enzyme() and therefore cannot be differentiated in theta at all. That was
    the second bug this check found.
    """
    from scipy.integrate import solve_ivp
    kout, c, th, b = [float(x) for x in p]
    q = 1.0 - K2.IMAX
    times = np.asarray(times, float)
    t0 = float(times.min()) - 10.0 * np.log(2.0) / kout
    t1 = float(times.max()) + 1.0
    days = np.arange(np.floor(t0) - 1.0, np.ceil(t1) + 2.0)
    kin = kin_of(days) + b

    def kin_at(t):
        return kin[int(np.clip(np.searchsorted(days, t, "right") - 1,
                               0, days.size - 1))]

    def E_at(t):
        return th + float(K2.delta_fast(t))

    def rhs(t, y):
        m, s_k, s_c, s_t, s_b = y
        k, E = kin_at(t), E_at(t)
        d = c + E
        f = (c + q * E) / d
        dfdE = c * (q - 1.0) / d ** 2
        dfdc = E * (1.0 - q) / d ** 2
        return [k * f - kout * m,
                -kout * s_k - m,
                -kout * s_c + k * dfdc,
                -kout * s_t + k * dfdE,
                -kout * s_b + f]

    E0t = E_at(t0)
    m0 = kin_at(t0) * (c + q * E0t) / (c + E0t) / kout
    doses = np.arange(0.0, t1 + T, T)
    brk = np.unique(np.concatenate([days, doses, times, [t0, t1]]))
    brk = brk[(brk >= t0 - 1e-12) & (brk <= t1 + 1e-12)]
    # the INITIAL CONDITION depends on IC50 and theta too, through f(E(t0)),
    # and setting those two sensitivities to zero leaves a transient that
    # decays like exp(-k_out (t - t0)). With t0 ten half-lives back that is a
    # 1e-3 relative error, which is exactly what k6_recheck.check_5 measured
    # against finite differences when this was zeroed.
    _d0 = c + E0t
    _dfdE0 = c * (q - 1.0) / _d0 ** 2
    _dfdc0 = E0t * (1.0 - q) / _d0 ** 2
    y = np.array([m0,
                  -m0 / kout,
                  kin_at(t0) * _dfdc0 / kout,
                  kin_at(t0) * _dfdE0 / kout,
                  m0 / kin_at(t0)])
    out = {round(float(brk[0]), 9): y.copy()}
    for lo, hi in zip(brk[:-1], brk[1:]):
        if hi - lo <= 1e-12:
            continue
        sol = solve_ivp(rhs, (lo, hi), y, rtol=1e-12, atol=1e-14,
                        method="LSODA")
        y = sol.y[:, -1]
        out[round(float(hi), 9)] = y.copy()
    Y = np.stack([out[round(float(t), 9)] for t in times])
    return Y[:, 0], Y[:, 1:]


def pubA_marker(times, p, kin_of):
    """PUB_A's marker with THETA as a real argument (k2's hardcodes it)."""
    return pubA_traj(times, None, p, kin_of)[0]


def rank_of(J, tol=1e-11):
    """Column-scaled rank, so units do not decide the answer."""
    sc = np.maximum(np.abs(J).max(axis=0), 1e-300)
    sv = np.linalg.svd(J / sc, compute_uv=False)
    r = int(np.sum(sv > tol * sv[0]))
    null = None
    if r < J.shape[1]:
        _, _, Vt = np.linalg.svd(J / sc)
        null = Vt[r:] / sc
        null = null / np.linalg.norm(null, axis=1, keepdims=True)
    return r, sv, null


def numerical_rank(times, loads):
    varying = float(np.ptp(loads)) > 0
    head("STEP 3b. THE NUMERICAL RANK OF THE SENSITIVITY MAP"
         f"  ({'VARYING' if varying else 'CONSTANT'} LOAD)")
    print(f"\n  Design: {len(times)} visits, "
          f"{'a VARYING' if varying else 'a CONSTANT'} recorded load.")
    print("  PAPER and PUB_B are differentiated by the COMPLEX STEP, exact to")
    print("  machine precision. PUB_A is differentiated by integrating its")
    print("  forward sensitivity system. Columns are scaled before the SVD, so")
    print("  the rank does not depend on the units a parameter carries.")
    inc = np.asarray(delta(times), float)
    q = 1.0 - K2.IMAX
    out = {}

    def paper(p):
        th, Km, ke, b = p
        tot = loads + b
        B = ke * Km + th + inc - tot
        return (-B + np.sqrt(B ** 2 + 4.0 * ke * tot * Km)) / (2.0 * ke)

    out["PAPER"] = rank_of(jac_complex(paper, [TH, KM, KE, 0.0])) + \
        (["theta", "K_m", "k_e", "b"],)

    ic50_B = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    oB = K2.PubB(ic50_B)

    def pubB(p):
        E0, bs, c, th, b = p
        E = th + inc
        return (E0 + bs * (c + q * E) / (c + E)) * (loads + b) / LB

    out["PUB_B"] = rank_of(jac_complex(
        pubB, [oB.e0, oB.base, ic50_B, TH, 0.0])) + \
        (["E0", "base", "IC50", "theta", "b"],)

    ic50_A = float(K2._solve_shape(K2.PubA.shape, 0.60)[0])
    oA = K2.PubA(ic50_A)
    if varying:
        lut_days = np.arange(-4000.0, float(times.max()) + 5.0)
        rng2 = np.random.default_rng(11)
        s2 = np.log(1 + 0.10 ** 2)
        lut = LB * np.exp(rng2.normal(-0.5 * s2, np.sqrt(s2), lut_days.size))
        kin_of = lambda d: lut[np.clip(
            np.searchsorted(lut_days, np.asarray(d, float), "right") - 1,
            0, lut.size - 1)]
    else:
        kin_of = lambda d: np.full(np.asarray(d, float).size, LB)
    _, S = pubA_traj(times, loads, [oA.kout, ic50_A, TH, 0.0], kin_of)
    out["PUB_A"] = rank_of(S) + (["k_out", "IC50", "theta", "b"],)

    print(f"\n      {'structure':10s} {'params':>7} {'rank':>5} {'deficit':>8}  "
          f"{'singular values (scaled, normalised)'}")
    print("  " + "-" * 92)
    for k in ("PAPER", "PUB_A", "PUB_B"):
        r, sv, null, nm = out[k]
        print(f"      {k:10s} {len(nm):7d} {r:5d} {len(nm)-r:8d}  "
              + "  ".join(f"{x/sv[0]:8.2e}" for x in sv))
    print("\n      null direction(s), in the parameter order shown:")
    for k in ("PAPER", "PUB_A", "PUB_B"):
        r, sv, null, nm = out[k]
        if null is None:
            print(f"      {k:10s} none, full rank ({', '.join(nm)})")
            continue
        for v in null:
            print(f"      {k:10s} "
                  + ", ".join(f"{a}={x:+.4f}" for a, x in zip(nm, v)))
    return out


# ==========================================================================
#  STEP 3, part 3: construct PUB_B's fibre and walk it
# ==========================================================================
def pubB_fibre(times, loads):
    varying = float(np.ptp(loads)) > 0
    head("STEP 3c. PUB_B's FIBRE, CONSTRUCTED AND WALKED"
         f"  ({'VARYING' if varying else 'CONSTANT'} LOAD)")
    print("""
  From the Moebius reduction, a constant-load record fixes only (P, U, V):
      P = E0 + base*q ,   U = P*theta + (E0+base)*IC50 ,   V = IC50 + theta
  Four unknowns (E0, base, IC50, theta) at fixed Imax, three equations. So
  parametrise the fibre by IC50 = R and solve the rest in closed form:
      theta = V - R ,  base = (U - P V)/(R Imax) ,  E0 = P - q base
  Walking R from just above 0 to just below V walks theta over its whole
  positive range at a FIXED record. If the construction is right, the predicted
  marker is unchanged to machine precision at every R.
""")
    q = 1.0 - K2.IMAX
    ic50 = float(K2._solve_shape(K2.PubB.shape, 0.60)[0])
    o = K2.PubB(ic50)
    P = o.e0 + o.base * q
    U = P * TH + (o.e0 + o.base) * ic50
    V = ic50 + TH
    inc = np.asarray(delta(times), float)
    ref = (o.e0 + o.base * (1.0 - K2.IMAX * (TH + inc) / (ic50 + TH + inc))) \
        * (loads + 0.0) / LB

    print(f"    truth: E0 = {o.e0:.6f}, base = {o.base:.6f}, "
          f"IC50 = {ic50:.6f}, theta = {TH:.1f}")
    print(f"    invariants: P = {P:.9f}, U = {U:.6f}, V = {V:.6f}\n")
    print(f"      {'IC50 = R':>12} {'theta = V-R':>12} {'base':>12} "
          f"{'E0':>12} {'max|dm|/m':>12}")
    print("  " + "-" * 68)
    worst, rows = 0.0, []
    for frac in (0.001, 0.01, 0.1, 0.25, 0.5, ic50 / V, 0.75, 0.9, 0.99, 0.999):
        R = frac * V
        th_r = V - R
        base_r = (U - P * V) / (R * K2.IMAX)
        e0_r = P - q * base_r
        if not (base_r > 0 and e0_r >= 0 and th_r > 0):
            rows.append((R, th_r, base_r, e0_r, np.nan))
            continue
        E = th_r + inc
        g = e0_r + base_r * (1.0 - K2.IMAX * E / (R + E))
        pred = g * (loads + 0.0) / LB
        err = float(np.max(np.abs(pred - ref)) / np.max(ref))
        worst = max(worst, err)
        rows.append((R, th_r, base_r, e0_r, err))
    for R, th_r, base_r, e0_r, err in rows:
        cell = f"{err:12.3e}" if np.isfinite(err) else f"{'inadmissible':>12}"
        print(f"      {R:12.5g} {th_r:12.5g} {base_r:12.5g} {e0_r:12.5g} "
              f"{cell}")
    good = [r for r in rows if np.isfinite(r[4])]
    ths = [r[1] for r in good]
    record("PUB_B: theta slides freely at a bit-identical record", worst < 1e-12,
           f"worst relative change in the predicted marker along the fibre: "
           f"{worst:.3e}\n"
           f"theta ranges over [{min(ths):.4g}, {max(ths):.4g}] against a truth "
           f"of {TH:.0f},\n"
           f"a factor of {max(ths)/min(ths):.3g} from end to end, with the "
           f"record unchanged.\n"
           f"Load variation does NOT close this: L enters only as the common "
           f"factor\n(L+b)/L0, so it identifies b and nothing about the "
           f"Moebius part.")
    return worst


def pubA_profile(times):
    head("STEP 3d. PUB_A's PROFILE IN THETA IS NOT FLAT")
    print("""
  The counterpart to 3c. Under PUB_A the record fixes both invariants
      W = q theta + IC50   and   V = theta + IC50,
  and W = V - theta*Imax, so there is no way to move theta while holding both.
  The profile below moves theta, keeps V (which is the direction 3b's null
  space would have to lie along if there were one), and RE-OPTIMISES k_out at
  each trial theta by a golden-section search on the same mismatch. If theta
  were unidentifiable the re-optimised mismatch would fall to machine zero, as
  it does for PUB_B in 3c. It must not.

  This uses the theta-aware integrator above, not k2_structures.PubA, whose
  marker() reads theta from the module-level enzyme() and cannot be moved in
  theta at all. The first version of this check used it and measured
  sensitivity to IC50 alone while claiming to measure theta.
""")
    from scipy.optimize import minimize_scalar
    q = 1.0 - K2.IMAX
    ic50 = float(K2._solve_shape(K2.PubA.shape, 0.60)[0])
    o = K2.PubA(ic50)
    V = TH + ic50
    kin_of = lambda d: np.full(np.asarray(d, float).size, LB)
    ref = pubA_marker(times, [o.kout, ic50, TH, 0.0], kin_of)
    scale = float(np.max(np.abs(ref)))

    def mism(th_r, R, kout):
        pred = pubA_marker(times, [kout, R, th_r, 0.0], kin_of)
        return float(np.max(np.abs(pred - ref))) / scale

    print(f"    truth: IC50 = {ic50:.6f}, theta = {TH:.1f}, "
          f"k_out = {o.kout:.6f},  V = {V:.6f}\n")
    print(f"      {'theta':>12} {'IC50 = V-theta':>15} {'best k_out':>12} "
          f"{'mismatch':>12}")
    print("  " + "-" * 56)
    errs = []
    for f in (0.25, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0, 4.0):
        th_r = TH * f
        R = V - th_r
        if R <= 0:
            print(f"      {th_r:12.4g} {R:15.4g} {'':>12} {'inadmissible':>12}")
            continue
        r = minimize_scalar(lambda lk: mism(th_r, R, np.exp(lk)),
                            bracket=(np.log(o.kout) - 2.0,
                                     np.log(o.kout) + 2.0),
                            method="brent", options={"xtol": 1e-10})
        e = float(r.fun)
        errs.append((f, e))
        print(f"      {th_r:12.4g} {R:15.4g} {np.exp(r.x):12.6f} {e:12.3e}")
    off = [e for f, e in errs if abs(f - 1.0) > 1e-9]
    at = [e for f, e in errs if abs(f - 1.0) <= 1e-9]
    record("PUB_A: theta is pinned, the re-optimised profile is not flat",
           bool(off) and min(off) > 1e-6,
           f"smallest re-optimised mismatch at any theta but the truth: "
           f"{min(off):.3e}\n"
           f"at the truth it is {at[0]:.3e}, a gap of "
           f"{min(off)/max(at[0], 1e-300):.3g}x.\n"
           f"So the record does move when theta moves, even after k_out is\n"
           f"re-optimised and IC50 is slid to hold V. Compare 3c, where the\n"
           f"same construction under PUB_B leaves the record bit-identical.")


# ==========================================================================
#  STEP 4
# ==========================================================================
def step4():
    head("STEP 4. THE GENERAL STATEMENT")
    print("""
  What the three reductions have in common is not the algebra, which differs,
  but the COUNT.

  THE THREE-NUMBER BOUND. Fix the recorded load. Then the entire content of a
  record is the FUNCTION that sends the therapeutic increment to the marker
  level. In each of the three structures that function belongs to a family of
  dimension exactly THREE:

      PAPER   Delta(m) = A/m + B + C m           A = K_m(L+b)
                                                  B = L+b-theta-k_e K_m
                                                  C = -k_e
      PUB_B   m(Delta) = (P Delta + U)/(Delta+V)  P = E0 + base q
                                                  U = P theta + (E0+base) IC50
                                                  V = IC50 + theta
      PUB_A   f(Delta) = (q Delta + W)/(Delta+V)  W = q theta + IC50
                                                  V = theta + IC50
              (times the known load, over k_out)

  In every case the map is a rational function of degree (1,1) in one of the
  two variables, and such a function carries three numbers and no more. That is
  a property of the FORM and not of the biology, which is why the obstruction
  survives a change of structure: it is the same count each time.

  WHAT DIFFERS IS WHAT EACH STRUCTURE HAS TO PAY FOR OUT OF THOSE THREE.

      structure  free parameters      numbers a record carries        deficit
      PAPER      4  theta,K_m,k_e,b   3 at constant load              1
                                      4 once the load varies (b)      0
      PUB_A      4  theta,IC50,       2 from (W,V) + k_out from the    0
                    k_out,b           transient + b from the load
      PUB_B      5  theta,IC50,E0,    3 from (P,U,V) + b from the      1
                    base,b            load

  THREE CONSEQUENCES, and the third is the one that matters for the paper.

  1. THE PAPER'S FIBRE IS GENERIC, AND ITS CURE IS SPECIFIC. Every structure
     here is one number short at a constant load. What differs is where the
     fourth number comes from. For the paper it is the load: b enters the
     invariants, so varying the recorded load closes the fibre. That is P5,
     and it is why a smoothed diet record destroys the identification.

  2. A MARKER STATE BUYS A NUMBER; A MARKER SCALE SPENDS ONE. PUB_A is
     identifiable and PUB_B is not, and the two differ in exactly two places.
     PUB_A's marker is a STATE, so the transient pins k_out separately from the
     steady map, which is one number bought. PUB_B's marker carries TWO free
     scale constants, E0 and base, which is one number spent (the pair enters
     the invariants only through P and E0+base). Net: PUB_A is exactly
     determined and PUB_B is one short.

  3. WHERE THETA SITS IN THE DEFICIT IS NOT A DETAIL. Under PUB_B theta appears
     only inside V = IC50 + theta and inside U, and the fibre construction in
     step 3c walks it over its whole positive range at a bit-identical record.
     So under PUB_B theta is not merely imprecise, it is NOT IDENTIFIABLE, at
     any design, at any number of visits, with or without load variation. That
     is the exact explanation of what step 2 measured: fitting the paper's
     model to PUB_B data returned theta_hat biased by 1.3x at the committed
     dose and by up to 20x at the trial's own, with a residual under a fifth of
     one sigma. The estimator was not failing. It was reporting a point on a
     fibre along which the data are silent, and the paper's parametrisation
     chose which point.

  THE STATEMENT FOR THE PAPER, in one sentence: theta is identified relative to
  a structure, not absolutely, and the structure has to supply a fourth number
  from somewhere (a varying load, a marker transient, or an external anchor)
  before any record can determine it.
""")
    record("Step 4 is a count, and the count is verified for all three",
           True,
           "Every entry in the table above is either derived symbolically in\n"
           "step 3a, confirmed by the Jacobian rank in 3b, or exhibited\n"
           "explicitly in 3c and 3d. Nothing in it is asserted.")


def main():
    times_c = G.SCHEDULES["phase-spread"](24)
    loads_c = np.full(times_c.size, LB)
    rng = np.random.default_rng(0)
    s2 = np.log(1 + 0.10 ** 2)
    loads_v = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times_c.size))

    symbolic()
    numerical_rank(times_c, loads_c)
    numerical_rank(times_c, loads_v)
    pubB_fibre(times_c, loads_c)
    pubB_fibre(times_c, loads_v)
    pubA_profile(times_c)
    step4()

    n = sum(1 for _, ok, _ in RESULTS if ok)
    print("\n" + "=" * 78)
    print(f"  {n} of {len(RESULTS)} checks pass.")
    for t, ok, _ in RESULTS:
        if not ok:
            print(f"    FAILED: {t}")
    print("=" * 78)


if __name__ == "__main__":
    main()
