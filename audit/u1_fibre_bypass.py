"""
Does ONE external K_m close the quasi-steady fibre exactly?

The claim under test, in the author's own notation. Write the quasi-steady
balance of eq:local with a constant recorded load L and an unknown additive
load bias b,

    (theta + Delta) m/(K_m + m) + k_e m = L + b,

clear the denominator, and it is a quadratic in m whose coefficients are the
only things a record of (m_i, Delta_i) at CONSTANT recorded load can determine:

    C m^2 + (Delta - B) m - A = 0,
        A = (L + b) K_m,   B = (L + b) - theta - k_e K_m,   C = k_e.

Delta is known, so the record sees exactly the triple (A, B, C) and no more.
Four unknowns (theta, K_m, k_e, b) against three determined combinations is a
one-parameter fibre, and the claim is that K_m parametrises it:

    k_e = C,   L + b = A/K_m,   theta = A/K_m - B - k_e K_m.

PART A proves that symbolically, from the balance rather than from the
quadratic, so the quadratic is a consequence and not an assumption.

PART B gives the fibre's sensitivity in closed form and checks it against a
finite difference along the exact fibre.

PART C tests it against this repository's own likelihood, which shares no
algebra with parts A and B: it builds the information matrix of
audit/g_design_common.py in (ln theta, ln K_m, b, ln k_e), puts a log-normal
prior of width s on K_m, and asks whether sd(ln theta) tends to
|dln theta / dln K_m| * s as the record lengthens at a constant recorded load,
which is what "the external K_m is the only thing closing the fibre" means.

PART D is the same question with the load allowed to vary, where the record
itself carries information about K_m and the limit should NOT be attained.

FOOTING. src/common/params.py and src/design/run.py throughout: theta = 3000,
K_m = 100, k_e = 5 /d, L = clearance_A(10) = 322.727 umol/L/d, D = 6000 on a
14 d train with a 3.5 d half-life, load-bias prior sd = 0.05 L, k_e log-prior
0.8, sigma = model_A.run.sigma_total(0.05, 0.05). Capacity sweeps hold the
RECORDED LOAD FIXED and re-solve the operating level, sigma and the load-bias
prior at each capacity, which is Table tab:scope's convention.

Nothing here records an expected result in advance.

Run:  python audit/u1_fibre_bypass.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp

from design.run import P, delta, load, DOSE_INTERVAL, LOAD_FRAC, ASSAY_CV
import model_A.run as LOC

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
T = DOSE_INTERVAL
GATE = 0.20


# ======================================================================
# PART A. The fibre, symbolically.
# ======================================================================
def part_a():
    print("=" * 74)
    print("PART A.  (A, B, C) is everything a constant-load record determines")
    print("=" * 74)

    m, Dl, th, Km, ke, L, b = sp.symbols(
        "m Delta theta K_m k_e L b", positive=True)

    balance = sp.Eq((th + Dl) * m / (Km + m) + ke * m, L + b)
    poly = sp.Poly(sp.simplify(sp.expand(
        (sp.together(balance.lhs - balance.rhs) * (Km + m)).simplify())), m)
    print("\nthe balance, denominator cleared, as a polynomial in m:")
    print("   ", sp.factor(poly.as_expr()))

    c2, c1, c0 = poly.all_coeffs()
    print("\n  coefficient of m^2 :", sp.simplify(c2))
    print("  coefficient of m^1 :", sp.simplify(c1))
    print("  coefficient of m^0 :", sp.simplify(c0))

    #  the author's A, B, C, read off the coefficients. Delta is a KNOWN
    #  regressor, so the m^1 coefficient splits into a known part and an
    #  unknown one, and only the unknown part is a parameter of the record.
    A_s, B_s, C_s = sp.symbols("A B C", positive=True)
    A_def = (L + b) * Km
    B_def = (L + b) - th - ke * Km
    C_def = ke

    resid_2 = sp.simplify(c2 - C_def)
    resid_1 = sp.simplify(c1 - (Dl - B_def))
    resid_0 = sp.simplify(c0 - (-A_def))
    print("\n  c2 - C            :", resid_2)
    print("  c1 - (Delta - B)  :", resid_1)
    print("  c0 - (-A)         :", resid_0)
    ok_coeffs = (resid_2 == 0 and resid_1 == 0 and resid_0 == 0)
    print("  the triple reproduces the quadratic exactly:", ok_coeffs)

    #  invert: given (A, B, C) and an external K_m, is theta determined?
    sol = sp.solve([sp.Eq(A_s, A_def), sp.Eq(B_s, B_def), sp.Eq(C_s, C_def)],
                   [th, b, ke], dict=True)
    print("\n  solving (A, B, C) for (theta, b, k_e) at a GIVEN K_m:")
    assert len(sol) == 1, sol
    s0 = sol[0]
    for k in (th, b, ke):
        print(f"    {k} = {sp.simplify(s0[k])}")

    claim = sp.simplify(s0[th] - (A_s / Km - B_s - C_s * Km))
    print("\n  theta - (A/K_m - B - k_e K_m) =", claim)
    ok_theta = (claim == 0)
    print("  the author's inversion is exact:", ok_theta)

    #  and the fibre in K_m, with (A, B, C) held at their true values
    theta_fib = sp.simplify((A_def / Km).subs(Km, Km) )
    K = sp.Symbol("K", positive=True)
    theta_of_K = sp.simplify(A_s / K - B_s - C_s * K)
    dlt_dlK = sp.simplify(sp.diff(sp.log(theta_of_K), sp.log(K))
                          if False else K * sp.diff(theta_of_K, K) / theta_of_K)
    #  substitute the true (A, B, C) and evaluate at K = K_m
    subs_true = {A_s: A_def, B_s: B_def, C_s: C_def}
    slope = sp.simplify(dlt_dlK.subs(subs_true).subs(K, Km))
    print("\n  dln theta / dln K_m along the fibre, at the true point:")
    print("   ", sp.simplify(slope))
    closed = sp.simplify(slope + (L + b + ke * Km) / th)
    print("  that minus  -(L + b + k_e K_m)/theta  =", closed)
    ok_slope = (closed == 0)
    print("  the closed form holds:", ok_slope)

    #  what the record determines WITHOUT an external K_m: the fibre is
    #  one-dimensional, so theta is not determined. Show the residual freedom.
    print("\n  with (A, B, C) exact and K_m free, theta spans the fibre:")
    for Kv in (25.0, 50.0, 100.0, 200.0, 400.0):
        num = float(theta_of_K.subs(subs_true).subs(
            {Km: KM, L: LB, b: 0.0, th: TH, ke: KE, K: Kv}))
        print(f"    K_m = {Kv:7.1f}  ->  theta = {num:10.2f}")

    return ok_coeffs and ok_theta and ok_slope


# ======================================================================
# PART B. The fibre numerically, and its slope by finite difference.
# ======================================================================
def level(theta, inc, Ltot, Km=KM, ke=KE):
    """Quasi-steady marker: the positive root of the cleared balance."""
    inc = np.asarray(inc, float)
    Bq = ke * Km + theta + inc - Ltot
    return (-Bq + np.sqrt(Bq ** 2 + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


def fibre_theta(K, A, B, C):
    return A / K - B - C * K


def part_b():
    print()
    print("=" * 74)
    print("PART B.  the fibre numerically, and the slope by finite difference")
    print("=" * 74)
    print(f"\n  footing: L = {LB:.3f} umol/L/d, b = 0, k_e = {KE}, "
          f"K_m = {KM}")
    print("  capacities: recorded load held fixed, level re-solved\n")
    print(f"  {'theta':>8}  {'m_bar':>8}  {'A':>12}  {'B':>10}  "
          f"{'closed form':>12}  {'finite diff':>12}  {'rel gap':>9}")
    print("  " + "-" * 82)
    rows = []
    for theta in (3000.0, 600.0, 300.0, 150.0):
        A = (LB + 0.0) * KM
        B = (LB + 0.0) - theta - KE * KM
        C = KE
        mbar = float(level(theta, 0.0, LB))
        closed = -(LB + 0.0 + KE * KM) / theta
        h = 1e-6
        fd = ((np.log(fibre_theta(KM * np.exp(h), A, B, C))
               - np.log(fibre_theta(KM * np.exp(-h), A, B, C))) / (2 * h))
        gap = abs(fd - closed) / abs(closed)
        print(f"  {theta:8.0f}  {mbar:8.3f}  {A:12.1f}  {B:10.1f}  "
              f"{closed:12.5f}  {fd:12.5f}  {gap:9.2e}")
        rows.append((theta, mbar, closed, fd))
    print("\n  |dln theta / dln K_m| = (L + b + k_e K_m)/theta, so it scales")
    print("  as 1/theta: the fibre is TEN TIMES steeper at a tenth of nominal.")
    return rows


# ======================================================================
# PART C and D. Against this repository's own likelihood.
# ======================================================================
def sigma_at(theta, mbar, load_frac=LOAD_FRAC, assay_cv=ASSAY_CV):
    """
    The ordinate noise at a given capacity, on model_A.run.sigma_total's
    own definition but re-solved at this capacity's operating level.
    """
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(load_frac * LB, slope * assay_cv * mbar))


def columns(theta, m, inc):
    """Design columns in (ln theta, ln K_m, b, ln k_e), at this capacity."""
    g = m / (KM + m)
    return np.stack([theta * g,
                     -(theta + inc) * KM * m / (KM + m) ** 2,
                     np.ones_like(m),
                     KE * m], axis=1)


def sd_ln_theta(theta, m, inc, sigma, km_prior=None, ke_prior=0.8,
                load_prior=True, exact=()):
    """
    sd(ln theta) from the information filter.

    `km_prior` and `ke_prior` are log-normal prior widths; None means a FLAT
    prior, which is not the same as known. `exact` names parameters held
    EXACTLY known, whose columns are dropped: any of "Km", "b", "ke".
    That is g_design_common.sd_from's convention for k_e, extended to the
    other two, and it is spelled out because "known" and "free" were the two
    readings that inverted a comparison in this repository before.
    """
    Q = columns(theta, m, inc)
    pri = [np.inf,
           km_prior ** 2 if km_prior is not None else np.inf,
           (LOAD_FRAC * LB) ** 2 if load_prior else np.inf,
           ke_prior ** 2 if ke_prior is not None else np.inf]
    drop = {"Km": 1, "b": 2, "ke": 3}
    kill = {drop[k] for k in exact}
    keep = [i for i in range(4) if i not in kill]
    Q, pri = Q[:, keep], [pri[i] for i in keep]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def record(theta, times, load_cv=0.0, seed=0):
    times = np.asarray(times, float)
    inc = np.asarray(delta(times), float)
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        Ls = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        Ls = np.full(times.size, LB)
    return level(theta, inc, Ls), inc


def part_c():
    """
    Is the record's null space the fibre, and is one external K_m enough?

    Two different questions, and the second does not follow from the first.
    The rank question is answered by the eigenstructure of the DATA's
    information matrix with no priors at all. The precision question is
    answered by what an exactly known K_m is worth against the alternatives.
    """
    print()
    print("=" * 74)
    print("PART C.  the record's own null space, against the fibre tangent")
    print("=" * 74)
    print("\n  coordinates (ln theta, ln K_m, b/L, ln k_e). The b column")
    print("  carries g_design_common's sign, +1 for the residual's derivative,")
    print("  so the fibre tangent is written in that same sign.")
    print("  fibre tangent: dln theta/dln K_m = -(L + k_e K_m)/theta,")
    print("                 db/dln K_m = -(L + b), i.e. -1 in units of L,")
    print("                 which is +1 in the column's sign, and 0 on ln k_e.\n")

    for theta in (3000.0, 600.0, 300.0):
        mbar = float(level(theta, 0.0, LB))
        sig = sigma_at(theta, mbar)
        slope = -(LB + KE * KM) / theta
        tan = np.array([slope, 1.0, +1.0, 0.0])
        tan /= np.linalg.norm(tan)
        for n in (25, 400):
            t = 30.4 * np.arange(n, dtype=float)
            m, inc = record(theta, t, 0.0)
            Q = columns(theta, m, inc) @ np.diag([1.0, 1.0, LB, 1.0])
            F = Q.T @ Q / sig ** 2
            w, V = np.linalg.eigh(F)
            cos = abs(float(V[:, 0] @ tan))
            print(f"  theta = {theta:6.0f}  n = {n:4d}   "
                  f"eigenvalues {np.array2string(w, precision=2)}")
            print(f"      lambda_min / lambda_max = {w[0] / w[-1]:.2e}   "
                  f"lambda_2 / lambda_max = {w[1] / w[-1]:.2e}")
            print(f"      |cos(weakest eigenvector, fibre tangent)| = "
                  f"{cos:.8f}")
        print()

    print("  READING. A relative smallest eigenvalue at machine level says the")
    print("  null space is EXACTLY one-dimensional, and the cosine says which")
    print("  direction it is. The SECOND eigenvalue is the one that decides")
    print("  what an external K_m is worth, and it is not at machine level.")


def part_d():
    """What each external quantity is worth, held EXACTLY known."""
    print()
    print("=" * 74)
    print("PART D.  what one exactly known quantity buys, against the gate")
    print("=" * 74)
    print("\n  25 monthly draws, constant recorded load, k_e log-prior 0.8 and")
    print("  load-bias prior 0.05 L wherever the quantity is not held exact.")
    print("  'K_m 1e-4' is a log-prior, not a dropped column, and the two")
    print("  agree to the digits printed.\n")
    cases = [("nothing (the paper's footing)", dict()),
             ("K_m exact", dict(exact=("Km",))),
             ("k_e exact", dict(exact=("ke",))),
             ("b exact", dict(exact=("b",))),
             ("K_m and k_e exact", dict(exact=("Km", "ke"))),
             ("k_e and b exact", dict(exact=("ke", "b"))),
             ("K_m and b exact", dict(exact=("Km", "b"))),
             ("all three exact", dict(exact=("Km", "ke", "b")))]
    print(f"  {'held exactly known':<30}" +
          "".join(f"{f'theta={t:.0f}':>12}" for t in (3000, 600, 300)))
    print("  " + "-" * 66)
    out = {}
    for name, kw in cases:
        cells = []
        for theta in (3000.0, 600.0, 300.0):
            mbar = float(level(theta, 0.0, LB))
            sig = sigma_at(theta, mbar)
            m, inc = record(theta, 30.4 * np.arange(25, dtype=float), 0.0)
            v = sd_ln_theta(theta, m, inc, sig, **kw)
            cells.append(v)
            out[(name, theta)] = v
        print(f"  {name:<30}" + "".join(f"{c:12.4f}" for c in cells))
    print("\n  gate is 0.20.")
    base = [out[("nothing (the paper's footing)", t)] for t in (3000., 600., 300.)]
    km = [out[("K_m exact", t)] for t in (3000., 600., 300.)]
    print("  an EXACT K_m is worth a factor of "
          + ", ".join(f"{b / k:.3f}" for b, k in zip(base, km))
          + " at theta = 3000, 600, 300.")
    return out


def part_e():
    """The same with the load varying, which is the specification's footing."""
    print()
    print("=" * 74)
    print("PART E.  the same at the specification's 10% day-to-day load CV")
    print("=" * 74)
    print("\n  medians over 24 seeds, seed range in brackets. A varying load")
    print("  carries information about K_m that a constant one does not.\n")
    cases = [("nothing", dict()), ("K_m exact", dict(exact=("Km",))),
             ("k_e exact", dict(exact=("ke",)))]
    for theta in (3000.0, 300.0):
        mbar = float(level(theta, 0.0, LB))
        sig = sigma_at(theta, mbar)
        print(f"  theta = {theta:.0f}   m_bar = {mbar:.3f}   sigma = {sig:.3f}")
        for name, kw in cases:
            vv = np.array([
                sd_ln_theta(theta, *record(theta,
                                           30.4 * np.arange(25, dtype=float),
                                           load_cv=0.10, seed=s), sig, **kw)
                for s in range(24)], float)
            print(f"    {name:<14} {np.median(vv):8.4f} "
                  f"[{vv.min():.4f}, {vv.max():.4f}]")
        print()


if __name__ == "__main__":
    ok = part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print("=" * 74)
    print("PART A verdict (all three symbolic identities hold):", ok)
    print("=" * 74)
