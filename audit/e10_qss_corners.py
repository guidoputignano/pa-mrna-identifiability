"""
Does the quasi-steady reduction of the carnitine state hold away from the
illustrative centre, at the corners of the cohort's sampled region?

WHAT WAS ALREADY CHECKED.  tests/test_models.py checks the reduction at one
parameter set: rhs_C settles onto the reduced attractor from below the
separatrix, an integrated steady state matches the closed-form equilibrium,
and halving the step does not move the answer.  All three are run at PC, the
illustrative centre.

WHAT THIS ADDS, AND WHAT IT FINDS.  The cohort of audit/e8_cohort_unblocked.py
spans theta and K_m over nine-fold ranges and the managed load over four-fold.
Running the same comparison at the corners of that box returns a relative
discrepancy at the level of machine epsilon everywhere -- and that result is
NOT evidence that the timescales separate.  It is an algebraic identity, and
it is worth being explicit about why, because the identity is easy to mistake
for a numerical confirmation of something it does not touch.

alpha(m) is DEFINED by solving dc/dt = 0.  At any fixed point of the two-state
system dc/dt = 0 holds exactly.  So the parent's fixed-point condition and the
reduced one, s = theta g(m) + alpha(m), are the same equation rearranged, and
they agree to machine precision at every parameter set, at every load, and for
every value of the conjugation constant k -- including values at which the
carnitine state is SLOWER than the metabolite and the reduction's premise is
false.  An equilibrium comparison cannot fail, so it cannot report degradation
either, and the existing tests inherit that blindness.

Where the premise does bite is off equilibrium, and this script measures it in
the two places the paper relies on:

    transient    integrating the parent and the reduced system from the same
                 state, the trajectories differ by 3.4 to 25.4 per cent at the
                 default k = 1.  fig_timeseries in model_C/run.py drives an
                 illness pulse through the reduced system; that is this error.

    separatrix   the basin boundary is the paper's central claim about Model C.
                 Started ON the slow manifold it is reproduced exactly, for the
                 same fixed-point reason.  Started OFF it the displacement is
                 SIGNED, and the sign belongs to the carnitine pool rather than
                 to the reduction.  From a fuller pool -- carnitine still at its
                 pre-perturbation value, which is what a sudden load spike
                 leaves behind -- the true boundary sits 2.9 to 49.6 per cent
                 ABOVE the reduced separatrix and the reduction calls
                 decompensation early, which is the safe direction.  From a pool
                 depleted to half of c*(m) it sits 0.5 to 31.9 per cent BELOW
                 it, and there the reduction calls safe a state the parent
                 escapes from.  Reporting only the first would be true and
                 misleading, because secondary carnitine depletion is the
                 mechanism this model is built on, so the unsafe side is the one
                 the disease produces.  (Every figure quoted in this docstring
                 is printed by the run, and the run computes its own closing
                 summary rather than repeating these.)

The quantity all of this tracks is the timescale ratio

    eps = tau_c / tau_m,    tau_c = 1/[k (u + rho)],   tau_m = 1/|C'(m*)|

and eps is proportional to 1/k.  k is the conjugation rate constant.  It
cancels out of alpha exactly -- common/params.py records it as DERIVED, "drops
out of the quasi-steady reduction exactly" -- and that is true of the reduced
CURVE.  It is not true of the reduction's VALIDITY, which k alone sets.  At the
k = 1.0 that rhs_C defaults to and the tests run at, eps is between 0.52 and 55
across the sampled region: the carnitine state is comparable to or slower than
the metabolite, and there is no separation of timescales to appeal to.

None of this disturbs a conclusion drawn at equilibrium, which is what the
estimator and every identifiability result in the note use.  It bounds the
dynamical claims, and it says that "carnitine is fast" is at present a choice
of k rather than a measurement.

Run:  python audit/e10_qss_corners.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dataclasses import replace

import numpy as np
from scipy.integrate import solve_ivp

from common.models import (rhs_C, alt_flux_C, clearance_C_quasistatic,
                           equilibria_nonmonotone)
from model_C.run import PC, S_BASE, R0_for_peak

# The solver configuration of the existing parent-integration test,
# tests/test_models.py::test_the_two_state_system_settles_on_the_reduced_equilibria:
# adaptive, stiff-capable, and tighter than the rtol 1e-8 asked for.
SOLVER = dict(method="LSODA", rtol=1e-10, atol=1e-12)
K_CONJ = 1.0            # rhs_C's own default, and what the tests run at
TOL = 1e-6              # the tolerance the existing tests assert at
GRID = np.linspace(1e-3, 3e4, 120000)   # as in the existing tests

# Horizons, in multiples of the slower of the two relaxation times.  The
# separatrix is a SADDLE: a trajectory launched near it leaves slowly, so the
# escape verdict needs a longer run than the settling test does.  At 400x the
# bisection reports a spurious 7e-5 on the widest-shape set purely because the
# horizon ran out; at 2000x that set agrees to 4e-12 like every other, and
# 10000x does not move it further.
T_SETTLE, T_SADDLE = 400.0, 2000.0


# --------------------------------------------------------------- the box
#
# The cohort of e8_cohort_unblocked.py is log-uniform over theta (9-fold),
# K_m (9-fold), k_e (2.2-fold) and the managed load (4-fold).  Three of those
# four axes carry over to the parent unchanged.  The fourth does not: k_e is
# the LOCAL model's linear alternative channel, and the parent has no such
# parameter -- its alternative route IS the carnitine mechanism, carried by
# (R0, Kc, Mi, rho0).  The last two sets below stand in for that axis by
# putting the parent's own shape constants at the extremes of the admissible
# ranges that model_C.run.collapse_check samples, with R0 re-derived through
# R0_for_peak so the calibration is held rather than silently broken.
FT, FK, FL = 3.0, 3.0, 2.0      # the cohort's own half-widths

SETS = [
    ("centre (control)",        1.0,   1.0,   1.0,  {}),
    ("theta lo, Km lo, load lo", 1/FT, 1/FK, 1/FL,  {}),
    ("theta lo, Km hi, load lo", 1/FT,  FK,  1/FL,  {}),
    ("theta hi, Km lo, load lo",  FT,  1/FK, 1/FL,  {}),
    ("theta hi, Km hi, load lo",  FT,   FK,  1/FL,  {}),
    ("theta hi, Km lo, load hi",  FT,  1/FK,  FL,   {}),
    ("theta hi, Km hi, load hi",  FT,   FK,   FL,   {}),
    ("theta lo, Km lo, load hi", 1/FT, 1/FK,  FL,   {}),
    ("shape lo (Kc,Mi,rho0)",    1.0,   1.0,  1.0,
     dict(Kc=10 ** 1.4, Mi=10 ** 2.0, rho0=0.35)),
    ("shape hi (Kc,Mi,rho0)",    1.0,   1.0,  1.0,
     dict(Kc=10 ** 2.9, Mi=10 ** 3.1, rho0=1.00)),
]


def build(ft, fk, fl, shape):
    """One parameter set: the parent at a corner, calibration held."""
    p = replace(PC, theta=PC.theta * ft, Km=PC.Km * fk)
    if shape:
        p = replace(p, **shape)
        p = replace(p, R0=R0_for_peak(p))    # re-derive, as collapse_check does
    return p, S_BASE * fl


# --------------------------------------------------------- the two states
def c_star(m, p, k=K_CONJ):
    """The slow manifold: the carnitine level the reduction assumes."""
    u = m / (p.Kc + m)
    rho = p.rho0 * (1.0 + m / p.Mi)
    return p.R0 / (k * (u + rho))


def timescales(m, p, k=K_CONJ):
    """
    (tau_c, tau_m, eps).  tau_c is the relaxation of the carnitine state at
    frozen m; tau_m is the relaxation of the metabolite ON the reduced
    manifold, 1/|C'(m)|.  eps = tau_c/tau_m is the small parameter the
    reduction needs, and it is proportional to 1/k.
    """
    u = m / (p.Kc + m)
    rho = p.rho0 * (1.0 + m / p.Mi)
    tau_c = 1.0 / (k * (u + rho))
    h = max(m * 1e-6, 1e-9)
    Cp = (clearance_C_quasistatic(m + h, p)
          - clearance_C_quasistatic(m - h, p)) / (2.0 * h)
    tau_m = np.inf if Cp == 0 else 1.0 / abs(Cp)
    return tau_c, tau_m, tau_c / tau_m


def eig_ratio(m, p, k=K_CONJ):
    """|slow eigenvalue| / |fast eigenvalue| of the parent's Jacobian at (m, c*)."""
    c = c_star(m, p, k)
    u, du = m / (p.Kc + m), p.Kc / (p.Kc + m) ** 2
    gp = p.Km / (p.Km + m) ** 2
    rho, drho = p.rho0 * (1.0 + m / p.Mi), p.rho0 / p.Mi
    J = np.array([[-p.theta * gp - k * c * du, -k * u],
                  [-k * c * du - k * drho * c, -k * (u + rho)]])
    w = np.sort(np.abs(np.linalg.eigvals(J)))
    return float(w[0] / w[1])


def reduced_equilibria(p, s):
    """(attractor, separatrix) of the closed form; either may be None."""
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, p, s, GRID)
    if len(roots) == 0:
        return None, None
    st = np.asarray(st, dtype=bool)
    attr = float(roots[st][0]) if st.any() else None
    sep = float(roots[~st][0]) if (~st).any() else None
    return attr, sep


def parent_steady_state(p, s, m0, c0, T, k=K_CONJ):
    """Integrate the two-state parent to steady state, with a convergence guard."""
    def run(horizon):
        sol = solve_ivp(lambda t, y: rhs_C(t, y, p, lambda _: s, k),
                        (0.0, horizon), [m0, c0], **SOLVER)
        return float(sol.y[0, -1])
    a, b = run(T), run(2.0 * T)
    return b, abs(b - a) / max(abs(b), 1e-300)


def trajectories(p, s, m0, T, k=K_CONJ, n=4001):
    """Parent and reduced m(t) from the same state, on a common grid."""
    tt = np.linspace(0.0, T, n)
    par = solve_ivp(lambda t, y: rhs_C(t, y, p, lambda _: s, k),
                    (0.0, T), [m0, c_star(m0, p, k)], dense_output=True,
                    **SOLVER).sol(tt)[0]
    red = solve_ivp(lambda t, y: [s - clearance_C_quasistatic(y[0], p)],
                    (0.0, T), [m0], dense_output=True, **SOLVER).sol(tt)[0]
    return np.max(np.abs(par - red) / np.maximum(np.abs(red), 1e-30))


def critical_level(p, s, sep, attr, T, start, k=K_CONJ):
    """
    The parent's true basin boundary: bisect on the initial metabolite level
    for the escape/no-escape verdict.  `start` chooses the carnitine the
    trajectory begins with, and it is the whole question, because the basin
    boundary is a CURVE in (m, c) and the reduced separatrix is one point on
    it rather than a level of m:

        "manifold"  c*(m0), the value the reduction assumes
        "pre"       c*(attractor), the pre-perturbation pool, which is what a
                    sudden load spike leaves behind.  c* falls with m, so this
                    is FULLER than the manifold and the parent is forgiving
        a float f   f * c*(m0), a pool depleted below the manifold, which is
                    what this model's own mechanism produces.  Here the parent
                    is LESS forgiving than the reduction and the error changes
                    sign
    """
    thresh = 10.0 * sep

    def escapes(m0):
        if start == "manifold":
            c0 = c_star(m0, p, k)
        elif start == "pre":
            c0 = c_star(attr, p, k)
        else:
            c0 = float(start) * c_star(m0, p, k)
        sol = solve_ivp(lambda t, y: rhs_C(t, y, p, lambda _: s, k),
                        (0.0, T), [m0, c0], **SOLVER)
        return float(sol.y[0, -1]) > thresh

    lo, hi = 0.3 * sep, 8.0 * sep
    if escapes(lo) or not escapes(hi):
        return np.nan
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if escapes(mid):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def escape_pair(p, s, m0=1.0, T=4000.0):
    """m(T) for parent and reduced where no equilibrium exists: both must run away."""
    par = solve_ivp(lambda t, y: rhs_C(t, y, p, lambda _: s, K_CONJ),
                    (0.0, T), [m0, c_star(m0, p)], **SOLVER)
    red = solve_ivp(lambda t, y: [s - clearance_C_quasistatic(y[0], p)],
                    (0.0, T), [m0], **SOLVER)
    return float(par.y[0, -1]), float(red.y[0, -1])


# ------------------------------------------------------------------ report
def flag(x):
    return " " if (not np.isfinite(x) or x <= TOL) else "*"


def main():
    print(__doc__)
    rows = []
    print("=" * 100)
    print("  THE COMPARISON THE REDUCTION IS USUALLY ASKED FOR: STEADY STATE")
    print("=" * 100)
    print(f"  solver: {SOLVER['method']}, rtol={SOLVER['rtol']:.0e}, "
          f"atol={SOLVER['atol']:.0e}; conjugation constant k = {K_CONJ}")
    print(f"  flagged with * if the relative discrepancy exceeds {TOL:.0e}\n")
    print(f"  {'parameter set':<26}{'eps=tau_c/tau_m':>16}{'|lam_s/lam_f|':>15}"
          f"{'m* (closed form)':>18}{'rel discrepancy':>17}")
    print("  " + "-" * 90)

    for name, ft, fk, fl, shape in SETS:
        p, s = build(ft, fk, fl, shape)
        attr, sep = reduced_equilibria(p, s)
        if attr is None:
            rows.append((name, p, s, None, sep, np.nan, np.nan, np.nan))
            print(f"  {name:<26}{'--':>16}{'--':>15}{'no equilibrium':>18}"
                  f"{'undefined':>17}")
            continue
        tau_c, tau_m, eps = timescales(attr, p)
        er = eig_ratio(attr, p)
        T = T_SETTLE * max(tau_m, tau_c)
        m_end, guard = parent_steady_state(p, s, 0.5 * attr, c_star(0.5 * attr, p), T)
        rel = abs(m_end - attr) / attr
        rows.append((name, p, s, attr, sep, eps, rel,
                     T_SADDLE * max(tau_m, tau_c)))
        print(f"  {name:<26}{eps:>16.3e}{er:>15.3e}{attr:>18.6f}"
              f"{rel:>16.3e}{flag(rel)}")
        if guard > 1e-10:
            print(f"  {'':<26}(warning: horizon guard {guard:.2e})")

    print("\n  Every entry is at machine epsilon.  That is an identity, not a "
          "measurement:\n  alpha is defined by dc/dt = 0, and dc/dt = 0 holds "
          "exactly at any fixed point, so\n  the two equilibrium conditions are "
          "one equation.  It holds at every eps above,\n  including eps > 1, "
          "where carnitine is SLOWER than the metabolite.")

    # ---------------------------------------------------------------- off equilibrium
    print("\n" + "=" * 100)
    print("  THE COMPARISON THAT CAN FAIL: OFF EQUILIBRIUM")
    print("=" * 100)
    print(f"  {'parameter set':<26}{'eps':>12}{'transient max rel':>19}"
          f"{'separatrix, on-mfld':>21}{'sep, full pool':>17}{'sep, depleted':>16}")
    print("  " + "-" * 112)
    off_eq = []
    for name, p, s, attr, sep, eps, _, T in rows:
        if attr is None:
            # No equilibrium exists, so there is nothing to compare a level
            # against.  What can still be checked is that the reduction gets
            # the QUALITATIVE verdict right: both systems must run away.
            a, b = escape_pair(p, s)
            verdict = ("both escape, m(T) agreeing to "
                       f"{abs(a - b) / b:.1e}" if (a > 1e3 and b > 1e3)
                       else f"VERDICTS DIFFER: parent {a:.3e}, reduced {b:.3e}")
            print(f"  {name:<26}{'--':>12}{verdict:>77}")
            continue
        tau_m = timescales(attr, p)[1]
        tr = trajectories(p, s, 0.5 * attr, 60.0 * tau_m)
        if sep is None:
            on = off = signed = depl = np.nan
        else:
            a = critical_level(p, s, sep, attr, T, "manifold")
            b = critical_level(p, s, sep, attr, T, "pre")
            d = critical_level(p, s, sep, attr, T, 0.5)
            on = abs(a - sep) / sep if np.isfinite(a) else np.nan
            off = abs(b - sep) / sep if np.isfinite(b) else np.nan
            signed = (b - sep) / sep if np.isfinite(b) else np.nan
            depl = (d - sep) / sep if np.isfinite(d) else np.nan
        f = lambda x: f"{x:.3e}{flag(x)}" if np.isfinite(x) else "no separatrix"
        g = lambda x: f"{100*x:+.2f}%" if np.isfinite(x) else "--"
        print(f"  {name:<26}{eps:>12.3e}{f(tr):>19}{f(on):>21}"
              f"{g(signed):>17}{g(depl):>16}")
        off_eq.append((tr, on, off, signed, depl))

    # ---------------------------------------------------------------- k sweep
    print("\n" + "=" * 100)
    print("  WHAT THE DISCREPANCY ACTUALLY TRACKS: THE CONJUGATION CONSTANT k")
    print("=" * 100)
    print("  k cancels out of alpha exactly, so the reduced curve does not depend on")
    print("  it.  The reduction's VALIDITY depends on nothing else: eps is "
          "proportional to 1/k.\n")
    p, s = build(1.0, 1.0, 1.0, {})
    attr, sep = reduced_equilibria(p, s)
    print(f"  {'k':>8}{'eps':>12}{'steady state':>16}{'transient':>13}"
          f"{'separatrix, off-mfld':>23}")
    print("  " + "-" * 72)
    for k in (0.1, 1.0, 10.0, 100.0):
        tau_c, tau_m, eps = timescales(attr, p, k)
        T = T_SADDLE * max(tau_m, tau_c)
        m_end, _ = parent_steady_state(p, s, 0.5 * attr, c_star(0.5 * attr, k=k, p=p),
                                       T, k=k)
        rel = abs(m_end - attr) / attr
        tr = trajectories(p, s, 0.5 * attr, 60.0 * tau_m, k=k)
        b = critical_level(p, s, sep, attr, T, False, k=k)
        off = abs(b - sep) / sep if np.isfinite(b) else np.nan
        print(f"  {k:>8.1f}{eps:>12.3e}{rel:>15.3e}{flag(rel)}"
              f"{tr:>12.3e}{flag(tr)}{off:>22.3e}{flag(off)}")
    print("\n" + "=" * 100)
    print("  WHAT THIS MEANS")
    print("=" * 100)

    ss = [r[6] for r in rows if np.isfinite(r[6])]
    tr_ = [a for a, _, _, _, _ in off_eq if np.isfinite(a)]
    on_ = [b for _, b, _, _, _ in off_eq if np.isfinite(b)]
    off_ = [c for _, _, c, _, _ in off_eq if np.isfinite(c)]
    sgn_ = [d for _, _, _, d, _ in off_eq if np.isfinite(d)]
    dep_ = [e for _, _, _, _, e in off_eq if np.isfinite(e)]
    # The sign is NOT a property of the reduction, it is a property of which
    # side of the slow manifold the carnitine pool starts on.  Reporting only
    # the pre-perturbation start would say "conservative" and be wrong by
    # omission, because this model's own mechanism depletes the pool.
    up = all(d > 0 for d in sgn_)
    down = all(e < 0 for e in dep_)
    conservative = (
        "signed, and the sign is the carnitine pool's, not the reduction's. From a\n"
        "      FULLER pool, the pre-perturbation start a load spike leaves behind, the\n"
        f"      parent tolerates MORE than the reduced model predicts ({100*min(sgn_):+.1f}% to\n"
        f"      {100*max(sgn_):+.1f}%), so there the error is conservative. From a pool DEPLETED to\n"
        f"      half of c*(m) it tolerates LESS ({100*min(dep_):+.1f}% to {100*max(dep_):+.1f}%), and there the\n"
        "      reduced model calls safe a state the parent escapes from. The reduced\n"
        "      separatrix is one point on a curve in (m, c), not a bound in either\n"
        "      direction, and the unsafe side is the one this model's own mechanism\n"
        "      produces"
        if (up and down) else
        f"      NOT cleanly signed: {sum(d > 0 for d in sgn_)}/{len(sgn_)} up from a full pool, "
        f"{sum(e < 0 for e in dep_)}/{len(dep_)} down from a depleted one")
    pct = lambda v: f"{100 * min(v):.1f} to {100 * max(v):.1f} per cent"

    print(f"""
  1.  At equilibrium the reduction is exact everywhere in the sampled region:
      the relative discrepancy runs {min(ss):.1e} to {max(ss):.1e} over {len(ss)} sets, which is
      machine epsilon and not a measurement of anything.  Nothing read off an
      equilibrium -- the estimator, the fibre, the identifiability results,
      theta_crit -- is exposed to the reduction at all.  The same holds for the
      separatrix approached ON the slow manifold ({min(on_):.1e} to {max(on_):.1e}), which is
      a fixed point for the same reason.

  2.  So the existing tests assert 1e-6 on a quantity whose value is ~1e-16 and
      which cannot move.  They are regression tests on the algebra, worth
      keeping as that, but they are not evidence of a separation of timescales.
      Neither is this table's first block.

  3.  Off equilibrium the reduction is wrong by {pct(tr_)} along the
      transient, and the basin boundary is displaced by {pct(off_)}
      against a perturbation that leaves carnitine behind.  The displacement is
      {conservative}.

  4.  All of that is set by k, the conjugation constant, which is currently a
      free scale.  eps is proportional to 1/k, and at the k = 1.0 the code
      defaults to, eps spans {min(r[5] for r in rows if np.isfinite(r[5])):.2f} to {max(r[5] for r in rows if np.isfinite(r[5])):.1f} across these sets -- carnitine is
      comparable to or SLOWER than the metabolite, which is the opposite of the
      reduction's premise.  Justifying "carnitine is fast" needs k large enough
      that eps << 1 at the corners, of order 100 from the sweep above, and
      nothing in the repository measures it.
""")


if __name__ == "__main__":
    main()
