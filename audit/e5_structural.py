"""
3.7 (completion) -- reconciling the machine's structural verdict with the
note's numerical findings.

StructuralIdentifiability.jl returns GLOBALLY IDENTIFIABLE for every
configuration in symbolic/julia/run_identifiability.jl except B. The note
reports obstructions. This script asks whether that is a real disagreement.

Three questions, each settled by computation rather than by argument:

  Q1  Does the Julia encoding "A2 = constant input" represent the note's
      "patient observed at steady state"?
  Q2  If not, what does the steady-state case actually give?
  Q3  What kind of object is the paper's rank-one deficiency, given that the
      model is globally identifiable?
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp

from common.models import ParamsA, g

P = ParamsA()
TH, KM, KE = P.theta, P.Km, P.ke


def banner(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


# ---------------------------------------------------------------------------
def q1_constant_input_is_not_steady_state():
    banner("Q1  a constant INPUT is not a constant STATE")

    # steady state under a constant load
    L = TH * g(20.0, KM) + KE * 20.0          # a load holding m at 20
    f = lambda t, y: [L - TH * g(y[0], KM) - KE * y[0]]

    print(f"  constant load L = {L:.1f}, steady state m* = 20.0")
    print("  StructuralIdentifiability.jl assumes GENERIC initial conditions,")
    print("  so the trajectory it analyses starts away from m* and relaxes.")
    print()
    print("  levels visited, from a generic start:")
    for m0 in (5.0, 20.0, 60.0):
        sol = solve_ivp(f, (0, 5.0), [m0], dense_output=True, rtol=1e-10, atol=1e-12)
        mm = sol.sol(np.linspace(0, 5.0, 4001))[0]
        print(f"    m(0) = {m0:5.1f}  ->  m spans [{mm.min():.3f}, {mm.max():.3f}]"
              f"   ({'an interval' if mm.max()-mm.min() > 1e-6 else 'a point'})")

    # how long the transient lasts
    sol = solve_ivp(f, (0, 5.0), [5.0], dense_output=True, rtol=1e-10, atol=1e-12)
    tg = np.linspace(0, 5.0, 200001)
    mm = sol.sol(tg)[0]
    gap = np.abs(mm - 20.0)
    tau = tg[np.argmax(gap < 0.01 * gap[0])]
    print(f"\n  time to close 99% of the gap: {tau*24:.2f} hours")
    print("  => the interval of Proposition 'global' IS visited, in hours.")
    print("  The Julia A2 therefore asks a DIFFERENT question from the note's")
    print("  steady-state patient. Its 'globally identifiable' is a statement")
    print("  about that transient. The comment in run_identifiability.jl which")
    print("  said A2 'encodes a patient held at steady state' was wrong.")
    return tau


# ---------------------------------------------------------------------------
def q2_the_steady_state_case():
    banner("Q2  what the steady-state case actually gives")

    print("  A patient AT steady state supplies one algebraic equation")
    print("      s = theta g(m; K_m) + k_e m")
    print("  in three unknowns. Deficiency 2, exactly as the note states.")
    print()
    print("  Exhibiting the two-dimensional fibre explicitly: fix (s, m) and")
    print("  sweep (K_m, k_e); theta is then determined and the residual is 0.")
    m_obs = 20.0
    s_obs = TH * g(m_obs, KM) + KE * m_obs
    print(f"    observed (s, m) = ({s_obs:.2f}, {m_obs:.1f})")
    print(f"    {'K_m':>10} {'k_e':>10} {'theta':>12} {'residual':>12}")
    for km in (KM / 3, KM, KM * 3):
        for ke in (KE / 2, KE, KE * 2):
            th = (s_obs - ke * m_obs) / g(m_obs, km)
            if th <= 0:
                continue
            r = s_obs - th * g(m_obs, km) - ke * m_obs
            print(f"    {km:10.2f} {ke:10.3f} {th:12.2f} {r:12.2e}")
    print("\n  Every row reproduces the observation exactly. The fibre is real,")
    print("  and it is NOT what Julia's A2 was asked about.")


# ---------------------------------------------------------------------------
def q3_what_kind_of_object_is_the_rank_deficiency(lam=np.log(2) / 3.5,
                                                  interval=14.0):
    banner("Q3  the rank-one deficiency is a FIRST-ORDER degeneracy, "
           "not a structural one")

    # the quasi-steady locus the clinic actually samples: m follows the dose
    # train algebraically, because the marker relaxes in hours and the dose
    # decays over days.
    L0 = 0.0
    L = TH * g(10.0, KM) + KE * 10.0

    def m_of(delta, th=TH, km=KM, ke=KE, load=None):
        load = L if load is None else load
        lo, hi = 1e-9, 1e6
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (th + delta) * g(mid, km) + ke * mid < load:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    # a month of biweekly dosing, sampled daily, increment amplitude 20% of theta
    t = np.arange(0.0, 28.0, 1.0)
    amp = 0.20 * TH
    dose = amp * np.exp(-lam * (t % interval))

    def predict(p):
        th, km, ke, b = p
        return np.array([m_of(d, th, km, ke, load=L + b) for d in dose])

    p0 = np.array([TH, KM, KE, L0])
    base = predict(p0)

    # column-scaled Jacobian in log parameters (b is scaled by the load)
    steps = np.array([1e-5 * TH, 1e-5 * KM, 1e-5 * KE, 1e-5 * L])
    J = np.zeros((len(t), 4))
    for j in range(4):
        pp = p0.copy(); pp[j] += steps[j]
        pm = p0.copy(); pm[j] -= steps[j]
        J[:, j] = (predict(pp) - predict(pm)) / (2 * steps[j]) * steps[j]
    sv = np.linalg.svd(J, compute_uv=False)
    print("  singular values of the 4-column design Jacobian, AT THE TRUTH:")
    print("   ", np.array2string(sv, formatter={'float_kind': lambda x: f'{x:.3e}'}))
    print(f"    s_min / s_max = {sv[-1]/sv[0]:.3e}")

    # the flat direction
    _, _, Vt = np.linalg.svd(J)
    v = Vt[-1]
    print("    flat direction (theta, K_m, k_e, b), column-scaled:")
    print("   ", np.array2string(v, formatter={'float_kind': lambda x: f'{x:+.4f}'}))

    # move along it and watch the noiseless misfit
    print("\n  noiseless chi^2 along the flat direction, in units of eps:")
    print(f"    {'eps':>10} {'rms misfit':>14} {'ratio to eps^2':>16}")
    for eps in (1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
        pp = p0 + eps * (v * steps)
        r = predict(pp) - base
        rms = np.sqrt(np.mean(r ** 2))
        print(f"    {eps:10.4f} {rms:14.6e} {rms/eps**2:16.4e}")
    # ---- the flat direction, derived by hand ---------------------------
    print()
    print("  THE SAME DIRECTION, DERIVED BY HAND FROM THE QUASI-STEADY BALANCE")
    print("  ----------------------------------------------------------------")
    print("  On the locus (theta + Delta) u + k_e m = L + b, with u = m/(K_m+m),")
    print("  the sensitivities of the balance are")
    print("      dG/dtheta = u,   dG/dK_m = -(theta+Delta) u(1-u)/K_m,")
    print("      dG/dk_e   = m,   dG/db   = -1,")
    print("  and (theta+Delta) u = L - k_e m eliminates the increment. Using")
    print("  m(1-u) = K_m u, a combination a,b,c,e vanishes for every Delta iff")
    print("      c = 0                     (m = K_m u/(1-u) is not polynomial in u)")
    print("      a + bL/K_m + b k_e = 0    (coefficient of u)")
    print("      e + bL/K_m = 0            (constant)")
    print("  so, taking b = -K_m,")
    print("      n = (L + k_e K_m,  -K_m,  0,  L).")
    n = np.array([L + KE * KM, -KM, 0.0, L])
    v_pred = n / steps
    v_pred = v_pred / np.linalg.norm(v_pred)
    if np.dot(v_pred, v) < 0:
        v_pred = -v_pred
    print(f"\n    predicted (column-scaled, normalised):")
    print("   ", np.array2string(v_pred, formatter={'float_kind': lambda x: f'{x:+.4f}'}))
    print(f"    computed:")
    print("   ", np.array2string(v, formatter={'float_kind': lambda x: f'{x:+.4f}'}))
    print(f"    max abs difference: {np.max(np.abs(v_pred - v)):.2e}")
    print(f"    ||J n|| / (||J|| ||n||) = "
          f"{np.linalg.norm(J @ (n/steps)) / (np.linalg.norm(J) * np.linalg.norm(n/steps)):.3e}")
    print()
    print("  Dividing through by theta L puts it in the note's normalisation:")
    print("      c_theta = (1 + K_m k_e / L)/theta,   c_Km = -1/L,")
    print("  which is the identity the note states. It is confirmed here by a")
    print("  route that shares no code with the design module.")

    print("\n  The misfit is quadratic in eps, not linear: the direction is flat")
    print("  to FIRST order and curved at second order. A singular Fisher matrix")
    print("  at a point is not a structural non-identifiability -- the map stays")
    print("  injective, and the second-order term does the identifying.")
    print()
    print("  Consequence for precision: a first-order-flat direction is resolved")
    print("  only when the second-order term clears the noise, so the recoverable")
    print("  displacement scales as the SQUARE ROOT of the noise level, not")
    print("  linearly. That is the sense in which the note's obstruction is real")
    print("  while the model is globally identifiable.")


if __name__ == "__main__":
    print(__doc__)
    q1_constant_input_is_not_steady_state()
    q2_the_steady_state_case()
    q3_what_kind_of_object_is_the_rank_deficiency()
