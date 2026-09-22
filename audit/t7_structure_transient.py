"""
6. THE TRANSIENT AS A STRUCTURE TEST, ON THE RIGHT CLEARANCE.

Scoping only. Nothing here reaches the paper.

Section 6.3 says the one measurement that would settle which structure a real
record came from is a transient: the published alternative's marker is
ALGEBRAIC in the enzyme, so it has no relaxation of its own, while a marker
that is a state pins its own rate constant separately from the steady map.
That is asserted and never measured. This measures it.

THE CAVEAT THAT DECIDES THE DRAW TIMES, AND IT COMES FIRST. Over a one-hour
window carnitine barely moves: its own time constant is tau_c = 0.92 d, so it
travels about 4% of the way to its new value in an hour. It is FROZEN on the
scale of the transient, and the clearance that governs the transient is
therefore

    k c_0 m / (K_c + m)          with c_0 the PRE-DOSE carnitine,

and NOT the equilibrium-reduced alpha(m) that every steady-state number in the
paper uses. The two differ in exactly the way that matters here: eliminating c
at equilibrium lets the carnitine pool fall as the marker rises, which is a
positive feedback that FLATTENS the clearance curve; freezing c removes it, so
the frozen curve is steeper and the relaxation faster. Part 1 computes both
and compares them against the local model's 0.81 h, because if the true
relaxation is faster the informative window is earlier than the local model
says and a draw schedule designed on the local model samples too late.

Part 2 then asks the discrimination question on whichever relaxation Part 1
says is right: how many post-dose draws, at what assay sigma, separate the two
structures at 3 sigma.

WHAT IS HELD FIXED. Both structures are matched at the steady state, which is
what makes this a test of the transient and not of the level: the alternative
is calibrated by the same routine Section 6.3's comparison uses, so its
untreated marker sits at the operating point and its treated fall matches the
paper's model's own.

Nothing here records an expected result in advance.

Run:  python audit/t7_structure_transient.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from common.models import ParamsC, g
from design.run import D_AMP, DOSE_INTERVAL, LAMBDA, load

TH_LOC, KM, KE = 3000.0, 100.0, 5.0
LB = load(0.0)
HOUR = 1.0 / 24.0
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
TROUGH = PEAK - D_AMP
M_BAR = 10.0
K_CONJ = 1.0                      # the registered conjugation rate constant


# ------------------------------------------------------- part 1, relaxation
def local_rate(theta, inc, m):
    """The local model's relaxation rate at level m: (theta+Delta) g'(m) + k_e."""
    return (theta + inc) * KM / (KM + m) ** 2 + KE


def parent_frozen_rate(p, theta, inc, m, c):
    """
    d/dm [ theta g(m) + k c m/(K_c+m) ] with c HELD, which is the transient.

    The carnitine state does not move on this timescale, so its own equation
    contributes nothing to the Jacobian the marker relaxes on.
    """
    return ((theta + inc) * KM / (KM + m) ** 2
            + K_CONJ * c * p.Kc / (p.Kc + m) ** 2)


def parent_reduced_rate(p, theta, inc, m, h=1e-6):
    """
    The same derivative with c ELIMINATED, which is what the paper's
    steady-state work uses. Differenced rather than differentiated by hand so
    a slip cannot agree with itself.
    """
    def cl(x):
        u = x / (p.Kc + x)
        c = p.R0 / (K_CONJ * (u + p.rho0 * (1.0 + x / p.Mi)))
        return (theta + inc) * g(x, KM) + K_CONJ * c * u
    return float((cl(m + h) - cl(m - h)) / (2.0 * h))


def parent_state(p, theta):
    """The parent's pre-dose attractor and the carnitine sitting with it."""
    def cl(x):
        u = x / (p.Kc + x)
        c = p.R0 / (K_CONJ * (u + p.rho0 * (1.0 + x / p.Mi)))
        return theta * g(x, KM) + K_CONJ * c * u
    grid = np.geomspace(1e-4, 1e5, 40000)
    peak = float(grid[np.argmax(cl(grid))])
    if cl(peak) <= LB:
        return None, None
    m0 = brentq(lambda x: cl(x) - LB, 1e-9, peak, xtol=1e-12)
    u0 = m0 / (p.Kc + m0)
    c0 = p.R0 / (K_CONJ * (u0 + p.rho0 * (1.0 + m0 / p.Mi)))
    return m0, c0


def carnitine_travel(p, m0, c0, hours=1.0):
    """How far c moves in `hours`, as a fraction of its own journey."""
    tau_c = 1.0 / (K_CONJ * (m0 / (p.Kc + m0)
                             + p.rho0 * (1.0 + m0 / p.Mi)))
    return 1.0 - np.exp(-hours * HOUR / tau_c), tau_c


# ------------------------------------------ part 2, the two marker structures
def flux_traj(te, theta=TH_LOC):
    """The paper's model through one dose, from the pre-dose trough."""
    B = KE * KM + theta + TROUGH - LB
    m0 = (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)
    inc = lambda t: PEAK * np.exp(-LAMBDA * t)
    sol = solve_ivp(lambda t, y: [LB - (theta + inc(t)) * y[0] / (KM + y[0])
                                  - KE * y[0]],
                    (0.0, float(te[-1])), [m0], t_eval=te,
                    method="LSODA", rtol=1e-12, atol=1e-14)
    return sol.y[0]


def algebraic_traj(te, theta=TH_LOC):
    """
    The published readout: the marker is an algebraic function of the enzyme,
    so it tracks the dose instantly and has no relaxation at all.

    Calibrated against the paper's own model at BOTH ends of the transient --
    the pre-dose trough and the post-dose settled level -- so the two agree on
    the steady map and differ only in how they travel between them. That is
    the most conservative calibration available: any looser one would let the
    test succeed on a level difference rather than on a shape.
    """
    B = KE * KM + theta + TROUGH - LB
    m_lo = (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)
    inc = PEAK * np.exp(-LAMBDA * np.asarray(te, float))
    B = KE * KM + theta + inc - LB
    m_qs = (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)
    return m_qs, m_lo


def separation(te, sigma_rel):
    """
    The chi between the two structures over the draws at `te`.

    sigma is a RELATIVE assay error applied at each draw's own level, which is
    how the paper carries the assay term, so a draw taken where the marker is
    low is not credited with impossible precision.
    """
    a = flux_traj(te)
    b, _ = algebraic_traj(te)
    s = sigma_rel * 0.5 * (a + b)
    return float(np.sqrt(np.sum(((a - b) / s) ** 2)))


def structure_numbers():
    """
    What the structure test costs, as generated numbers.

    WHAT IT RULES OUT, which is narrower than "settles which structure a
    record came from" and is the claim the paper makes. The test separates a
    marker that is a STATE, relaxing towards the balance with its own rate
    constant, from one that is ALGEBRAIC in the enzyme and has none. That
    rules out the published alternative and every structure sharing its
    defining feature. It does not separate this paper's structure from a third
    one whose marker is also a state but whose steady map differs, and
    Section 6.1 already records that of four plausible ways to introduce a
    substrate load into the alternative, three leave theta unidentifiable and
    one restores it. A test that rules out a family is worth having; one
    reported as settling the question is not the same thing.

    The two structures are matched at BOTH ends of the transient, the pre-dose
    trough and the post-dose settled level, so they share the steady map
    exactly and are separated by the travel between levels and by nothing
    else. Any looser calibration would let the test succeed on a level
    difference, which the steady record already sees.
    """
    from common.numbers import write_numbers
    p = ParamsC(theta=300.0, Km=KM)
    m0, c0 = parent_state(p, TH_LOC)
    trav, tau_c = carnitine_travel(p, *parent_state(p, p.theta))
    N = {"struct_sigma_3min": separation(np.array([0.05 * HOUR]), 0.05),
         "struct_sigma_15min": separation(np.array([0.25 * HOUR]), 0.05),
         "struct_sigma_30min": separation(np.array([0.50 * HOUR]), 0.05),
         "struct_sigma_1h": separation(np.array([1.00 * HOUR]), 0.05),
         "struct_sigma_2h": separation(np.array([2.00 * HOUR]), 0.05),
         "carn_travel_pct": 100.0 * trav,
         "tau_frozen_h": 24.0 / parent_frozen_rate(p, TH_LOC, 0.0, m0, c0),
         "tau_reduced_h": 24.0 / parent_reduced_rate(p, TH_LOC, 0.0, m0)}
    #  the cheapest schedule reaching 3 sigma at a deliberately poor assay,
    #  because the interesting claim is that it survives one
    for n in range(1, 65):
        if separation(np.geomspace(0.05, 1.0, n) * HOUR, 0.40) >= 3.0:
            N["struct_draws_40pct"] = n
            break
    else:
        N["struct_draws_40pct"] = float("inf")
    fmt = {k: "{:.1f}" for k in N}
    fmt.update({"struct_sigma_1h": "{:.2f}", "struct_sigma_2h": "{:.2f}",
                "struct_draws_40pct": "{:.0f}",
                "tau_frozen_h": "{:.3f}", "tau_reduced_h": "{:.3f}"})
    write_numbers("structure", N, fmt)
    return N


def main():
    print(__doc__)
    bar = "=" * 78
    p = ParamsC(theta=300.0, Km=KM)

    print(bar)
    print("1. THE RELAXATION THE TRANSIENT ACTUALLY SEES")
    print(bar)
    m0, c0 = parent_state(p, p.theta)
    trav, tau_c = carnitine_travel(p, m0, c0)
    print(f"  the parent's pre-dose attractor is m = {m0:.3f} with "
          f"c = {c0:.3f}")
    print(f"  carnitine's own time constant is {tau_c:.3f} d, so in one hour "
          f"it moves {100*trav:.1f}% of its journey: frozen\n")
    print(f"  {'structure':<46} {'rate /d':>9} {'tau (h)':>9}")
    #  the parent block runs at theta = 300 and the local block at 3000, so a
    #  bare comparison of the two confounds the structure with the capacity.
    #  The parent is therefore also evaluated at the LOCAL block's capacity,
    #  at its own attractor there, and that is the row that isolates the
    #  structure. R0 is not re-derived: the calibration that fixes it is a
    #  statement about the alternative channel and does not move with theta.
    m0h, c0h = parent_state(p, TH_LOC)
    rows = [
        ("local model at its own theta = 3000, Delta = 0",
         local_rate(TH_LOC, 0.0, M_BAR)),
        ("local model at the post-dose peak",
         local_rate(TH_LOC, PEAK, 3.3638)),
        ("parent at its own theta = 300, c ELIMINATED (the reduction)",
         parent_reduced_rate(p, p.theta, 0.0, m0)),
        ("parent at its own theta = 300, c FROZEN (the transient)",
         parent_frozen_rate(p, p.theta, 0.0, m0, c0)),
    ]
    if m0h is not None:
        rows += [
            ("parent at theta = 3000, c ELIMINATED (matched capacity)",
             parent_reduced_rate(p, TH_LOC, 0.0, m0h)),
            ("parent at theta = 3000, c FROZEN (matched capacity)",
             parent_frozen_rate(p, TH_LOC, 0.0, m0h, c0h)),
        ]
    for nm, r in rows:
        print(f"  {nm:<46} {r:9.3f} {24.0 / r:9.4f}")
    red = rows[2][1]
    fro = rows[3][1]
    print(f"\n  freezing c makes the parent's relaxation a factor of "
          f"{fro / red:.2f} {'faster' if fro > red else 'slower'}, which is "
          f"the positive")
    print(f"  feedback removed, and at its own capacity the frozen parent "
          f"against the local")
    print(f"  model is a factor of {fro / rows[0][1]:.2f}.")
    if len(rows) > 4:
        mred, mfro = rows[4][1], rows[5][1]
        print(f"  AT MATCHED CAPACITY, which is the comparison that isolates "
              f"the structure:")
        print(f"    parent frozen {24.0/mfro:.4f} h against the local model's "
              f"{24.0/rows[0][1]:.4f} h, a factor of "
              f"{rows[0][1]/mfro:.2f}")
        print(f"    and freezing c is worth {mfro/mred:.2f} there")
        print(f"  So the informative window is "
              f"{'EARLIER' if mfro > rows[0][1] else 'LATER'} than the local "
              f"model implies, by that factor.")

    print("\n" + bar)
    print("2. DISCRIMINATING THE TWO STRUCTURES, AT 3 SIGMA")
    print(bar)
    print("  The published readout is algebraic in the enzyme, so it has no")
    print("  relaxation; the paper's marker is a state. Both are matched at")
    print("  the steady map, so only the travel between levels differs.\n")
    print(f"  {'assay CV':>9} " +
          "".join(f"{f'n={n}':>9}" for n in (2, 4, 8, 16, 32)))
    for cv in (0.01, 0.02, 0.05, 0.10, 0.20):
        row = [f"  {cv:9.0%} "]
        for n in (2, 4, 8, 16, 32):
            te = np.geomspace(0.05, 1.0, n) * HOUR
            row.append(f"{separation(te, cv):9.1f}")
        print("".join(row))
    print("\n  the same, with the draws where the local model would put them")
    print(f"  (0.25 to 6 h, which is what audit/t1 coded):")
    print(f"  {'assay CV':>9} " +
          "".join(f"{f'n={n}':>9}" for n in (2, 4, 8, 16, 32)))
    for cv in (0.05, 0.10, 0.20):
        row = [f"  {cv:9.0%} "]
        for n in (2, 4, 8, 16, 32):
            te = np.geomspace(0.25, 6.0, n) * HOUR
            row.append(f"{separation(te, cv):9.1f}")
        print("".join(row))

    print("\n" + bar)
    print("3. THE CHEAPEST SCHEDULE THAT REACHES 3 SIGMA")
    print(bar)
    for cv in (0.05, 0.10, 0.20, 0.40):
        found = None
        for n in range(1, 65):
            te = np.geomspace(0.05, 1.0, n) * HOUR
            if separation(te, cv) >= 3.0:
                found = n
                break
        print(f"  assay CV {cv:4.0%}: "
              + (f"{found} draw(s) over the first hour"
                 if found else "not reached at 64 draws over the first hour"))
    print()
    single = [(t, separation(np.array([t * HOUR]), 0.05))
              for t in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 6.0)]
    print("  and what ONE draw is worth at a 5% assay, by when it is taken:")
    for t, s in single:
        print(f"    {t:5.2f} h after the dose: {s:6.2f} sigma")

    print("\n" + bar)
    print("R1. RE-CHECK: the two structures agree where they must")
    print(bar)
    te = np.geomspace(0.05, 240.0, 40) * HOUR
    a = flux_traj(te)
    b, m_lo = algebraic_traj(te)
    print(f"  at ten days after the dose the flux balance reads {a[-1]:.6f} "
          f"and the algebraic readout {b[-1]:.6f}")
    print(f"  relative difference {abs(a[-1]-b[-1])/b[-1]:.2e}: they share the "
          f"steady map, so the test is on the transient alone")
    print(f"  both start from the same pre-dose trough, {m_lo:.6f}")

    print("\nR2. RE-CHECK: the frozen-c rate against a direct integration")
    inc0 = 0.0
    def two_state(t, y):
        m, c = y
        conj = K_CONJ * c * m / (p.Kc + m)
        return [LB - (p.theta + inc0) * g(m, KM) - conj,
                p.R0 - conj - p.rho0 * (1.0 + m / p.Mi) * c]
    eps = 1e-4 * m0
    sol = solve_ivp(two_state, (0.0, 0.02), [m0 + eps, c0],
                    dense_output=True, rtol=1e-12, atol=1e-14, method="Radau")
    t1, t2 = 0.002, 0.004
    d1 = sol.sol(t1)[0] - m0
    d2 = sol.sol(t2)[0] - m0
    rate = -np.log(d2 / d1) / (t2 - t1)
    print(f"  decay of a small displacement over the first minutes gives "
          f"{rate:.3f} /d")
    print(f"  against the frozen-c Jacobian's {fro:.3f} /d "
          f"({abs(rate-fro)/fro:.2%} apart)")


if __name__ == "__main__":
    main()
    structure_numbers()
