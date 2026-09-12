"""
The dynamic model, its sensitivity system, and the observability Gramian.

The note's identifiability results are all algebraic: they are statements
about the equilibrium condition s = C(m), read at two or more steady states.
This module is the same question asked of the trajectory, which is what a
record actually contains once the sampling is fast enough to catch the patient
away from steady state.

    dm/dt = L(t) - (theta + Delta(t)) g(m) - (ke + kappa u(t)) m

with Delta the therapeutic increment on the enzymatic channel, driven by a
known dosing schedule, and kappa*u the carnitine increment on the alternative
channel, driven by the recorded dose. L is the propiogenic load.

Relative sensitivities q_p = p * dm/dp obey

    d q_theta / dt = -theta g(m)                 - J(t) q_theta
    d q_Km    / dt = -(theta + Delta(t)) Km g_Km - J(t) q_Km
    d q_ke    / dt = -ke m                       - J(t) q_ke
    J(t) = (theta + Delta(t)) g'(m) + ke + kappa u(t)

and the Gramian is G_ij = int q_i q_j dt, dimensionless because the
sensitivities are relative.

WHAT THE GRAMIAN SAYS THAT THE DETERMINANT DOES NOT.

Proposition 2 says the 2x2 system for (theta, ke) has the same coefficient
matrix whichever channel was perturbed, so its determinant depends on the
excursion and not on the channel. That is exact, and it generalises: at a
matched trajectory the two sensitivity FORCINGS above, -theta g(m) and -ke m,
are functions of the state alone and are therefore channel-blind identically.
The only channel dependence anywhere in the system is the shared decay rate,
and at a matched trajectory

    J_therapy - J_carnitine = -kappa u g(m) <= 0

so the therapeutic input decays the sensitivities more slowly and returns a
LARGER Gramian. Posing the question on carnitine understates what the therapy
would deliver, by a few per cent over every excursion tried here.

It stops generalising the moment K_m joins the parameter vector. The K_m
forcing carries theta + Delta(t), and that factor is modulated by the
therapeutic input and constant under carnitine. So the ratio of the K_m
forcing to the theta forcing is -Km/(Km+m) under carnitine, nearly constant
while the patient sits far below K_m, and the two sensitivities are then
collinear to four decimal places. Under the therapy the same ratio carries
1 + Delta(t)/theta, which the schedule moves at will.

The consequence is the paper's own thesis in another form. Carnitine buys the
excursion; it does not buy the curvature. Curvature is bought either by a
perturbation that multiplies the enzymatic term -- which only the therapy does
-- or by carrying the patient into the neighbourhood of K_m, which only the
load does, and which is the excursion the clinic exists to prevent.
"""

import numpy as np
from scipy.integrate import solve_ivp


ZERO = staticmethod(lambda t: 0.0)


def _zero(t):
    return 0.0


def simulate(theta, Km, ke, load, delta=None, kappa_u=None, m0=None,
             T=40.0, npts=8001, rtol=1e-10, atol=1e-13):
    """
    Integrate the state and its three relative sensitivities.

    `load`, `delta` and `kappa_u` are callables of time. `m0` defaults to the
    steady state at load(0) with no perturbation, so the trajectory starts at
    rest and everything that follows is the response to the input.

    Returns (t, m, Q) with Q the 3 x npts array of relative sensitivities in
    the order (theta, Km, ke).
    """
    delta = _zero if delta is None else delta
    kappa_u = _zero if kappa_u is None else kappa_u

    def g(m):
        return m / (Km + m)

    def gp(m):
        return Km / (Km + m) ** 2

    def dg_dKm(m):
        return -m / (Km + m) ** 2

    if m0 is None:
        from scipy.optimize import brentq
        L0 = load(0.0)
        m0 = brentq(lambda x: theta * g(x) + ke * x - L0, 1e-12, 1e9)

    def rhs(t, y):
        m, qt, qK, qk = y
        m = max(m, 1e-14)
        th = theta + delta(t)
        k = ke + kappa_u(t)
        J = th * gp(m) + k
        return [load(t) - th * g(m) - k * m,
                -theta * g(m) - J * qt,
                -th * Km * dg_dKm(m) - J * qK,
                -ke * m - J * qk]

    sol = solve_ivp(rhs, (0.0, T), [m0, 0.0, 0.0, 0.0], dense_output=True,
                    rtol=rtol, atol=atol, method="LSODA")
    tt = np.linspace(0.0, T, npts)
    Y = sol.sol(tt)
    return tt, Y[0], Y[1:]


def gramian(tt, Q, idx=(0, 1, 2)):
    """G_ij = int q_i q_j dt over the coordinates in `idx`."""
    S = Q[list(idx)]
    n = len(idx)
    return np.array([[np.trapezoid(S[i] * S[j], tt) for j in range(n)]
                     for i in range(n)])


def collinearity(G, i=0, j=1):
    """|cos| between two sensitivity directions, as the Gramian sees them."""
    return float(abs(G[i, j]) / np.sqrt(G[i, i] * G[j, j]))


def sd_log_theta(G):
    """sd(ln theta_hat) per unit output noise density, from the full Gramian."""
    return float(np.sqrt(np.linalg.inv(G)[0, 0]))


def step(amp, t_on=5.0, t_off=25.0):
    """A single sustained increment, the shape the note's two steady states assume."""
    return lambda t: amp if t_on <= t <= t_off else 0.0


def pulses(amp, period=14.0, duty=0.4, t0=2.0):
    """
    A repeating dose, the shape the mRNA-3927 schedule has.

    The default period is the biweekly interval selected for that trial, so
    the inter-dose recovery this produces is the excursion a real schedule
    makes available without anyone being asked to decompensate.
    """
    return lambda t: amp if t >= t0 and ((t - t0) % period) < duty * period else 0.0
