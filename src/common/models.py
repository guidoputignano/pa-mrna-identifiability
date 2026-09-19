"""
Shared model definitions for the propionate dynamics.

Three alternative hypotheses about the non-enzymatic clearance term:

    A   dm/dt = s - theta*g(m) - ke*m                   linear, unbounded
    B   dm/dt = s - theta*g(m) - ke*h(m)                both saturating
    C   dm/dt = s - theta*g(m) - k*c*m                  carnitine depleted
        dc/dt = R - k*c*m

with g(m) = m/(Km+m) and h(m) = m/(Ke+m).

Both g and h are MICHAELIS-MENTEN SATURATION FUNCTIONS. In the standard
biochemical form of the Michaelis-Menten law the velocity of a single-substrate
enzymatic reaction is

    v = V_max * [S] / (Km + [S])

with V_max the maximal velocity and Km the Michaelis constant, the substrate
concentration at which the reaction runs at half of V_max. The dictionary used
throughout this repository and the paper is therefore

    theta  ->  V_max of PCC, the maximal enzymatic velocity (patient-specific)
    Km     ->  the Michaelis constant of PCC
    g(m)   ->  the fractional saturation, g in [0, 1)
    ke     ->  V_max of the alternative route (Model B; a rate constant in A)
    Ke     ->  the Michaelis constant of the alternative route (Model B)
    h(m)   ->  the fractional saturation of the alternative route

The same functional form appears under other names in adjacent literatures: the
Monod function in microbial growth kinetics, and the Hill function with
cooperativity coefficient n = 1.

All functions take metabolite level `m` first and parameters as keywords,
so that they can be passed to scipy solvers uniformly.
"""

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------- kinetics
def g(m, Km):
    """
    Michaelis-Menten saturation of the enzymatic channel.

    Fractional saturation of PCC: the fraction of the maximal velocity V_max
    (here theta) currently in use, so that the enzymatic flux is theta*g(m).
    Km is the Michaelis constant, g(Km) = 1/2, and g saturates at 1.
    """
    return m / (Km + m)


def h(m, Ke):
    """
    Michaelis-Menten saturation of the alternative channel, Model B.

    Same functional form as g, with Ke the Michaelis constant of the
    alternative route and ke its maximal velocity, so that the alternative
    flux is ke*h(m). Saturates at 1.
    """
    return m / (Ke + m)


# ---------------------------------------------------------------- parameters
@dataclass
class ParamsA:
    theta: float = 3000.0   # V_max of PCC, maximal enzymatic velocity
                            # (patient-specific)
    Km: float = 100.0       # Michaelis constant of PCC (population)
    ke: float = 5.0         # alternative clearance rate constant (population);
                            # in Model A this channel is linear, not
                            # Michaelis-Menten


@dataclass
class ParamsB:
    theta: float = 3000.0   # V_max of PCC
    Km: float = 100.0       # Michaelis constant of PCC
    ke: float = 400.0       # NOTE: now the V_max of the alternative route,
                            # a maximal flux and not a rate constant
    Ke: float = 300.0       # Michaelis constant of the alternative route


@dataclass
class ParamsC:
    """
    Patient-scale parameters. Note theta is small: the enzyme is defective,
    so the alternative route carries a large share of the flux and its
    behaviour dominates the geometry.
    """
    theta: float = 300.0    # residual V_max of PCC, ~10% of normal
    Km: float = 100.0       # Michaelis constant of PCC
    Kc: float = 100.0       # half-saturation of the conjugation step.
                            # ARBITRARY. There is no measurement behind it and
                            # no reason it should equal Km beyond both being
                            # half-saturations at the plasma level. Its cost is
                            # bounded rather than argued away: over Kc from 25
                            # to 800, a 32-fold range, theta_crit spans 186 to
                            # 120, a factor of 1.56. Section 5.4 reports that.
    Mi: float = 300.0       # metabolite level at which acylcarnitine
                            # competition doubles the non-conjugative loss
                            # rate: delta(m) = delta_0 (1 + m/Mi)
    rho0: float = 0.5       # delta_0 / k, the basal non-conjugative loss
                            # relative to the conjugation rate constant.
                            # Bounded above by Di Donato: the drain in PA is
                            # dominated by propionylcarnitine, so rho0 < 1.
                            # Bounded below only by a DIRECTION, not a number:
                            # free-carnitine excretion is measured falling while
                            # short-chain acylcarnitine excretion rises. A figure
                            # of "about 0.35" was quoted here and is RETRACTED;
                            # re-derived it is 0.195 (REVIEW_LOG.md 0.1). See the
                            # rho_0 row of common/params.py.
    R0: float = 1122.0085   # carnitine replenishment, CONSTANT. Not free: set
                            # so the alternative route alone peaks at 500
                            # against a baseline load of 600, so a patient with
                            # theta = 0 cannot be maintained. That is the
                            # severity of null genotypes, and it is the same
                            # constraint that previously forced R0 down from
                            # 1200. Changing Kc, Mi or rho0 moves the peak and
                            # R0 has to be re-derived; test_models.py checks it.


# ---------------------------------------------------------------- clearance
def clearance_A(m, p: ParamsA):
    return p.theta * g(m, p.Km) + p.ke * m


def clearance_B(m, p: ParamsB):
    return p.theta * g(m, p.Km) + p.ke * h(m, p.Ke)


def alt_flux_C(m, p: ParamsC):
    """
    Carnitine-mediated elimination after adiabatic elimination of the fast
    carnitine state. This is alpha(m) of the paper, Section 5.1.

    One mechanism supplies both features the curve needs.

    Acylcarnitines compete with free carnitine at the renal transporter that
    recovers it from the filtrate, so the fraction recovered falls as the
    metabolite accumulates. That gives the free pool an exit other than
    conjugation, and that exit grows with m:

        delta(m) = delta_0 (1 + m/Mi),    rho(m) = delta(m)/k = rho0 (1 + m/Mi)

    Having a second exit at all is what makes the flux vanish at m = 0: with
    conjugation as the only route out, dc/dt = 0 would force the flux to equal
    the replenishment at every m, so it would not go to zero with the
    metabolite and total clearance would not pass through the origin. The exit
    growing with m is what makes the flux decline rather than merely saturate.

    Solving dc/dt = 0 in rhs_C and substituting gives what is returned here.
    The reduction preserves the constants: R0 and Kc appear in alpha exactly as
    they appear in the two-state system, with no rescaling.
    """
    u = m / (p.Kc + m)
    return p.R0 * u / (u + p.rho0 * (1.0 + m / p.Mi))


def clearance_C_quasistatic(m, p: ParamsC):
    """Total clearance under Model C. Rises, peaks, then declines."""
    return p.theta * g(m, p.Km) + alt_flux_C(m, p)


# ---------------------------------------------------------------- dynamics
def rhs_A(t, y, p: ParamsA, s):
    m = y[0]
    return [s(t) - clearance_A(m, p)]


def rhs_B(t, y, p: ParamsB, s):
    m = y[0]
    return [s(t) - clearance_B(m, p)]


def rhs_C(t, y, p: ParamsC, s, k=1.0):
    """
    The two-state system of the paper, Section 5.1.

    Kept in step with alt_flux_C, which is its quasi-steady reduction:
    setting dc/dt = 0 here returns alt_flux_C(m) for every m, with the same
    R0 and the same Kc. tests/test_models.py checks that directly.

    k is the conjugation rate constant. It cancels out of the reduction, so
    the reduced curve does not depend on it and it is a free scale here.
    """
    m, c = y
    conj = k * c * m / (p.Kc + m)
    delta = p.rho0 * k * (1.0 + m / p.Mi)
    return [s(t) - p.theta * g(m, p.Km) - conj,
            p.R0 - conj - delta * c]


# ---------------------------------------------------------------- equilibria
def steady_state(clearance_fn, p, s, lo=1e-12, hi=1e7):
    """Metabolite level at which clearance balances the load, if it exists."""
    f = lambda m: clearance_fn(m, p) - s
    if f(hi) < 0:
        return np.inf          # load exceeds the ceiling: no equilibrium
    try:
        return brentq(f, lo, hi)
    except ValueError:
        return np.nan


def equilibria_nonmonotone(clearance_fn, p, s, grid=None):
    """
    All equilibria of a system whose clearance is non-monotone.
    Returns (values, stability) with stability True for stable.
    """
    if grid is None:
        grid = np.linspace(1e-6, 2000.0, 20000)
    f = np.array([clearance_fn(m, p) for m in grid]) - s
    roots, stab = [], []
    for i in range(len(grid) - 1):
        if f[i] * f[i + 1] < 0:
            r = brentq(lambda m: clearance_fn(m, p) - s, grid[i], grid[i + 1])
            # dm/dt = s - C(m); stable iff d/dm (s - C) < 0, i.e. C increasing
            eps = max(r * 1e-4, 1e-6)
            slope = (clearance_fn(r + eps, p) - clearance_fn(r - eps, p)) / (2 * eps)
            roots.append(r)
            stab.append(slope > 0)
    return np.array(roots), np.array(stab)


# ---------------------------------------------------------------- sensitivity
def sensitivity_A(m, p: ParamsA):
    """d(clearance)/d(V_max, Michaelis constant, ke) at metabolite level m."""
    return np.array([
        g(m, p.Km),
        -p.theta * m / (p.Km + m) ** 2,
        m,
    ])


def sensitivity_B(m, p: ParamsB):
    """d(clearance)/d(theta, Km, ke, Ke) at metabolite level m."""
    return np.array([
        g(m, p.Km),
        -p.theta * m / (p.Km + m) ** 2,
        h(m, p.Ke),
        -p.ke * m / (p.Ke + m) ** 2,
    ])


def condition_number(sens_fn, p, levels):
    """Condition number of the column-normalised sensitivity matrix."""
    S = np.array([sens_fn(m, p) for m in levels])
    S = S / np.linalg.norm(S, axis=0)
    return np.linalg.cond(S)


# ---------------------------------------------------------------- threshold
def theta_crit(s_peak, m_tol, p: ParamsA):
    """
    Critical capacity, MODEL A: the value of V_max (theta) at which, under
    peak load, the equilibrium sits exactly at the tolerable metabolite level.

    The argument is m_tol, not m_bar. They were the same name and they are
    different quantities: m_bar is the point the affine fit of Proposition 1
    is taken at, and the sample mean of the observed band; m_tol is a clinical
    convention about what level is tolerable. Equation (17) mixed them.

    Model B's alternative channel saturates, so its critical capacity is
    (s_peak - ke*h(m_tol))/g(m_tol) and not this. That is what
    model_B/run.py::fig_thetacrit computes.

    This definition imports two quantities that are not measured, s_peak and
    m_tol. See the note for why that matters.
    """
    return (s_peak - p.ke * m_tol) / g(m_tol, p.Km)
