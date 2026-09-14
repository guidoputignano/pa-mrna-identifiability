"""
The excitation condition, and the design quantities that follow from it.

PARAMETRISATION, PINNED BEFORE ANYTHING IS PROVED.

    dm/dt = L(t) + L_0 - (theta + Delta(t)) g(m) - k_e m ,   g(m) = m/(K_m + m)

theta IS the residual maximal velocity of the enzymatic channel, absorbed: the
flux is theta*g(m) and not V_max*theta*g(m). Carrying a separate V_max adds a
parameter that the marker cannot separate from theta at any excursion, since
only the product multiplies g. Splitting them needs an external anchor, and the
1-13C-propionate breath test is what such an anchor would look like.

L(t) is the recorded component of the propiogenic load and L_0 the unrecorded
one, carried as an unknown constant. Delta(t) is the therapeutic increment,
known from the dosing schedule.

THE CONDITION.

Both ends of the saturation curve are degenerate, and they fail differently.

  m << K_m.  g(m) -> m/K_m and total clearance collapses to
             (theta/K_m + k_e) m. Three of the four sensitivity directions
             become proportional to m, so only the lumped ratio is identifiable
             and K_m and k_e are lost separately. The clearance is still
             constrained; what is lost is its decomposition.

  m >> K_m.  g(m) -> 1 and clearance becomes theta + k_e m, so the dynamics are
             dm/dt = L + L_0 - theta - k_e m. theta enters only through
             (L_0 - theta). WITH L_0 UNKNOWN, THETA AND THE LOAD ARE
             CONFOUNDED and cannot be separated from the marker trajectory
             alone. This is the worse failure: at the low end a nuisance
             parameter is lost, at the high end the target is.

Identifiability therefore comes from CURVATURE, and the requirement is to
TRAVERSE THE NEIGHBOURHOOD OF K_m, not to exceed it.

Two consequences, and the second is the one a design turns on.

  1. The theta/K_m excitation has a closed form and an optimum. The ratio of
     the K_m sensitivity to the theta sensitivity is -(1 - g(m)) = -phi(m) with
     phi(m) = K_m/(K_m + m), so what separates them over a band [lo, hi] is

         E = phi(lo) - phi(hi) = K_m (hi - lo) / ((K_m + lo)(K_m + hi)) .

     At a fixed fold range r = hi/lo, E is maximised at lo = K_m/sqrt(r),
     hi = K_m sqrt(r), that is when sqrt(lo*hi) = K_m exactly, and there

         E_max = (sqrt(r) - 1) / (sqrt(r) + 1) ,

     so a required E needs r >= ((1+E)/(1-E))^2. The band must STRADDLE K_m in
     geometric mean. A band of the same fold range sitting a decade below K_m
     delivers a quarter of the excitation.

  2. Because L_0 is unknown, the full problem's optimum sits a little above
     K_m rather than on it -- the load-bias degeneracy pulls up where the K_m
     degeneracy pulls down -- and, crucially, it is INTERIOR. sd(ln theta) has
     a genuine minimum in the operating point, which it does not have when the
     load is known. So load-recording precision is not a mechanical detail: it
     is what creates the upper obstruction, and without it saturation would be
     the best place to measure rather than the worst.

The excursion this asks for can be bought two ways, and they are not equally
available: a wide band from the operating point, or a narrow band held at a
much higher level. The second is far cheaper in width and far more dangerous
in level, and the clinical constraint is on the level.
"""

import numpy as np
from scipy.optimize import brentq


def phi(m, Km):
    """K_m/(K_m + m) = 1 - g(m). The shape the K_m direction is carried by."""
    return Km / (Km + m)


def excitation(lo, hi, Km):
    """E over a band. Zero in both asymptotes, maximal straddling K_m."""
    return float(phi(lo, Km) - phi(hi, Km))


def excitation_max(r):
    """The most any band of fold range r can deliver, at its optimal placement."""
    s = np.sqrt(r)
    return float((s - 1.0) / (s + 1.0))


def fold_range_for(E):
    """The fold range a required excitation needs, at optimal placement."""
    return float(((1.0 + E) / (1.0 - E)) ** 2)


def optimal_band(r, Km):
    """The band of fold range r that maximises E: geometric mean exactly K_m."""
    s = np.sqrt(r)
    return Km / s, Km * s


def sensitivities(m, theta, Km, ke, load_bias=True):
    """
    d(clearance)/d(ln theta, ln K_m, ln k_e[, load bias]) at level m.

    The last column is the sensitivity of the ORDINATE to an unrecorded
    constant in the load, which is 1 at every level. Including it is what makes
    the saturated end degenerate; dropping it makes saturation look like the
    best place to measure.
    """
    m = np.asarray(m, dtype=float)
    cols = [theta * m / (Km + m),
            -theta * Km * m / (Km + m) ** 2,
            ke * m]
    if load_bias:
        cols.append(np.ones_like(m))
    return np.stack(cols, axis=-1)


def information(levels, theta, Km, ke, sigma=1.0, load_bias=True):
    """Fisher information of a steady-state design, in log-parameters."""
    Q = sensitivities(np.asarray(levels, dtype=float), theta, Km, ke, load_bias)
    return (Q.T @ Q) / sigma ** 2


def sd_log_theta(levels, theta, Km, ke, sigma=1.0, load_bias=True):
    """sd(ln theta_hat) for a steady-state design. inf when the design is singular."""
    F = information(levels, theta, Km, ke, sigma, load_bias)
    try:
        v = np.linalg.inv(F)[0, 0]
    except np.linalg.LinAlgError:
        return np.inf
    return float(np.sqrt(v)) if v > 0 else np.inf


def band_levels(centre, r, n=9):
    """n levels spread geometrically over a band of fold range r about `centre`."""
    s = np.sqrt(r)
    return np.geomspace(centre / s, centre * s, int(n))


def fold_range_for_precision(target, centre, theta, Km, ke, sigma=1.0, n=9,
                             r_max=1e7):
    """
    The fold range a band about `centre` needs to reach `target` on sd(ln theta).

    Returns None when no width about that centre reaches it, which is the
    common case once the centre is far from K_m: the band runs into one
    degeneracy or the other before it gets wide enough.
    """
    f = lambda lr: sd_log_theta(band_levels(centre, np.exp(lr), n),
                                theta, Km, ke, sigma) - target
    lo = 1e-3
    while lo < np.log(r_max) and not np.isfinite(f(lo)):
        lo *= 1.5
    if not np.isfinite(f(lo)) or f(lo) < 0:
        return 1.0 if np.isfinite(f(lo)) and f(lo) < 0 else None
    hi = np.log(r_max)
    if f(hi) > 0:
        return None
    return float(np.exp(brentq(f, lo, hi)))


def top_of_excursion_for_precision(target, m_bar, theta, Km, ke, sigma=1.0,
                                   n=9, cap=1e6):
    """
    How high the marker must go, starting the band at the operating point.

    This is the form the clinic can be asked for: the band runs from m_bar
    upward, and what is reported is its top. Returns None when unreachable
    below `cap`.
    """
    f = lambda hi: sd_log_theta(np.geomspace(m_bar, hi, n),
                                theta, Km, ke, sigma) - target
    lo = m_bar * 1.01
    while lo < cap and not np.isfinite(f(lo)):
        lo *= 1.5
    if not np.isfinite(f(lo)):
        return None
    if f(lo) < 0:
        return float(lo)
    if f(cap) > 0:
        return None
    return float(brentq(f, lo, cap))
