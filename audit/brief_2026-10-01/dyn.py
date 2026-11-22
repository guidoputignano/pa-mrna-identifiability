"""
The dynamic model for the washout and drug-forcing experiments (Tasks 6, 7).

    dm/dt = u(t) + L0 - (theta + Delta(t)) g(m) - k_e m,   g = m/(K_m + m)

u(t)      the KNOWN load (diet, or an infused load), piecewise constant
L0        the UNRECORDED load (gut, catabolism, odd-chain fatty acids): unknown
Delta(t)  a known increment on the enzyme channel (Task 7's drug forcing),
          kappa * delta(t) with delta known and kappa optionally unknown
m(0-)     the steady state under u(0-), Delta(0-): the experiment starts at rest

Parameters p = (ln theta, ln k_e, ln K_m, L0[, ln kappa]). Observation
ln y_k = ln m(t_k) + e_k, e_k ~ N(0, cv^2), cv = 5% (the repository's assay).
Time in DAYS, as in the repository (k_e = 5 /d, so 1/k_e = 4.8 h).

ASSUMPTIONS THE REPOSITORY DOES NOT SUPPLY, stated once:
  * the split of the 5% patient's load L = 322.73 into a known part u and an
    unrecorded part L0. Default L0 = 0.5 L. Swept in Task 7.
  * that L0 stays constant through the experiment. The repository withdraws a
    fast for exactly this reason (u15: fasting triggers catabolism, so the
    unrecorded load RISES); a glucose-protected diet withdrawal is its
    candidate substitute (u17), and constancy of L0 under it is unverified.
"""
import numpy as np
from scipy.integrate import solve_ivp

import core as C

TH, KM, KE = C.TH_5, C.KM0, C.KE
LTOT = float(C.N4.anchor(TH, KM, KE, "A")[0])
CV = 0.05
NP = 4


class Experiment:
    """A load protocol u(t) (steps) plus an optional increment kappa*delta(t)."""

    def __init__(self, steps, L0=0.5 * LTOT, delta=None, kappa_free=False):
        # steps: list of (t_start, u) with the first at t <= 0 giving u(0-)
        self.steps = sorted(steps)
        self.L0 = L0
        self.delta = delta if delta is not None else (lambda t: 0.0)
        self.kappa_free = kappa_free

    def u(self, t):
        v = self.steps[0][1]
        for ts, us in self.steps:
            if t >= ts:
                v = us
        return v

    def breaks(self, T):
        return sorted({ts for ts, _ in self.steps if 0 < ts < T})


def params(theta=TH, ke=KE, Km=KM, L0=None, kappa=1.0, ex=None):
    L0 = ex.L0 if L0 is None else L0
    p = [np.log(theta), np.log(ke), np.log(Km), L0]
    if ex is not None and ex.kappa_free:
        p.append(np.log(kappa))
    return np.array(p, float)


def unpack(p):
    th, ke, K, L0 = np.exp(p[0]), np.exp(p[1]), np.exp(p[2]), p[3]
    kap = np.exp(p[4]) if len(p) > 4 else 1.0
    return th, ke, K, L0, kap


def steady(p, ex, t=0.0):
    th, ke, K, L0, kap = unpack(p)
    D = kap * ex.delta(t)
    return float(C.m_ss(th + D, K, ke, ex.u(-1e-9 if t == 0 else t) + L0))


def simulate(p, ex, times, sens=True, rtol=1e-9):
    """m(t_k) and dm/dp at the sample times, by forward sensitivities."""
    th, ke, K, L0, kap = unpack(p)
    npar = len(p)
    m0 = steady(p, ex)
    # initial sensitivities: the steady state moves with every parameter
    D0 = kap * ex.delta(0.0)
    g0 = m0 / (K + m0)
    Cp0 = (th + D0) * K / (K + m0) ** 2 + ke
    dC = [th * g0, ke * m0, -(th + D0) * K * m0 / (K + m0) ** 2, -1.0]
    if npar > 4:
        dC.append(D0 * g0)
    s0 = [-d / Cp0 for d in dC]

    def rhs(t, y):
        m = y[0]
        D = kap * ex.delta(t)
        g = m / (K + m)
        f = ex.u(t) + L0 - (th + D) * g - ke * m
        if not sens:
            return [f]
        Cp = (th + D) * K / (K + m) ** 2 + ke
        dCdp = [th * g, ke * m, -(th + D) * K * m / (K + m) ** 2, -1.0]
        if npar > 4:
            dCdp.append(D * g)
        return [f] + [-Cp * y[1 + j] - dCdp[j] for j in range(npar)]

    times = np.asarray(times, float)
    T = times[-1]
    edges = [0.0] + ex.breaks(T) + [T]
    y = np.array([m0] + (s0 if sens else []), float)
    cols = {}
    for a, b in zip(edges[:-1], edges[1:]):
        last = b == T
        sel = times[(times >= a) & ((times <= b) if last else (times < b))]
        te = np.unique(np.r_[sel, b])
        sol = solve_ivp(rhs, (a, b), y, t_eval=te, method="LSODA",
                        rtol=rtol, atol=1e-12)
        for k, t in enumerate(sol.t):
            if np.any(np.isclose(sel, t, rtol=0, atol=1e-12)):
                cols[t] = sol.y[:, k]
        y = sol.y[:, -1]
    if 0.0 in times and 0.0 not in cols:
        cols[0.0] = np.array([m0] + (s0 if sens else []))
    yy = np.array([cols[t] for t in sorted(cols)]).T
    m = yy[0]
    S = yy[1:].T if sens else None
    return m, S


def fim(p, ex, times, cv=CV, priors=None):
    """Information in p from ln y_k = ln m_k + N(0, cv^2), plus diagonal priors."""
    m, S = simulate(p, ex, times)
    R = S / m[:, None]
    F = R.T @ R / cv ** 2
    if priors:
        for j, sd in priors.items():
            if sd:
                F[j, j] += sd ** -2
    return F


def sd_lnth(F):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def lam(p, ex, m, t=0.0):
    """The relaxation eigenvalue C'(m) = (theta + Delta) g'(m) + k_e."""
    th, ke, K, L0, kap = unpack(p)
    return (th + kap * ex.delta(t)) * K / (K + m) ** 2 + ke
