"""
Independent re-implementation of the paper's local model and its information
matrix. Written from the equations in the paper and the constants in
src/common/params.py; imports nothing from the repository.

Balance (local model, quasi-steady):
    s = (theta + kappa*Delta(t)) * m/(Km+m) + ke*m + b
Ordinate s is the recorded load; noise sigma in load units combines a 5%
load-record error and a 5% assay CV propagated through the clearance slope.
Parameters for the Fisher matrix: (ln theta, ln Km, b, ln ke) [+ ln kappa].
Priors: b ~ N(0, (0.05 L)^2), ln ke ~ N(ln ke0, 0.8^2); none on theta, Km.
"""
import numpy as np

TH_NOM, KM, KE, MBAR = 3000.0, 100.0, 5.0, 10.0
D_AMP, TAU, HALF = 6000.0, 14.0, 3.5
LOAD_FRAC, ASSAY_CV, LOAD_CV = 0.05, 0.05, 0.10

# convention A: the diet is the one a normal-capacity patient balances at m=10
L_A = TH_NOM * MBAR / (KM + MBAR) + KE * MBAR


def train(t, D=D_AMP, half=HALF, tau=TAU, t0=0.0):
    """Single-exponential dose train, doses at t0, t0+tau, ... (closed form)."""
    t = np.asarray(t, float)
    lam = np.log(2.0) / half
    k = np.floor((t - t0) / tau)
    k = np.where(t >= t0, k, -1)
    n = k + 1                                  # doses given so far
    since = (t - t0) - k * tau                 # time since last dose
    r = np.exp(-lam * tau)
    geo = np.where(n > 0, (1 - r ** n) / (1 - r), 0.0)
    return np.where(n > 0, D * np.exp(-lam * since) * geo, 0.0)


def train_biexp(t, D=1.0, tau=TAU, kdeg=0.0075, kq=0.00474):
    """
    Two-compartment protein (Attarwala Suppl. S1(B)), impulse input into PCC:
        PCC' = -(kq+kdeg) PCC + kq PCP,  PCP' = kq PCC - kq PCP   (per hour)
    Returns the PCC response to a train of unit impulses, scaled by D.
    """
    A = np.array([[-(kq + kdeg), kq], [kq, -kq]]) * 24.0     # per day
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.array([1.0, 0.0]))
    coef = V[0, :] * c                                        # PCC(t)=sum coef e^{w t}
    t = np.asarray(t, float)
    out = np.zeros_like(t)
    kmax = int(np.ceil(t.max() / tau)) + 1
    for k in range(kmax):
        on = t >= k * tau
        dt = t[on] - k * tau
        out[on] += np.real(coef[0] * np.exp(w[0] * dt) + coef[1] * np.exp(w[1] * dt))
    return D * out


def level(theta, Km, ke, L, inc):
    """Positive root of (theta+inc) m/(Km+m) + ke m = L."""
    B = ke * Km + theta + inc - L
    return (-B + np.sqrt(B * B + 4.0 * ke * L * Km)) / (2.0 * ke)


def sigma_for(theta, Km=KM, ke=KE, L=L_A):
    m = level(theta, Km, ke, L, 0.0)
    slope = theta * Km / (Km + m) ** 2 + ke
    return float(np.hypot(LOAD_FRAC * L, slope * ASSAY_CV * m)), float(m)


def loads(n, seed, L=L_A, cv=LOAD_CV):
    if cv == 0:
        return np.full(n, L)
    rng = np.random.default_rng(10_000 + seed)
    s2 = np.log1p(cv * cv)
    return L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))


def fisher_sd(theta, times, seed, ke_prior=0.8, b_prior=None, D=D_AMP,
              kappa_cv=None, sigma_scale=1.0, Km=KM, ke=KE, L=L_A, load_cv=LOAD_CV):
    """sd(ln theta_hat) from the Fisher matrix at the truth, with the priors."""
    sig, _ = sigma_for(theta, Km, ke, L)
    sig *= sigma_scale
    Ls = loads(times.size, seed, L, load_cv)
    inc = train(times, D)
    m = level(theta, Km, ke, Ls, inc)
    g = m / (Km + m)
    cols = [theta * g, -(theta + inc) * Km * m / (Km + m) ** 2,
            np.ones_like(m), ke * m]
    # a prior of 0 (or b_prior == "known") means the parameter is known exactly:
    # drop its column rather than add an infinite prior
    bp = (1.0 / (LOAD_FRAC * L) ** 2) if b_prior is None else b_prior
    pri = [0.0, 0.0, bp, (0.0 if ke_prior is None else
                          (np.inf if ke_prior == 0 else 1.0 / ke_prior ** 2))]
    if kappa_cv is not None:
        cols.append(inc * g)
        pri.append(np.inf if kappa_cv == 0 else 1.0 / kappa_cv ** 2)
    keep = [i for i, p in enumerate(pri) if not np.isinf(p)]
    Q = np.stack([cols[i] for i in keep], axis=1)
    F = np.diag([pri[i] for i in keep]) + Q.T @ Q / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def med_sd(theta, times, seeds=range(24), **kw):
    v = np.array([fisher_sd(theta, times, s, **kw) for s in seeds])
    v = v[np.isfinite(v)]
    return (float(np.median(v)), float(v.min()), float(v.max())) if v.size else (np.inf,) * 3
