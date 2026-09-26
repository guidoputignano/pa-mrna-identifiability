"""MAP fitter for the local model, independent of the repository's estimators."""
import numpy as np
from scipy.optimize import least_squares
from core import *


def fit(m, inc, y, sigma, L, ke_prior=0.8, b_sd=None, kappa_cv=None,
        theta0=TH_NOM, starts=None, ke0=KE):
    """
    Minimise sum((y - pred)/sigma)^2 + priors over
        p = (ln theta, ln Km, b, ln ke [, ln kappa])
    b_sd: None -> 5% of L; 0 -> b fixed at 0; ke_prior 0 -> ke fixed.
    Returns dict with theta_hat etc and the chi2 of the data term.
    """
    b_sd = LOAD_FRAC * L if b_sd is None else b_sd
    free_b, free_ke = b_sd > 0, ke_prior > 0
    use_k = kappa_cv is not None and kappa_cv > 0

    def unpack(p):
        i = 0
        lt, lk = p[0], p[1]; i = 2
        b = p[i] if free_b else 0.0; i += free_b
        lke = p[i] if free_ke else np.log(ke0); i += free_ke
        lka = p[i] if use_k else 0.0
        return lt, lk, b, lke, lka

    def res(p):
        lt, lk, b, lke, lka = unpack(p)
        th, K, ke, ka = np.exp(lt), np.exp(lk), np.exp(lke), np.exp(lka)
        pred = (th + ka * inc) * m / (K + m) + ke * m + b
        r = [(y - pred) / sigma]
        if free_b: r.append([b / b_sd])
        if free_ke: r.append([(lke - np.log(ke0)) / ke_prior])
        if use_k: r.append([lka / kappa_cv])
        return np.concatenate(r)

    if starts is None:
        starts = []
        for dt in (0.0, -1.0, 1.0, -2.5, 2.0):
            for dk in (0.0, -1.0, 1.0):
                starts.append((np.log(theta0) + dt, np.log(KM) + dk))
    best = None
    for lt0, lk0 in starts:
        p0 = [lt0, lk0] + ([0.0] if free_b else []) + ([np.log(ke0)] if free_ke else []) + ([0.0] if use_k else [])
        try:
            r = least_squares(res, p0, method="trf", x_scale="jac", max_nfev=4000,
                              xtol=1e-12, ftol=1e-12, gtol=1e-12)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    lt, lk, b, lke, lka = unpack(best.x)
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    pred = (th + np.exp(lka) * inc) * m / (K + m) + ke * m + b
    rms = float(np.sqrt(np.mean(((y - pred) / sigma) ** 2)))
    return dict(theta=th, Km=K, ke=ke, b=b, kappa=np.exp(lka), rms=rms, cost=best.cost)


def record(theta, times, seed, inc_true=None, noisy=True, L=L_A, Km=KM, ke=KE,
           b_true=None, theta_t=None, load_cv=LOAD_CV):
    """
    Simulate a record. inc_true: the TRUE increment at the draw times (defaults
    to the single-exponential train). theta_t: optional time-varying theta.
    b_true: optional time-varying unrecorded load (enters the true balance but
    not the recorded ordinate).
    Returns m (marker), y (recorded load), sigma.
    """
    sig, _ = sigma_for(theta if theta_t is None else float(np.mean(theta_t)), Km, ke, L)
    Ls = loads(times.size, seed, L, load_cv)
    inc = train(times) if inc_true is None else inc_true
    th = theta if theta_t is None else theta_t
    btrue = 0.0 if b_true is None else b_true
    m = level(th, Km, ke, Ls + btrue, inc)       # the body sees recorded + unrecorded
    y = Ls.copy()                                 # the analyst sees only the recorded load
    if noisy:
        rng = np.random.default_rng(777 + seed)
        y = y + rng.normal(0.0, sig, times.size)
    return m, y, sig, Ls
