"""Can the record carry the half-life itself? Fit with ln(half-life) free."""
import numpy as np
from scipy.optimize import least_squares
from core import *
from fit import record

def fit_hl(m, t, y, sig, L, hl_prior=None, hl0=HALF):
    bsd = LOAD_FRAC * L
    def res(p):
        lt, lk, b, lke, lhl = p
        th, K, ke, hl = np.exp([lt, lk, lke, lhl])
        pred = (th + train(t, half=hl)) * m / (K + m) + ke * m + b
        r = [(y - pred) / sig, [b / bsd], [(lke - np.log(KE)) / 0.8]]
        if hl_prior: r.append([(lhl - np.log(hl0)) / hl_prior])
        return np.concatenate(r)
    best = None
    for dt in (0, -1, 1):
        for dh in (0, 0.5, -0.5, 1.0):
            r = least_squares(res, [np.log(TH_NOM) + dt, np.log(KM), 0, np.log(KE), np.log(hl0) + dh],
                              x_scale="jac", max_nfev=5000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
            if best is None or r.cost < best.cost: best = r
    lt, lk, b, lke, lhl = best.x
    return np.exp(lt), np.exp(lhl), best

fine = np.arange(0.0, 730.0, 0.01); ss = (fine > 600) & (fine < 614)
a_pk = train(fine)[ss].max()
truths = {"7 d": lambda t: train(t, half=7.0) * (a_pk / train(fine, half=7.0)[ss].max()),
          "biexp": lambda t: train_biexp(t) * (a_pk / train_biexp(fine)[ss].max()),
          "3.5 d (control)": lambda t: train(t)}
# note: 7 d truth is peak-matched, so its dose amplitude differs from 6000; the
# analyst's train(t, half) keeps D=6000, so a free half-life alone cannot match
# amplitude and shape together. Report both with and without that freedom.
monthly = np.arange(0.0, 730.0, 30.4)
print("monthly design, noiseless, half-life FREE (flat prior); amplitude fixed at D=6000")
for name, f in truths.items():
    out = []
    for s in range(6):
        m, y, sig, Ls = record(TH_NOM, monthly, s, inc_true=f(monthly), noisy=False)
        th, hl, r = fit_hl(m, monthly, y, sig, L_A)
        out.append((th / TH_NOM, hl, np.sqrt(np.mean(r.fun[:monthly.size] ** 2))))
    o = np.median(np.array(out), axis=0)
    print(f"  truth {name:16s} theta_hat/theta {o[0]:.2f}  half-life_hat {o[1]:.2f} d  rms {o[2]:.2f}")

# precision cost of carrying the half-life, CRB at nominal, monthly
def crb_with_hl(theta, times, seed, hl_prior=None):
    sig, _ = sigma_for(theta)
    Ls = loads(times.size, seed)
    inc = train(times)
    m = level(theta, KM, KE, Ls, inc)
    g = m / (KM + m)
    h = 1e-5
    dinc = (train(times, half=HALF * np.exp(h)) - train(times, half=HALF * np.exp(-h))) / (2 * h)
    Q = np.stack([theta * g, -(theta + inc) * KM * m / (KM + m) ** 2, np.ones_like(m), KE * m, dinc * g], 1)
    pri = [0, 0, 1 / (LOAD_FRAC * L_A) ** 2, 1 / 0.8 ** 2, 0 if hl_prior is None else 1 / hl_prior ** 2]
    F = np.diag(pri) + Q.T @ Q / sig ** 2
    return np.sqrt(np.linalg.inv(F)[0, 0])
for sched_name, t in (("monthly", monthly), ("14.5 d", np.arange(0.5, 730, 14.5)), ("2-daily", np.arange(1.0, 730, 2.0))):
    base = np.median([fisher_sd(TH_NOM, t, s) for s in range(12)])
    free = np.median([crb_with_hl(TH_NOM, t, s) for s in range(12)])
    p50 = np.median([crb_with_hl(TH_NOM, t, s, hl_prior=0.2) for s in range(12)])
    print(f"  CRB sd(ln theta) {sched_name:8s}: half-life known {base:.3f} | free {free:.3f} | 20% prior {p50:.3f}")
