import numpy as np
from scipy.optimize import least_squares
from core import *
from fit import record
fine = np.arange(0.0, 730.0, 0.01); ss = (fine > 600) & (fine < 614)
a_pk = train(fine)[ss].max()
monthly = np.arange(0.0, 730.0, 30.4)
truths = {"7 d": lambda t: train(t, half=7.0) * (a_pk / train(fine, half=7.0)[ss].max()),
          "biexp": lambda t: train_biexp(t) * (a_pk / train_biexp(fine)[ss].max())}
def fit2(m, t, y, sig):
    bsd = LOAD_FRAC * L_A
    def res(p):
        lt, lk, b, lke, lhl, lka = p
        th, K, ke, hl, ka = np.exp([lt, lk, lke, lhl, lka])
        pred = (th + ka * train(t, half=hl)) * m / (K + m) + ke * m + b
        return np.concatenate([(y - pred) / sig, [b / bsd], [(lke - np.log(KE)) / 0.8]])
    best = None
    for dt in (0, -1, 1):
        for dh in (0, 0.5, -0.5):
            r = least_squares(res, [np.log(3000) + dt, np.log(KM), 0, np.log(KE), np.log(3.5) + dh, 0],
                              x_scale="jac", max_nfev=6000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
            if best is None or r.cost < best.cost: best = r
    return best
print("half-life AND amplitude free (flat), monthly, noiseless")
for name, f in truths.items():
    out = []
    for s in range(6):
        m, y, sig, Ls = record(TH_NOM, monthly, s, inc_true=f(monthly), noisy=False)
        r = fit2(m, monthly, y, sig)
        out.append([np.exp(r.x[0]) / 3000, np.exp(r.x[4]), np.exp(r.x[5]), np.sqrt(np.mean(r.fun[:25] ** 2))])
    o = np.median(np.array(out), 0)
    print(f"  truth {name:6s} theta_hat/theta {o[0]:.2f}  hl {o[1]:.2f} d  kappa {o[2]:.2f}  rms {o[3]:.2f}")
