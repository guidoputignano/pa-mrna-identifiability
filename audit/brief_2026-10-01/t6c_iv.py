"""
Item 4. One fix for the finite-difference (FD) estimator before writing it off.

t6's FD regression (t6_estimators.est_fd), row k = interval [t_k, t_{k+1}]:
    z_k = psi_k m_k   on   X_k = [1, m_k, -m_k^2, -psi_k],   psi_k = (y_{k+1} - y_k)/tau - u_k
    lumped p = (L0 K_m, L0 - theta - k_e K_m, k_e, K_m);  L0 = p1/p4, theta = p1/p4 - p2 - p3 p4.
m_k is y_k (left point) or (y_k + y_{k+1})/2 (midpoint). Row k's error involves the noise of
samples k and k+1 only.

THE NOISE MODEL (t6's benches): y_k = m_k exp(e_k), e_k ~ N(0, cv^2), independent. Exactly,
    E[y_k^i y_{k+1}^j] = m_k^i m_{k+1}^j exp((i^2 + j^2) cv^2 / 2).
Row by row, every entry of X'X and X'z, and every entry of X_k and z_k, is a polynomial of
degree <= 4 in (y_k, y_{k+1}). Replacing each monomial y_k^i y_{k+1}^j by
y_k^i y_{k+1}^j exp(-(i^2 + j^2) cv^2 / 2) gives an unbiased estimate of its noiseless value.
That is the corrected-score construction, and it is exact to all orders for lognormal noise.
The script builds the polynomials symbolically, so nothing is truncated. To second order,
with m = m_k, m' = m_{k+1} and the left point:
    E[psi^2]     - psi0^2     = cv^2 (2 m^2 - 2 m m' + 2 m'^2) / tau^2 - cv^2 u (m' - m) / tau
    E[psi^2 m]   - psi0^2 m   = cv^2 (5/2 m m'^2 - 5 m^2 m' + 9/2 m^3) / tau^2 + O(cv^2 u m^2 / tau)
    E[psi m]     - psi0 m     = cv^2 (m m' - 2 m^2) / tau - cv^2 u m / 2
The psi-psi entry of X'X and the psi-z (psi * psi m) entry of X'z carry the large terms,
~2 cv^2 m^2/tau^2 and ~2 cv^2 m^3/tau^2. Their ratio is ~m. So where the -psi column is
dominated by its own noise, OLS returns K_m -> -(a noise-weighted mean of m). The script
computes that limit exactly and compares it with the Monte Carlo.

ESTIMATORS, each with the left-point and the midpoint form:
  OLS      t6's est_fd (checked equal).
  IV-JI    just-identified IV. The instruments come from samples k-1 and k+2, which are
           independent of row k's noise: W = [1, m~, -m~^2, -psi~], with m~ = (y_{k-1} + y_{k+2})/2
           and psi~ = (y_{k+2} - y_{k-1})/(3 tau) - mean u. The regressors and response are the
           CORRECTED rows, because psi m and m^2 have noise means of their own. Then
           E[w_k (z^_k - x^_k'p)] = E[w_k] (z0_k - x0_k'p), which vanishes at the noiseless
           solution (the truth, up to discretisation).
  2SLS     over-identified: W = [1, y_{k-1}, y_{k+2}, y_{k-1}^2, y_{k+2}^2, y_{k-1} y_{k+2}], plus u_k
           where the known load varies over the rows.
  GMM      two-step, from 2SLS. Bartlett HAC weight over lags 1-3: g_k = w_k e_k shares samples
           with g_{k+l} only for l <= 3, so this is the MA(1)-robust weight for these moments.
  wIV      IV-JI with the instruments scaled by 1/[m~^2 (m~ + K_w)^2], K_w = median(y). The
           equation error's sd is ~ m (m + K_m)/tau, so unweighted rows are badly heteroscedastic.
  BCLS     bias-compensated LS: the corrected X'X and X'z, solved.
  BCLS-2   the same with the second-order factor 1 - (i^2 + j^2) cv^2/2, as a check.
  wBCLS    BCLS with rows weighted as in wIV. The weights come from samples k-1 and k+2, so
           the correction stays unbiased, apart from K_w's weak dependence on all samples.
  GTLS-d   generalised TLS, the approximate-ML (Sampson) cost  sum e_k(p)^2 / var_k(p).
           Here e = z - X p (raw), var_k = (de_k/dy_k)^2 cv^2 y_k^2 + (de_k/dy_{k+1})^2 cv^2 y_{k+1}^2.
  GTLS-t   the same with the full tridiagonal covariance Omega = G V G' (G = de/dy, bidiagonal),
           i.e. e' Omega(p)^-1 e. This is the MA(1). Neither GTLS argmin depends on the cv
           level, only on the noise STRUCTURE. Both are constrained to K_m > 0 and k_e > 0
           (log-parametrised), and the script prints what happens without that constraint.
           Starts use no truth: K_m on a 25-point log grid from 1 to 1e4, with (p1, p2, p3) from
           weighted LS of z + K psi on [1, m, -m^2]. The best 3 by the diagonal cost, plus OLS,
           BCLS, wBCLS and wIV when physical. The lowest cost wins.

BENCHES: A and B are t6's (T.bench). C is an intermediate bench of this script's choosing:
theta = 150 (5%) on B's protocol (load x6 for 24 h then withdrawn, 15-min sampling, 1% CV).
NLS references are QUOTED from out/t6_estimators.txt for A and B. NLS on C is computed here
with t6's est_nls and its truth-centred starts, on fewer replicates (see output).
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import re
import time
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky_banded, solve_banded
from scipy.optimize import least_squares

import core as C
import dyn as D
import t6_estimators as T
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

HERE = Path(__file__).resolve().parent
NREP = int(os.environ.get("T6C_NREP", 300))           # replicates per bench and CV, every estimator
NREP_G_SWEEP = int(os.environ.get("T6C_NREP_G", 100))  # GTLS replicates per CV in the sweep (runtime)
NREP_NLS_C = int(os.environ.get("T6C_NREP_NLS", 100))  # NLS replicates on bench C (~1 s per fit)
SWEEP = (1e-2, 1e-3, 1e-4, 1e-5)
WORKERS = 3
KGRID = np.exp(np.linspace(np.log(1.0), np.log(1e4), 25))


# ================================================================ benches
def bench(name):
    if name in ("A", "B"):
        return T.bench(name)
    th = C.TH_5
    L = float(C.N4.anchor(th, C.KM0, C.KE, "A")[0])
    ex = D.Experiment([(-1, 0.5 * L), (0, 5.5 * L), (1, 0.0)], L0=0.5 * L)
    t = np.arange(0, 2 + 1e-9, 1 / 96)
    p = D.params(theta=th, ex=ex)
    m, _ = D.simulate(p, ex, t, sens=False)
    return dict(th=th, L=L, ex=ex, t=t, cv=0.01, p=p, m=m)


_BENCH = {}


def getb(name):
    if name not in _BENCH:
        _BENCH[name] = bench(name)
    return _BENCH[name]


# ================================================================ polynomials in (y_k, y_{k+1})
def padd(P, Q, s=1.0):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0.0) + s * v
    return R


def pmul(P, Q):
    R = {}
    for (i, j), v in P.items():
        for (k, l), w in Q.items():
            R[(i + k, j + l)] = R.get((i + k, j + l), 0.0) + v * w
    return R


def pscale(P, c):
    return {k: c * v for k, v in P.items()}


def pdiff(P, var):
    R = {}
    for (i, j), v in P.items():
        if var == 0 and i > 0:
            R[(i - 1, j)] = R.get((i - 1, j), 0.0) + i * v
        if var == 1 and j > 0:
            R[(i, j - 1)] = R.get((i, j - 1), 0.0) + j * v
    return R


class FD:
    """The FD regression for one bench and one discretisation, as polynomials per row."""

    def __init__(self, b, mid):
        self.b, self.mid = b, mid
        t, ex = b["t"], b["ex"]
        tau = np.diff(t)
        self.u = T.u_bar(ex, t)
        n1 = tau.size
        one = {(0, 0): np.ones(n1)}
        A, B = {(1, 0): np.ones(n1)}, {(0, 1): np.ones(n1)}
        psi = {(0, 1): 1 / tau, (1, 0): -1 / tau, (0, 0): -self.u}
        mm = pscale(padd(A, B), 0.5) if mid else A
        self.X = [one, mm, pscale(pmul(mm, mm), -1.0), pscale(psi, -1.0)]
        self.z = pmul(psi, mm)
        self.MM = [[pmul(self.X[i], self.X[j]) for j in range(4)] for i in range(4)]
        self.Mz = [pmul(self.X[i], self.z) for i in range(4)]
        self.dX = [[pdiff(c, v) for v in (0, 1)] for c in self.X]
        self.dz = [pdiff(self.z, v) for v in (0, 1)]
        self.n = t.size
        self.k = np.arange(1, self.n - 2)          # rows with both neighbours k-1, k+2

    @staticmethod
    def ev(P, pa, pb, fac):
        out = 0.0
        for (i, j), v in P.items():
            out = out + v * pa[i] * pb[j] * fac[i * i + j * j]
        return out * np.ones(pa[0].size)

    def prep(self, y, cv, s, order="exact"):
        a, b = y[:-1], y[1:]
        pa = [a ** i for i in range(5)]
        pb = [b ** i for i in range(5)]
        q = np.arange(33)
        fac = (np.exp(s * q * cv ** 2 / 2) if order == "exact"
               else 1 + s * q * cv ** 2 / 2)
        return pa, pb, fac

    def rows(self, y, cv, s):
        pa, pb, fac = self.prep(y, cv, s)
        Xr = np.stack([self.ev(c, pa, pb, fac) for c in self.X], axis=1)
        return Xr, self.ev(self.z, pa, pb, fac)

    def moments(self, y, cv, s, w=None, rows=None, order="exact"):
        pa, pb, fac = self.prep(y, cv, s, order)
        sel = slice(None) if rows is None else rows
        ww = 1.0 if w is None else w
        M = np.empty((4, 4))
        for i in range(4):
            for j in range(i, 4):
                M[i, j] = M[j, i] = np.sum(ww * self.ev(self.MM[i][j], pa, pb, fac)[sel])
        r = np.array([np.sum(ww * self.ev(self.Mz[i], pa, pb, fac)[sel]) for i in range(4)])
        return M, r


# ================================================================ estimators
def unlump(p):
    p = np.asarray(p, float)
    with np.errstate(all="ignore"):
        Km = p[..., 3]
        L0 = p[..., 0] / Km
        th = L0 - p[..., 1] - p[..., 2] * Km
    return th, p[..., 2], Km, L0


def ncond(X):
    return T.ncond(X)


def neighbours(fd, y):
    k = fd.k
    return k, y[k - 1], y[k + 2]


def weights(fd, y, Kw):
    k, lo, hi = neighbours(fd, y)
    mt = 0.5 * (lo + hi)
    return k, 1.0 / (mt ** 2 * (mt + Kw) ** 2)


def inst_ji(fd, y):
    k, lo, hi = neighbours(fd, y)
    t = fd.b["t"]
    u3 = (fd.u[k - 1] + fd.u[k] + fd.u[k + 1]) / 3
    psit = (hi - lo) / (t[k + 2] - t[k - 1]) - u3
    mt = 0.5 * (lo + hi)
    return np.stack([np.ones(k.size), mt, -mt ** 2, -psit], 1)


def inst_oi(fd, y):
    k, lo, hi = neighbours(fd, y)
    cols = [np.ones(k.size), lo, hi, lo ** 2, hi ** 2, lo * hi]
    if np.ptp(fd.u[k]) > 0:
        cols.append(fd.u[k])
    return np.stack(cols, 1)


def bartlett(g, L=3):
    n = g.shape[0]
    S = g.T @ g / n
    for l in range(1, L + 1):
        G = g[l:].T @ g[:-l] / n
        S += (1 - l / (L + 1)) * (G + G.T)
    return S


def proj(W, X):
    Q, _ = np.linalg.qr(W)
    return Q @ (Q.T @ X)


def partial_r2(W, x):
    r_full = x - proj(W, x[:, None])[:, 0]
    r_one = x - x.mean()
    return 1 - np.sum(r_full ** 2) / np.sum(r_one ** 2)


def linear_all(fd, y, cv):
    """All the closed-form estimators on one record. Returns {name: (p, diag)}."""
    with np.errstate(all="ignore"):
        return _linear_all(fd, y, cv)


def _linear_all(fd, y, cv):
    out = {}
    Xr, zr = fd.rows(y, cv, 0)
    Xc, zc = fd.rows(y, cv, -1)
    k = fd.k
    out["OLS"] = (np.linalg.lstsq(Xr, zr, rcond=None)[0], ncond(Xr))
    # IV, just identified
    W = inst_ji(fd, y)
    try:
        p = np.linalg.solve(W.T @ Xc[k], W.T @ zc[k])
    except np.linalg.LinAlgError:
        p = np.full(4, np.nan)
    out["IV-JI"] = (p, ncond(proj(W, Xc[k])))
    r2_ji = partial_r2(W, Xc[k, 3])
    # 2SLS, over identified
    Wo = inst_oi(fd, y)
    PX = proj(Wo, Xc[k])
    p2 = np.linalg.lstsq(PX, zc[k], rcond=None)[0]
    out["2SLS"] = (p2, ncond(PX))
    r2_oi = partial_r2(Wo, Xc[k, 3])
    # conditional (Sanderson-Windmeijer) first stage: -psi^ net of its IV projection on the
    # other regressors [1, m^, -m^2^], regressed on the instruments
    X3 = Xc[k, :3]
    P3 = proj(Wo, X3)
    dlt = np.linalg.lstsq(P3, Xc[k, 3], rcond=None)[0]
    r2_sw = partial_r2(Wo, Xc[k, 3] - X3 @ dlt)
    # GMM, two step, Bartlett HAC over lags 1-3
    e = zc[k] - Xc[k] @ p2
    try:
        Si = np.linalg.inv(bartlett(Wo * e[:, None]))
        A = Xc[k].T @ Wo @ Si @ Wo.T @ Xc[k]
        pg = np.linalg.solve(A, Xc[k].T @ Wo @ Si @ Wo.T @ zc[k])
    except np.linalg.LinAlgError:
        pg = p2
    out["GMM"] = (pg, ncond(PX))
    # weighted IV
    Kw = float(np.median(y))
    _, w = weights(fd, y, Kw)
    Ww = W * w[:, None]
    try:
        pw = np.linalg.solve(Ww.T @ Xc[k], Ww.T @ zc[k])
    except np.linalg.LinAlgError:
        pw = np.full(4, np.nan)
    out["wIV"] = (pw, ncond(proj(Ww, Xc[k] * np.sqrt(w)[:, None])))
    # BCLS, exact and second order, and weighted
    for nm, kw in (("BCLS", dict()), ("BCLS-2", dict(order=2)),
                   ("wBCLS", dict(w=w, rows=k))):
        M, r = fd.moments(y, cv, -1, **kw)
        try:
            p = np.linalg.solve(M, r)
        except np.linalg.LinAlgError:
            p = np.full(4, np.nan)
        dd = np.sqrt(np.abs(np.diag(M)))
        ev = np.linalg.eigvalsh(M / np.outer(dd, dd))
        out[nm] = (p, float(np.sqrt(np.abs(ev).max() / np.abs(ev).min())), bool(ev.min() > 0))
    out["_r2"] = (r2_ji, r2_oi, r2_sw)
    return out


# ---------------------------------------------------------------- GTLS
class Gtls:
    def __init__(self, fd, y, cv):
        a, b = y[:-1], y[1:]
        pa, pb, fac = fd.prep(y, cv, 0)
        self.X = np.stack([fd.ev(c, pa, pb, fac) for c in fd.X], axis=1)
        self.z = fd.ev(fd.z, pa, pb, fac)
        self.DXa = np.stack([fd.ev(fd.dX[i][0], pa, pb, fac) for i in range(4)], 1)
        self.DXb = np.stack([fd.ev(fd.dX[i][1], pa, pb, fac) for i in range(4)], 1)
        self.dza = fd.ev(fd.dz[0], pa, pb, fac)
        self.dzb = fd.ev(fd.dz[1], pa, pb, fac)
        self.vy = cv ** 2 * y ** 2
        self.n1 = a.size

    def parts(self, p):
        e = self.z - self.X @ p
        ga = self.dza - self.DXa @ p
        gb = self.dzb - self.DXb @ p
        d = ga ** 2 * self.vy[:-1] + gb ** 2 * self.vy[1:]
        return e, ga, gb, d

    def white(self, p, tri, M):
        """Omega(p)^-1/2 applied to M (a vector or the columns of a matrix)."""
        e, ga, gb, d = self.parts(p)
        if not tri:
            return (M.T / np.sqrt(d)).T
        off = gb[:-1] * ga[1:] * self.vy[1:-1]
        ab = np.zeros((2, self.n1))
        ab[0, 1:] = off
        ab[1] = d
        U = cholesky_banded(ab, lower=False)
        Ut = np.zeros((2, self.n1))
        Ut[0] = U[1]
        Ut[1, :-1] = U[0, 1:]
        return solve_banded((1, 0), Ut, M)

    def res(self, p, tri):
        with np.errstate(all="ignore"):
            try:
                r = self.white(p, tri, self.z - self.X @ p)
            except (np.linalg.LinAlgError, ValueError):
                return np.full(self.n1, 1e8)
        if not np.all(np.isfinite(r)):
            return np.full(self.n1, 1e8)
        return r

    def cost(self, p, tri):
        return 0.5 * float(np.sum(self.res(p, tri) ** 2))

    def kstarts(self):
        m = self.X[:, 1]
        out = []
        for K in KGRID:
            lhs = self.z - K * self.X[:, 3]
            sw = 1.0 / (np.abs(m) * (np.abs(m) + K))
            q = np.linalg.lstsq(self.X[:, :3] * sw[:, None], lhs * sw, rcond=None)[0]
            out.append(np.r_[q[0], q[1], max(q[2], 1e-2), K])
        c = [self.cost(p, False) for p in out]
        return [out[i] for i in np.argsort(c)[:3]]

    def fit(self, starts, tri, constrained=True, nfev=200):
        best = None
        for p0 in starts:
            p0 = np.asarray(p0, float)
            if not np.all(np.isfinite(p0)):
                continue
            if constrained and not (p0[2] > 0 and p0[3] > 0):
                continue
            sc = np.maximum(np.abs(p0), 1e-6)
            if constrained:
                sc[2:] = 1.0
                q0 = np.r_[p0[:2] / sc[:2], np.log(p0[2]), np.log(p0[3])]

                def back(q):
                    return np.r_[q[:2] * sc[:2], np.exp(np.clip(q[2:], -40, 40))]
            else:
                q0 = p0 / sc

                def back(q):
                    return q * sc
            with np.errstate(all="ignore"):
                try:
                    s = least_squares(lambda q: self.res(back(q), tri), q0, method="lm",
                                      xtol=1e-12, ftol=1e-12, max_nfev=nfev)
                except (np.linalg.LinAlgError, ValueError):
                    continue
            if not np.isfinite(s.cost):
                continue
            if best is None or s.cost < best[1]:
                best = (back(s.x), s.cost)
        return best

    def wcond(self, p, tri):
        with np.errstate(all="ignore"):
            try:
                return ncond(self.white(p, tri, self.X))
            except (np.linalg.LinAlgError, ValueError):
                return np.nan


# ================================================================ Monte Carlo jobs
LIN = ("OLS", "IV-JI", "2SLS", "GMM", "wIV", "BCLS", "BCLS-2", "wBCLS")
GT = ("GTLS-d", "GTLS-t")
NAMES = LIN + GT


def job(args):
    bname, mid, cv, reps, reps_g, gt_d, unc_n, seed = args
    b = getb(bname)
    fd = FD(b, mid)
    rng = np.random.default_rng(seed)
    P = {k: np.full((reps, 4), np.nan) for k in NAMES}
    cond = {k: np.full(reps, np.nan) for k in NAMES}
    pd_frac = {"BCLS": [], "wBCLS": []}
    r2 = []
    Mc, Mr = [], []
    unc = []
    t0 = time.time()
    for r in range(reps):
        y = b["m"] * np.exp(rng.normal(0, cv, b["m"].size))
        L = linear_all(fd, y, cv)
        for k in LIN:
            P[k][r] = L[k][0]
            cond[k][r] = L[k][1]
        for k in pd_frac:
            pd_frac[k].append(L[k][2])
        r2.append(L["_r2"])
        if r < 300 and cv == b["cv"]:
            Mc.append(fd.moments(y, cv, -1)[0])
            Mr.append(fd.moments(y, cv, 0)[0])
        if r < reps_g:
            g = Gtls(fd, y, cv)
            starts = g.kstarts() + [L[k][0] for k in ("OLS", "BCLS", "wBCLS", "wIV")]
            for nm, tri in (("GTLS-d", False), ("GTLS-t", True)):
                if nm == "GTLS-d" and not gt_d:
                    continue
                f = g.fit(starts, tri)
                if f is not None:
                    P[nm][r] = f[0]
                    cond[nm][r] = g.wcond(f[0], tri)
            if r < unc_n:
                f = g.fit([L[k][0] for k in ("OLS", "BCLS", "wBCLS", "wIV")], False,
                          constrained=False)
                unc.append(f[0][3] if f is not None else np.nan)
    return dict(key=(bname, mid, cv), P=P, cond=cond, pd={k: np.mean(v) for k, v in pd_frac.items()},
                r2=np.array(r2), Mc=np.array(Mc), Mr=np.array(Mr), unc=np.array(unc),
                secs=time.time() - t0, reps=reps, reps_g=reps_g)


def nls_job(args):
    bname, seed, reps = args
    b = getb(bname)
    th0, ke0, K0, L00, _ = D.unpack(b["p"])
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(reps):
        y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
        ps = [b["p"], b["p"] + np.array([0.7, -0.3, 0.5, -0.2 * L00]),
              b["p"] + np.array([-0.7, 0.3, -0.5, 0.2 * L00])]
        out.append(T.est_nls(y, b["t"], b["ex"], ps))
    return np.array(out)


# ================================================================ reporting helpers
def summ(P, b):
    th0, ke0, K0, L00, _ = D.unpack(b["p"])
    ok = np.all(np.isfinite(P), axis=1)
    ok[ok] = np.isfinite(unlump(P[ok])[0])
    th, ke, Km, L0 = unlump(P[ok])
    a = T.summarise(th, th0)
    c = T.summarise(ke, ke0)
    return a, c, int(ok.sum())


def fmt_row(lab, cn, a, c, n):
    def f(v, s):
        return f"{v:{s}}" if np.isfinite(v) else f"{'n/a':>{s.lstrip('+').split('.')[0]}}"
    return (f"     {lab:16} {cn:9.3g} {f(a['bias'], '+8.3f')} {f(a['rmse'], '9.3f')} "
            f"{f(a['robust'], '8.3f')} {a['neg']:6.1%} | {f(c['bias'], '+8.3f')} "
            f"{f(c['rmse'], '8.3f')} {f(c['robust'], '8.3f')} {c['neg']:6.1%} {n:5d}")


def fim_sd(b, cv):
    F = D.fim(b["p"], b["ex"], b["t"], cv=cv)
    Ci = np.linalg.inv(F)
    return float(np.sqrt(Ci[0, 0])), float(np.sqrt(Ci[1, 1]))


def noiseless_limits(fd, b):
    m = b["m"]
    L = linear_all(fd, m, 0.0)
    out = {k: L[k][0] for k in LIN}
    g = Gtls(fd, m, 1.0)            # the GTLS argmin does not depend on the cv level
    starts = g.kstarts() + [L[k][0] for k in ("OLS", "BCLS", "wBCLS", "wIV")]
    for nm, tri in (("GTLS-d", False), ("GTLS-t", True)):
        f = g.fit(starts, tri)
        out[nm] = f[0] if f is not None else np.full(4, np.nan)
    return out


def quote_t6():
    txt = (HERE / "out" / "t6_estimators.txt").read_text().splitlines()
    out, cur = {}, None
    for ln in txt:
        mb = re.match(r"\s+BENCH (\w):", ln)
        if mb:
            cur = mb.group(1)
        if cur and re.match(r"\s+\((a|c)\) (FD, left point|FD, midpoint|NLS on the ODE)", ln):
            out.setdefault(cur, []).append(ln.rstrip())
    return out


# ================================================================ main
def main():
    T0 = time.time()
    print("=" * 112)
    print("  ITEM 4. ONE FIX FOR THE FINITE-DIFFERENCE ESTIMATOR: INSTRUMENTS, BIAS COMPENSATION, GENERALISED TLS")
    print("=" * 112)
    names = ("A", "B", "C")
    B = {nm: getb(nm) for nm in names}
    desc = {"A": "t6's TB-A: theta 150, diet withdrawn, hourly 48 h",
            "B": "t6's TB-B: theta 600, load x6 for 24 h then out, 15 min",
            "C": "intermediate: theta 150 (5%), TB-B's x6 protocol, 15 min"}
    print("\n  0. BENCHES (no priors anywhere in this script)")
    for nm in names:
        b = B[nm]
        th0, ke0, K0, L00, _ = D.unpack(b["p"])
        s1, k1 = fim_sd(b, b["cv"])
        su, _ = fim_sd(b, 1.0)
        print(f"     {nm}: {desc[nm]}, CV {100*b['cv']:.0f}%; {b['t'].size} samples, m {b['m'].min():.1f}"
              f"-{b['m'].max():.1f}; FIM sd ln theta {s1:.4g}, ln k_e {k1:.4g} (sd ln theta per unit"
              f" CV {su:.5g}, so {su*1e-3:.4g} at 0.1%)")

    # ------------------------------------------------ Monte Carlo jobs, in parallel
    jobs = []
    for bi, nm in enumerate(names):
        b = B[nm]
        cvs = [b["cv"]] + [c for c in SWEEP if not np.isclose(c, b["cv"])]
        for mid in (False, True):
            for ci, cv in enumerate(cvs):
                nominal = np.isclose(cv, b["cv"])
                jobs.append((nm, mid, cv, NREP, NREP if nominal else NREP_G_SWEEP,
                             nominal, 50 if nominal else 0, 7000 + 100 * bi + 10 * mid + ci))
    nls_jobs = [("C", 900 + i, NREP_NLS_C // 4) for i in range(4)]
    ctx = mp.get_context("fork")
    with ctx.Pool(WORKERS) as pool:
        nres = pool.map_async(nls_job, nls_jobs)
        res = pool.map(job, sorted(jobs, key=lambda j: -j[3] * (3 if j[5] else 1)), chunksize=1)
        nls_c = np.concatenate(nres.get())
    R = {r["key"]: r for r in res}

    def get(nm, mid, cv):
        for k, v in R.items():
            if k[0] == nm and k[1] == mid and np.isclose(k[2], cv):
                return v
        raise KeyError((nm, mid, cv))

    # ------------------------------------------------ 1. checks
    print("\n  1. CHECKS")
    for nm in names:
        b = B[nm]
        for mid in (False, True):
            fd = FD(b, mid)
            rng = np.random.default_rng(5)
            y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
            p_ols = linear_all(fd, y, b["cv"])["OLS"][0]
            p_t6 = T.est_fd(y, b["t"], b["ex"], mid=mid)["p"]
            err = float(np.max(np.abs(p_ols / p_t6 - 1)))
            assert err < 1e-8, err
            print(f"     bench {nm}, {'midpoint' if mid else 'left point'}: OLS through the polynomial "
                  f"route = t6's est_fd to {err:.1e}")
    # exactness of the correction: mean corrected moments against the noiseless ones
    for nm in names:
        b = B[nm]
        fd = FD(b, True)
        r = get(nm, True, b["cv"])
        M0, _ = fd.moments(b["m"], 0.0, 0)
        Mc, Mr = r["Mc"], r["Mr"]
        se = Mc.std(0, ddof=1) / np.sqrt(Mc.shape[0])
        iu = np.triu_indices(4)
        zrf = (Mr.mean(0) - M0) / np.where(se > 0, se, 1)
        zc = ((Mc.mean(0) - M0) / np.where(se > 0, se, 1))[iu]
        zr = zrf[iu]
        assert np.max(np.abs(zc)) < 5, zc
        print(f"     bench {nm}, midpoint, {Mc.shape[0]} records: mean CORRECTED X'X against the noiseless X'X,"
              f" max |z| over the 10 entries {np.max(np.abs(zc)):.2f}; RAW X'X max |z| "
              f"{np.max(np.abs(zr)):.1f} (psi-psi entry z = {zrf[3, 3]:.1f})")
    # IV exogeneity: E[W'(z^ - X^ p)] = E[W]'(z0 - X0 p), corrected against raw rows
    for nm in ("A", "B"):
        b = B[nm]
        fd = FD(b, True)
        lt = T.lumped_true(b)
        X0, z0 = fd.rows(b["m"], 0.0, 0)
        d0 = (z0 - X0 @ lt)[fd.k]
        rng = np.random.default_rng(17)
        gc, gr = [], []
        for _ in range(2000):
            y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
            W = inst_ji(fd, y)
            Xc, zc_ = fd.rows(y, b["cv"], -1)
            Xr, zr_ = fd.rows(y, b["cv"], 0)
            gc.append(W.T @ (zc_ - Xc @ lt)[fd.k] - W.T @ d0)
            gr.append(W.T @ (zr_ - Xr @ lt)[fd.k] - W.T @ d0)
        gc, gr = np.array(gc), np.array(gr)
        zc = gc.mean(0) / (gc.std(0, ddof=1) / np.sqrt(len(gc)))
        zr = gr.mean(0) / (gr.std(0, ddof=1) / np.sqrt(len(gr)))
        assert np.max(np.abs(zc)) < 5, zc
        print(f"     bench {nm}, midpoint, 2000 records: IV-JI moments at the true lumped p, minus "
              f"E[W]'(discretisation residual): corrected rows max |z| {np.max(np.abs(zc)):.2f}, "
              f"raw rows max |z| {np.max(np.abs(zr)):.1f}")
    for nm in names:
        b = B[nm]
        r = get(nm, True, b["cv"])
        d = np.abs(unlump(r["P"]["BCLS"])[0] - unlump(r["P"]["BCLS-2"])[0]) / np.abs(unlump(r["P"]["BCLS"])[0])
        print(f"     bench {nm}, midpoint: exact against second-order correction, median |relative "
              f"difference| in theta_hat {np.median(d):.1e}")

    # ------------------------------------------------ 2. mechanism
    print("\n  2. MECHANISM, QUANTIFIED (exact expectations from the lognormal moments; no Monte Carlo"
          " except where marked MC)")
    for nm in names:
        b = B[nm]
        lt = T.lumped_true(b)
        for mid in (False, True):
            fd = FD(b, mid)
            r = get(nm, mid, b["cv"])
            X0, z0 = fd.rows(b["m"], 0.0, 0)
            M0, r0 = fd.moments(b["m"], 0.0, 0)
            ME, rE = fd.moments(b["m"], b["cv"], +1)            # E[X'X], E[X'z] exactly
            dd = np.sqrt(np.diag(M0))
            ev = np.linalg.eigvalsh(M0 / np.outer(dd, dd))
            share = (np.diag(ME) - np.diag(M0)) / np.diag(ME)
            p_lim = np.linalg.solve(ME, rE)
            Kn = (rE[3] - r0[3]) / (ME[3, 3] - M0[3, 3])         # noise-only coefficient on -psi
            Kmc = np.median(r["P"]["OLS"][:, 3])
            # partial signal of the psi column: psi0 off [1, m, m^2], against psi's noise variance
            psi0 = -X0[:, 3]
            Q = proj(X0[:, :3], psi0[:, None])[:, 0]
            rs = psi0 - Q
            pa, pb, fac = fd.prep(b["m"], b["cv"], +1)
            psiP = pscale(fd.X[3], -1.0)
            Epsi2 = fd.ev(pmul(psiP, psiP), pa, pb, fac)
            Epsi = fd.ev(psiP, pa, pb, fac)
            vpsi = Epsi2 - Epsi ** 2
            rho = np.sum(rs ** 2) / np.sum(vpsi)
            w = 1.0 / vpsi
            sw = np.sqrt(w)
            rsw = psi0 * sw - proj(X0[:, :3] * sw[:, None], (psi0 * sw)[:, None])[:, 0]
            rho_w = np.sum(rsw ** 2) / np.sum(w * vpsi)
            # the same with psi's full (tridiagonal) noise covariance: consecutive differences
            # share a sample, cov(psi_k, psi_k+1) = -Var(y_k+1)/(tau_k tau_k+1)
            tau = np.diff(b["t"])
            vy_ = b["m"] ** 2 * np.exp(b["cv"] ** 2) * np.expm1(b["cv"] ** 2)
            Om = np.diag(vpsi)
            cv_ = -vy_[1:-1] / (tau[:-1] * tau[1:])
            Om += np.diag(cv_, 1) + np.diag(cv_, -1)
            Lc = np.linalg.cholesky(Om)
            Wx = np.linalg.solve(Lc, X0[:, :3])
            Wp = np.linalg.solve(Lc, psi0)
            rho_t = np.sum((Wp - proj(Wx, Wp[:, None])[:, 0]) ** 2) / psi0.size
            Mc = r["Mc"]
            if Mc.size == 0:
                Mc = np.array([fd.moments(b["m"] * np.exp(np.random.default_rng(i).normal(
                    0, b["cv"], b["m"].size)), b["cv"], -1)[0] for i in range(100)])
            sdpp = (Mc[:, 3, 3].std(ddof=1)) / M0[3, 3] if mid else np.nan
            lab = "midpoint" if mid else "left point"
            print(f"   bench {nm}, {lab}:")
            print(f"     noiseless design: cond(X0) {ncond(X0):.4g}; normalised eigenvalues of X0'X0 "
                  f"{ev[0]:.2e} ... {ev[-1]:.3f}. With noise (MC median) cond(X) {np.median(r['cond']['OLS']):.3g}")
            print(f"     noise share of E[X'X] diagonal (1, m, m^2, psi): {share[0]:.1e}, {share[1]:.1e}, "
                  f"{share[2]:.1e}, {share[3]:.4f}; of the psi entry of E[X'z]: "
                  f"{(rE[3]-r0[3])/rE[3]:.4f}")
            print(f"     partial signal-to-noise of the psi column (psi0 off [1, m, m^2] / noise variance): "
                  f"unweighted {rho:.3g}, noise-precision weighted {rho_w:.3g}, whitened with the MA(1) "
                  f"covariance {rho_t:.3g} (per row); each scales as CV^-2")
            th_lim = unlump(p_lim)[0]
            ratio = p_lim / lt
            mcr = np.median(r["P"]["OLS"], axis=0) / lt
            print(f"     OLS limit E[X'X]^-1 E[X'z]: coefficients / truth {ratio[0]:+.3f} {ratio[1]:+.3f} "
                  f"{ratio[2]:+.3f} {ratio[3]:+.3f}, theta {th_lim:.4g}; MC median {mcr[0]:+.3f} {mcr[1]:+.3f} "
                  f"{mcr[2]:+.3f} {mcr[3]:+.3f}")
            print(f"     K_m coefficient if the -psi column fitted only its own noise: {Kn:.1f}; OLS limit "
                  f"{p_lim[3]:.1f}; MC median {Kmc:.1f} (truth {lt[3]:.0f})")
            if mid:
                print(f"     corrected psi-psi entry: MC sd / noiseless value {sdpp:.2e}, against the smallest "
                      f"normalised eigenvalue {ev[0]:.2e} (ratio {sdpp/ev[0]:.0f})")

    # ------------------------------------------------ 3. Monte Carlo at the nominal CV
    print(f"\n  3. MONTE CARLO AT EACH BENCH'S OWN CV, {NREP} records each (same records for every "
          "estimator)")
    print("     cond: OLS cond(X); IV cond(P_W X^); BCLS sqrt cond of the normalised corrected X'X;"
          " GTLS cond of Omega^-1/2 X at the fit (medians).")
    print("     bias = median error, RMSE over positive estimates, robust = half the 16-84% range"
          " (t6's summarise)")
    t6q = quote_t6()
    NOM = {}
    for nm in names:
        b = B[nm]
        th0, ke0, K0, L00, _ = D.unpack(b["p"])
        s1, k1 = fim_sd(b, b["cv"])
        print(f"\n   BENCH {nm} (CV {100*b['cv']:.0f}%), FIM sd ln theta {s1:.4g}, ln k_e {k1:.4g}")
        print(f"     {'estimator':16} {'cond':>9} {'ln th':>8} {'RMSE':>9} {'robust':>8} {'th<=0':>6} |"
              f" {'ln ke':>8} {'RMSE':>8} {'robust':>8} {'ke<=0':>6} {'n':>5}")
        for mid in (False, True):
            r = get(nm, mid, b["cv"])
            print(f"     -- {'midpoint' if mid else 'left point'}")
            for k in NAMES:
                a, c, n = summ(r["P"][k], b)
                NOM[(nm, mid, k)] = (a, c)
                print(fmt_row(k, float(np.nanmedian(r["cond"][k])), a, c, n))
            print(f"     BCLS corrected X'X positive definite in {r['pd']['BCLS']:.1%} of records, "
                  f"wBCLS {r['pd']['wBCLS']:.1%}; first-stage partial R^2 of -psi^ on the instruments, "
                  f"median: JI {np.median(r['r2'][:, 0]):.3f}, OI {np.median(r['r2'][:, 1]):.3f}; "
                  f"conditional on the other regressors (Sanderson-Windmeijer), OI "
                  f"{np.median(r['r2'][:, 2]):.4f}")
            u = r["unc"]
            print(f"     GTLS-d WITHOUT the K_m, k_e > 0 constraint, from the four linear starts, "
                  f"{u.size} records: K_m <= 0 in {np.mean(u <= 0):.0%}, K_m > 1e4 in "
                  f"{np.mean(u > 1e4):.0%}")
            print(f"     ({r['secs']:.0f} s for this job)")
        if nm in t6q:
            print("     quoted from out/t6_estimators.txt (t6's own records: 400 linear, 80 NLS):")
            for ln in t6q[nm]:
                print("       " + ln.strip())
        if nm == "C":
            a = T.summarise(np.exp(nls_c[:, 0]), th0)
            c = T.summarise(np.exp(nls_c[:, 1]), ke0)
            NOM[("C", None, "NLS")] = (a, c)
            print(f"     NLS on the ODE (t6's est_nls, truth-centred starts), {nls_c.shape[0]} records "
                  f"(not 300: ~1 s a fit):")
            print(fmt_row("NLS", np.nan, a, c, nls_c.shape[0]))

    # ------------------------------------------------ 4. consistency sweep
    print("\n  4. CONSISTENCY: CV -> 0. Cells: median ln(theta_hat/theta) / robust spread / fraction"
          f" theta_hat <= 0. Linear estimators {NREP} records per CV; GTLS {NREP} at the nominal CV,"
          f" {NREP_G_SWEEP} elsewhere.")
    print("     'limit' = the estimator applied to NOISELESS samples (its CV -> 0 limit, discretisation only).")
    print("     converged at a CV: theta_hat <= 0 below 1%, robust spread <= 0.10 and |median - limit|"
          " <= 0.05, there and at every smaller CV in the sweep")
    CONV = {}
    for nm in names:
        b = B[nm]
        th0 = b["th"]
        cvs = sorted({b["cv"], *SWEEP}, reverse=True)
        fims = [fim_sd(b, c)[0] for c in cvs]
        for mid in (False, True):
            fd = FD(b, mid)
            lim = noiseless_limits(fd, b)
            print(f"\n   BENCH {nm}, {'midpoint' if mid else 'left point'}  (FIM sd ln theta: "
                  + ", ".join(f"{100*c:g}% {f:.3g}" for c, f in zip(cvs, fims)) + ")")
            print(f"     {'estimator':10} {'limit':>8} " + " ".join(f"{'CV ' + format(100*c, 'g') + '%':>22}"
                                                            for c in cvs) + "   converged from")
            for k in NAMES:
                if k == "GTLS-d":
                    continue
                with np.errstate(all="ignore"):
                    tl = unlump(lim[k])[0]
                el = np.log(tl / th0) if np.isfinite(tl) and tl > 0 else np.nan
                cells, okc = [], []
                for c in cvs:
                    r = get(nm, mid, c)
                    a, _, n = summ(r["P"][k], b)
                    cells.append(f"{a['bias']:+7.3f}/{a['robust']:6.3f}/{a['neg']:5.1%}"
                                 if np.isfinite(a['bias']) else f"{'all <= 0':>22}")
                    okc.append(bool(np.isfinite(a["bias"]) and a["neg"] < 0.01 and a["robust"] <= 0.10
                                    and np.isfinite(el) and abs(a["bias"] - el) <= 0.05))
                conv = None
                for c, o in zip(cvs[::-1], okc[::-1]):
                    if o:
                        conv = c
                    else:
                        break
                CONV[(nm, mid, k)] = conv
                cs = f"{100*conv:g}%" if conv else "not by 0.001%"
                ls = f"{el:+8.3f}" if np.isfinite(el) else f"{'th<=0':>8}"
                print(f"     {k:10} {ls} " + " ".join(f"{s:>22}" for s in cells) + f"   {cs}")
            if mid:
                for k in ("OLS", "BCLS", "IV-JI", "wIV", "GTLS-t"):
                    with np.errstate(all="ignore"):
                        tl = unlump(lim[k])[0]
                    assert tl > 0 and abs(np.log(tl / th0)) < 0.05, (nm, k, tl)

    # ------------------------------------------------ 5. verdict lines
    print("\n  5. SUMMARY AT THE NOMINAL CV, midpoint (median ln-theta error / robust spread / theta_hat <= 0)")
    for nm in names:
        b = B[nm]
        s1, _ = fim_sd(b, b["cv"])
        parts = []
        for k in ("OLS", "IV-JI", "GMM", "wIV", "BCLS", "wBCLS", "GTLS-d", "GTLS-t"):
            a, _ = NOM[(nm, True, k)]
            parts.append(f"{k} {a['bias']:+.2f}/{a['robust']:.2f}/{a['neg']:.0%}" if np.isfinite(a["bias"])
                         else f"{k} all<=0")
        print(f"     {nm} (FIM {s1:.3g}): " + "; ".join(parts))
        a, _ = NOM[(nm, True, "GTLS-t")]
        print(f"     {nm}: GTLS-t robust spread / FIM sd = {a['robust'] / s1:.3g}")
    for nm in names:
        print(f"     {nm} converged-from CV (midpoint): " + ", ".join(
            f"{k} {100*CONV[(nm, True, k)]:g}%" if CONV[(nm, True, k)] else f"{k} none"
            for k in ("OLS", "IV-JI", "2SLS", "GMM", "wIV", "BCLS", "wBCLS", "GTLS-t")))
        print(f"     {nm} converged-from CV (left point): " + ", ".join(
            f"{k} {100*CONV[(nm, False, k)]:g}%" if CONV[(nm, False, k)] else f"{k} none"
            for k in ("OLS", "IV-JI", "GMM", "wIV", "BCLS", "wBCLS", "GTLS-t")))

    # ------------------------------------------------ figure
    fig, axs = plt.subplots(1, 3, figsize=(13.0, 3.9), sharey=True)
    series = (("OLS", "OLS-FD (t6)", C1, "o"), ("GMM", "GMM-HAC (IV)", C2, "s"),
              ("wBCLS", "weighted BCLS", C3, "^"), ("GTLS-t", "GTLS, MA(1) covariance", C4, "D"))
    for ax, nm in zip(axs, names):
        b = B[nm]
        cvs = np.array(sorted({b["cv"], *SWEEP}))
        ax.plot(100 * cvs, [fim_sd(b, c)[0] for c in cvs], color=INK, lw=1.2, ls="--",
                label="FIM sd (efficient)")
        for k, lab, col, mk in series:
            rob, neg = [], []
            for c in cvs:
                a, _, _ = summ(get(nm, True, c)["P"][k], b)
                rob.append(a["robust"] if np.isfinite(a["robust"]) else np.nan)
                neg.append(a["neg"])
            rob, neg = np.array(rob), np.array(neg)
            rob[~(rob > 0)] = np.nan
            ax.plot(100 * cvs, rob, color=col, lw=2, label=lab)
            good = neg < 0.5
            ax.scatter(100 * cvs[good], rob[good], s=40, color=col, marker=mk, zorder=3,
                       edgecolors="white", linewidths=1.0)
            ax.scatter(100 * cvs[~good], rob[~good], s=40, facecolors="white", edgecolors=col,
                       marker=mk, zorder=3, linewidths=1.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(100 * b["cv"], color=MUTED, lw=0.8, ls=":")
        ax.text(100 * b["cv"], 12, "own CV ", color=MUTED, fontsize=8, va="bottom", ha="right")
        ax.set_xlabel("assay CV (%)")
        ax.set_title({"A": "TB-A: 5%, washout, hourly", "B": "TB-B: 20%, ×6 step, 15 min",
                      "C": "C: 5%, ×6 step, 15 min"}[nm])
        plain(ax)
        from matplotlib.ticker import NullFormatter
        ax.xaxis.set_minor_formatter(NullFormatter())
        finish(ax)
    axs[0].set_ylabel("robust spread of ln θ̂ (midpoint form)")
    axs[0].set_ylim(1e-5, 30)
    hs, ls = axs[0].get_legend_handles_labels()
    fig.legend(hs, ls, loc="lower center", ncol=5, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.07))
    fig.text(0.5, -0.12, "hollow marker: θ̂ ≤ 0 in at least half the records (spread over the positive"
             " ones only); a gap: every θ̂ ≤ 0", ha="center", fontsize=8, color=MUTED)
    save(fig, "t6c_iv")
    print(f"  runtime {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
