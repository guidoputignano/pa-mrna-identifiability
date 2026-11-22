"""
Adversarial recheck of t6c_iv.py (item 4: instruments, bias compensation and generalised TLS
for the finite-difference regression).

Shares NO code with t6c_iv.py, core.py, dyn.py or any t*.py: nothing from this directory is
imported. Every route is rebuilt from the model statement:

  dm/dt = u(t) + L0 - theta m/(K_m + m) - k_e m,  m(0-) at rest,  y_k = m(t_k) exp(e_k),
  e_k ~ N(0, cv^2) independent;  L = 3000*10/110 + 5*10 (convention A), L0 = L/2.

  bench A  theta 150, u = L/2 before t = 0 and 0 after, hourly for 48 h, cv 5%
  bench B  theta 600, u = L/2, then 5.5 L on [0, 1) d, then 0, every 15 min for 48 h, cv 1%
  bench C  theta 150 on B's protocol, cv 1%

Routes, each different from the investigator's:
  R1  the ODE solved with DOP853 (not LSODA with forward sensitivities); the information
      matrix from CENTRAL FINITE DIFFERENCES of ln m with Richardson extrapolation,
      two step sizes printed, instead of the sensitivity equations.
  R2  the noiseless FD design from my own simulation: cond(X0) and its smallest eigenvalue.
  R3  the OLS limit E[X'X]^-1 E[X'z] by DIRECT Monte Carlo averaging of X'X and X'z over many
      records (the investigator used exact lognormal moments).
  R4  the bias-corrected rows and moments built SYMBOLICALLY with sympy (sp.Poly, each
      monomial a^i b^j scaled by exp(-(i^2+j^2) cv^2/2)), not the investigator's dict
      polynomials; IV-JI, wIV, 2SLS and BCLS re-implemented from the definitions in its
      docstring, on NEW seeds. Also an IV with wider instruments (k-2, k+3), which the request
      allows ("k-1 and k+2 or further") and the investigator did not try.
  R5  GTLS-t re-implemented in PHYSICAL coordinates (theta, ln k_e, ln K_m, L0; theta and L0
      free, so the same constraint set as the investigator's K_m, k_e > 0), with a dense
      Cholesky whitening, scipy's trf optimiser and a different start grid in (k_e, K_m).
      A start AT THE TRUTH is run as a diagnostic only, to see whether the no-truth starts
      find the lowest cost.
  R6  the lag structure of the IV moments g_k = w_k e_k, by Monte Carlo at the true
      lumped coefficients.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

T0 = time.time()
TH_N, KM, KE, M_N = 3000.0, 100.0, 5.0, 10.0
LT = TH_N * M_N / (KM + M_N) + KE * M_N


def mss(th, ke, K, Ltot):
    """Rest state of th m/(K+m) + ke m = Ltot (positive root)."""
    q = th + ke * K - Ltot
    return (-q + np.sqrt(q * q + 4 * ke * Ltot * K)) / (2 * ke)


# =========================================================================== benches
class Bench:
    def __init__(self, name):
        self.name = name
        self.L0 = 0.5 * LT
        if name == "A":
            self.th, self.cv = 150.0, 0.05
            self.segs = [(0.0, 2.0, 0.0)]
            self.u_pre = 0.5 * LT
            self.t = np.arange(49) / 24.0
        else:
            self.th = 600.0 if name == "B" else 150.0
            self.cv = 0.01
            self.segs = [(0.0, 1.0, 5.5 * LT), (1.0, 2.0, 0.0)]
            self.u_pre = 0.5 * LT
            self.t = np.arange(193) / 96.0
        self.ptrue = np.array([np.log(self.th), np.log(KE), np.log(KM), self.L0])
        self.m = self.sim(self.ptrue)
        tm = 0.5 * (self.t[:-1] + self.t[1:])
        self.ubar = np.array([self.useg(x) for x in tm])
        self.tau = np.diff(self.t)

    def useg(self, x):
        for lo, hi, u in self.segs:
            if lo <= x < hi:
                return u
        return self.segs[-1][2]

    def sim(self, prm, rtol=1e-12):
        th, ke, K, L0 = np.exp(prm[0]), np.exp(prm[1]), np.exp(prm[2]), prm[3]
        m0 = mss(th, ke, K, self.u_pre + L0)
        out = np.empty(self.t.size)
        y0 = m0
        for lo, hi, u in self.segs:
            sel = (self.t >= lo - 1e-12) & (self.t <= hi + 1e-12)
            te = self.t[sel]
            sol = solve_ivp(lambda tt, mm: [u + L0 - th * mm[0] / (K + mm[0]) - ke * mm[0]],
                            (lo, hi), [y0], t_eval=np.clip(te, lo, hi), method="DOP853",
                            rtol=rtol, atol=1e-12)
            out[sel] = sol.y[0]
            y0 = sol.y[0, -1]
        return out

    def lumped(self):
        th, ke, K, L0 = self.th, KE, KM, self.L0
        return np.array([L0 * K, L0 - th - ke * K, ke, K])


def fim_fd(b, cv, h):
    """F = S'S/cv^2, S = d ln m / d(ln th, ln ke, ln K, L0), central differences + Richardson."""
    steps = np.array([h, h, h, h * b.L0])

    def dcol(j, hh):
        e = np.zeros(4); e[j] = hh
        return (np.log(b.sim(b.ptrue + e)) - np.log(b.sim(b.ptrue - e))) / (2 * hh)
    S = np.empty((b.t.size, 4))
    for j in range(4):
        d1, d2 = dcol(j, steps[j]), dcol(j, steps[j] / 2)
        S[:, j] = (4 * d2 - d1) / 3
    F = S.T @ S / cv ** 2
    Ci = np.linalg.inv(F)
    return np.sqrt(Ci[0, 0]), np.sqrt(Ci[1, 1]), S


# =========================================================================== FD design
def design(b, y, mid):
    psi = np.diff(y) / b.tau - b.ubar
    mm = 0.5 * (y[:-1] + y[1:]) if mid else y[:-1]
    X = np.column_stack([np.ones_like(mm), mm, -mm ** 2, -psi])
    return X, psi * mm


def unlump(p):
    K = p[3]
    L0 = p[0] / K
    return L0 - p[1] - p[2] * K, p[2], K, L0


def ncond(X):
    Xn = X / np.linalg.norm(X, axis=0)
    s = np.linalg.svd(Xn, compute_uv=False)
    return s[0] / s[-1], s[-1] ** 2


# =========================================================================== sympy corrections
A_, B_, TT, UU, SS = sp.symbols("a b T U s")
PP = sp.symbols("p1 p2 p3 p4")


def corrected(expr):
    P = sp.Poly(sp.expand(expr), A_, B_)
    out = 0
    for (i, j), c in P.terms():
        out += c * A_ ** i * B_ ** j * sp.exp(-sp.Rational(i * i + j * j, 2) * SS ** 2)
    return out


class Sym:
    def __init__(self, mid):
        m = (A_ + B_) / 2 if mid else A_
        psi = (B_ - A_) / TT - UU
        X = [sp.Integer(1), m, -m ** 2, -psi]
        z = psi * m
        args = (A_, B_, TT, UU, SS)
        self.row = [sp.lambdify(args, corrected(x) + 0 * A_, "numpy") for x in X]
        self.zrow = sp.lambdify(args, corrected(z), "numpy")
        self.MM = {}
        for i in range(4):
            for j in range(i, 4):
                self.MM[(i, j)] = sp.lambdify(args, corrected(X[i] * X[j]) + 0 * A_, "numpy")
        self.Mz = [sp.lambdify(args, corrected(X[i] * z), "numpy") for i in range(4)]
        e = z - sum(PP[i] * X[i] for i in range(4))
        ga = (A_, B_, TT, UU) + PP
        self.e = sp.lambdify(ga, e, "numpy")
        self.ea = sp.lambdify(ga, sp.diff(e, A_), "numpy")
        self.eb = sp.lambdify(ga, sp.diff(e, B_), "numpy")

    def rows(self, b, y, s):
        a, bb = y[:-1], y[1:]
        Xc = np.column_stack([np.broadcast_to(f(a, bb, b.tau, b.ubar, s), a.shape) for f in self.row])
        return Xc, self.zrow(a, bb, b.tau, b.ubar, s)

    def moments(self, b, y, s):
        a, bb = y[:-1], y[1:]
        M = np.empty((4, 4))
        for (i, j), f in self.MM.items():
            M[i, j] = M[j, i] = np.sum(np.broadcast_to(f(a, bb, b.tau, b.ubar, s), a.shape))
        r = np.array([np.sum(f(a, bb, b.tau, b.ubar, s)) for f in self.Mz])
        return M, r


SYM = {True: Sym(True), False: Sym(False)}


# =========================================================================== IV estimators
def instruments(b, y, span):
    """Rows k with k-span+1 >= 0 ... ; lo = y_{k-(span-1)}... span=1: (k-1, k+2); span=2: (k-2, k+3)."""
    n = y.size
    k = np.arange(span, n - 1 - span)
    lo, hi = y[k - span], y[k + 1 + span]
    dt = b.t[k + 1 + span] - b.t[k - span]
    ub = np.array([b.ubar[kk - span:kk + 1 + span].mean() for kk in k])
    psit = (hi - lo) / dt - ub
    mt = 0.5 * (lo + hi)
    return k, lo, hi, np.column_stack([np.ones(k.size), mt, -mt ** 2, -psit])


def linear_estimates(b, y, s, mid):
    sym = SYM[mid]
    Xr, zr = design(b, y, mid)
    Xc, zc = sym.rows(b, y, s)
    out = {"OLS": np.linalg.lstsq(Xr, zr, rcond=None)[0]}
    k, lo, hi, W = instruments(b, y, 1)
    out["IV-JI"] = np.linalg.solve(W.T @ Xc[k], W.T @ zc[k])
    mt = 0.5 * (lo + hi)
    Kw = np.median(y)
    w = 1.0 / (mt ** 2 * (mt + Kw) ** 2)
    out["wIV"] = np.linalg.solve((W * w[:, None]).T @ Xc[k], (W * w[:, None]).T @ zc[k])
    cols = [np.ones(k.size), lo, hi, lo ** 2, hi ** 2, lo * hi]
    if np.ptp(b.ubar[k]) > 0:
        cols.append(b.ubar[k])
    Wo = np.column_stack(cols)
    Q, _ = np.linalg.qr(Wo)
    PX = Q @ (Q.T @ Xc[k])
    out["2SLS"] = np.linalg.lstsq(PX, zc[k], rcond=None)[0]
    # conditional first stage, the investigator's definition re-coded
    P3 = PX[:, :3]
    dl = np.linalg.lstsq(P3, Xc[k, 3], rcond=None)[0]
    v = Xc[k, 3] - Xc[k, :3] @ dl
    res = v - Q @ (Q.T @ v)
    out["_sw"] = 1 - np.sum(res ** 2) / np.sum((v - v.mean()) ** 2)
    k2, _, _, W2 = instruments(b, y, 2)
    out["IV-wide"] = np.linalg.solve(W2.T @ Xc[k2], W2.T @ zc[k2])
    M, r = sym.moments(b, y, s)
    out["BCLS"] = np.linalg.solve(M, r)
    return out


def summary(ths, th0):
    v = np.asarray(ths, float)
    v = v[np.isfinite(v)]
    pos = v > 0
    neg = 1 - pos.mean()
    if pos.sum() < 3:
        return np.nan, np.nan, np.nan, neg, v.size
    e = np.log(v[pos] / th0)
    return (float(np.median(e)), float(np.sqrt(np.mean(e ** 2))),
            float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16))), float(neg), v.size)


# =========================================================================== GTLS-t, physical coordinates
class GT:
    def __init__(self, b, y, mid):
        self.b, self.y, self.sym = b, y, SYM[mid]
        self.a, self.bb = y[:-1], y[1:]
        self.v = y ** 2              # cv^2 drops out of the argmin
        X, z = design(b, y, mid)
        self.X, self.z = X, z

    @staticmethod
    def lump(x):
        th, ke, K, L0 = 100 * x[0], np.exp(x[1]), np.exp(x[2]), 100 * x[3]
        return np.array([L0 * K, L0 - th - ke * K, ke, K])

    def wres(self, x, tri=True):
        return self.wres_p(self.lump(x), tri)

    def wres_p(self, p, tri=True):
        """Whitened residual at LUMPED p (no constraint on any coordinate)."""
        b = self.b
        args = (self.a, self.bb, b.tau, b.ubar, *p)
        e = self.sym.e(*args)
        ga = np.broadcast_to(self.sym.ea(*args), e.shape)
        gb = np.broadcast_to(self.sym.eb(*args), e.shape)
        d = ga ** 2 * self.v[:-1] + gb ** 2 * self.v[1:]
        if not tri:
            r = e / np.sqrt(d)
            return r if np.all(np.isfinite(r)) else np.full(e.size, 1e6)
        off = gb[:-1] * ga[1:] * self.v[1:-1]
        Om = np.diag(d) + np.diag(off, 1) + np.diag(off, -1)
        try:
            Lc = np.linalg.cholesky(Om)
        except np.linalg.LinAlgError:
            return np.full(e.size, 1e6)
        r = solve_triangular(Lc, e, lower=True)
        if not np.all(np.isfinite(r)):
            return np.full(e.size, 1e6)
        return r

    def cost(self, x, tri=True):
        r = self.wres(x, tri)
        return float(r @ r)

    def grid_starts(self, nbest=3, tri=True):
        m, psi = self.X[:, 1], -self.X[:, 3]
        st = []
        for ke in np.exp(np.linspace(np.log(0.3), np.log(30), 5)):
            for K in np.exp(np.linspace(np.log(2), np.log(5000), 8)):
                lhs = self.z + ke * m ** 2 + K * psi + ke * K * m
                Z = np.column_stack([K + m, -m])
                sw = 1 / (np.abs(m) * (np.abs(m) + K))
                L0, th = np.linalg.lstsq(Z * sw[:, None], lhs * sw, rcond=None)[0]
                x = np.array([th / 100, np.log(ke), np.log(K), L0 / 100])
                st.append((self.cost(x, tri), x))
        st.sort(key=lambda q: q[0])
        return [x for _, x in st[:nbest]]

    def fit(self, starts, tri=True):
        best = None
        lb = [-np.inf, -30, -30, -np.inf]
        ub = [np.inf, 30, 30, np.inf]
        for x0 in starts:
            try:
                s = least_squares(lambda x: self.wres(x, tri), x0, method="trf", bounds=(lb, ub),
                                  xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=400)
            except (ValueError, np.linalg.LinAlgError):
                continue
            c = float(s.fun @ s.fun)
            if best is None or c < best[0]:
                best = (c, s.x)
        return best


    def fit_free(self, starts_p, tri=True):
        """Unconstrained fit in LUMPED coordinates (K_m and k_e may go negative)."""
        best = None
        for p0 in starts_p:
            p0 = np.asarray(p0, float)
            if not np.all(np.isfinite(p0)):
                continue
            sc = np.maximum(np.abs(p0), 1e-6)
            try:
                s = least_squares(lambda q: self.wres_p(q * sc, tri), p0 / sc, method="trf",
                                  xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=200)
            except (ValueError, np.linalg.LinAlgError):
                continue
            c = float(s.fun @ s.fun)
            if np.isfinite(c) and (best is None or c < best[0]):
                best = (c, s.x * sc)
        return best


def truth_x(b):
    return np.array([b.th / 100, np.log(KE), np.log(KM), b.L0 / 100])


# =========================================================================== main
def main():
    print("=" * 100)
    print("  RECHECK OF t6c_iv.py (item 4), independent code: DOP853 + finite-difference FIM, sympy")
    print("  corrections, own IV/BCLS/GTLS-t, new seeds")
    print("=" * 100)
    print(f"  L (convention A) = {LT:.4f}; L0 = L/2 = {0.5*LT:.4f}")
    B = {nm: Bench(nm) for nm in "ABC"}

    # ---------------------------------------------------------------- R1 FIM
    print("\n  R1. FIM sd ln theta, no priors, from central differences of ln m (Richardson), two steps")
    FIMSD = {}
    for nm, b in B.items():
        r1 = fim_fd(b, b.cv, 1e-3)
        r2 = fim_fd(b, b.cv, 2e-4)
        FIMSD[nm] = r2[0]
        print(f"     bench {nm}: m {b.m.min():.1f}-{b.m.max():.1f}, {b.t.size} samples, cv {b.cv:g}; "
              f"sd ln theta {r1[0]:.4g} (h=1e-3) / {r2[0]:.4g} (h=2e-4); sd ln k_e {r2[1]:.4g}; "
              f"at cv 0.1%: {r2[0]*1e-3/b.cv:.4g}")
    assert abs(FIMSD["B"] / 0.08927 - 1) < 0.02, FIMSD["B"]
    assert abs(FIMSD["C"] / 0.3218 - 1) < 0.02, FIMSD["C"]

    # ---------------------------------------------------------------- R2 noiseless design
    print("\n  R2. Noiseless FD design X0 = [1, m, -m^2, -psi0] from my own simulation")
    for nm, b in B.items():
        for mid in (False, True):
            X0, z0 = design(b, b.m, mid)
            c, lmin = ncond(X0)
            p0 = np.linalg.lstsq(X0, z0, rcond=None)[0]
            th = unlump(p0)[0]
            lab = "mid " if mid else "left"
            th_s = f"{np.log(th/b.th):+.3f}" if th > 0 else f"theta {th:.3g} <= 0"
            print(f"     {nm} {lab}: cond(X0) {c:.4g}, smallest normalised eigenvalue {lmin:.3g}; "
                  f"OLS on noiseless samples: ln(theta/truth) {th_s}")

    # ---------------------------------------------------------------- R3 OLS limit, direct MC
    print("\n  R3. OLS limit by DIRECT Monte Carlo averaging of X'X and X'z (midpoint)")
    NMC = {"A": 40000, "B": 20000, "C": 20000}
    for nm, b in B.items():
        rng = np.random.default_rng(91000 + ord(nm))
        X0, z0 = design(b, b.m, True)
        M0, r0 = X0.T @ X0, X0.T @ z0
        SM, Sr = np.zeros((4, 4)), np.zeros(4)
        Ks = []
        for i in range(NMC[nm]):
            y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
            X, z = design(b, y, True)
            SM += X.T @ X; Sr += X.T @ z
            if i < 2000:
                Ks.append(np.linalg.lstsq(X, z, rcond=None)[0][3])
        ME, rE = SM / NMC[nm], Sr / NMC[nm]
        plim = np.linalg.solve(ME, rE)
        Kn = (rE[3] - r0[3]) / (ME[3, 3] - M0[3, 3])
        share = (ME[3, 3] - M0[3, 3]) / ME[3, 3]
        thl = unlump(plim)[0]
        print(f"     {nm}: {NMC[nm]} records: K_m noise-only {Kn:.1f}, OLS limit {plim[3]:.1f}, median of "
              f"2000 OLS fits {np.median(Ks):.1f} (truth 100); limit coefficients/truth "
              + " ".join(f"{q:+.3f}" for q in plim / b.lumped())
              + f"; theta limit {thl:.4g}; noise share of psi-psi entry {share:.4f}")

    # ---------------------------------------------------------------- R5-pre partial SNR, closed form
    print("\n  partial SNR of the psi column, closed-form lognormal variances (midpoint)")
    for nm, b in B.items():
        X0, _ = design(b, b.m, True)
        psi0 = -X0[:, 3]
        s2 = b.cv ** 2
        vy = b.m ** 2 * np.exp(s2) * np.expm1(s2)
        vpsi = (vy[:-1] + vy[1:]) / b.tau ** 2
        Z = X0[:, :3]
        rs = psi0 - Z @ np.linalg.lstsq(Z, psi0, rcond=None)[0]
        rho = rs @ rs / vpsi.sum()
        Om = np.diag(vpsi) + np.diag(-vy[1:-1] / (b.tau[:-1] * b.tau[1:]), 1) \
            + np.diag(-vy[1:-1] / (b.tau[:-1] * b.tau[1:]), -1)
        Lc = np.linalg.cholesky(Om)
        Zw, pw = solve_triangular(Lc, Z, lower=True), solve_triangular(Lc, psi0, lower=True)
        rw = pw - Zw @ np.linalg.lstsq(Zw, pw, rcond=None)[0]
        print(f"     {nm}: unweighted {rho:.3g}; MA(1)-whitened {rw @ rw / psi0.size:.3g} per row")

    # ---------------------------------------------------------------- R4 linear corrected estimators
    print("\n  R4. Linear corrected estimators, own sympy corrections, NEW seeds, 300 records each")
    print("     cells: median ln(theta/truth) over theta>0 / RMSE / robust (16-84 half range) / frac theta<=0")
    # exactness check of the sympy correction, direct
    b = B["B"]
    rng = np.random.default_rng(5150)
    X0, z0 = design(b, b.m, True)
    M0 = X0.T @ X0
    Mc = []
    for _ in range(400):
        y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
        Mc.append(SYM[True].moments(b, y, b.cv)[0])
    Mc = np.array(Mc)
    se = Mc.std(0, ddof=1) / np.sqrt(len(Mc))
    zz = np.where(se > 0, (Mc.mean(0) - M0) / np.where(se > 0, se, 1), 0.0)
    iu = np.triu_indices(4)
    print(f"     sympy-corrected X'X vs noiseless, B mid, 400 records: max |z| {np.max(np.abs(zz[iu])):.2f};"
          f" corrected psi-psi sd/value {Mc[:, 3, 3].std(ddof=1)/M0[3, 3]:.3g}")
    assert np.max(np.abs(zz[iu])) < 5
    # at cv = 0, every corrected estimator must equal OLS on the noiseless samples
    L0e = linear_estimates(b, b.m, 0.0, True)
    print("     cv = 0 limits, B mid (ln theta/truth): " + ", ".join(
        f"{k} {np.log(unlump(L0e[k])[0]/b.th):+.4f}" for k in ("OLS", "IV-JI", "wIV", "2SLS", "BCLS", "IV-wide")))
    LIN = ("OLS", "IV-JI", "wIV", "2SLS", "IV-wide", "BCLS")
    LINRES = {}
    for nm, cv, mid in (("B", 0.01, True), ("C", 0.01, True), ("A", 0.05, True), ("B", 0.01, False),
                        ("B", 0.001, True), ("C", 0.001, True), ("B", 0.0001, True)):
        b = B[nm]
        rng = np.random.default_rng(int(31337 + 1000 * ord(nm) + 1e6 * cv + 7 * mid))
        th = {k: [] for k in LIN}
        ke = {k: [] for k in LIN}
        sw = []
        for _ in range(300):
            y = b.m * np.exp(rng.normal(0, cv, b.m.size))
            with np.errstate(all="ignore"):
                L = linear_estimates(b, y, cv, mid)
            for k in LIN:
                t_, k_, _, _ = unlump(L[k])
                th[k].append(t_); ke[k].append(k_)
            sw.append(L["_sw"])
        LINRES[(nm, cv, mid)] = {k: summary(th[k], b.th) for k in LIN}
        print(f"     bench {nm}, cv {100*cv:g}%, {'midpoint' if mid else 'left point'}: (SW conditional R^2 "
              f"median {np.median(sw):.4f})")
        for k in LIN:
            md, rm, rb, ng, n = LINRES[(nm, cv, mid)][k]
            kk = np.asarray(ke[k]); kk = kk[kk > 0]
            print(f"        {k:8} {md:+8.3f} / {rm:7.3f} / {rb:7.3f} / {ng:6.1%}   ln k_e median "
                  f"{np.median(np.log(kk/KE)):+.3f}  (n {n})")

    # ---------------------------------------------------------------- R6 lag structure of IV moments
    print("\n  R6. Lag structure of the IV-JI moments g_k = w_k (z^_k - x^_k'p_true), B mid, 3000 records")
    b = B["B"]
    rng = np.random.default_rng(777)
    lt = b.lumped()
    Xc0, zc0 = SYM[True].rows(b, b.m, 0.0)
    acc = {l: [] for l in range(1, 6)}
    for _ in range(3000):
        y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
        Xc, zc = SYM[True].rows(b, y, b.cv)
        k, lo, hi, W = instruments(b, y, 1)
        e = (zc - Xc @ lt)[k] - (zc0 - Xc0 @ lt)[k]   # remove the deterministic discretisation part
        g = W[:, 3] * e                                 # the -psi~ instrument's moment
        g = g / g.std()
        for l in acc:
            acc[l].append(np.mean(g[l:] * g[:-l]))
    print("     mean lag-l product of the standardised psi~ moment (z-score): " + ", ".join(
        f"l={l}: {np.mean(v):+.4f} (z {np.mean(v)/(np.std(v, ddof=1)/np.sqrt(len(v))):+.1f})"
        for l, v in acc.items()))

    # ---------------------------------------------------------------- R5 GTLS-t
    print("\n  R5. GTLS-t (tridiagonal covariance), physical coordinates, own starts; truth start as a")
    print("      DIAGNOSTIC only (reported separately, never used for the estimate)")
    for nm, mid in (("B", False), ("C", False), ("B", True), ("C", True)):
        b = B[nm]
        g = GT(b, b.m, mid)
        f = g.fit(g.grid_starts() + [truth_x(b)])
        th = 100 * f[1][0]
        print(f"     noiseless limit, {nm} {'mid ' if mid else 'left'}: ln(theta/truth) {np.log(th/b.th):+.3f}, "
              f"ln k_e {f[1][1]-np.log(KE):+.3f}, ln K_m {f[1][2]-np.log(KM):+.3f}")
    # the investigator's C left-point GTLS-t 'noiseless limit' (ln theta +1.829). Its parameters,
    # read from t6c_iv.noiseless_limits in a scratch diagnostic, are used ONLY as a start point here:
    # is that point a minimum of the GTLS-t cost?
    b = B["C"]
    g = GT(b, b.m, False)
    xI = np.array([933.8388 / 100, np.log(4.548065), np.log(966.2467), 141.3539 / 100])
    cI = g.cost(xI)
    gr = np.array([(g.cost(xI + 1e-6 * e) - g.cost(xI - 1e-6 * e)) / 2e-6 for e in np.eye(4)])
    fI = g.fit([xI])
    print(f"     C left, at the investigator's limit (theta 934, k_e 4.55, K_m 966): cost {cI:.6g}, gradient "
          + " ".join(f"{q:+.2e}" for q in gr) + f"; descending from it ends at cost {fI[0]:.6g}, theta "
          f"{100*fI[1][0]:.4g} (ln {np.log(fI[1][0]*100/b.th):+.3f}), k_e {np.exp(fI[1][1]):.3g}, K_m {np.exp(fI[1][2]):.4g}")
    ends = []
    for ke0 in (1.0, 5.0):
        for K0 in (100.0, 300.0, 1000.0):
            f = g.fit([np.array([9.34, np.log(ke0), np.log(K0), b.L0 / 100])])
            ends.append((f[0], 100 * f[1][0], np.exp(f[1][1])))
    ends = np.array(ends)
    print(f"     C left, six starts at theta 934: end costs {ends[:, 0].min():.6g}-{ends[:, 0].max():.6g}, end theta "
          f"{ends[:, 1].min():.3g}-{ends[:, 1].max():.3g}, end k_e {ends[:, 2].min():.3g}-{ends[:, 2].max():.3g}"
          " (a flat valley, not an isolated minimum)")
    assert fI[0] < cI
    NG = {"B": 300, "C": 200}
    for nm in ("B", "C"):
        b = B[nm]
        rng = np.random.default_rng(424242 + ord(nm))
        ths, kes, gap, tt = [], [], [], time.time()
        for r in range(NG[nm]):
            y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
            g = GT(b, y, True)
            with np.errstate(all="ignore"):
                f = g.fit(g.grid_starts())
                ft = g.fit([truth_x(b)])
            ths.append(100 * f[1][0]); kes.append(np.exp(f[1][1]))
            gap.append((f[0] - ft[0]) / ft[0])
        md, rm, rb, ng, n = summary(ths, b.th)
        mk, rmk, rbk, _, _ = summary(kes, KE)
        gap = np.array(gap)
        print(f"     bench {nm} mid, cv 1%, {NG[nm]} records: ln theta {md:+.3f} / RMSE {rm:.3f} / robust {rb:.3f} / "
              f"theta<=0 {ng:.1%}; ln k_e {mk:+.3f} / {rmk:.3f} / {rbk:.3f}; FIM sd {FIMSD[nm]:.4g}, "
              f"robust/FIM {rb/FIMSD[nm]:.3f}")
        print(f"        truth-start diagnostic: no-truth fit's cost exceeds the truth-start fit's by >1e-6 "
              f"(relative) in {np.mean(gap > 1e-6):.1%} of records, by >1e-2 in {np.mean(gap > 1e-2):.1%}; "
              f"({time.time()-tt:.0f} s)")
    # ---------------------------------------------------------------- R7 docstring formulas
    print("\n  R7. The investigator's hand-derived second-order noise terms (left point), checked symbolically")
    m1, m2, s_ = sp.symbols("m m2 s", positive=True)

    def expect(expr):
        P = sp.Poly(sp.expand(expr), A_, B_)
        return sum(c * m1 ** i * m2 ** j * sp.exp(sp.Rational(i * i + j * j, 2) * s_ ** 2)
                   for (i, j), c in P.terms())

    def noiseless(expr):
        return sp.expand(expr).subs({A_: m1, B_: m2})

    psi = (B_ - A_) / TT - UU
    checks = (("E[psi^2] - psi0^2", psi ** 2,
               s_ ** 2 * (2 * m1 ** 2 - 2 * m1 * m2 + 2 * m2 ** 2) / TT ** 2 - s_ ** 2 * UU * (m2 - m1) / TT, False),
              ("E[psi m] - psi0 m", psi * A_,
               s_ ** 2 * (m1 * m2 - 2 * m1 ** 2) / TT - s_ ** 2 * UU * m1 / 2, False),
              ("E[psi^2 m] - psi0^2 m (1/tau^2 part)", psi ** 2 * A_,
               s_ ** 2 * (sp.Rational(5, 2) * m1 * m2 ** 2 - 5 * m1 ** 2 * m2 + sp.Rational(9, 2) * m1 ** 3) / TT ** 2,
               True))
    for lab, ex, claimed, u0 in checks:
        d = expect(ex) - noiseless(ex)
        if u0:
            d = d.subs(UU, 0)
        ser = sp.series(d, s_, 0, 3).removeO()
        diff = sp.simplify(sp.expand(ser - claimed))
        print(f"     {lab}: second-order term minus the docstring's = {diff}")
        assert diff == 0

    # ---------------------------------------------------------------- R8 does GTLS need the constraint?
    print("\n  R8. Is the K_m, k_e > 0 constraint needed, or did the unconstrained run start in K_m < 0?")
    print("      per record: constrained fit (physical coordinates, grid starts); unconstrained fits in")
    print("      lumped coordinates from (i) linear starts OLS, BCLS, wIV, IV-JI (the investigator used OLS,")
    print("      BCLS, wBCLS, wIV), (ii) the positive grid starts, (iii) all of these plus the constrained optimum")
    for nm in ("B", "C"):
        b = B[nm]
        rng = np.random.default_rng(55000 + ord(nm))
        NR = 20
        for tri in (False, True):
            k_lin, k_grid, k_all, lower, kstart = [], [], [], [], []
            for _ in range(NR):
                y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
                with np.errstate(all="ignore"):
                    L = linear_estimates(b, y, b.cv, True)
                    g = GT(b, y, True)
                    gs = g.grid_starts(tri=tri)
                    fc = g.fit(gs, tri)
                    lin = [L[k] for k in ("OLS", "BCLS", "wIV", "IV-JI")]
                    kstart.append(np.mean([q[3] <= 0 for q in lin]))
                    fl = g.fit_free(lin, tri)
                    fg = g.fit_free([g.lump(x) for x in gs], tri)
                    fa = g.fit_free(lin + [g.lump(x) for x in gs] + [g.lump(fc[1])], tri)
                k_lin.append(fl[1][3]); k_grid.append(fg[1][3]); k_all.append(fa[1][3])
                lower.append(fa[0] < fc[0] * (1 - 1e-6) and fa[1][3] <= 0)
            print(f"     bench {nm} mid, {'GTLS-t' if tri else 'GTLS-d'}, {NR} records: linear starts with K_m <= 0 "
                  f"{np.mean(kstart):.0%}; unconstrained K_m <= 0 from (i) {np.mean(np.array(k_lin) <= 0):.0%}, "
                  f"(ii) {np.mean(np.array(k_grid) <= 0):.0%}, (iii) {np.mean(np.array(k_all) <= 0):.0%}; a K_m <= 0 "
                  f"minimum BELOW the constrained one in {np.mean(lower):.0%}")

    # ---------------------------------------------------------------- R9 seed sensitivity, GMM kernel
    print("\n  R9a. Sampling spread of the theta_hat <= 0 fraction, B mid, cv 1%, 300 records per seed set.")
    print("       Seed 7110 regenerates the investigator's own B-midpoint records (its job seed; my m equals")
    print("       its m to ~1e-9, and y = m exp(N(0, cv^2)) is drawn the same way), the others are new")
    b = B["B"]
    for seed in (7110, 12345, 2026):
        rng = np.random.default_rng(seed)
        tb, ti = [], []
        for _ in range(300):
            y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
            with np.errstate(all="ignore"):
                L = linear_estimates(b, y, b.cv, True)
            tb.append(unlump(L["BCLS"])[0]); ti.append(unlump(L["IV-JI"])[0])
        print(f"       seed {seed}: BCLS theta<=0 {np.mean(np.array(tb) <= 0):.1%}, IV-JI {np.mean(np.array(ti) <= 0):.1%}"
              f" (binomial sd at p = 0.5 and n = 300: {np.sqrt(0.25/300):.1%})")

    print("\n  R9b. GMM weight: Bartlett lags 1-3 (the investigator's) against the uniform kernel over lags 1")
    print("       (exact MA(1)) and 1-2, from 2SLS residuals; pseudo-inverse when the uniform S is not PD")

    def kern(g, w):
        n = g.shape[0]
        S = g.T @ g / n
        for l, wl in enumerate(w, 1):
            G = g[l:].T @ g[:-l] / n
            S = S + wl * (G + G.T)
        return S
    for nm in ("B", "C"):
        b = B[nm]
        rng = np.random.default_rng(8080 + ord(nm))
        res = {k: [] for k in ("2SLS", "Bartlett 1-3", "uniform 1", "uniform 1-2")}
        pdk = {k: [] for k in res}
        for _ in range(300):
            y = b.m * np.exp(rng.normal(0, b.cv, b.m.size))
            Xc, zc = SYM[True].rows(b, y, b.cv)
            k, lo, hi, _ = instruments(b, y, 1)
            cols = [np.ones(k.size), lo, hi, lo ** 2, hi ** 2, lo * hi]
            if np.ptp(b.ubar[k]) > 0:
                cols.append(b.ubar[k])
            Wo = np.column_stack(cols)
            Q, _ = np.linalg.qr(Wo)
            X, z = Xc[k], zc[k]
            p2 = np.linalg.lstsq(Q @ (Q.T @ X), z, rcond=None)[0]
            res["2SLS"].append(unlump(p2)[0]); pdk["2SLS"].append(True)
            g = Wo * (z - X @ p2)[:, None]
            for nmk, w in (("Bartlett 1-3", (0.75, 0.5, 0.25)), ("uniform 1", (1.0,)), ("uniform 1-2", (1.0, 1.0))):
                S = kern(g, w)
                okpd = np.linalg.eigvalsh(S).min() > 0
                Si = np.linalg.inv(S) if okpd else np.linalg.pinv(S)
                A = X.T @ Wo @ Si @ Wo.T @ X
                with np.errstate(all="ignore"):
                    pg = np.linalg.solve(A, X.T @ Wo @ Si @ Wo.T @ z)
                res[nmk].append(unlump(pg)[0]); pdk[nmk].append(okpd)
        for kk, v in res.items():
            md, rm, rb, ng, n = summary(v, b.th)
            print(f"       {nm} mid {kk:13}: {md:+7.3f} / {rm:6.3f} / {rb:6.3f} / theta<=0 {ng:6.1%}  (S PD in "
                  f"{np.mean(pdk[kk]):.0%})")

    print(f"\n  runtime {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
