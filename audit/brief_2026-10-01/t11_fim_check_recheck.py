"""
Adversarial recheck of t11_fim_check.py (follow-up item 6).

Shares NO code with t11_fim_check.py, core.py, dyn.py or any t*.py module: none
of them is imported. The footing (convention A, the monthly record, the dose
train, sigma, the priors) is re-derived here from the repository's definitions
(audit/n4_theta_scope.py, src/design/run.py, src/model_A/run.py, read, not
imported), and checked against numbers those scripts PRINTED (n4's 3.0673,
t5's sd5 values, t2's cohort FIM), by assertions that can fail.

Routes that differ from the investigator's:
  * information matrix built in LINEAR coordinates (theta, K_m, b, k_e), the
    cohort one as the full (2 + 2N)-square matrix inverted directly (no Schur
    complement), and the ln-k_e prior carried into k_e by the delta method;
  * the MAP minimised by bounded Nelder-Mead (and L-BFGS-B) on the projected
    2-D cost written directly from residual arrays, with a coarse grid start;
    and, on a subset, by a FULL-dimensional L-BFGS-B over all 2 + 2N
    parameters (theta linear, bounded) from truth-free starts;
  * an independent seed scheme (numpy SeedSequence keyed on (part, N, k)),
    with fresh intake draws, so the Monte Carlo shares no random numbers with
    the investigator's; and, separately, the investigator's own draws
    re-created here so that its printed numbers can be reproduced by a
    different estimator;
  * a direct test of the claimed mechanism (positivity of the OTHER patients,
    acting through the shared k_e), which the investigator inferred but did
    not test.

Run from this directory:
  python3 -W error t11_fim_check_recheck.py > out/t11_fim_check_recheck.txt 2>&1
"""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import time

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

T_START = time.time()


def el():
    return f"[{time.time() - T_START:.0f} s]"


# =====================================================================
#  Footing, re-derived (not imported)
# =====================================================================
TH0, KM, KE, MB0 = 3000.0, 100.0, 5.0, 10.0
TH = 0.05 * TH0
LN_TH = np.log(TH)
L = TH0 * MB0 / (KM + MB0) + KE * MB0          # convention A: held as capacity falls
LF, CV = 0.05, 0.05                             # load-record scatter, assay CV
KE_SD = 0.8
S2 = np.log1p(0.10 ** 2)                        # lognormal day-to-day load, CV 10%
T = 30.4 * np.arange(25)
LAM = np.log(2.0) / 3.5
DOSE = 6000.0
TH_HI = 1e5
LK_BOX = (np.log(5.0), np.log(5000.0))
LKE_BOX = (np.log(1e-3), np.log(500.0))
CHI_HALF = 0.5 * 3.841458820694124


def increment(t):
    """Dose train: 6000 every 14 d, 3.5 d half-life, 200 doses, summed directly."""
    tk = 14.0 * np.arange(200)
    dt = np.asarray(t, float)[:, None] - tk[None, :]
    on = dt >= 0
    return np.where(on, DOSE * np.exp(-LAM * np.where(on, dt, 0.0)), 0.0).sum(1)


def mss(th, K, ke, load, inc):
    """Positive root of (th + inc) m/(K + m) + ke m = load (quadratic in m)."""
    B = ke * K + th + inc - load
    return (-B + np.sqrt(B * B + 4.0 * ke * load * K)) / (2.0 * ke)


INC = increment(T)
MBAR = float(mss(TH, KM, KE, L, 0.0))
SLOPE = TH * KM / (KM + MBAR) ** 2 + KE


def sig(lf=LF):
    """n4's flux-space sigma: load-record term in quadrature with the assay term at m_bar."""
    return float(np.hypot(lf * L, SLOPE * CV * MBAR))


SIG = sig()


def check_footing():
    print("FOOTING (re-derived)")
    m_norm = float(mss(TH0, KM, KE, L, 0.0))
    s = TH * MBAR / (KM + MBAR) / L
    resid = (TH + INC) * mss(TH, KM, KE, L, INC) / (KM + mss(TH, KM, KE, L, INC)) + KE * mss(TH, KM, KE, L, INC) - L
    print(f"  L = {L:.4f}; m at normal {m_norm:.4f}; m at 5% untreated {MBAR:.4f}; s = {s:.4f}; "
          f"sigma = {SIG:.4f}; max |balance residual| {np.abs(resid).max():.1e}")
    assert abs(L - 322.73) < 5e-3 and abs(m_norm - 10.0) < 1e-9 and abs(MBAR - 54.02) < 5e-3
    assert abs(s - 0.163) < 5e-4 and np.abs(resid).max() < 1e-9


# ------------------------------------------------------------- records
def loads_repo(seed):
    """The repository's intake draw for seed s (same generator call as the footing)."""
    return L * np.exp(np.random.default_rng(seed).normal(-0.5 * S2, np.sqrt(S2), T.size))


def loads_from(u):
    return L * np.exp(-0.5 * S2 + np.sqrt(S2) * u)


def marker(loads):
    return mss(TH, KM, KE, loads, INC)


def ytrue(m, inc):
    return (TH + inc) * m / (KM + m) + KE * m


# =====================================================================
#  Information matrices in LINEAR coordinates
# =====================================================================
def fim_single(m, inc, sg, ske, sb, flux_Y=None):
    """(theta, K_m, b, k_e); ln k_e prior -> k_e prior sd ske * k_e. Returns (F, sd ln theta)."""
    g = m / (KM + m)
    J = np.stack([g, -(TH + inc) * m / (KM + m) ** 2, np.ones_like(m), m], 1) / np.asarray(sg)[..., None]
    F = J.T @ J
    F[2, 2] += sb ** -2
    F[3, 3] += (ske * KE) ** -2
    if flux_Y is not None:
        g0 = MBAR / (KM + MBAR)
        F0 = TH * g0
        sdf = F0 * np.hypot(flux_Y, KM / (KM + MBAR) * CV)
        row = np.array([g0, -TH * MBAR / (KM + MBAR) ** 2, 0.0, 0.0])
        F += np.outer(row, row) / sdf ** 2
    P = np.linalg.inv(F)
    return F, float(np.sqrt(P[0, 0]) / TH)


def fim_cohort(Ms, sg=SIG):
    """Full (2 + 2N) information in (K_m, k_e, theta_1..N, b_1..N), priors included."""
    N, nd = Ms.shape
    J = np.zeros((N * nd, 2 + 2 * N))
    for i in range(N):
        m = Ms[i]
        rows = slice(i * nd, (i + 1) * nd)
        J[rows, 0] = -(TH + INC) * m / (KM + m) ** 2 / sg
        J[rows, 1] = m / sg
        J[rows, 2 + i] = m / (KM + m) / sg
        J[rows, 2 + N + i] = 1.0 / sg
    F = J.T @ J
    Pp = np.zeros(2 + 2 * N)
    Pp[1] = (KE_SD * KE) ** -2
    Pp[2 + N:] = (LF * L) ** -2
    F = F + np.diag(Pp)
    return F, J, Pp


def cohort_sd(Ms):
    F, J, Pp = fim_cohort(Ms)
    P = np.linalg.inv(F)
    Sand = P @ (F - np.diag(Pp)) @ P                # centres held at the truth
    return (float(np.sqrt(P[2, 2]) / TH), float(np.sqrt(P[1, 1]) / KE),
            float(np.sqrt(Sand[2, 2]) / TH), np.sqrt(np.diag(P)[2:2 + Ms.shape[0]]) / TH, P, J)


# =====================================================================
#  The MAP, own implementation
# =====================================================================
class MAP:
    """
    Cost = 0.5 sum_ij ((y_ij - (theta_i + inc_ij) g_ij - k_e m_ij - b_i)/sg)^2
         + 0.5 ((ln k_e - kec)/ske)^2 + 0.5 sum_i ((b_i - bc_i)/sb)^2
         [+ 0.5 sum_i ((z_i - theta_i g(m0))/sdf)^2]
    theta_i in [lo_i, 1e5]; K_m and k_e in t2's box. lo may be per patient.
    For given (K_m, k_e) the (theta_i, b_i) minimiser is closed form (convex
    2-variable quadratic, box on theta only: clip, then b given theta).
    """

    def __init__(self, M, INCm, Y, sg, sb, bc, kec, ske, lo=1e-3, flux=None, fix0=None):
        self.M, self.I, self.Y = (np.atleast_2d(np.asarray(a, float)) for a in (M, INCm, Y))
        self.N = self.M.shape[0]
        self.W = np.broadcast_to(np.asarray(sg, float) ** -2.0, self.M.shape)
        self.sb = np.broadcast_to(np.asarray(sb, float), (self.N,)).copy()
        self.bc = np.broadcast_to(np.asarray(bc, float), (self.N,)).copy()
        self.kec, self.ske = float(kec), float(ske)
        self.lo = np.broadcast_to(np.asarray(lo, float), (self.N,)).copy()
        self.flux = flux
        self.fix0 = fix0
        self.ub = ((LKE_BOX[0] - self.kec) / self.ske, (LKE_BOX[1] - self.kec) / self.ske)

    def _solve(self, K, ke):
        """K scalar; ke array (n,). Returns cost (n,), theta (n, N)."""
        ke = np.atleast_1d(ke)
        G = self.M / (K + self.M)
        Z = (self.Y - self.I * G)[None] - ke[:, None, None] * self.M[None]      # (n, N, nd)
        W = self.W[None]
        A11 = (self.W * G * G).sum(-1)
        A12 = (self.W * G).sum(-1)
        A22 = self.W.sum(-1) + self.sb ** -2
        r1 = (W * G[None] * Z).sum(-1)
        r2 = (W * Z).sum(-1) + self.bc / self.sb ** 2
        if self.flux is not None:
            zf, sdf, m0 = self.flux
            g0 = m0 / (K + m0)
            A11 = A11 + g0 ** 2 / sdf ** 2
            r1 = r1 + g0 * zf / sdf ** 2
        th = (A22 * r1 - A12 * r2) / (A11 * A22 - A12 ** 2)
        th = np.minimum(np.maximum(th, self.lo), TH_HI)
        if self.fix0 is not None:
            th[:, 0] = self.fix0
        b = (r2 - A12 * th) / A22
        R = Z - th[..., None] * G[None] - b[..., None]
        c = 0.5 * (W * R * R).sum((-1, -2)) + 0.5 * (((b - self.bc) / self.sb) ** 2).sum(-1)
        if self.flux is not None:
            c = c + 0.5 * (((zf - th * g0) / sdf) ** 2).sum(-1)
        c = c + 0.5 * ((np.log(ke) - self.kec) / self.ske) ** 2
        return c, th, b

    def f(self, x):
        lk = min(max(x[0], LK_BOX[0]), LK_BOX[1])
        u = min(max(x[1], self.ub[0]), self.ub[1])
        return float(self._solve(np.exp(lk), np.exp(self.kec + self.ske * u))[0][0])

    def fit(self, starts=(), nk=26, nu=41, urange=7.0, top=4):
        lks = np.linspace(*LK_BOX, nk)
        us = np.linspace(max(-urange, self.ub[0]), min(urange, self.ub[1]), nu)
        kes = np.exp(self.kec + self.ske * us)
        Cg = np.array([self._solve(np.exp(lk), kes)[0] for lk in lks])
        idx = np.argsort(Cg.ravel())[:top]
        cand = [np.array([lks[i // nu], us[i % nu]]) for i in idx] + [np.asarray(s, float) for s in starts]
        bnds = [LK_BOX, self.ub]
        best = None
        for x0 in cand:
            x0 = np.array([np.clip(x0[0], *LK_BOX), np.clip(x0[1], *self.ub)])
            simplex = np.array([x0, x0 + [0.15, 0.0], x0 + [0.0, 0.4]])
            simplex[:, 0] = np.clip(simplex[:, 0], *LK_BOX)
            simplex[:, 1] = np.clip(simplex[:, 1], *self.ub)
            if np.linalg.matrix_rank(simplex[1:] - simplex[0]) < 2:
                simplex = np.array([x0, x0 - [0.15, 0.0], x0 - [0.0, 0.4]])
            r = minimize(self.f, x0, method="Nelder-Mead", bounds=bnds,
                         options=dict(initial_simplex=simplex, xatol=1e-9, fatol=1e-13,
                                      maxiter=4000, maxfev=8000))
            if best is None or r.fun < best.fun:
                best = r
        r2 = minimize(self.f, best.x, method="L-BFGS-B", bounds=bnds,
                      options=dict(ftol=1e-15, gtol=1e-11, maxiter=500))
        if r2.fun < best.fun:
            best = r2
        lk, u = best.x
        ke = float(np.exp(self.kec + self.ske * u))
        c, th, b = self._solve(np.exp(lk), np.array([ke]))
        return dict(cost=float(c[0]), th=th[0], b=b[0], lk=float(lk), lke=float(np.log(ke)),
                    x=np.array(best.x))


def full_dim_lbfgs(p, starts):
    """
    A different optimiser: L-BFGS-B over ALL 2 + 2N parameters, theta linear and
    box-bounded, analytic gradient. Variables (ln K, u, theta/TH, b/sb).
    Returns the best cost over the starts.
    """
    N, M, I, Y, W = p.N, p.M, p.I, p.Y, p.W
    sgi = np.sqrt(W)

    def fg(x):
        lk, u = x[0], x[1]
        t, c = x[2:2 + N], x[2 + N:]
        K, ke = np.exp(lk), np.exp(p.kec + p.ske * u)
        th, b = TH * t, p.sb * c
        G = M / (K + M)
        rd = (Y - (th[:, None] + I) * G - ke * M - b[:, None]) * sgi
        cost = 0.5 * (rd ** 2).sum() + 0.5 * u * u + 0.5 * ((c - p.bc / p.sb) ** 2).sum()
        g = np.zeros_like(x)
        g[0] = (rd * sgi * (th[:, None] + I) * K * M / (K + M) ** 2).sum()
        g[1] = -(rd * sgi * ke * M).sum() * p.ske + u
        g[2:2 + N] = -(rd * sgi * G).sum(1) * TH
        g[2 + N:] = -(rd * sgi).sum(1) * p.sb + (c - p.bc / p.sb)
        return cost, g

    bnds = ([LK_BOX, p.ub] + [(lo / TH if np.isfinite(lo) else None, TH_HI / TH) for lo in p.lo]
            + [(None, None)] * N)
    # gradient check on the first start
    x0 = np.asarray(starts[0], float)
    c0, g0 = fg(x0)
    h = 1e-6
    gn = np.array([(fg(x0 + h * e)[0] - fg(x0 - h * e)[0]) / (2 * h) for e in np.eye(x0.size)[:4]])
    assert np.allclose(gn, g0[:4], rtol=1e-4, atol=1e-4 * (1 + np.abs(g0[:4]).max())), (gn, g0[:4])
    best = np.inf
    for x0 in starts:
        r = minimize(fg, np.asarray(x0, float), jac=True, method="L-BFGS-B", bounds=bnds,
                     options=dict(maxiter=30000, maxfun=60000, ftol=1e-16, gtol=1e-10, maxcor=30))
        best = min(best, float(r.fun))
    return best


# =====================================================================
#  Draws
# =====================================================================
def draw_own(n, k, part=1):
    """Independent stream: SeedSequence keyed on (part, n, k). Fresh intake, noise, centres."""
    rng = np.random.default_rng(np.random.SeedSequence([20261001, part, n, k]))
    U = rng.standard_normal((n, T.size))
    E = rng.standard_normal((n, T.size))
    ekc = rng.standard_normal()
    eb = rng.standard_normal(n)
    Ms = marker(loads_from(U))
    Ys = ytrue(Ms, INC) + SIG * E
    return Ms, Ys, np.log(KE) + KE_SD * ekc, LF * L * eb


def draw_theirs(n, k, rng):
    """The investigator's (and t2's) replicate k, re-created from the stated seed scheme."""
    Ms, Ys = [], []
    for i in range(n):
        s = 1000 * k + i
        m = marker(loads_repo(s))
        Ms.append(m)
        Ys.append(ytrue(m, INC) + np.random.default_rng(10_000 + s).normal(0.0, np.full(T.size, SIG)))
    kec = np.log(KE) + rng.normal(0, KE_SD)
    bc = np.array([rng.normal(0, LF * L) for _ in range(n)])
    return np.array(Ms), np.array(Ys), kec, bc


def problem(Ms, Ys, kec, bc, **kw):
    n = Ms.shape[0]
    return MAP(Ms, np.broadcast_to(INC, Ms.shape), Ys, SIG, LF * L, bc, kec, KE_SD, **kw)


# =====================================================================
#  Statistics
# =====================================================================
def robust(e):
    e = np.asarray(e).ravel()
    return float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16)))


def rmse(e):
    e = np.asarray(e).ravel()
    return float(np.sqrt(np.mean(e ** 2)))


def boot(e, fn, B=2000, seed=1):
    """Percentile bootstrap over rows (replicates or whole cohorts)."""
    e = np.asarray(e)
    rng = np.random.default_rng(seed)
    v = np.array([fn(e[rng.integers(0, e.shape[0], e.shape[0])]) for _ in range(B)])
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def lnerr(r):
    """ln(theta_hat/theta) from the RATIO r = theta_hat/theta (constrained fits only: r > 0)."""
    r = np.asarray(r, float)
    assert np.all(r > 0)
    return np.log(r)


def bs(e, fn, d=3):
    lo, hi = boot(e, fn)
    return f"[{lo:.{d}f}, {hi:.{d}f}]"


# =====================================================================
#  CHECK 1. FIM re-derived against printed values
# =====================================================================
SCEN = [
    ("current design", dict(ske=0.8, lf=0.05, Y=None), 3.0673),
    ("flux 10%", dict(ske=0.8, lf=0.05, Y=0.10), 0.1120),
    ("flux 20%", dict(ske=0.8, lf=0.05, Y=0.20), 0.2064),
    ("k_e exact", dict(ske=1e-4, lf=0.05, Y=None), 0.3879),
    ("k_e to 5%", dict(ske=np.log(1.05), lf=0.05, Y=None), 0.4447),
    ("k_e exact + load exact", dict(ske=1e-4, lf=1e-4, Y=None), 0.2044),
    ("flux 20% + k_e 1.2 + load 2%", dict(ske=np.log(1.2), lf=0.02, Y=0.20), 0.1989),
]


def scen_arrays(m, sc):
    """Add t5's untreated baseline draw (load L, no increment, m = m_bar) for flux scenarios."""
    if sc["Y"] is None:
        return m, INC
    return np.r_[MBAR, m], np.r_[0.0, INC]


def check1():
    print("\nCHECK 1. Information matrices re-derived in linear coordinates, against printed values")
    print("  (a) Task 5 scenarios, median over the repository's intake seeds 0..23")
    print(f"     {'scenario':32} {'here':>8} {'printed (t5/t11)':>17}")
    fims = {}
    for lab, sc, ref in SCEN:
        v = []
        for s in range(24):
            m, inc = scen_arrays(marker(loads_repo(s)), sc)
            v.append(fim_single(m, inc, np.full(m.size, sig(sc["lf"])), sc["ske"], sc["lf"] * L, sc["Y"])[1])
        med = float(np.median(v))
        fims[lab] = med
        print(f"     {lab:32} {med:8.4f} {ref:17.4f}")
        assert abs(med - ref) < 6e-5, (lab, med, ref)
    print("  (b) cohort FIM (full matrix, direct inverse), t2's seeds 0..N-1")
    for n, ref in ((1, 3.082), (16, 1.331), (64, 0.806)):
        Ms = np.array([marker(loads_repo(s)) for s in range(n)])
        sd = cohort_sd(Ms)
        print(f"     N = {n:2d}: sd ln theta_1 {sd[0]:.4f} (t2 printed {ref}); sd ln k_e {sd[1]:.4f}; "
              f"sandwich (centres held) {sd[2]:.4f}")
        assert abs(sd[0] - ref) < 6e-4
    return fims


# =====================================================================
#  CHECK 2. The cohort MC: investigator's draws, own estimator
# =====================================================================
def run_mc(n, reps, source, variants):
    """Returns dict of arrays. source 'theirs' or 'own'."""
    rng = np.random.default_rng(0)
    out = {v: [] for v in variants}
    out.update(fim=[], fimke=[], sand=[], kes={v: [] for v in variants}, nb=[], dfl=[], gn=[], surr=[])
    for k in range(reps):
        if source == "theirs":
            Ms, Ys, kec, bc = draw_theirs(n, k, rng)
        else:
            Ms, Ys, kec, bc = draw_own(n, k)
        sd0, sdk, sand, sds, P0, J0 = cohort_sd(Ms)
        out["fim"].append(sds)
        out["fimke"].append(sdk)
        out["sand"].append(sand)
        # linear-Gaussian surrogate on the same draws: one GN step from the truth
        w = np.r_[((Ys - ytrue(Ms, INC)) / SIG).ravel()]
        rhs = J0.T @ w
        rhs[1] += (kec - np.log(KE)) * KE / (KE_SD * KE) ** 2          # ln-k_e prior centre -> k_e
        rhs[2 + n:] += bc / (LF * L) ** 2
        d = P0 @ rhs
        out["surr"].append(d[2:2 + n] / TH)
        base = None
        for v in variants:
            kw = {}
            if v == "c3":
                kw = dict(lo=1e-3)
            elif v == "c2":
                kw = dict(lo=1e-2)
            elif v == "c5":
                kw = dict(lo=1e-5)
            elif v == "un":
                kw = dict(lo=-np.inf)
            elif v == "others":            # patient 0 free, the other N-1 constrained
                kw = dict(lo=np.r_[-np.inf, np.full(n - 1, 1e-3)])
            elif v == "self":              # only patient 0 constrained
                kw = dict(lo=np.r_[1e-3, np.full(n - 1, -np.inf)])
            elif v == "held":
                kw = dict(lo=1e-3)
            elif v == "floor":
                kw = dict(lo=1e-3, fix0=TH * np.exp(-11.5))
            if v == "held":
                p = problem(Ms, Ys, np.log(KE), np.zeros(n), **kw)
            else:
                p = problem(Ms, Ys, kec, bc, **kw)
            starts = () if base is None else (base["x"],)
            f = p.fit(starts=starts)
            if v == "c3":
                base = f
                out["nb"].append(int(np.sum(f["th"] <= 1e-3 * (1 + 1e-9))))
            if v == "floor":
                out["dfl"].append(f["cost"] - base["cost"])
                continue
            out[v].append(f["th"] / TH)
            out["kes"][v].append(f["lke"] - np.log(KE))
        # Gauss-Newton information at the constrained MAP, active patients removed
        if "c3" in variants and base["th"][0] > 1e-3 * (1 + 1e-9):
            K, ke, th = np.exp(base["lk"]), np.exp(base["lke"]), base["th"]
            act = th <= 1e-3 * (1 + 1e-9)
            nd = T.size
            J = np.zeros((n * nd, 2 + 2 * n))
            for i in range(n):
                m = Ms[i]
                rows = slice(i * nd, (i + 1) * nd)
                J[rows, 0] = -(th[i] + INC) * m / (K + m) ** 2 / SIG
                J[rows, 1] = m / SIG
                J[rows, 2 + i] = m / (K + m) / SIG
                J[rows, 2 + n + i] = 1.0 / SIG
            F = J.T @ J
            F[1, 1] += (KE_SD * ke) ** -2
            F[2 + n:, 2 + n:] += np.eye(n) * (LF * L) ** -2
            keep = [j for j in range(2 + 2 * n) if not (2 <= j < 2 + n and act[j - 2])]
            Pg = np.linalg.inv(F[np.ix_(keep, keep)])
            out["gn"].append(np.sqrt(Pg[keep.index(2), keep.index(2)]) / th[0])
    res = {}
    for kk, v in out.items():
        if kk == "kes":
            res[kk] = {a: np.array(b) for a, b in v.items()}
        else:
            res[kk] = np.array(v)
    return res


def summarise_cohort(n, R, label):
    fo = float(np.sqrt(np.mean(R["fim"][:, 0] ** 2)))
    fp = float(np.sqrt(np.mean(R["fim"] ** 2)))
    e0 = lnerr(R["c3"][:, 0])
    ep = lnerr(R["c3"])
    reps = e0.size
    print(f"\n  {label}: N = {n}, {reps} replicates")
    print(f"     FIM sd ln theta_1 (rms over own records) {fo:.3f}; pooled {fp:.3f}; "
          f"sandwich (centres held) {np.sqrt(np.mean(R['sand'] ** 2)):.3f}; FIM sd ln k_e {np.sqrt(np.mean(R['fimke'] ** 2)):.3f}")
    atb = np.mean(R["c3"][:, 0] <= 1e-3 / TH * (1 + 1e-9))
    print(f"     patient 0, bound 1e-3: RMSE {rmse(e0):.3f} {bs(e0, rmse)}  robust {robust(e0):.3f} "
          f"{bs(e0, robust)}  median {np.median(e0):+.3f}  at bound {atb:.1%}")
    for v, bd in (("c2", 1e-2), ("c3", 1e-3), ("c5", 1e-5)):
        if v in R:
            e = lnerr(R[v][:, 0])
            eP = lnerr(R[v])
            print(f"     bound {bd:.0e}: patient 0 RMSE {rmse(e):.3f}, robust {robust(e):.3f}, P84 {np.percentile(e, 84):+.3f};"
                  f" pooled RMSE {rmse(eP):.3f}, pooled robust {robust(eP):.3f};"
                  f" 0.5 (P84 - ln(bound/theta)) = {0.5 * (np.percentile(e, 84) - np.log(bd / TH)):.3f}")
    lin3 = R["c3"] - 1.0
    print(f"     pooled constrained: RMSE ln {rmse(ep):.3f} [cohort-bootstrap {bs(ep, rmse)}], robust ln "
          f"{robust(ep):.3f} {bs(ep, robust)}, median ln {np.median(ep):+.3f}; linear RMSE {rmse(lin3):.3f}, "
          f"mean linear err {np.mean(lin3):+.3f} {bs(lin3, np.mean)}")
    print(f"     cohorts with a patient at the bound {np.mean(R['nb'] > 0):.1%}; patient-fits at the bound "
          f"{np.mean(R['c3'] <= 1e-3 / TH * (1 + 1e-9)):.2%}")
    if "un" in R:
        lu = R["un"] - 1.0
        print(f"     unconstrained (pooled, linear): sd {np.std(lu):.3f} {bs(lu, np.std)}, robust {robust(lu):.3f}, "
              f"mean {np.mean(lu):+.3f}, theta_hat < 0 in {np.mean(lu < -1):.1%} "
              f"(Gaussian prediction Phi(-1/FIM) = {norm.cdf(-1 / fp):.1%})")
        print(f"     linear surrogate on the same draws (pooled): sd {np.std(R['surr']):.3f}, robust {robust(R['surr']):.3f}")
    ke = R["kes"]
    fk = float(np.sqrt(np.mean(R["fimke"] ** 2)))
    line = f"     ln k_e error: FIM +-1.96 sd = +-{1.96 * fk:.3f};"
    for v in ("c3", "un", "others", "self"):
        if v in ke and ke[v].size:
            p = np.percentile(ke[v], [2.5, 50, 97.5])
            line += f" {v} P2.5/P50/P97.5 {p[0]:+.3f}/{p[1]:+.3f}/{p[2]:+.3f};"
    print(line)
    return dict(fo=fo, fp=fp, e0=e0, ep=ep)


def check2():
    print("\nCHECK 2. COHORT MONTE CARLO ON THE INVESTIGATOR'S OWN DRAWS, WITH AN INDEPENDENT ESTIMATOR")
    print("  (draws re-created from the seed scheme stated in t11/t2: intake seed 1000k+i, noise seed")
    print("  10000+1000k+i, prior centres from default_rng(0) in t2's order; estimator: own code, Nelder-Mead)")
    out = {}
    for n, reps in ((1, 200), (64, 200)):
        R = run_mc(n, reps, "theirs", ("c3", "c2", "c5", "un"))
        out[n] = R
        S = summarise_cohort(n, R, "investigator's draws")
        if n == 64:
            e60 = lnerr(R["c3"][:60, 0])
            print(f"     first 60 replicates (= t2's table E): RMSE {rmse(e60):.3f}, robust {robust(e60):.3f}, "
                  f"median {np.median(e60):+.3f}")
            assert abs(rmse(e60) - 0.494) < 1.5e-3 and abs(robust(e60) - 0.455) < 1.5e-3
        print(f"  {el()}")
    return out


# =====================================================================
#  CHECK 3. The cohort MC on independent draws, with the mechanism test
# =====================================================================
def check3():
    print("\nCHECK 3. COHORT MONTE CARLO ON INDEPENDENT DRAWS (own seed stream, fresh intake)")
    out = {}
    plan = ((1, 400, ("c3", "c2", "c5", "un")),
            (16, 200, ("c3", "un", "others")),
            (64, 200, ("c3", "c2", "c5", "un", "others", "self", "held", "floor")))
    for n, reps, var in plan:
        R = run_mc(n, reps, "own", var)
        out[n] = R
        S = summarise_cohort(n, R, "own draws")
        if "others" in R:
            l3, lo_, lu = R["c3"][:, 0] - 1, R["others"][:, 0] - 1, R["un"][:, 0] - 1
            print(f"     MECHANISM TEST, patient 0 in linear theta (theta_hat/theta - 1):")
            print(f"       all constrained:            robust {robust(l3):.3f} {bs(l3, robust)}, sd {np.std(l3):.3f}, mean {np.mean(l3):+.3f}")
            print(f"       others constrained, 0 free: robust {robust(lo_):.3f} {bs(lo_, robust)}, sd {np.std(lo_):.3f}, mean {np.mean(lo_):+.3f},"
                  f" theta_0 < 0 in {np.mean(lo_ < -1):.1%}")
            if "self" in R:
                ls = R["self"][:, 0] - 1
                print(f"       only patient 0 constrained: robust {robust(ls):.3f} {bs(ls, robust)}, sd {np.std(ls):.3f}, mean {np.mean(ls):+.3f}")
            print(f"       none constrained:           robust {robust(lu):.3f} {bs(lu, robust)}, sd {np.std(lu):.3f}, mean {np.mean(lu):+.3f}")
            print(f"       FIM sd (linear = ln at the truth): {S['fo']:.3f}")
            if "self" in R:
                print(f"       ln theta_1 robust: all constrained {robust(lnerr(R['c3'][:, 0])):.3f}; only patient 0 "
                      f"constrained {robust(lnerr(R['self'][:, 0])):.3f}")
        if "held" in R:
            eh = lnerr(R["held"][:, 0])
            print(f"     centres held at the truth: patient 0 robust ln {robust(eh):.3f}, median {np.median(eh):+.3f}")
        if "floor" in R:
            d = R["dfl"]
            print(f"     profile at theta_1 = 150 e^-11.5 (the investigator's grid floor): -log L rise over the MAP "
                  f"below chi2/2 = {CHI_HALF:.3f} (interval open towards 0) in {np.mean(d < CHI_HALF):.1%} of {d.size} cohorts "
                  f"{bs((d < CHI_HALF).astype(float), np.mean)}; linear-Gaussian prediction "
                  f"P(theta_hat/theta < 1.96 FIM) = {norm.cdf(1.96 - 1 / S['fo']):.1%}")
        if R["gn"].size:
            print(f"     GN sd at the MAP, patient 0 (not at bound): median {np.median(R['gn']):.3f}, P10 "
                  f"{np.percentile(R['gn'], 10):.3f}, P90 {np.percentile(R['gn'], 90):.3f} (n = {R['gn'].size}); FIM {S['fo']:.3f}")
        print(f"  {el()}")
    return out


# =====================================================================
#  CHECK 4. Global optimum by a different optimiser
# =====================================================================
def check4():
    print("\nCHECK 4. IS THE PROJECTED MAP THE GLOBAL OPTIMUM? Full-dimensional L-BFGS-B, truth-free starts")
    print(f"     {'N':>3} {'reps':>5} {'max (proj - full) cost':>23} {'full lower by > 1e-6':>21}")
    for n, reps in ((1, 100), (16, 30), (64, 30)):
        worst, cnt = -np.inf, 0
        for k in range(reps):
            Ms, Ys, kec, bc = draw_own(n, k, part=4)
            p = problem(Ms, Ys, kec, bc)
            f = p.fit()
            starts = []
            for lk in (np.log(30.0), np.log(1000.0)):
                for u in (-1.0, 0.0, 1.0):
                    starts.append(np.r_[lk, u, np.full(n, 0.5), np.zeros(n)])
            starts.append(np.r_[np.log(300.0), 0.5, np.full(n, 3.0), bc / (LF * L)])
            cf = full_dim_lbfgs(p, starts)
            worst = max(worst, f["cost"] - cf)
            cnt += (f["cost"] - cf) > 1e-6
        print(f"     {n:3d} {reps:5d} {worst:23.2e} {cnt:21d}")
    print(f"  {el()}")


# =====================================================================
#  CHECK 5. Linear-Gaussian surrogate: investigator's seeds and independent seeds
# =====================================================================
def check5():
    print("\nCHECK 5. LINEAR-GAUSSIAN SURROGATE, N = 64: whitened error of theta_1 must be N(0, 1) exactly")
    print("  (conditional on the records, w is N(0, I) and the GN step is N(0, inv F): no approximation)")
    for label, reps, src in (("investigator's draws k = 0..1999", 2000, "theirs"),
                             ("own independent draws", 6000, "own")):
        rng = np.random.default_rng(0)
        z = []
        for k in range(reps):
            if src == "theirs":
                Ms, Ys, kec, bc = draw_theirs(64, k, rng)
            else:
                Ms, Ys, kec, bc = draw_own(64, k, part=5)
            F, J, Pp = fim_cohort(Ms)
            P = np.linalg.inv(F)
            rhs = J.T @ ((Ys - ytrue(Ms, INC)) / SIG).ravel()
            rhs[1] += (kec - np.log(KE)) * KE / (KE_SD * KE) ** 2
            rhs[2 + 64:] += bc / (LF * L) ** 2
            d = P @ rhs
            z.append(d[2] / np.sqrt(P[2, 2]))
        z = np.array(z)
        half = 1.96 / np.sqrt(2 * reps)
        print(f"     {label:34} reps {reps:5d}: whitened sd {np.std(z):.3f}  95% band [{1 - half:.3f}, {1 + half:.3f}]"
              f"  mean {np.mean(z):+.3f}")
    print(f"  {el()}")


# =====================================================================
#  CHECK 6. Task 5 scenarios on independent draws
# =====================================================================
def check6(fims):
    print("\nCHECK 6. TASK 5 SCENARIOS AT 5%: MC MAP on independent draws (fresh intake each replicate),")
    print("  400 replicates per scenario, common random numbers across scenarios, ratios to the FIM rms")
    print("  over the replicates' own records. Measurements simulated about the truth (as t11 and t2 do).")
    reps = 400
    res = {}
    for lab, sc, ref in SCEN[1:]:
        e, fs, dfl = [], [], []
        for j in range(reps):
            rng = np.random.default_rng(np.random.SeedSequence([20261001, 6, j]))
            u = rng.standard_normal(T.size)
            en = rng.standard_normal(T.size + 1)
            ef, eke, eb = rng.standard_normal(3)
            m25 = marker(loads_from(u))
            m, inc = scen_arrays(m25, sc)
            sg = sig(sc["lf"])
            y = ytrue(m, inc) + sg * (en if m.size == T.size + 1 else en[1:])
            sb = sc["lf"] * L
            flux = None
            if sc["Y"] is not None:
                F0 = TH * MBAR / (KM + MBAR)
                sdf = F0 * np.hypot(sc["Y"], KM / (KM + MBAR) * CV)
                flux = (np.array([F0 + sdf * ef]), np.array([sdf]), np.array([MBAR]))
            kec = np.log(KE) + sc["ske"] * eke
            p = MAP(m, inc, y, sg, sb, sb * eb, kec, sc["ske"], lo=1e-3, flux=flux)
            f = p.fit()
            e.append(f["th"][0] / TH)
            fs.append(fim_single(m, inc, np.full(m.size, sg), sc["ske"], sb, sc["Y"])[1])
            if lab in ("k_e exact", "k_e to 5%"):
                pf = MAP(m, inc, y, sg, sb, sb * eb, kec, sc["ske"], lo=1e-3, flux=flux, fix0=TH * np.exp(-11.5))
                dfl.append(pf.fit(starts=(f["x"],))["cost"] - f["cost"])
        r = np.array(e)
        le = lnerr(r)
        fr = float(np.sqrt(np.mean(np.square(fs))))
        res[lab] = dict(le=le, lin=r - 1, fr=fr)
        atb = np.mean(r <= 1e-3 / TH * (1 + 1e-9))
        ci = boot(le, rmse)
        cb = boot(le, robust)
        cl = boot(r - 1, rmse)
        print(f"     {lab:30} FIM median(0..23) {fims[lab]:.4f}, rms own {fr:.4f} | ln RMSE {rmse(le):.4f} "
              f"[{ci[0]:.4f}, {ci[1]:.4f}] ({rmse(le) / fr:.3f}) | robust {robust(le):.4f} [{cb[0]:.4f}, {cb[1]:.4f}] "
              f"({robust(le) / fr:.3f}) | lin RMSE {rmse(r - 1):.4f} [{cl[0]:.4f}, {cl[1]:.4f}] ({rmse(r - 1) / fr:.3f}) "
              f"| median {np.median(le):+.4f} | at bound {atb:.1%}")
        if dfl:
            d = np.array(dfl)
            sl = fims[lab]
            print(f"       profile open towards theta -> 0 (rise at the floor < {CHI_HALF:.3f}) in {np.mean(d < CHI_HALF):.1%} "
                  f"{bs((d < CHI_HALF).astype(float), np.mean)}; linear-Gaussian prediction Phi(1.96 - 1/FIM) = "
                  f"{norm.cdf(1.96 - 1 / sl):.1%}; investigator: 6/24 = 25.0%")
    print("  paired comparison (same random numbers): combo minus 'k_e exact + load exact'")
    a, b = res["flux 20% + k_e 1.2 + load 2%"], res["k_e exact + load exact"]
    both = np.stack([a["le"], b["le"]], 1)
    dr = lambda x: rmse(x[:, 0]) - rmse(x[:, 1])
    print(f"     ln RMSE difference {dr(both):+.4f} {bs(both, dr)}; FIM difference {fims['flux 20% + k_e 1.2 + load 2%'] - fims['k_e exact + load exact']:+.4f}")
    lc = res["flux 20% + k_e 1.2 + load 2%"]["le"]
    print(f"     combo ln RMSE {rmse(lc):.4f}, P(bootstrap RMSE <= 0.20) = "
          f"{np.mean([rmse(lc[np.random.default_rng(3 + i).integers(0, lc.size, lc.size)]) <= 0.20 for i in range(2000)]):.3f}")
    ratios = [rmse(res[l]["le"]) / res[l]["fr"] for l, _, _ in SCEN[1:] if l not in ("k_e exact", "k_e to 5%")]
    lr = [rmse(res[l]["lin"]) / res[l]["fr"] for l, _, _ in SCEN[1:]]
    print(f"     ln RMSE / FIM over the four scenarios without bound hits: {', '.join(f'{x:.3f}' for x in ratios)}")
    print(f"     linear RMSE / FIM, all six: {', '.join(f'{x:.3f}' for x in lr)}")
    order_f = sorted(fims, key=fims.get)
    order_m = sorted(res, key=lambda l: rmse(res[l]["le"]))
    print(f"     rank by FIM:     {[l for l in order_f if l != 'current design']}")
    print(f"     rank by MC RMSE: {order_m}")
    print(f"  {el()}")


def scen_map(m, inc, y, sc, eb, eke, ef, fix0=None):
    sg = sig(sc["lf"])
    sb = sc["lf"] * L
    flux = None
    if sc["Y"] is not None:
        F0 = TH * MBAR / (KM + MBAR)
        sdf = F0 * np.hypot(sc["Y"], KM / (KM + MBAR) * CV)
        flux = (np.array([F0 + sdf * ef]), np.array([sdf]), np.array([MBAR]))
    return MAP(m, inc, y, sg, sb, sb * eb, np.log(KE) + sc["ske"] * eke, sc["ske"], lo=1e-3,
               flux=flux, fix0=fix0)


def pl_ends(make, fmap):
    """95% profile interval in ln theta (absolute), by stepping out then bisection (own route)."""
    from scipy.optimize import brentq
    lt0 = float(np.log(fmap["th"][0]))
    c0 = fmap["cost"]
    lowest = [c0]

    def h(lt):
        c = make(np.exp(lt)).fit(starts=(fmap["x"],), nk=14, nu=21, top=2)["cost"]
        lowest[0] = min(lowest[0], c)
        return c - c0 - CHI_HALF

    ends = []
    for sgn in (+1, -1):
        a, step = lt0, 0.05
        floor = LN_TH - 11.5
        hit = None
        while True:
            b = a + sgn * step
            if sgn < 0 and b <= floor:
                hit = None if h(floor) <= 0 else (a, floor)
                break
            if h(b) > 0:
                hit = (a, b)
                break
            a, step = b, step * 1.6
            if sgn > 0 and b > LN_TH + 4:
                break
        if hit is None:
            ends.append(floor if sgn < 0 else np.inf)
        else:
            ends.append(brentq(h, min(hit), max(hit), xtol=1e-4))
    return ends[1], ends[0], lowest[0] < c0 - 1e-7


def check7():
    print("\nCHECK 7. TASK 5 PROFILE WIDTHS IN ln theta: the investigator's 24 draws against independent draws")
    print("  (investigator's draws re-created: intake seed s, then default_rng(90000 + s): 26 noise normals,")
    print("  then the flux, k_e and load-bias errors). Own profile: step out, then bisection.")
    efs = []
    for s in range(24):
        rng = np.random.default_rng(90_000 + s)
        rng.standard_normal(26)
        efs.append(rng.standard_normal(3)[0])
    efs = np.array(efs)
    print(f"     the 24 flux-reading errors (standard units) on those draws: median {np.median(efs):+.3f}, "
          f"mean {np.mean(efs):+.3f}, {int(np.sum(efs > 0))} of 24 positive; P(median of 24 N(0,1) >= this) = "
          f"{norm.sf(np.median(efs) / (1.2533 / np.sqrt(24))):.4f} (normal approximation)")
    for lab, sc, ref in [x for x in SCEN if x[0] in ("flux 10%", "flux 20%", "flux 20% + k_e 1.2 + load 2%")]:
        s_lin = ref
        at_truth = np.log((1 + 1.96 * s_lin) / (1 - 1.96 * s_lin)) / (3.92 * s_lin)
        for src, nd in (("investigator's draws", 24), ("own independent draws", 40)):
            rat, maps = [], []
            for j in range(nd):
                if src.startswith("inv"):
                    m, inc = scen_arrays(marker(loads_repo(j)), sc)
                    rng = np.random.default_rng(90_000 + j)
                    e = rng.standard_normal(26)
                    ef, eke, eb = rng.standard_normal(3)
                else:
                    rng = np.random.default_rng(np.random.SeedSequence([20261001, 7, j]))
                    m, inc = scen_arrays(marker(loads_from(rng.standard_normal(T.size))), sc)
                    e = rng.standard_normal(26)
                    ef, eke, eb = rng.standard_normal(3)
                e = e if m.size == 26 else e[1:]
                y = ytrue(m, inc) + sig(sc["lf"]) * e
                fm = scen_map(m, inc, y, sc, eb, eke, ef).fit()
                lo, hi, lower_found = pl_ends(lambda fx: scen_map(m, inc, y, sc, eb, eke, ef, fix0=fx), fm)
                fsd = fim_single(m, inc, np.full(m.size, sig(sc["lf"])), sc["ske"], sc["lf"] * L, sc["Y"])[1]
                rat.append((hi - lo) / (3.92 * fsd))
                maps.append(np.log(fm["th"][0] / TH))
                assert not lower_found
            print(f"     {lab:30} {src:22} n = {nd:2d}: median PL width / FIM width {np.median(rat):.3f} "
                  f"(range {np.min(rat):.3f}-{np.max(rat):.3f}); median MAP {np.median(maps):+.3f}")
        print(f"     {'':30} linear-Gaussian ln-width ratio with theta_hat at the truth: {at_truth:.3f}")
    print(f"  {el()}")


def main():
    print("=" * 100)
    print("  RECHECK OF t11_fim_check.py (item 6). No import of core, dyn, t*.py. Simulated patients only.")
    print("=" * 100)
    check_footing()
    fims = check1()
    print(f"  {el()}")
    check2()
    check3()
    check4()
    check5()
    check6(fims)
    check7()
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
