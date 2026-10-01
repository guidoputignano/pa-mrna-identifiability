"""
Adversarial recheck of t6d_ekf.py (follow-up item 5: the one-pass EKF + RTS).

Shares no code with t6d_ekf.py, t6_estimators.py, dyn.py or core.py: nothing
from this directory is imported. The model, the benches, the draws, the
filters and the information matrix are rebuilt from scratch here.

Independent routes used
  * the transition is the EXACT closed-form flow of
        dm/dt = Lam - theta m/(K+m) - k_e m      (Lam = u + L0, constant per interval)
    obtained by separation of variables. With r1 > r2 the roots of
    k_e m^2 + (k_e K + theta - Lam) m - Lam K and d = r1 - r2,
        A ln(w1/w0) + B ln((w1 + d)/(w0 + d)) = -k_e tau,   w = m - r1,
    A = (K + r1)/d, B = 1 - A, solved for w1 by monotone Newton. t6 integrates
    the ODE with LSODA and forward sensitivities; here there is no ODE solver
    in the filter at all;
  * every Jacobian (flow map, steady state, trajectory sensitivities) is a
    Richardson-extrapolated central difference of that closed form, not a
    variational equation;
  * the EKF update is in Joseph form, the RTS smoother written out anew;
  * a bootstrap of the 80 records gives Monte Carlo intervals, a second seed
    set (default_rng(20261001)) checks that the draws are not special, and
    single records are removed to test outlier dependence;
  * two estimators the investigator did not run: a scaled unscented filter
    WITH a centre point, and a one-pass particle filter (Liu-West kernel
    shrinkage), to test the general verdict "a one-pass filter cannot be made
    a fair comparator"; and NLS started from the drawn prior centre only, to
    test the claim that t6's NLS figure relies on a start at the truth.

The model constants are restated from the repository's footing (normal
theta = 3000, K_m = 100, k_e = 5 /d, m = 10 at normal, convention A) and the
benches from t6_estimators.bench as read, not imported.
"""
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

T_START = time.time()

# ------------------------------------------------------------------ footing
TH_N, KM0, KE0, M_N = 3000.0, 100.0, 5.0, 10.0
LTOT = TH_N * M_N / (KM0 + M_N) + KE0 * M_N          # convention A, held as theta falls
L00 = 0.5 * LTOT                                     # unrecorded share (t6's default)
SD0 = np.array([0.3, 1.5, 0.8, 1.0, 0.5 * L00])      # t6.run_bench's prior sds
Q = 1e-10                                            # t6's process noise
REPS = 80
LOG = np.log


def bench(name):
    """t6's two benches, restated: steps of the KNOWN load u(t), then samples."""
    if name == "A":
        th = 0.05 * TH_N
        steps = [(-np.inf, 0.5 * LTOT), (0.0, 0.0)]
        t = np.arange(0, 2 + 1e-9, 1 / 24)
        cv = 0.05
    else:
        th = 0.20 * TH_N
        steps = [(-np.inf, 0.5 * LTOT), (0.0, 5.5 * LTOT), (1.0, 0.0)]
        t = np.arange(0, 2 + 1e-9, 1 / 96)
        cv = 0.01

    def u(tt):
        v = steps[0][1]
        for ts, us in steps:
            if tt >= ts:
                v = us
        return v
    tau = np.diff(t)
    umid = np.array([u(0.5 * (a + b)) for a, b in zip(t[:-1], t[1:])])
    # every break must fall on a sample time, so u is constant on each interval
    for ts, _ in steps[1:]:
        assert np.min(np.abs(t - ts)) < 1e-12
    p = np.array([LOG(th), LOG(KE0), LOG(KM0), L00])  # (ln th, ln ke, ln K, L0), t6's order
    return dict(name=name, th=th, t=t, tau=tau, umid=umid, upre=0.5 * LTOT, cv=cv, p=p)


# ------------------------------------------------------------------ closed form
def roots(th, ke, K, Lam):
    b = ke * K + th - Lam
    disc = np.sqrt(b * b + 4.0 * ke * Lam * K)
    # numerically stable pair of roots (np.where evaluates both branches)
    with np.errstate(divide="ignore", invalid="ignore"):
        r1 = np.where(b > 0, 2.0 * Lam * K / (b + disc), (-b + disc) / (2.0 * ke))
        r2 = np.where(b > 0, (-b - disc) / (2.0 * ke), -2.0 * Lam * K / (-b + disc))
    return r1, r2


def m_ss(th, ke, K, Lam):
    return roots(th, ke, K, Lam)[0]


def flow(m0, th, ke, K, Lam, tau):
    """Exact m(tau) from m(0) = m0 under constant total load Lam (vectorised)."""
    m0, th, ke, K, Lam = np.broadcast_arrays(*(np.asarray(v, float) for v in (m0, th, ke, K, Lam)))
    r1, r2 = roots(th, ke, K, Lam)
    d = r1 - r2
    A = (K + r1) / d
    B = 1.0 - A
    w0 = m0 - r1
    kt = ke * tau
    # start right of the root when w0 > 0 (f convex), left of it when w0 < 0 (concave)
    s_lo = (-kt - B * np.log(d / np.maximum(w0 + d, 1e-300))) / A
    s = np.where(w0 >= 0, 0.0, s_lo)

    def newton(s):
        w1 = w0 * np.exp(s)
        f = A * s + B * np.log((w1 + d) / (w0 + d)) + kt
        fp = A + B * w1 / (w1 + d)
        return s - f / fp, f
    for _ in range(200):
        s_new, _ = newton(s)
        conv = np.all(np.abs(s_new - s) < 1e-13 * np.maximum(1.0, np.abs(s)))
        s = s_new
        if conv:
            break
    s, f = newton(s)               # one more step at rounding level
    scale = kt + np.abs(A * s) + 1.0
    if not np.all(np.abs(f) < 1e-12 * scale):
        raise RuntimeError("flow: Newton did not converge")
    m1 = r1 + w0 * np.exp(s)
    return np.maximum(m1, 1e-9)


def lnflow(x, Lam_u, tau):
    """ln m1 from x = (ln m0, ln th, ln ke, ln K, L0) stacked as columns (n, 5)."""
    x = np.atleast_2d(x)
    return np.log(flow(np.exp(x[:, 0]), np.exp(x[:, 1]), np.exp(x[:, 2]), np.exp(x[:, 3]),
                       Lam_u + x[:, 4], tau))


def lnss(p, upre):
    p = np.atleast_2d(p)
    return np.log(np.maximum(m_ss(np.exp(p[:, 0]), np.exp(p[:, 1]), np.exp(p[:, 2]),
                                  upre + p[:, 3]), 1e-9))


def steps_for(x):
    return np.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-3 * max(abs(x[-1]), 1.0)])[-len(x):]


def jac_rich(fun, x):
    """d fun / d x by Richardson-extrapolated central differences (one vectorised call)."""
    n = len(x)
    h = steps_for(x)
    pts = [x]
    for j in range(n):
        for c in (1.0, -1.0, 0.5, -0.5):
            e = np.zeros(n); e[j] = c * h[j]
            pts.append(x + e)
    v = fun(np.array(pts))
    f0 = v[0]
    J = np.zeros((v.shape[1] if v.ndim > 1 else 1, n)) if v.ndim > 1 else np.zeros(n)
    for j in range(n):
        a, b_, c, d_ = v[1 + 4 * j: 5 + 4 * j]
        D1 = (a - b_) / (2 * h[j])
        D2 = (c - d_) / h[j]
        if v.ndim > 1:
            J[:, j] = (4 * D2 - D1) / 3
        else:
            J[j] = (4 * D2 - D1) / 3
    return f0, J


# ------------------------------------------------------------------ trajectories
def traj(P, b, lnm0=None):
    """ln m at the sample times for parameter rows P (n, 4); m0 = steady state unless given."""
    P = np.atleast_2d(P)
    lm = lnss(P, b["upre"]) if lnm0 is None else np.asarray(lnm0, float)
    out = [lm]
    th, ke, K = np.exp(P[:, 0]), np.exp(P[:, 1]), np.exp(P[:, 2])
    for k in range(len(b["tau"])):
        lm = np.log(flow(np.exp(lm), th, ke, K, b["umid"][k] + P[:, 3], b["tau"][k]))
        out.append(lm)
    return np.array(out).T          # (n_rows, n_samples)


def fim(b, free_m0, prior=None):
    """Information in (ln th, ln ke, ln K, L0) [m0 tied] or (ln m0, ...) [m0 free]."""
    x = b["p"] if not free_m0 else np.r_[lnss(b["p"], b["upre"])[0], b["p"]]
    if free_m0:
        fun = lambda X: traj(X[:, 1:], b, lnm0=X[:, 0])
    else:
        fun = lambda X: traj(X, b)
    _, S = jac_rich(fun, x)
    F = S.T @ S / b["cv"] ** 2
    if prior is not None:
        F = F + np.diag(prior ** -2)
    return F, S


def sd_th(F, idx):
    return float(np.sqrt(np.linalg.inv(F)[idx, idx]))


# ------------------------------------------------------------------ draws
def draws(b, m_true, seed, reps):
    """t6.run_bench's order: per replicate the record y, then five normals for x0."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(reps):
        y = m_true * np.exp(rng.normal(0, b["cv"], m_true.size))
        x0 = np.r_[np.log(m_true[0]), b["p"]] + SD0 * rng.normal(size=5)
        x0[0] = lnss(x0[1:], b["upre"])[0]
        out.append((y, x0))
    return out


# ------------------------------------------------------------------ EKF + RTS
def ekf_rts(y, b, x0, P0, q=Q, rts=True):
    ly = np.log(y)
    n = len(y)
    I5 = np.eye(5)
    xf, Pf, xpred, Ppred, Fl = [], [], [], [], []
    x, P = x0.astype(float).copy(), P0.astype(float).copy()
    for k in range(n):
        if k > 0:
            fun = lambda X, kk=k: lnflow(X, b["umid"][kk - 1], b["tau"][kk - 1])
            f0, J = jac_rich(fun, x)
            F = I5.copy(); F[0, :] = J
            xn = x.copy(); xn[0] = f0
            P = F @ P @ F.T + q * I5
            x = xn
            xpred.append(x.copy()); Ppred.append(P.copy()); Fl.append(F)
        # Joseph-form update on ln y_k = x[0] + e
        S = P[0, 0] + b["cv"] ** 2
        Kg = P[:, 0] / S
        x = x + Kg * (ly[k] - x[0])
        IKH = I5 - np.outer(Kg, I5[0])
        P = IKH @ P @ IKH.T + np.outer(Kg, Kg) * b["cv"] ** 2
        xf.append(x.copy()); Pf.append(P.copy())
    xs, Ps = xf[-1].copy(), Pf[-1].copy()
    for k in (range(n - 2, -1, -1) if rts else ()):
        G = np.linalg.solve(Ppred[k].T, (Pf[k] @ Fl[k].T).T).T
        xs = xf[k] + G @ (xs - xpred[k])
        Ps = Pf[k] + G @ (Ps - Ppred[k]) @ G.T
    return xs, Ps, xf[-1], Pf[-1]


def tied_P0(x0, b):
    """ln m0 = h(p): P0 = [J; I] S [J; I]^T with J = d ln m_ss / dp at the centre."""
    _, J = jac_rich(lambda X: lnss(X, b["upre"]), x0[1:])
    S = np.diag(SD0[1:] ** 2)
    T = np.vstack([J, np.eye(4)])
    return T @ S @ T.T


# ------------------------------------------------------------------ UKF (scaled, centre point)
def ukf(y, b, x0, P0, alpha, beta=2.0, kappa=0.0, q=Q, tied=True):
    """Scaled unscented filter (Julier/Wan-van der Merwe) on (ln m, p). With tied=True the
    initial points are drawn on p alone and ln m0 = h(p) per point."""
    ly = np.log(y)

    def pts(mu, P):
        n = len(mu)
        lam = alpha ** 2 * (n + kappa) - n
        w, V = np.linalg.eigh(0.5 * (P + P.T))
        L = V * np.sqrt(np.clip(w, 0, None) * (n + lam))
        X = np.vstack([mu] + [mu + L[:, j] for j in range(n)] + [mu - L[:, j] for j in range(n)])
        wm = np.full(2 * n + 1, 0.5 / (n + lam)); wc = wm.copy()
        wm[0] = lam / (n + lam); wc[0] = wm[0] + 1 - alpha ** 2 + beta
        return X, wm, wc

    def mom(X, wm, wc):
        mu = wm @ X
        dX = X - mu
        return mu, (dX * wc[:, None]).T @ dX

    if tied:
        Xp, wm, wc = pts(x0[1:], P0[1:, 1:])
        X = np.c_[lnss(Xp, b["upre"]), Xp]
        x, P = mom(X, wm, wc)
    else:
        x, P = x0.copy(), P0.copy()
    for k in range(len(y)):
        if k > 0:
            X, wm, wc = pts(x, P)
            X[:, 0] = lnflow(X, b["umid"][k - 1], b["tau"][k - 1])
            x, P = mom(X, wm, wc)
            P = P + q * np.eye(5)
        S = P[0, 0] + b["cv"] ** 2
        Kg = P[:, 0] / S
        x = x + Kg * (ly[k] - x[0])
        P = P - np.outer(Kg, P[0, :])
        P = 0.5 * (P + P.T)
    return x, P


# ------------------------------------------------------------------ particle filter
def pf(y, b, x0, N, rng, a_shrink=0.98):
    """One-pass bootstrap particle filter for (ln m, p), p static, with Liu-West kernel
    shrinkage (a = a_shrink) at each resampling. Each y_k is used once, in time order."""
    ly = np.log(y)
    cv = b["cv"]
    Pp = x0[1:] + SD0[1:] * rng.normal(size=(N, 4))
    lm = lnss(Pp, b["upre"])
    logw = np.zeros(N)
    h2 = 1.0 - a_shrink ** 2
    nres = 0
    for k in range(len(y)):
        if k > 0:
            lm = np.log(flow(np.exp(lm), np.exp(Pp[:, 0]), np.exp(Pp[:, 1]), np.exp(Pp[:, 2]),
                             b["umid"][k - 1] + Pp[:, 3], b["tau"][k - 1]))
        logw = logw - 0.5 * ((ly[k] - lm) / cv) ** 2
        w = np.exp(logw - logw.max()); w /= w.sum()
        ess = 1.0 / np.sum(w ** 2)
        if ess < 0.5 * N:
            nres += 1
            mu = w @ Pp
            C = (Pp - mu).T @ ((Pp - mu) * w[:, None])
            # systematic resampling
            u0 = (rng.random() + np.arange(N)) / N
            idx = np.minimum(np.searchsorted(np.cumsum(w), u0), N - 1)
            Pp, lm = Pp[idx], lm[idx]
            # Liu-West: shrink to the mean and jitter, keeping the first two moments
            Lc = np.linalg.cholesky(C + 1e-12 * np.eye(4))
            Pp = a_shrink * Pp + (1 - a_shrink) * mu + np.sqrt(h2) * rng.normal(size=(N, 4)) @ Lc.T
            logw = np.zeros(N)
    w = np.exp(logw - logw.max()); w /= w.sum()
    mu = w @ Pp
    return mu, float(np.sqrt(w @ (Pp[:, 0] - mu[0]) ** 2)), nres


# ------------------------------------------------------------------ NLS from the centre
def nls_from(y, b, p0):
    ly = np.log(y)

    def res(p):
        return traj(p, b)[0] - ly

    def jac(p):
        _, J = jac_rich(lambda X: traj(X, b), p)
        return J
    lo = np.array([np.log(1e-2), np.log(1e-3), np.log(1.0), -1e4])
    hi = np.array([np.log(1e5), np.log(1e3), np.log(1e5), 1e5])
    s = least_squares(res, np.clip(p0, lo + 1e-6, hi - 1e-6), jac=jac, bounds=(lo, hi),
                      method="trf", x_scale="jac", xtol=1e-10, ftol=1e-10, max_nfev=200)
    return s.x, s.cost


# ------------------------------------------------------------------ statistics
def stats(e):
    e = np.asarray(e, float)
    return dict(bias=float(np.median(e)), rmse=float(np.sqrt(np.mean(e ** 2))),
                robust=float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16))),
                gt1=float(np.mean(np.abs(e) > 1)))


def boot(e, B=4000, seed=99):
    rng = np.random.default_rng(seed)
    e = np.asarray(e, float)
    idx = rng.integers(0, e.size, size=(B, e.size))
    E = e[idx]
    rm = np.sqrt(np.mean(E ** 2, axis=1))
    rb = 0.5 * (np.percentile(E, 84, axis=1) - np.percentile(E, 16, axis=1))
    return np.percentile(rm, [2.5, 97.5]), np.percentile(rb, [2.5, 97.5])


def paired(e1, e2, B=4000, seed=98):
    rng = np.random.default_rng(seed)
    e1, e2 = np.asarray(e1, float), np.asarray(e2, float)
    idx = rng.integers(0, e1.size, size=(B, e1.size))
    A, C = e1[idx], e2[idx]
    rm = np.sqrt(np.mean(A ** 2, 1)) - np.sqrt(np.mean(C ** 2, 1))
    rb = (0.5 * (np.percentile(A, 84, 1) - np.percentile(A, 16, 1))
          - 0.5 * (np.percentile(C, 84, 1) - np.percentile(C, 16, 1)))
    return np.percentile(rm, [2.5, 97.5]), np.percentile(rb, [2.5, 97.5])


def line(lab, e, extra=""):
    s = stats(e)
    return (f"   {lab:46} {s['bias']:+7.3f} {s['rmse']:7.3f} {s['robust']:7.3f} {s['gt1']:6.0%}"
            + (f"  {extra}" if extra else ""))


HDR = f"   {'estimator':46} {'bias':>7} {'RMSE':>7} {'robust':>7} {'|e|>1':>6}"


# ------------------------------------------------------------------ main
def main():
    print("=" * 100)
    print("  RECHECK OF t6d_ekf.py (item 5), on a closed-form transition, no shared code")
    print("=" * 100)
    print(f"  footing: L = {LTOT:.4f}, L0 = {L00:.4f}, prior sd0 = {np.round(SD0, 3)}")
    B = {nm: bench(nm) for nm in ("A", "B")}

    # ------------------------------------------------ 0. closed form against an ODE solver
    print("\n  0. CLOSED-FORM FLOW against scipy Radau (rtol 1e-12) on the true trajectories")
    M = {}
    for nm, b in B.items():
        lm = traj(b["p"], b)[0]
        M[nm] = np.exp(lm)
        th, ke, K, L0 = np.exp(b["p"][0]), np.exp(b["p"][1]), np.exp(b["p"][2]), b["p"][3]
        worst, m = 0.0, M[nm][0]
        for k in range(len(b["tau"])):
            Lam = b["umid"][k] + L0
            s = solve_ivp(lambda t_, z: [Lam - th * z[0] / (K + z[0]) - ke * z[0]],
                          (0, b["tau"][k]), [m], method="Radau", rtol=1e-12, atol=1e-12)
            worst = max(worst, abs(s.y[0, -1] / M[nm][k + 1] - 1))
            m = M[nm][k + 1]
        print(f"   bench {nm}: m0 = {M[nm][0]:.4f}, m from {M[nm].min():.2f} to {M[nm].max():.2f};"
              f" max rel. one-step difference {worst:.1e}")
        assert worst < 1e-8

    # ------------------------------------------------ 1. information bounds
    print("\n  1. INFORMATION BOUNDS for ln theta, from Richardson differences of the closed form")
    BOUND = {}
    for nm, b in B.items():
        F0, _ = fim(b, free_m0=False)
        Ft, _ = fim(b, free_m0=False, prior=SD0[1:])
        Ff, _ = fim(b, free_m0=True, prior=SD0)
        BOUND[nm] = dict(none=sd_th(F0, 0), tied=sd_th(Ft, 0), free=sd_th(Ff, 1))
        print(f"   bench {nm}: no prior {BOUND[nm]['none']:.3f};  FIM+prior, m0 tied "
              f"{BOUND[nm]['tied']:.5f};  FIM+prior, m0 free {BOUND[nm]['free']:.5f};  "
              f"cond(FIM, no prior) {np.linalg.cond(F0):.1e}")

    # ------------------------------------------------ 2. draws
    DR = {nm: draws(B[nm], M[nm], 1, REPS) for nm in B}
    DR2 = {"B": draws(B["B"], M["B"], 20261001, REPS)}

    # ------------------------------------------------ 3. t6's one-pass EKF + RTS, rebuilt
    print(f"\n  2. t6's CONFIGURATION (diag P0, q = 1e-10, RTS time-0 estimate), rebuilt here, {REPS}"
          f" records; t6's draws (default_rng(1), y then x0)")
    print(HDR)
    E = {}
    for nm, b in B.items():
        t0 = time.time()
        e, sdr, ekm = [], [], []
        for y, x0 in DR[nm]:
            xs, Ps, _, _ = ekf_rts(y, b, x0, np.diag(SD0 ** 2))
            e.append(xs[1] - b["p"][0]); sdr.append(np.sqrt(Ps[1, 1])); ekm.append(xs[1:] - b["p"])
        E[(nm, "ekf")] = np.array(e)
        E[(nm, "ekf_all")] = np.array(ekm)
        E[(nm, "prior")] = np.array([x0[1] - b["p"][0] for _, x0 in DR[nm]])
        ci_rm, ci_rb = boot(E[(nm, "ekf")])
        print(line(f"bench {nm}: EKF+RTS, t6 configuration", e,
                   f"median rep. sd {np.median(sdr):.3f}; boot 95% RMSE [{ci_rm[0]:.2f}, {ci_rm[1]:.2f}],"
                   f" robust [{ci_rb[0]:.2f}, {ci_rb[1]:.2f}]; {time.time() - t0:.0f} s"))
        ci_rm, ci_rb = boot(E[(nm, "prior")])
        print(line(f"bench {nm}: prior centre alone", E[(nm, "prior")],
                   f"boot 95% RMSE [{ci_rm[0]:.2f}, {ci_rm[1]:.2f}], robust [{ci_rb[0]:.2f}, {ci_rb[1]:.2f}]"))
    eB = E[("B", "ekf")]
    o = np.argsort(-np.abs(eB))
    print(f"   bench B, the largest |error| records: " + ", ".join(f"#{i}: {eB[i]:+.2f}" for i in o[:4]))
    share = eB[o[0]] ** 2 / np.sum(eB ** 2)
    print(f"   share of the bench-B EKF sum of squared errors from ONE record: {share:.1%}")
    print(line("bench B: EKF+RTS without that one record", np.delete(eB, o[0])))
    print(line("bench B: prior centre alone without that record", np.delete(E[("B", "prior")], o[0])))
    xe = E[("B", "ekf_all")]
    print(f"   bench B: corr(error ln th, error ln Km) {np.corrcoef(xe[:, 0], xe[:, 2])[0, 1]:+.3f};"
          f" records with |e| > 2: {int(np.sum(np.abs(eB) > 2))}")

    # second seed set
    b = B["B"]
    e2 = []
    for y, x0 in DR2["B"]:
        xs, _, _, _ = ekf_rts(y, b, x0, np.diag(SD0 ** 2))
        e2.append(xs[1] - b["p"][0])
    E[("B2", "ekf")] = np.array(e2)
    E[("B2", "prior")] = np.array([x0[1] - b["p"][0] for _, x0 in DR2["B"]])
    print(line("bench B, NEW seeds (20261001): EKF+RTS", e2))
    print(line("bench B, NEW seeds: prior centre alone", E[("B2", "prior")]))
    d_rm, d_rb = paired(E[("B", "ekf")], E[("B", "prior")])
    print(f"   paired bootstrap, t6's draws, EKF minus prior centre: RMSE difference 95% "
          f"[{d_rm[0]:+.2f}, {d_rm[1]:+.2f}], robust difference 95% [{d_rb[0]:+.2f}, {d_rb[1]:+.2f}]")
    # the RTS smoother on a static parameter: time-0 smoothed vs filtered end
    dmax = 0.0
    for y, x0 in DR["B"][:20]:
        xs, _, xfn, _ = ekf_rts(y, b, x0, np.diag(SD0 ** 2))
        dmax = max(dmax, float(np.max(np.abs(xs[1:] - xfn[1:]))))
    print(f"   RTS: max |smoothed p(t0) - filtered p(t_n)| over 20 records {dmax:.1e}")

    # ------------------------------------------------ 4. at the truth, prior scaled
    print("\n  3. t6's EKF STARTED AT THE TRUTH on bench B (t6's draws), parameter prior sds scaled;")
    print("     reference: sd of the LINEARISED MAP with the prior centre fixed at the truth,")
    print("     (I_d + I_p)^-1 I_d (I_d + I_p)^-1, m0 free (what a correct estimator would give)")
    print(HDR + f" {'lin. MAP sd':>12}")
    F_d, _ = fim(b, free_m0=True)
    x_true = np.r_[np.log(M["B"][0]), b["p"]]
    for sc in (0.1, 0.3, 0.5, 1.0):
        s = SD0.copy(); s[1:] *= sc
        Ip = np.diag(s ** -2)
        Wi = np.linalg.inv(F_d + Ip)
        lin = float(np.sqrt((Wi @ F_d @ Wi)[1, 1]))
        e = []
        for y, _ in DR["B"]:
            xs, _, _, _ = ekf_rts(y, b, x_true, np.diag(s ** 2))
            e.append(xs[1] - b["p"][0])
        E[("Btruth", sc)] = np.array(e)
        print(line(f"EKF at the truth, prior sd x {sc}", e) + f" {lin:12.4f}")

    # ------------------------------------------------ 5. tied P0 (variant a)
    print("\n  4. VARIANT (a), TIED P0, filtered end point (as t6d reports it), t6's draws, bench B")
    print(HDR)
    e = []
    for y, x0 in DR["B"]:
        _, _, xfn, _ = ekf_rts(y, b, x0, tied_P0(x0, b))
        e.append(xfn[1] - b["p"][0])
    E[("B", "tied")] = np.array(e)
    ci_rm, ci_rb = boot(e)
    print(line("tied P0 (rebuilt)", e, f"boot 95% RMSE [{ci_rm[0]:.2f}, {ci_rm[1]:.2f}], "
                                       f"robust [{ci_rb[0]:.2f}, {ci_rb[1]:.2f}]"))
    d_rm, d_rb = paired(E[("B", "tied")], E[("B", "prior")])
    print(f"   paired bootstrap, tied minus prior centre: RMSE difference 95% "
          f"[{d_rm[0]:+.2f}, {d_rm[1]:+.2f}], robust difference 95% [{d_rb[0]:+.2f}, {d_rb[1]:+.2f}]")
    print(f"   seconds so far {time.time() - T_START:.0f}")

    # ------------------------------------------------ 6. unscented filter with a centre point
    print("\n  5. SCALED UNSCENTED FILTER WITH A CENTRE POINT (not run by t6d), tied prior, bench B")
    print(HDR)
    for alpha, kappa in ((1e-2, 0.0), (0.5, 0.0), (1.0, 1.0)):
        e = []
        for y, x0 in DR["B"]:
            x, P = ukf(y, b, x0, np.diag(SD0 ** 2), alpha=alpha, kappa=kappa)
            e.append(x[1] - b["p"][0])
        E[("B", "ukf", alpha)] = np.array(e)
        print(line(f"UKF alpha {alpha}, beta 2, kappa {kappa}", e))

    # ------------------------------------------------ 7. particle filter
    print("\n  6. ONE-PASS PARTICLE FILTER (Liu-West), tied prior, posterior mean of ln th, bench B")
    print("     shrinkage a tuned on SEPARATE draws (default_rng(7), 12 records, N = 20000) by RMSE")
    print(HDR)
    DT = draws(b, M["B"], 7, 12)
    tune = {}
    for a in (0.98, 0.995, 0.998):
        rng = np.random.default_rng(5)
        e = [pf(y, b, x0, 20000, rng, a_shrink=a)[0][0] - b["p"][0] for y, x0 in DT]
        tune[a] = stats(e)["rmse"]
        print(line(f"tuning: a = {a}", e))
    a_best = min(tune, key=tune.get)
    print(f"   chosen a = {a_best}")
    for N, reps in ((20000, REPS), (50000, 40)):
        rng = np.random.default_rng(5)
        t0 = time.time()
        e, nr = [], []
        for y, x0 in DR["B"][:reps]:
            mu, sdp, n_r = pf(y, b, x0, N, rng, a_shrink=a_best)
            e.append(mu[0] - b["p"][0]); nr.append(n_r)
        E[("B", "pf", N)] = np.array(e)
        ci_rm, ci_rb = boot(e)
        print(line(f"PF, N = {N}, a = {a_best}, first {reps} of t6's draws", e,
                   f"boot 95% RMSE [{ci_rm[0]:.2f}, {ci_rm[1]:.2f}], robust [{ci_rb[0]:.2f}, "
                   f"{ci_rb[1]:.2f}]; median resamplings {np.median(nr):.0f}; {time.time() - t0:.0f} s"))
        if reps < REPS:
            print(line(f"   t6's EKF on the same {reps} records", E[("B", "ekf")][:reps]))
            print(line(f"   prior centre on the same {reps} records", E[("B", "prior")][:reps]))
    print(f"   seconds so far {time.time() - T_START:.0f}")

    # ------------------------------------------------ 8. NLS started at the drawn centre
    print("\n  7. NLS (no priors, as t6) STARTED ONLY AT THE DRAWN PRIOR CENTRE, bench B, t6's draws")
    print(HDR)
    nre = 40
    t0 = time.time()
    e, et = [], []
    for y, x0 in DR["B"][:nre]:
        pc, cc = nls_from(y, b, x0[1:])
        pt, ct = nls_from(y, b, b["p"])
        e.append(pc[0] - b["p"][0]); et.append(pt[0] - b["p"][0])
    e, et = np.array(e), np.array(et)
    print(line(f"NLS from the prior centre, first {nre} records", e, f"{time.time() - t0:.0f} s (both rows)"))
    print(line(f"NLS from the truth, same {nre} records", et))
    print(f"   records where the centre start ends > 0.05 from the truth start in ln th: "
          f"{int(np.sum(np.abs(e - et) > 0.05))} of {nre}")

    # ------------------------------------------------ verdict line, on this file's numbers
    nls_rm = stats(e)["rmse"]
    print("\n  8. FAIRNESS against this file's own bound and NLS-from-centre, bench B")
    cands = {"EKF+RTS, t6 configuration": E[("B", "ekf")], "tied P0": E[("B", "tied")],
             f"PF N = 20000": E[("B", "pf", 20000)], f"PF N = 50000 (40 records)": E[("B", "pf", 50000)]}
    for a_ in (1e-2, 0.5, 1.0):
        cands[f"UKF alpha {a_}"] = E[("B", "ukf", a_)]
    for lab, v in cands.items():
        st = stats(v)
        ok = st["robust"] <= 1.5 * BOUND["B"]["tied"] and st["rmse"] <= 2 * nls_rm
        print(f"   {lab:32} robust {st['robust']:.3f} ({st['robust'] / BOUND['B']['tied']:.1f} x bound), "
              f"RMSE {st['rmse']:.3f} ({st['rmse'] / nls_rm:.1f} x NLS-from-centre) -> "
              f"{'fair' if ok else 'not fair'}")

    # ------------------------------------------------ 9. smaller checks of t6d's component claims
    print("\n  9. COMPONENT CLAIMS, on this file's code")
    for nm, bb in B.items():
        dl = np.array([x0[0] - np.log(M[nm][0]) for _, x0 in DR[nm]])
        sdt = np.array([np.sqrt(tied_P0(x0, bb)[0, 0]) for _, x0 in DR[nm]])
        print(f"   bench {nm}: ln m_ss(drawn centre) - ln m0(true): sd {dl.std():.3f}, max |.| "
              f"{np.abs(dl).max():.3f}; tied prior sd(ln m0) median {np.median(sdt):.3f}")
        # information identity on noiseless data at the truth
        xt = np.r_[np.log(M[nm][0]), bb["p"]]
        xs, Ps, _, _ = ekf_rts(M[nm], bb, xt, np.diag(SD0 ** 2))
        xs2, Ps2, _, _ = ekf_rts(M[nm], bb, np.r_[lnss(bb["p"], bb["upre"])[0], bb["p"]],
                                 tied_P0(np.r_[0.0, bb["p"]], bb))
        print(f"   bench {nm}: noiseless at the truth, EKF+RTS sd(ln th) {np.sqrt(Ps[1, 1]):.5f} (diag P0),"
              f" {np.sqrt(Ps2[1, 1]):.5f} (tied); bounds {BOUND[nm]['free']:.5f}, {BOUND[nm]['tied']:.5f}")
    # IEKF measurement update on h(x) = x[0]: iterate x_i -> xp + K (y - h(x_i) - H (xp - x_i))
    rng_ = np.random.default_rng(3)
    worst = 0.0
    for _ in range(200):
        xp = rng_.normal(size=5); A_ = rng_.normal(size=(5, 5)); Pp = A_ @ A_.T + np.eye(5)
        yk = rng_.normal(); H = np.eye(5)[0]
        K_ = Pp @ H / (H @ Pp @ H + 0.01 ** 2)
        x1 = xp + K_ * (yk - xp[0])
        xi = xp.copy()
        for _ in range(5):
            xi = xp + K_ * (yk - xi[0] - H @ (xp - xi))
        worst = max(worst, float(np.max(np.abs(xi - x1))))
    print(f"   IEKF update with h(x) = x[0], 5 iterations vs one update, 200 random cases: max |dx| {worst:.1e}")
    xe = E[("B", "ekf_all")]
    big = np.abs(xe[:, 0]) > 2
    kb = xe[:, 2]
    print(f"   bench B, t6 EKF: of the {big.sum()} records with |e ln th| > 2, ln Km error > +2 in "
          f"{int((big & (kb > 2)).sum())}, < -2 in {int((big & (kb < -2)).sum())}; of all {REPS}: "
          f"{int((kb > 2).sum())} and {int((kb < -2).sum())}")
    print(f"   bench B samples {B['B']['t'].size}, transitions per pass {B['B']['t'].size - 1}")
    v = np.clip(E[("B", "ekf")], -8, 8)
    q1, q3 = np.percentile(v, [25, 75])
    lo_w = v[v >= q1 - 1.5 * (q3 - q1)].min(); hi_w = v[v <= q3 + 1.5 * (q3 - q1)].max()
    print(f"   bench B, t6 EKF: quartiles {q1:+.2f} / {q3:+.2f}, 1.5-IQR whiskers {lo_w:+.2f} to {hi_w:+.2f}")
    e0 = []
    for y, x0 in DR["B"]:
        _, _, xfn, _ = ekf_rts(y, B["B"], x0, np.diag(SD0 ** 2), q=0.0, rts=False)
        e0.append(xfn[1] - B["B"]["p"][0])
    print(line("bench B, t6 configuration with q = 0, filtered end", e0))
    print(f"   random-walk sd from q = 1e-10 over the record: B {np.sqrt(1e-10 * 192):.1e}, "
          f"A {np.sqrt(1e-10 * 48):.1e}")

    # ------------------------------------------------ summary assertions (can fail)
    assert BOUND["B"]["tied"] < 0.1 and BOUND["A"]["none"] > 100
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
