"""
Follow-up item 5. Is t6's one-pass EKF + RTS a bug, a misconfiguration, or an
estimator that cannot be a fair comparator here?

t6.est_ekf: augmented state x = (ln m, ln theta, ln k_e, ln K_m, L0), prior
centre drawn about the truth with sd0 = (0.3, 1.5, 0.8, 1.0, 0.5 L0), ln m0
then reset to the steady state of the drawn parameters, P0 = diag(sd0^2),
process noise q = 1e-10 on every state, the dynamics linearised about the
running filtered estimate, a linear measurement ln y = x[0] + e, and an RTS
smoother whose time-0 parameter estimate is returned.

This script
  0. reproduces t6's EKF on t6's own draws (default_rng(1); y, then the five
     prior-centre normals, for each of the first 80 replicates), and parses
     out/t6_estimators.txt so the comparison rows are t6's printed numbers;
  1. checks each component the request names: the initial covariance, q, the
     log-state Jacobian against finite differences, the measurement update and
     RTS gain (by the exact-information identity on noiseless data at the
     truth), and the integrator failure path;
  2. runs fixed one-pass variants on both benches, tuning any free setting on
     separate draws (default_rng(2)) and scoring on t6's draws:
       (a)  consistent initial covariance: ln m0 = h(p), P0 = [J; I] S [J; I]^T
       (a') (a) plus an iterated (Gauss-Newton) update on y_0, where the
            measurement IS nonlinear in the parameters
       (b)  (a') plus annealed parameter process noise
       (c)  (a') plus an iterated update that re-linearises the latest
            transition about its lag-1 smoothed estimate (one data pass)
       (c0) t6's filter with an iterated measurement update (IEKF in the
            Bell-Cathey sense), to show whether it changes anything
       (d)  a cubature (unscented-type, third-degree spherical-radial) filter
            with the tied prior built by pushing the parameter points through
            the steady state;
  3. runs references that are NOT fair comparators but locate the failure: the
     prior centre alone, t6's filter started at the truth, and t6's filter at
     the truth with its prior narrowed.

"One pass" here means: every y_k is processed once, in time order, and no
transition before the latest one is ever re-linearised. The iterated smoother
(t6.est_ieks) re-linearises every transition on every pass; no variant here
does. The number of transition evaluations per replicate is printed for each.
"""
import re
import time

import numpy as np
from scipy.integrate import solve_ivp

import core as C
import dyn as D
import t6_estimators as T6
from plotstyle import plt, C1, C2, C3, INK, MUTED, finish, save, plain

REPS = 80          # t6's reps_nl: the same 80 (y, prior centre) pairs as t6
TUNE_REPS = 40     # tuning draws, separate seed
TUNE_SEED = 2
EVAL_SEED = 1      # t6.run_bench's seed
STEPS = [0]        # transition evaluations, counted
E0 = np.eye(5)[0]


# ------------------------------------------------------------------ helpers
def sd_prior(b):
    L00 = D.unpack(b["p"])[3]
    return np.array([0.3, 1.5, 0.8, 1.0, 0.5 * L00])


def draws(b, reps, seed):
    """t6.run_bench's draw order: per replicate y, then 5 normals for x0."""
    rng = np.random.default_rng(seed)
    sd0 = sd_prior(b)
    out = []
    for _ in range(reps):
        y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
        x0 = np.r_[np.log(b["m"][0]), b["p"]] + sd0 * rng.normal(size=5)
        x0[0] = np.log(D.steady(x0[1:], b["ex"]))
        out.append((y, x0))
    return out


def stats(e, sd=None):
    e = np.asarray(e, float)
    fin = np.isfinite(e)
    ef = e[fin]
    d = dict(bias=float(np.median(ef)), rmse=float(np.sqrt(np.mean(ef ** 2))),
             robust=float(0.5 * (np.percentile(ef, 84) - np.percentile(ef, 16))),
             gt1=float(np.mean(np.abs(ef) > 1)), div=int((~fin).sum()))
    d["repsd"] = float(np.median(np.asarray(sd)[fin])) if sd is not None else np.nan
    return d


def step(xb, ex, t0, t1):
    """t6's transition with its Jacobian in log state; flags the failure path."""
    STEPS[0] += 1
    mb = np.exp(xb[0])
    with np.errstate(all="ignore"):
        m1, dm0, dp = T6._step(mb, xb[1:], ex, t0, t1)
    failed = (dm0 == 0.0) and not np.any(dp)
    m1 = max(m1, 1e-9)
    F = np.eye(5)
    F[0, 0] = dm0 * mb / m1
    F[0, 1:] = dp / m1
    fx = xb.copy()
    fx[0] = np.log(m1)
    return F, fx, failed


def ss_lnm(p, ex):
    """ln m_ss(p) and d ln m_ss / dp (implicit differentiation of the balance)."""
    with np.errstate(all="ignore"):
        th, ke, K, L0, _ = D.unpack(p)
        Lt = max(ex.u(-1e-9) + L0, 1e-6)
        m0 = float(C.m_ss(th, K, ke, Lt))
        g0 = m0 / (K + m0)
        Cp = th * K / (K + m0) ** 2 + ke
        dC = np.array([th * g0, ke * m0, -th * K * m0 / (K + m0) ** 2, -1.0])
    return np.log(m0), -dC / Cp / m0


def tied_cov(J, S):
    P = np.zeros((5, 5))
    P[1:, 1:] = S
    P[0, 1:] = J @ S
    P[1:, 0] = S @ J
    P[0, 0] = J @ S @ J
    return P


def kf_update(xp, Pp, yk, cv, iters=1):
    """Measurement h(x) = x[0]; iters > 1 is the Bell-Cathey iterated update."""
    xi = xp.copy()
    for _ in range(iters):
        H = E0                      # Jacobian of h at the iterate
        S = H @ Pp @ H + cv ** 2
        K = Pp @ H / S
        xi = xp + K * (np.log(yk) - xi[0] - H @ (xp - xi))
    return xi, Pp - np.outer(K, H @ Pp)


def init_gn(y0, x0, S0, ex, cv, iters=30, tol=1e-10):
    """Iterated update on ln y_0 = h(p), h = ln m_ss: nonlinear in p."""
    p0 = x0[1:]
    pb = p0.copy()
    for _ in range(iters):
        lm, J = ss_lnm(pb, ex)
        S = J @ S0 @ J + cv ** 2
        K = S0 @ J / S
        pn = p0 + K * (np.log(y0) - lm - J @ (p0 - pb))
        if not np.all(np.isfinite(pn)):
            break
        done = np.max(np.abs(pn - pb)) < tol
        pb = pn
        if done:
            break
    lm, J = ss_lnm(pb, ex)
    Pn = S0 - np.outer(K, J @ S0)
    return np.r_[lm, pb], tied_cov(J, Pn)


# ------------------------------------------------------------------ filters
def ekf(y, t, ex, x0, P0, cv, q=1e-10, qpar=None, meas_it=1, dyn_it=1,
        init=None, rts=True):
    """
    One-pass EKF with optional pieces. With every option at its default this is
    t6.est_ekf line for line (asserted below). init = (S0, iters) replaces the
    update on y_0 by init_gn on the tied prior.
    Returns (smoothed x_0, its P, filtered x_n, its P, number of failed steps).
    """
    n = len(y)
    x, P = x0.copy(), P0.copy()
    xs, Ps, xp_l, Pp_l, Fs = [], [], [], [], []
    nfail = 0
    for k in range(n):
        if k == 0:
            if init is not None:
                x, P = init_gn(y[0], x0, init[0], ex, cv, init[1])
            else:
                x, P = kf_update(x, P, y[0], cv, meas_it)
        else:
            xprev, Pprev = x, P
            Q = q * np.eye(5)
            if qpar is not None:
                Q[1:, 1:] += np.diag(qpar(k, n))
            xb = xprev.copy()
            for i in range(dyn_it):
                F, fx, failed = step(xb, ex, t[k - 1], t[k])
                xp = fx + F @ (xprev - xb)
                Pp = F @ Pprev @ F.T + Q
                xn, Pn = kf_update(xp, Pp, y[k], cv, meas_it)
                if i < dyn_it - 1:
                    G = np.linalg.solve(Pp.T, (Pprev @ F.T).T).T
                    xb = xprev + G @ (xn - xp)
            nfail += failed
            xp_l.append(xp); Pp_l.append(Pp); Fs.append(F)
            x, P = xn, Pn
        xs.append(x.copy()); Ps.append(P.copy())
    xsm, Psm = xs[-1].copy(), Ps[-1].copy()
    if rts:
        for k in range(n - 2, -1, -1):
            G = Ps[k] @ Fs[k].T @ np.linalg.inv(Pp_l[k])
            xsm = xs[k] + G @ (xsm - xp_l[k])
            Psm = Ps[k] + G @ (Psm - Pp_l[k]) @ G.T
    return xsm, Psm, xs[-1], Ps[-1], nfail


def _sqrtm(P):
    try:
        return np.linalg.cholesky(P + 1e-14 * np.eye(len(P)))
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh(0.5 * (P + P.T))
        return V * np.sqrt(np.clip(w, 0, None))


def cub_points(x, P):
    n = len(x)
    L = _sqrtm(P) * np.sqrt(n)
    return np.vstack([x + L[:, j] for j in range(n)] + [x - L[:, j] for j in range(n)])


def moments(X):
    mu = X.mean(0)
    d = X - mu
    return mu, d.T @ d / len(X)


def ss_lnm_vec(Pts, ex):
    with np.errstate(all="ignore"):
        th, ke, K, L0 = np.exp(Pts[:, 0]), np.exp(Pts[:, 1]), np.exp(Pts[:, 2]), Pts[:, 3]
        Lt = np.maximum(ex.u(-1e-9) + L0, 1e-6)
        return np.log(np.maximum(C.m_ss(th, K, ke, Lt), 1e-9))


def prop_vec(lm, Pts, ex, t0, t1):
    STEPS[0] += 1
    with np.errstate(all="ignore"):
        th, ke, K, L0 = np.exp(Pts[:, 0]), np.exp(Pts[:, 1]), np.exp(Pts[:, 2]), Pts[:, 3]

        def rhs(t, z):
            mm = np.clip(z, 1e-9, 1e6)
            return ex.u(t) + L0 - th * mm / (K + mm) - ke * mm
        s = solve_ivp(rhs, (t0, t1), np.exp(lm), method="LSODA", rtol=1e-8, atol=1e-12)
        z = s.y[:, -1]
        return np.log(np.clip(z, 1e-9, 1e6))


def ckf(y, t, ex, x0, S0, cv, q=1e-10):
    """Cubature Kalman filter, tied prior: points on p, ln m0 = h(p) per point."""
    X = cub_points(x0[1:], S0)
    x, P = moments(np.c_[ss_lnm_vec(X, ex), X])
    for k in range(len(y)):
        if k > 0:
            X = cub_points(x, P)
            X[:, 0] = prop_vec(X[:, 0], X[:, 1:], ex, t[k - 1], t[k])
            x, P = moments(X)
            P = P + q * np.eye(5)
        x, P = kf_update(x, P, y[k], cv)
        P = 0.5 * (P + P.T)
    return x, P


def anneal(c, f):
    """Parameter process noise totalling c x the prior variance, decaying
    linearly to zero over the first fraction f of the record."""
    def qp(k, n, S=None):
        kk = f * n
        return c * qp.S * 2.0 / kk * max(0.0, 1.0 - k / kk)
    return qp


# ---------------------------------------------------------------- variants
def variant(name, b, y, x0, opt):
    """Run one estimator on one record. Returns (ln th error, reported sd, extras)."""
    t, ex, cv = b["t"], b["ex"], b["cv"]
    sd0 = sd_prior(b)
    S0 = np.diag(sd0[1:] ** 2)
    kind = opt["kind"]
    if kind == "prior":
        return x0[1] - b["p"][0], sd0[1], {}
    if kind == "ckf":
        xc = x0.copy()
        if opt.get("at_truth"):
            xc = np.r_[0.0, b["p"]]
        x, P = ckf(y, t, ex, xc, S0, cv)
        return x[1] - b["p"][0], np.sqrt(max(P[1, 1], 0.0)), {"x": x}
    xc, P0 = x0.copy(), np.diag(sd0 ** 2)
    if opt.get("at_truth"):
        s = sd0.copy()
        s[1:] *= opt.get("scale", 1.0)
        xc, P0 = np.r_[np.log(b["m"][0]), b["p"]], np.diag(s ** 2)
    init = None
    if opt.get("tied"):
        lm, J = ss_lnm(xc[1:], ex)
        xc = np.r_[lm, xc[1:]]
        P0 = tied_cov(J, S0)
        if opt.get("init_it"):
            init = (S0, opt["init_it"])
    qpar = None
    if opt.get("anneal"):
        qpar = anneal(*opt["anneal"])
        qpar.S = sd0[1:] ** 2
    xsm, Psm, xf, Pf, nf = ekf(y, t, ex, xc, P0, cv, q=opt.get("q", 1e-10), qpar=qpar,
                               meas_it=opt.get("meas_it", 1), dyn_it=opt.get("dyn_it", 1),
                               init=init, rts=opt.get("rts", True))
    use_f = opt.get("report") == "filtered"
    x, P = (xf, Pf) if use_f else (xsm, Psm)
    return x[1] - b["p"][0], np.sqrt(max(P[1, 1], 0.0)), dict(
        xsm=xsm, xf=xf, nfail=nf, x=x)


def run(name, b, D_, opt):
    STEPS[0] = 0
    t0 = time.time()
    e, s, ex_ = [], [], []
    for y, x0 in D_:
        ei, si, xi = variant(name, b, y, x0, opt)
        e.append(ei); s.append(si); ex_.append(xi)
    return dict(e=np.array(e), sd=np.array(s), extra=ex_, sec=time.time() - t0,
                steps=STEPS[0] / len(D_), st=stats(e, s))


# ---------------------------------------------------------------- checks
def jac_check(b, xb, k):
    """F from t6's _step against central differences of ln m1 over (ln m0, p)."""
    t, ex = b["t"], b["ex"]
    F, fx, _ = step(xb, ex, t[k - 1], t[k])
    Ffd = np.zeros(5)
    for j in range(5):
        h = 1e-5 * (max(abs(xb[4]), 1.0) if j == 4 else 1.0)
        a, c = xb.copy(), xb.copy()
        a[j] += h; c[j] -= h
        Ffd[j] = (step(a, ex, t[k - 1], t[k])[1][0] - step(c, ex, t[k - 1], t[k])[1][0]) / (2 * h)
    return F[0], Ffd


def ss_jac_check(b, p):
    _, J = ss_lnm(p, b["ex"])
    Jfd = np.zeros(4)
    for j in range(4):
        h = 1e-6 * (max(abs(p[3]), 1.0) if j == 3 else 1.0)
        a, c = p.copy(), p.copy()
        a[j] += h; c[j] -= h
        Jfd[j] = (np.log(D.steady(a, b["ex"])) - np.log(D.steady(c, b["ex"]))) / (2 * h)
    return J, Jfd


def fim_post_sd(b, tied):
    """Posterior sd of ln theta from FIM + prior, m0 tied to p or free."""
    sd0 = sd_prior(b)
    m, S = D.simulate(b["p"], b["ex"], b["t"])
    if tied:
        R = S / m[:, None]
        F = R.T @ R / b["cv"] ** 2 + np.diag(sd0[1:] ** -2)
        return float(np.sqrt(np.linalg.inv(F)[0, 0]))
    # m0 free: d m_k / d m0 by chaining t6's one-step derivative on the truth
    c, x, dmdm0 = 1.0, m[0], [1.0]
    for k in range(1, b["t"].size):
        with np.errstate(all="ignore"):
            x1, d0, _ = T6._step(x, b["p"], b["ex"], b["t"][k - 1], b["t"][k])
        c *= d0; x = x1; dmdm0.append(c)
    dmdm0 = np.array(dmdm0)
    R = np.hstack([(dmdm0 * m[0] / m)[:, None], (S - np.outer(dmdm0, S[0])) / m[:, None]])
    F = R.T @ R / b["cv"] ** 2 + np.diag(sd0 ** -2)
    return float(np.sqrt(np.linalg.inv(F)[1, 1]))


def parse_t6():
    """t6's printed rows for the three ODE-based estimators, both benches."""
    txt = (C.HERE / "out" / "t6_estimators.txt").read_text()
    out = {}
    for bench, block in re.findall(r"BENCH (\w):(.*?)(?=BENCH|\Z)", txt, re.S):
        for lab, key in (("(c) NLS on the ODE", "nls"), ("(c) EKF + RTS, one pass", "ekf"),
                         ("(c) iterated EKF smoother", "ieks")):
            line = [ln for ln in block.splitlines() if lab in ln][0]
            nums = re.findall(r"[+-]?\d+\.\d+(?:e[+-]?\d+)?|nan", line.split(lab)[1])
            out[(bench, key)] = dict(bias=float(nums[1]), rmse=float(nums[2]),
                                     robust=float(nums[3]), line=line.strip())
    return out


def whiskers(e):
    """matplotlib boxplot whiskers (1.5 IQR) on errors clipped to +-8, as t6's figure."""
    v = np.clip(e, -8, 8)
    q1, q3 = np.percentile(v, [25, 75])
    iqr = q3 - q1
    lo = v[v >= q1 - 1.5 * iqr].min()
    hi = v[v <= q3 + 1.5 * iqr].max()
    return q1, q3, lo, hi


def row(lab, r, steps=True):
    s = r["st"]
    tail = f" {r['steps']:7.0f}" if steps else ""
    div = f"  [{s['div']} diverged]" if s["div"] else ""
    return (f"   {lab:44} {s['bias']:+7.3f} {s['rmse']:7.3f} {s['robust']:7.3f} {s['gt1']:6.0%} "
            f"{s['repsd']:8.3f} {s['rmse'] / s['repsd']:6.1f}{tail}{div}")


HDR = (f"   {'variant':44} {'bias':>7} {'RMSE':>7} {'robust':>7} {'|e|>1':>6} "
       f"{'rep. sd':>8} {'RMSE/sd':>6} {'steps':>7}")


# ---------------------------------------------------------------- main
def main():
    print("=" * 100)
    print("  ITEM 5. THE ONE-PASS EKF + RTS OF t6: DIAGNOSIS AND FIXED VARIANTS")
    print("=" * 100)
    B = {nm: T6.bench(nm) for nm in ("A", "B")}
    DR = {nm: draws(B[nm], REPS, EVAL_SEED) for nm in B}
    DT = {nm: draws(B[nm], TUNE_REPS, TUNE_SEED) for nm in B}
    T6P = parse_t6()

    # ---------------------------------------------------------------- 0
    print("\n  0. REPRODUCTION ON t6'S DRAWS (default_rng(1), first 80 replicates, y then x0)")
    base = {}
    for nm, b in B.items():
        sd0 = sd_prior(b)
        dmax, e_t6 = 0.0, []
        for y, x0 in DR[nm]:
            xs_t6, _ = T6.est_ekf(y, b["t"], b["ex"], x0, np.diag(sd0 ** 2), b["cv"])
            xs_me = ekf(y, b["t"], b["ex"], x0, np.diag(sd0 ** 2), b["cv"])[0]
            dmax = max(dmax, float(np.max(np.abs(xs_t6 - xs_me))))
            e_t6.append(xs_t6[1] - b["p"][0])
        assert dmax < 1e-9, dmax
        st = stats(e_t6)
        p = T6P[(nm, "ekf")]
        print(f"   bench {nm}: t6.est_ekf on these draws: bias {st['bias']:+.3f}, RMSE {st['rmse']:.3f}, "
              f"robust {st['robust']:.3f};  t6 printed: {p['bias']:+.3f}, {p['rmse']:.3f}, {p['robust']:.3f}")
        assert abs(st["bias"] - p["bias"]) < 6e-4 and abs(st["rmse"] - p["rmse"]) < 6e-4 \
            and abs(st["robust"] - p["robust"]) < 6e-4
        print(f"            this file's copy of the filter equals t6.est_ekf to {dmax:.1e} (max |dx|)")
        q1, q3, lo, hi = whiskers(np.array(e_t6))
        print(f"            box-plot of ln(th_hat/th) as t6 draws it: quartiles {q1:+.2f} / {q3:+.2f}, "
              f"whiskers {lo:+.2f} to {hi:+.2f}; max |e| {np.max(np.abs(e_t6)):.2f}")
        base[nm] = np.array(e_t6)
    print("   t6's printed rows for comparison:")
    for nm in B:
        for k in ("nls", "ieks"):
            print(f"     bench {nm}: {T6P[(nm, k)]['line']}")
    print("   NOTE t6.run_bench starts NLS from [truth, truth +- offsets]: its NLS sees the true")
    print("   parameters as a start. The Kalman estimators start from the drawn prior centre.")

    # ---------------------------------------------------------------- 1
    print("\n  1. COMPONENT CHECKS")
    bad = []
    for nm, b in B.items():
        sd0 = sd_prior(b)
        print(f"   bench {nm}: theta {b['th']:.0f}, {b['t'].size} samples, CV {b['cv']:.0%}")
        # (1) initialisation
        J, Jfd = ss_jac_check(b, b["p"])
        dl = np.array([x0[0] - np.log(b["m"][0]) for _, x0 in DR[nm]])
        Jc = [ss_lnm(x0[1:], b["ex"])[1] for _, x0 in DR[nm]]
        sd_tied = np.array([np.sqrt(j @ np.diag(sd0[1:] ** 2) @ j) for j in Jc])
        print(f"     (1) init: d ln m0/dp at truth {np.round(J, 4)} (central diff. max rel. err "
              f"{np.max(np.abs(J - Jfd) / np.abs(Jfd)):.1e})")
        print(f"         ln m_ss(drawn centre) - ln m0(true): sd over the draws {dl.std():.3f}, "
              f"max |.| {np.abs(dl).max():.3f}; t6's P0 gives it sd {sd0[0]:.2f}, independent of p")
        print(f"         the tied prior gives sd(ln m0) = sqrt(J S J') = {np.median(sd_tied):.3f} "
              f"(median over the drawn centres)")
        # (2) process noise
        print(f"     (2) q = 1e-10: random-walk sd accumulated over the record "
              f"{np.sqrt(1e-10 * (b['t'].size - 1)):.1e}")
        # (3) Jacobian
        x_true = np.r_[np.log(b["m"][0]), b["p"]]
        worst = 0.0
        for xb, k in ((x_true, 1), (np.r_[np.log(b["m"][b["t"].size // 2]), b["p"]], b["t"].size // 2 + 1),
                      (DR[nm][0][1], 1), (DR[nm][1][1], 5)):
            F, Ffd = jac_check(b, xb, k)
            worst = max(worst, float(np.max(np.abs(F - Ffd) / np.maximum(np.abs(Ffd), 1e-12))))
        print(f"     (3) log-state Jacobian F[0,:] against central differences of ln m1, "
              f"4 points: max rel. err {worst:.1e}")
        if worst > 1e-3:
            bad.append(f"{nm}: Jacobian")
        # (4) measurement update and RTS: exact-information identity
        y0 = b["m"].copy()
        xs, Ps, xf, Pf, _ = ekf(y0, b["t"], b["ex"], x_true, np.diag(sd0 ** 2), b["cv"])
        fr = fim_post_sd(b, tied=False)
        lm, Jt = ss_lnm(b["p"], b["ex"])
        xs2, Ps2, _, _, _ = ekf(y0, b["t"], b["ex"], np.r_[lm, b["p"]],
                                tied_cov(Jt, np.diag(sd0[1:] ** 2)), b["cv"])
        ft = fim_post_sd(b, tied=True)
        print(f"     (4) noiseless data at the truth: EKF+RTS sd(ln th) {np.sqrt(Psm_(Ps)):.5f} vs "
              f"FIM+prior, m0 free {fr:.5f};  tied P0 {np.sqrt(Psm_(Ps2)):.5f} vs FIM+prior, "
              f"m0 tied {ft:.5f}")
        assert abs(np.sqrt(Psm_(Ps)) / fr - 1) < 1e-3 and abs(np.sqrt(Psm_(Ps2)) / ft - 1) < 1e-3
        print(f"         estimate error there {xs[1] - b['p'][0]:+.1e} and {xs2[1] - b['p'][0]:+.1e}")
        dsm = max(float(np.max(np.abs(r["xsm"][1:] - r["xf"][1:])))
                  for r in (variant(nm, b, y, x0, dict(kind="ekf"))[2] for y, x0 in DR[nm]))
        print(f"         RTS: max over the 80 records of |smoothed p at t0 - filtered p at t_n| {dsm:.1e}")
        if dsm > 1e-2:
            bad.append(f"{nm}: RTS")
        # (5) failure path, and the IEKF update
        r0 = run(nm, b, DR[nm], dict(kind="ekf"))
        ri = run(nm, b, DR[nm], dict(kind="ekf", meas_it=5))
        nf = sum(x["nfail"] for x in r0["extra"])
        dd = float(np.max(np.abs(r0["e"] - ri["e"])))
        print(f"     (5) _step's silent failure branch (returns F = 0) taken {nf} times in 80 records")
        print(f"         iterated measurement update (5 iterations) vs t6: max |d ln th_hat| {dd:.1e}")
        if nf:
            bad.append(f"{nm}: integrator failures")
        if dd > 1e-5:
            bad.append(f"{nm}: IEKF differs")
    if bad:
        print("   RESULT: a check failed: " + "; ".join(bad))
    else:
        print("   RESULT: every check passed (Jacobian to 1e-3, RTS shift of p below 1e-2, no")
        print("   integrator failure, IEKF update equal to 1e-5, and the asserted information")
        print("   identity at the truth). With h(x) = x[0] linear in the state, an IEKF")
        print("   measurement update cannot differ from t6's update except by rounding.")

    # ---------------------------------------------------------------- mechanism
    print("\n  2. WHERE THE ERROR COMES FROM (bench B, t6's filter, t6's draws)")
    b = B["B"]
    r0 = run("B", b, DR["B"], dict(kind="ekf"))
    xe = np.array([x["x"][1:] - b["p"] for x in r0["extra"]])
    big = np.abs(xe[:, 0]) > 2
    print(f"   reported sd(ln th) median {r0['st']['repsd']:.3f}; actual RMSE {r0['st']['rmse']:.3f}; "
          f"median |error|/reported sd {np.median(np.abs(r0['e']) / r0['sd']):.1f}")
    print(f"   corr(error ln th, error ln Km) over the 80 records {np.corrcoef(xe[:, 0], xe[:, 2])[0, 1]:+.3f}")
    print(f"   records with |error ln th| > 2: {big.sum()}; their errors (ln th, ln ke, ln Km, L0), "
          f"and ln th - ln Km:")
    for v in xe[big]:
        print(f"     {v[0]:+7.2f} {v[1]:+6.2f} {v[2]:+7.2f} {v[3]:+8.1f}   {v[0] - v[2]:+6.2f}")
    kb = xe[:, 2]
    print(f"   of those {big.sum()}: ln Km error > +2 (towards the linear limit, K_m >> m) in "
          f"{int((big & (kb > 2)).sum())}, < -2 (towards the saturated limit) in "
          f"{int((big & (kb < -2)).sum())}; of all 80 records: {int((kb > 2).sum())} and "
          f"{int((kb < -2).sum())}")
    print("   started AT THE TRUTH, same data, prior width scaled (not fair comparators: the")
    print("   prior is centred on the truth):")
    print(HDR)
    scal = {}
    for sc in (0.1, 0.3, 0.5, 1.0):
        r = run("B", b, DR["B"], dict(kind="ekf", at_truth=True, scale=sc))
        scal[sc] = r
        print(row(f"t6 EKF at truth, prior sd x {sc}", r))
    rct = run("B", b, DR["B"], dict(kind="ckf", at_truth=True))
    print(row("(d) cubature filter at truth, prior sd x 1", rct))

    # ---------------------------------------------------------------- tuning
    print(f"\n  3. TUNING on separate draws (default_rng({TUNE_SEED}), {TUNE_REPS} records, bench B),"
          f" chosen by robust spread")
    print(HDR)
    tune = {}
    for c in (0.05, 0.2, 1.0):
        for f in (0.25, 0.5, 0.8):
            r = run("B", b, DT["B"], dict(kind="ekf", tied=True, init_it=30, anneal=(c, f),
                                          report="filtered"))
            tune[("anneal", c, f)] = r
            print(row(f"(b) anneal: total {c} x prior var over first {f:.0%}", r))
    for it in (3, 5):
        r = run("B", b, DT["B"], dict(kind="ekf", tied=True, init_it=30, dyn_it=it,
                                      report="filtered"))
        tune[("dyn", it)] = r
        print(row(f"(c) lag-1 iterated transition, {it} iterations", r))
    best_a = min((k for k in tune if k[0] == "anneal"), key=lambda k: tune[k]["st"]["robust"])
    best_d = min((k for k in tune if k[0] == "dyn"), key=lambda k: tune[k]["st"]["robust"])
    print(f"   chosen: (b) total {best_a[1]} x prior var over first {best_a[2]:.0%};"
          f" (c) {best_d[1]} iterations")

    # ---------------------------------------------------------------- evaluation
    V = [("t6 EKF+RTS (diag P0, q 1e-10)", dict(kind="ekf")),
         ("t6 EKF, q = 0, filtered end", dict(kind="ekf", q=0.0, rts=False, report="filtered")),
         ("(c0) t6 + iterated measurement update", dict(kind="ekf", meas_it=5)),
         ("(a) tied P0", dict(kind="ekf", tied=True, report="filtered")),
         ("(a') tied P0 + iterated update on y0", dict(kind="ekf", tied=True, init_it=30,
                                                       report="filtered")),
         ("(b) (a') + annealed parameter noise", dict(kind="ekf", tied=True, init_it=30,
                                                      anneal=best_a[1:], report="filtered")),
         ("(c) (a') + lag-1 iterated transition", dict(kind="ekf", tied=True, init_it=30,
                                                       dyn_it=best_d[1], report="filtered")),
         ("(b)+(c)", dict(kind="ekf", tied=True, init_it=30, dyn_it=best_d[1],
                          anneal=best_a[1:], report="filtered")),
         ("(d) cubature filter, tied prior", dict(kind="ckf"))]
    REF = [("prior centre alone (no data)", dict(kind="prior")),
           ("t6 EKF started at the truth", dict(kind="ekf", at_truth=True, scale=1.0))]
    print(f"\n  4. EVALUATION on t6's draws ({REPS} records per bench). ln(th_hat/th): median bias,")
    print("     RMSE, robust = half the 16-84% range, share with |error| > 1, median reported sd,")
    print("     RMSE / reported sd, transition evaluations per record")
    R = {}
    for nm, b in B.items():
        print(f"\n   BENCH {nm}: FIM+prior sd(ln th) {fim_post_sd(b, True):.3f} (m0 tied), "
              f"{fim_post_sd(b, False):.3f} (m0 free); FIM sd, no prior "
              f"{np.sqrt(np.linalg.inv(D.fim(b['p'], b['ex'], b['t'], cv=b['cv']))[0, 0]):.3f}")
        print(HDR)
        for lab, opt in V + REF:
            r = run(nm, b, DR[nm], opt)
            R[(nm, lab)] = r
            print(row(lab, r))
        for k in ("nls", "ieks"):
            p = T6P[(nm, k)]
            print(f"   {'t6 printed: ' + ('NLS on the ODE' if k == 'nls' else 'iterated EKF smoother'):44}"
                  f" {p['bias']:+7.3f} {p['rmse']:7.3f} {p['robust']:7.3f}")
        if nm == "B":
            n = b["t"].size
            print(f"   t6.est_ieks re-linearises all {n - 1} transitions per pass, up to 60 passes"
                  f" (up to {60 * (n - 1) + (n - 1)} evaluations)")

    # ---------------------------------------------------------------- verdict
    one_pass = [lab for lab, _ in V]
    bestB = min(one_pass, key=lambda l: R[("B", l)]["st"]["robust"])
    sb = R[("B", bestB)]["st"]
    nls, ieks = T6P[("B", "nls")], T6P[("B", "ieks")]
    fimB = fim_post_sd(B["B"], True)
    fair = sb["robust"] <= 1.5 * fimB and sb["rmse"] <= 2 * nls["rmse"]
    pr = R[("B", "prior centre alone (no data)")]["st"]
    print("\n  5. VERDICT")
    print(f"   best one-pass variant on bench B by robust spread: {bestB}")
    print(f"     robust {sb['robust']:.3f}, RMSE {sb['rmse']:.3f}, |e|>1 in {sb['gt1']:.0%}; "
          f"against FIM+prior {fimB:.3f}, t6 NLS RMSE {nls['rmse']:.3f}, t6 IEKS robust {ieks['robust']:.3f}")
    print(f"     prior centre alone: RMSE {pr['rmse']:.3f}, robust {pr['robust']:.3f}")
    print(f"   criterion for 'fair': robust <= 1.5 x FIM+prior sd AND RMSE <= 2 x NLS RMSE -> "
          f"{'MET' if fair else 'NOT MET'}")
    worse_than_prior = [l for l in one_pass if R[("B", l)]["st"]["rmse"] > pr["rmse"]]
    print(f"   one-pass variants with RMSE above the prior centre's on bench B: "
          f"{len(worse_than_prior)} of {len(one_pass)}")
    rt = R[("B", "t6 EKF started at the truth")]["st"]
    print(f"   t6's filter started at the truth: RMSE {rt['rmse']:.3f}, robust {rt['robust']:.3f}"
          f" (prior centre error removed entirely)")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw=dict(width_ratios=[1.5, 1]))
    a = ax[0]
    labs = one_pass + ["t6 EKF started at the truth", "prior centre alone (no data)"]
    data = [np.clip(R[("B", l)]["e"][np.isfinite(R[("B", l)]["e"])], -8, 8) for l in labs]
    bp = a.boxplot(data, orientation="horizontal", widths=0.55, patch_artist=True,
                   showfliers=False, medianprops=dict(color=INK, lw=1.5))
    for j, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(C1 if j < len(one_pass) else C3)
        patch.set_alpha(0.35); patch.set_edgecolor(C1 if j < len(one_pass) else C3)
    from matplotlib.patches import Patch
    a.legend(handles=[Patch(facecolor=C1, alpha=0.35, edgecolor=C1, label="one-pass filter"),
                      Patch(facecolor=C3, alpha=0.35, edgecolor=C3, label="reference")],
             loc="upper right", ncol=2, fontsize=7.5, frameon=False)
    a.axvspan(-fimB, fimB, color=C2, alpha=0.18, lw=0)
    a.text(fimB + 0.1, len(labs) + 0.35, "± FIM sd", color=INK, fontsize=7.5)
    a.set_yticks(range(1, len(labs) + 1), labs, fontsize=7.5)
    a.set_ylim(0.4, len(labs) + 1.4)
    a.axvline(0, color=INK, lw=0.8)
    a.set_xlabel(r"ln($\hat\theta$/θ), clipped at ±8")
    a.set_title("TB-B: one-pass filters (blue), references (green)")
    finish(a)
    a = ax[1]
    sc = np.array(sorted(scal))
    a.plot(sc, [scal[s]["st"]["rmse"] for s in sc], "o-", color=C1, lw=2, ms=6, label="actual RMSE")
    a.plot(sc, [scal[s]["st"]["repsd"] for s in sc], "s--", color=C2, lw=2, ms=6,
           label="filter's own sd")
    a.text(sc[-1] * 0.9, scal[sc[-1]]["st"]["rmse"], "actual RMSE", color=INK, fontsize=8,
           ha="right", va="center")
    a.text(sc[-1] * 0.9, scal[sc[-1]]["st"]["repsd"] * 0.88, "filter's own sd", color=INK,
           fontsize=8, ha="right", va="top")
    a.legend(loc="upper left", fontsize=7.5, frameon=False)
    a.axhline(fimB, color=MUTED, lw=1, ls=":")
    a.text(sc[0] * 1.05, fimB * 0.9, "FIM+prior sd at x 1", color=MUTED, fontsize=7.5, va="top")
    a.set_xscale("log"); a.set_yscale("log")
    plain(a)
    a.set_xlabel("prior sd, multiple of t6's (filter started at the truth)")
    a.set_ylabel("sd of ln θ̂")
    a.set_title("TB-B: t6's EKF centred on the truth")
    finish(a)
    fig.tight_layout()
    save(fig, "t6d_ekf")


def Psm_(P):
    return float(P[1, 1])


if __name__ == "__main__":
    main()
