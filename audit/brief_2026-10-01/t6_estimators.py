"""
Task 6. Three estimators for a transient, compared on synthetic noisy records.

  (a)   FINITE DIFFERENCE, linear in lumped parameters (the brief's). Multiply
        the ODE by (K_m + m), psi_k = (y_{k+1} - y_k)/tau - u_k:
            psi m = L0 K_m + (L0 - th - ke K_m) m - ke m^2 - K_m psi
        regress z = psi*m on X = [1, m, -m^2, -psi]; p = (L0 Km, L0 - th - ke Km,
        ke, Km); unlump L0 = p1/p4, theta = p1/p4 - p2 - p3 p4.
        (a-mid) the same with m read at the interval midpoint.
  (b)   INTEGRAL, saturated (the brief's): y_{k+1} - y_k - int u = (L0 - th) tau
        - ke int m.  Two unknowns, c = L0 - theta and k_e. theta is NOT
        separately identified by this form at all: a saturated enzyme is a
        constant sink, indistinguishable from a constant source.
  (b')  INTEGRAL, general, per interval: ... = L0 tau - th int g - ke int m,
        linear given K_m, K_m profiled (variable projection).
  (b'') the same CUMULATIVE from t_0, which is where integration actually
        averages noise rather than merely rescaling the difference.
  (c)   NONLINEAR least squares on ln y through the ODE (forward
        sensitivities), and an augmented-state EKF with an RTS smoother.

Two benches, because at the brief's target nothing is identifiable:
  TB-A  the brief's target: theta = 150, diet (half the load) withdrawn at t = 0,
        hourly for 48 h, 5% assay CV. The information matrix gives sd(ln theta)
        = 383 with no prior (Task 7).
  TB-B  an identifiable BENCH, not a design: theta = 600 (20%), load raised
        six-fold for 24 h then withdrawn, every 15 min for 48 h, 1% CV, so
        that bias and variance are measurable at all.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize_scalar

import core as C
import dyn as D
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain


def bench(name):
    if name == "A":
        th = C.TH_5
        L = float(C.N4.anchor(th, C.KM0, C.KE, "A")[0])
        ex = D.Experiment([(-1, 0.5 * L), (0, 0.0)], L0=0.5 * L)
        t = np.arange(0, 2 + 1e-9, 1 / 24)
        cv = 0.05
    else:
        th = 0.20 * C.TH0
        L = float(C.N4.anchor(th, C.KM0, C.KE, "A")[0])
        ex = D.Experiment([(-1, 0.5 * L), (0, 5.5 * L), (1, 0.0)], L0=0.5 * L)
        t = np.arange(0, 2 + 1e-9, 1 / 96)
        cv = 0.01
    p = D.params(theta=th, ex=ex)
    m, _ = D.simulate(p, ex, t, sens=False)
    return dict(th=th, L=L, ex=ex, t=t, cv=cv, p=p, m=m)


def u_bar(ex, t):
    """Mean known load over each interval; breaks fall on sample times."""
    return np.array([ex.u(0.5 * (a + b)) for a, b in zip(t[:-1], t[1:])])


# ------------------------------------------------------------- (a) FD
def est_fd(y, t, ex, mid=False):
    tau = np.diff(t)
    psi = np.diff(y) / tau - u_bar(ex, t)
    mm = 0.5 * (y[:-1] + y[1:]) if mid else y[:-1]
    X = np.stack([np.ones_like(mm), mm, -mm ** 2, -psi], axis=1)
    z = psi * mm
    pp, *_ = np.linalg.lstsq(X, z, rcond=None)
    res = z - X @ pp
    Km = pp[3]
    L0 = pp[0] / Km
    th = L0 - pp[1] - pp[2] * Km
    r1 = np.corrcoef(res[:-1], res[1:])[0, 1]
    return dict(th=th, ke=pp[2], Km=Km, L0=L0, p=pp, X=X, r1=r1)


def lumped_true(b):
    th, ke, K, L0, _ = D.unpack(b["p"])
    return np.array([L0 * K, L0 - th - ke * K, ke, K])


# ------------------------------------------------------------- (b) integral
def est_int_sat(y, t, ex):
    tau = np.diff(t)
    d = np.diff(y) - u_bar(ex, t) * tau
    I = tau * 0.5 * (y[:-1] + y[1:])
    X = np.stack([tau, -I], axis=1)
    pp, *_ = np.linalg.lstsq(X, d, rcond=None)
    return dict(c=pp[0], ke=pp[1], X=X)


def _vp_design(y, t, ex, Km, cumulative):
    tau = np.diff(t)
    g = y / (Km + y)
    G = tau * 0.5 * (g[:-1] + g[1:])
    I = tau * 0.5 * (y[:-1] + y[1:])
    d = np.diff(y) - u_bar(ex, t) * tau
    if cumulative:
        tau, G, I, d = np.cumsum(tau), np.cumsum(G), np.cumsum(I), np.cumsum(d)
    return np.stack([tau, -G, -I], axis=1), d


def est_int_gen(y, t, ex, cumulative=False):
    def rss(lk):
        X, d = _vp_design(y, t, ex, np.exp(lk), cumulative)
        pp, *_ = np.linalg.lstsq(X, d, rcond=None)
        return float(np.sum((d - X @ pp) ** 2))
    grid = np.linspace(np.log(5), np.log(5000), 61)
    j = int(np.argmin([rss(g) for g in grid]))
    lo, hi = grid[max(j - 1, 0)], grid[min(j + 1, grid.size - 1)]
    lk = minimize_scalar(rss, bounds=(lo, hi), method="bounded").x
    X, d = _vp_design(y, t, ex, np.exp(lk), cumulative)
    pp, *_ = np.linalg.lstsq(X, d, rcond=None)
    return dict(L0=pp[0], th=pp[1], ke=pp[2], Km=np.exp(lk), X=X)


# ------------------------------------------------------------- (c) NLS
def est_nls(y, t, ex, p0s):
    ly = np.log(y)

    def res(p):
        m, _ = D.simulate(p, ex, t, sens=False, rtol=1e-7)
        return np.log(np.maximum(m, 1e-12)) - ly

    def jac(p):
        m, S = D.simulate(p, ex, t, rtol=1e-7)
        return S / m[:, None]

    lo = [np.log(1e-2), np.log(1e-3), np.log(1.0), -1e4]
    hi = [np.log(1e5), np.log(1e3), np.log(1e5), 1e5]
    best = None
    for p0 in p0s:
        try:
            s = least_squares(res, np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6),
                              jac=jac, bounds=(lo, hi), method="trf", x_scale="jac",
                              xtol=1e-10, ftol=1e-10, max_nfev=200)
        except Exception:
            continue
        if best is None or s.cost < best.cost:
            best = s
    return best.x


# ------------------------------------------------------------- (c) EKF + RTS
def _step(m, p, ex, t0, t1):
    """Propagate m over [t0, t1] with dm1/dm0 and dm1/dp (zero initial dm/dp)."""
    th, ke, K, L0, kap = D.unpack(p)

    def rhs(t, z):
        mm = min(max(z[0], 1e-9), 1e6)
        g = mm / (K + mm)
        Cp = th * K / (K + mm) ** 2 + ke
        dCdp = [th * g, ke * mm, -th * K * mm / (K + mm) ** 2, -1.0]
        return [ex.u(t) + L0 - th * g - ke * mm, -Cp * z[1]] + \
               [-Cp * z[2 + j] - dCdp[j] for j in range(4)]
    with np.errstate(all="ignore"):
        s = solve_ivp(rhs, (t0, t1), [m, 1.0, 0, 0, 0, 0], method="LSODA",
                      rtol=1e-8, atol=1e-12)
    z = s.y[:, -1]
    if not (s.success and np.all(np.isfinite(z))):
        return float(np.clip(m, 1e-9, 1e6)), 0.0, np.zeros(4)
    return float(np.clip(z[0], 1e-9, 1e6)), z[1], z[2:]


def est_ekf(y, t, ex, x0, P0, cv, q=1e-10):
    """Augmented state x = (ln m, ln th, ln ke, ln Km, L0); RTS smoother."""
    n = len(y)
    xs, Ps, xp_l, Pp_l, Fs = [], [], [], [], []
    x, P = x0.copy(), P0.copy()
    # the first state is the pre-step steady state: start from the prior and
    # update on y_0 through ln m
    for k in range(n):
        if k > 0:
            m = np.exp(x[0])
            m1, dm0, dp = _step(m, x[1:], ex, t[k - 1], t[k])
            m1 = max(m1, 1e-9)
            F = np.eye(5)
            F[0, 0] = dm0 * m / m1
            F[0, 1:] = dp / m1
            xp = x.copy(); xp[0] = np.log(m1)
            Pp = F @ P @ F.T + q * np.eye(5)
            xp_l.append(xp); Pp_l.append(Pp); Fs.append(F)
            x, P = xp, Pp
        Hx = np.zeros(5); Hx[0] = 1.0
        S = Hx @ P @ Hx + cv ** 2
        Kg = P @ Hx / S
        x = x + Kg * (np.log(y[k]) - x[0])
        P = P - np.outer(Kg, Hx @ P)
        xs.append(x.copy()); Ps.append(P.copy())
    # RTS
    xsm, Psm = xs[-1].copy(), Ps[-1].copy()
    for k in range(n - 2, -1, -1):
        G = Ps[k] @ Fs[k].T @ np.linalg.inv(Pp_l[k])
        xsm = xs[k] + G @ (xsm - xp_l[k])
        Psm = Ps[k] + G @ (Psm - Pp_l[k]) @ G.T
    return xsm, Psm


def est_ieks(y, t, ex, x0, P0, cv, iters=60, q=1e-10, tol=1e-7):
    """
    Iterated extended Kalman smoother: the filter linearises the dynamics about
    the previous pass's SMOOTHED trajectory rather than about its own running
    mean (Bell 1994). With static parameters and deterministic dynamics this is
    Gauss-Newton on the MAP objective, so converged it is the MAP; a single pass
    of it is est_ekf.
    """
    n = len(y)
    xbar = [x0.copy()]
    for k in range(1, n):          # nominal trajectory from the prior mean
        m1, _, _ = _step(np.exp(xbar[-1][0]), xbar[-1][1:], ex, t[k - 1], t[k])
        xx = xbar[-1].copy(); xx[0] = np.log(max(m1, 1e-9)); xbar.append(xx)
    for it in range(iters):
        xs, Ps, xp_l, Pp_l, Fs = [], [], [], [], []
        x, P = x0.copy(), P0.copy()
        for k in range(n):
            if k > 0:
                xb = xbar[k - 1]
                mb = np.exp(xb[0])
                m1, dm0, dp = _step(mb, xb[1:], ex, t[k - 1], t[k])
                m1 = max(m1, 1e-9)
                F = np.eye(5); F[0, 0] = dm0 * mb / m1; F[0, 1:] = dp / m1
                fx = xb.copy(); fx[0] = np.log(m1)
                xp = fx + F @ (x - xb)
                Pp = F @ P @ F.T + q * np.eye(5)
                xp_l.append(xp); Pp_l.append(Pp); Fs.append(F)
                x, P = xp, Pp
            Hx = np.zeros(5); Hx[0] = 1.0
            S = Hx @ P @ Hx + cv ** 2
            Kg = P @ Hx / S
            x = x + Kg * (np.log(y[k]) - x[0])
            P = P - np.outer(Kg, Hx @ P)
            xs.append(x.copy()); Ps.append(P.copy())
        sm = [None] * n
        sm[-1] = xs[-1].copy()
        Psm = Ps[-1].copy()
        for k in range(n - 2, -1, -1):
            G = Ps[k] @ Fs[k].T @ np.linalg.inv(Pp_l[k])
            sm[k] = xs[k] + G @ (sm[k + 1] - xp_l[k])
        # keep the parameters inside a physical box between passes
        Lb = abs(x0[4]) * 2 + 1.0
        new = [np.r_[np.clip(v[0], -20, 14), np.clip(v[1:4], [-5, -7, 0], [12, 7, 12]),
                     np.clip(v[4], -20 * Lb, 20 * Lb)] for v in sm]
        #  damping: Gauss-Newton overshoots from a prior centre that can be 20x
        #  off in theta, so a pass may move the parameters by at most 0.5 in
        #  log units (a trust region), the state following the parameters
        d = new[0][1:4] - xbar[0][1:4]
        a = min(1.0, 0.5 / max(np.abs(d).max(), 1e-12))
        new = [b_ + a * (v - b_) for v, b_ in zip(new, xbar)]
        step = max(np.abs(a_ - b_).max() for a_, b_ in zip(new, xbar))
        xbar = new
        if step < tol:
            break
    return xbar[0], it + 1


def ncond(X):
    X = np.asarray(X, float)
    return float(np.linalg.cond(X / np.linalg.norm(X, axis=0)))


def summarise(vals, truth, log=True):
    v = np.asarray(vals, float)
    if log:
        pos = v > 0
        if not pos.any():
            return dict(bias=np.nan, rmse=np.nan, robust=np.nan, neg=1.0)
        e = np.log(v[pos]) - np.log(truth)
        return dict(bias=float(np.median(e)), rmse=float(np.sqrt(np.mean(e ** 2))),
                    robust=float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16))),
                    neg=float(1 - pos.mean()))
    e = v - truth
    return dict(bias=float(np.median(e)), rmse=float(np.sqrt(np.mean(e ** 2))),
                robust=float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16))), neg=0.0)


def run_bench(name, reps_lin=400, reps_nl=80, seed=1):
    b = bench(name)
    rng = np.random.default_rng(seed)
    th0, ke0, K0, L00, _ = D.unpack(b["p"])
    lt = lumped_true(b)
    out = {k: [] for k in ("fd", "fdm", "sat", "gen", "cum", "nls", "ekf", "ieks")}
    pl, r1s = [], []
    # noiseless (discretisation only) for (a): truncation bias, no EIV
    fd0 = est_fd(b["m"], b["t"], b["ex"])
    fd0m = est_fd(b["m"], b["t"], b["ex"], mid=True)
    for k in range(reps_lin):
        y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
        a = est_fd(y, b["t"], b["ex"]); out["fd"].append(a); pl.append(a["p"]); r1s.append(a["r1"])
        out["fdm"].append(est_fd(y, b["t"], b["ex"], mid=True))
        out["sat"].append(est_int_sat(y, b["t"], b["ex"]))
        out["gen"].append(est_int_gen(y, b["t"], b["ex"]))
        out["cum"].append(est_int_gen(y, b["t"], b["ex"], cumulative=True))
        if k < reps_nl:
            ps = [b["p"], b["p"] + np.array([0.7, -0.3, 0.5, -0.2 * L00]),
                  b["p"] + np.array([-0.7, 0.3, -0.5, 0.2 * L00])]
            out["nls"].append(est_nls(y, b["t"], b["ex"], ps))
            # EKF prior: centres drawn about the truth (pseudo-observations)
            sd0 = np.array([0.3, 1.5, 0.8, 1.0, 0.5 * L00])
            x0 = np.r_[np.log(b["m"][0]), b["p"]] + sd0 * rng.normal(size=5)
            x0[0] = np.log(D.steady(x0[1:], b["ex"]))
            xs, _ = est_ekf(y, b["t"], b["ex"], x0, np.diag(sd0 ** 2), b["cv"])
            out["ekf"].append(xs[1:])
            xi, _ = est_ieks(y, b["t"], b["ex"], x0, np.diag(sd0 ** 2), b["cv"])
            out["ieks"].append(xi[1:])
    pl = np.array(pl)
    # theoretical MA(1): consecutive differences share one sample
    res = dict(b=b, lt=lt, fd0=fd0, fd0m=fd0m, out=out, pl=pl, r1=np.array(r1s))
    F = D.fim(b["p"], b["ex"], b["t"], cv=b["cv"])
    res["fim_sd"] = np.sqrt(np.diag(np.linalg.inv(F)))
    m, S = D.simulate(b["p"], b["ex"], b["t"])
    res["cond_nls"] = ncond(S / m[:, None])
    return res


def report(name, R):
    b, out = R["b"], R["out"]
    th0, ke0, K0, L00, _ = D.unpack(b["p"])
    print(f"\n  BENCH {name}: theta = {th0:.0f}, m from {b['m'].min():.1f} to {b['m'].max():.1f}, "
          f"{b['t'].size} samples, tau = {24*(b['t'][1]-b['t'][0]):.2f} h, CV {100*b['cv']:.0f}%")
    lam = D.lam(b["p"], b["ex"], b["m"])
    print(f"     lambda along the path {lam.min():.2f} to {lam.max():.2f} /d (1/lambda "
          f"{24/lam.max():.2f} to {24/lam.min():.2f} h); theta g' share of lambda "
          f"{((lam-ke0)/lam).min():.1%} to {((lam-ke0)/lam).max():.1%}")
    print(f"     FIM sd, no priors: ln th {R['fim_sd'][0]:.3f}, ln ke {R['fim_sd'][1]:.3f}, "
          f"ln Km {R['fim_sd'][2]:.3f}, L0 {R['fim_sd'][3]:.1f} (truth L0 {L00:.1f})")
    lt = R["lt"]
    pl = R["pl"]
    print("     (a) FD lumped coefficients: truth | noiseless, left point | noiseless, midpoint"
          " | noisy left point, median (ratio to truth)")
    for j, nm in enumerate(("L0 Km", "L0-th-ke Km", "ke", "Km")):
        print(f"        {nm:12} {lt[j]:12.2f} | {R['fd0']['p'][j]:12.2f} | {R['fd0m']['p'][j]:12.2f}"
              f" | {np.median(pl[:, j]):12.2f} ({np.median(pl[:, j])/lt[j]:+.3f})")
    print(f"        unlumped from NOISELESS samples: left point theta {R['fd0']['th']:.1f}, midpoint "
          f"theta {R['fd0m']['th']:.1f} (truth {th0:.0f})")
    print(f"     (a) lag-1 residual autocorrelation: median {np.median(R['r1']):+.3f} "
          f"(pure shared-sample noise predicts -0.5)")
    print(f"     {'estimator':30} {'cond(X)':>9} {'ln th bias':>10} {'RMSE':>8} {'robust':>8} "
          f"{'th<=0':>6} {'ln ke bias':>10} {'RMSE':>8}")
    rows = (("(a) FD, left point", "fd", ncond(out["fd"][0]["X"])),
            ("(a) FD, midpoint", "fdm", ncond(out["fdm"][0]["X"])),
            ("(b) integral, saturated", "sat", ncond(out["sat"][0]["X"])),
            ("(b') integral, general", "gen", ncond(out["gen"][0]["X"])),
            ("(b'') integral, cumulative", "cum", ncond(out["cum"][0]["X"])),
            ("(c) NLS on the ODE", "nls", R["cond_nls"]),
            ("(c) EKF + RTS, one pass", "ekf", np.nan),
            ("(c) iterated EKF smoother", "ieks", np.nan))
    summ = {}
    for lab, k, cn in rows:
        if k == "sat":
            ke = summarise([o["ke"] for o in out[k]], ke0)
            c = np.array([o["c"] for o in out[k]])
            print(f"     {lab:30} {cn:9.2e} {'n/a: only L0 - th':>36} {ke['bias']:+10.3f} "
                  f"{ke['rmse']:8.3f}   [L0 - th: truth {L00-th0:.1f}, median {np.median(c):.1f}]")
            summ[k] = dict(ke=ke)
            continue
        if k in ("nls", "ekf", "ieks"):
            thv = np.exp([o[0] for o in out[k]]); kev = np.exp([o[1] for o in out[k]])
        else:
            thv = [o["th"] for o in out[k]]; kev = [o["ke"] for o in out[k]]
        a, c = summarise(thv, th0), summarise(kev, ke0)
        summ[k] = dict(th=a, ke=c)
        print(f"     {lab:30} {cn:9.2e} {a['bias']:+10.3f} {a['rmse']:8.3f} {a['robust']:8.3f} "
              f"{a['neg']:6.1%} {c['bias']:+10.3f} {c['rmse']:8.3f}")
    return summ


def main():
    print("=" * 92)
    print("  TASK 6. DYNAMIC ESTIMATORS")
    print("=" * 92)
    RA = run_bench("A")
    SA = report("A", RA)
    RB = run_bench("B")
    SB = report("B", RB)

    # errors-in-variables isolated: bench B, (a) with the regressors noiseless but
    # psi noisy, against fully noisy
    b = RB["b"]
    rng = np.random.default_rng(7)
    att = []
    for _ in range(300):
        y = b["m"] * np.exp(rng.normal(0, b["cv"], b["m"].size))
        tau = np.diff(b["t"]); psi = np.diff(y) / tau - u_bar(b["ex"], b["t"])
        for mm, tag in ((b["m"][:-1], "clean m"), (y[:-1], "noisy m")):
            X = np.stack([np.ones_like(mm), mm, -mm ** 2, -psi], axis=1)
            pp, *_ = np.linalg.lstsq(X, psi * mm, rcond=None)
            att.append((tag, pp))
    lt = RB["lt"]
    for tag in ("clean m", "noisy m"):
        P = np.array([p for t_, p in att if t_ == tag])
        print(f"     bench B, (a) with {tag}: median coefficient / truth = "
              + ", ".join(f"{np.median(P[:, j])/lt[j]:.3f}" for j in range(4)))
    print("     psi is noisy in both rows and sits on both sides of the regression;")
    print("     The two rows agree, so the bias is from psi used as a regressor, not from noise in m.")

    # ---- figure: bench B, distributions of ln(theta_hat/theta)
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 3.8))
    for a, (R, nm) in zip(ax, ((RA, "TB-A: 5%, washout, hourly, CV 5%"),
                               (RB, "TB-B bench: 20%, ×6 step, 15 min, CV 1%"))):
        th0 = R["b"]["th"]
        keys = (("fd", "(a) FD"), ("fdm", "(a) FD mid"), ("gen", "(b′) integral"),
                ("cum", "(b″) cumulative"), ("nls", "(c) NLS"), ("ekf", "(c) EKF+RTS"),
                ("ieks", "(c) iterated EKS"))
        data = []
        for k, lab in keys:
            v = (np.exp([o[0] for o in R["out"][k]]) if k in ("nls", "ekf", "ieks")
                 else np.array([o["th"] for o in R["out"][k]]))
            v = v[v > 0]
            data.append(np.clip(np.log(v / th0), -8, 8))
        bp = a.boxplot(data, orientation="horizontal", widths=0.55, patch_artist=True, showfliers=False,
                       medianprops=dict(color=INK, lw=1.5))
        for patch in bp["boxes"]:
            patch.set_facecolor(C1); patch.set_alpha(0.35); patch.set_edgecolor(C1)
        a.set_yticks(range(1, len(keys) + 1), [lab for _, lab in keys])
        a.axvline(0, color=INK, lw=0.8)
        a.set_xlabel(r"ln($\hat\theta$/θ), positive estimates only, clipped at ±8")
        a.set_title(nm)
        finish(a)
    save(fig, "t6_estimators")


if __name__ == "__main__":
    main()
