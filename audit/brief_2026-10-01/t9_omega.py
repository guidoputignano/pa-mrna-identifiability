"""
Follow-up item 2. N = 24 and the between-patient spread omega of ln k_e.

WHAT THE "24 RECORDS" ARE. design.run.SEEDS = range(24) are 24 draws of the
day-to-day intake sequence of ONE simulated patient, every one with k_e = 5.
The repository holds no patient data (README: "Parameter values are
illustrative, not fitted"). The spread of k_e across the 24 records is
therefore exactly zero by construction, and any omega "estimated" from them
measures nothing about patients. This script says so in its output and does
not use them as patients anywhere.

WHAT IS COMPUTED
  A  Task 2's curves read at N = 24 (t2_cohort.cohort), for the disease cohort
     (all at 5%) and q14's mixed cohort (5/8/10/15/20% cycled), at omega in
     {0, 0.1, 0.2, 0.33, 0.5, 0.8}. Then a Monte Carlo joint MAP at N = 24 for
     omega = 0 and 0.33, with per-patient random effects eta_i ~ N(0, omega^2)
     (t2_cohort.joint_map extended here; t2_cohort.py is not edited).
  B  The approximation "plateau ~ sqrt(0.39^2 + ((1-s)/s omega)^2)" against the
     large-N values, and the effective amplification the posterior implies.
  C  Could 24 records identify omega? Each patient's DATA-ONLY information on
     their own ln k_e_i (theta_i, b_i profiled; K_m shared, or per patient, or
     known) as a Gaussian summary z_i with covariance V; marginally
     z ~ N(mu 1, V + omega^2 I). (i) Fisher information for omega^2; (ii) many
     simulated cohorts of 24, each with its 95% profile-likelihood interval
     for omega. Designs: 24 at 5% monthly, q14's mix monthly, both every 2 d,
     and hypothetical per-patient k_e readings. Then a check of the Gaussian
     summary itself against the per-patient profile likelihood.
  D  The repository's only statements about omega, quoted from the files at
     run time.

Run from this directory:  python3 -W error t9_omega.py > out/t9_omega.txt
"""
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm

import core as C
import t2_cohort as T2

T_START = time.time()
POP = C.KE_PRIOR
OMEGAS_A = (0.0, 0.1, 0.2, 0.33, 0.5, 0.8)
OMEGAS_C = (0.0, 0.2, 0.33, 0.5)
N24 = 24
CHI95 = 3.841458820694124          # chi2_1 0.95 quantile
CHI90 = 2.705543454095404          # chi2_1 0.90: the 5% boundary test for omega = 0
L_REC = C.N4.anchor(C.TH0, C.KM0, C.KE, "A")[0]   # 322.73, the recorded load
ROOT = C.ROOT


def banner(s):
    print("\n" + "=" * 96)
    print("  " + s)
    print("=" * 96)


# ----------------------------------------------------------------- records
def rec_ke(theta, ke, seed, times=C.MONTHLY, load_cv=0.10, D=C.D0):
    """
    core.record with a patient-specific k_e. ASSUMPTION: the load is the
    patient's diet, so it stays at convention A's 322.73 for every patient;
    the marker m_i solves the balance at the patient's own k_e. (core.record's
    `ke` argument would instead move the anchor load with k_e.)
    """
    times = np.asarray(times, float)
    L = L_REC
    m_bar = C.N4.level(theta, C.KM0, ke, L)
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        loads = np.full(times.size, L)
    inc = C.N4.delta_at(times, D) if D else np.zeros(times.size)
    m = C.N4.level(theta, C.KM0, ke, loads, inc)
    return dict(theta=theta, Km=C.KM0, ke=ke, L=L, m_bar=m_bar, loads=loads,
                inc=inc, m=m, times=times)


# ============================================================ A. N = 24
def joint_map_re(thetas, recs, ys, sgs, rng, omega):
    """
    t2_cohort.joint_map with random effects. Shared (ln Km, ln ke_pop);
    per patient (ln th_i, b_i, eta_i), ln ke_i = ln ke_pop + eta_i.
    Priors: ln ke_pop ~ N(c, sp^2) with sp^2 = 0.8^2 - omega^2 and the centre
    c drawn about the truth (t2's pseudo-observation reading); eta_i ~
    N(0, omega^2), a genuine random effect; b_i ~ N(b_c,i, (5% L)^2).
    omega = 0 drops eta and is then t2's joint_map. Returns
    (ln th_hat, ln ke_pop_hat, ln ke_i_hat).
    """
    n = len(recs)
    re = omega > 0
    sp = np.sqrt(POP ** 2 - omega ** 2)
    ke_c = np.log(C.KE) + rng.normal(0, sp)
    sb = np.array([C.LOAD_FRAC * r["L"] for r in recs])
    b_c = rng.normal(0.0, sb)
    nd = recs[0]["m"].size
    npar = 2 + 2 * n + (n if re else 0)

    def unpack(x):
        eta = x[2 + 2 * n:] if re else np.zeros(n)
        return x[0], x[1], x[2:2 + n], x[2 + n:2 + 2 * n], eta

    def res(x):
        lk, lke, lth, b, eta = unpack(x)
        K = np.exp(lk)
        out = []
        for i, r in enumerate(recs):
            m, inc = r["m"], r["inc"]
            pred = (np.exp(lth[i]) + inc) * m / (K + m) + np.exp(lke + eta[i]) * m + b[i]
            out.append((ys[i] - pred) / sgs[i])
        out.append([(lke - ke_c) / sp])
        out.append((b - b_c) / sb)
        if re:
            out.append(eta / omega)
        return np.concatenate(out)

    def jac(x):
        lk, lke, lth, b, eta = unpack(x)
        K = np.exp(lk)
        J = np.zeros((n * nd + 1 + n + (n if re else 0), npar))
        for i, r in enumerate(recs):
            m, inc = r["m"], r["inc"]
            th, ke = np.exp(lth[i]), np.exp(lke + eta[i])
            rows = slice(i * nd, (i + 1) * nd)
            J[rows, 0] = (th + inc) * K * m / (K + m) ** 2 / sgs[i]
            J[rows, 1] = -ke * m / sgs[i]
            J[rows, 2 + i] = -th * m / (K + m) / sgs[i]
            J[rows, 2 + n + i] = -1.0 / sgs[i]
            if re:
                J[rows, 2 + 2 * n + i] = -ke * m / sgs[i]
        J[n * nd, 1] = 1.0 / sp
        J[n * nd + 1:n * nd + 1 + n, 2 + n:2 + 2 * n] = np.diag(1.0 / sb)
        if re:
            J[n * nd + 1 + n:, 2 + 2 * n:] = np.eye(n) / omega
        return J

    lo = np.r_[np.log(5.0), np.log(1e-3), np.full(n, np.log(1e-3)), np.full(n, -np.inf)]
    hi = np.r_[np.log(5000.0), np.log(500.0), np.full(n, np.log(1e5)), np.full(n, np.inf)]
    if re:
        lo = np.r_[lo, np.full(n, -8.0)]
        hi = np.r_[hi, np.full(n, 8.0)]
    best = None
    for th0, ke0 in ((np.log(thetas), ke_c), (np.full(n, np.log(600.0)), ke_c - 0.5),
                     (np.full(n, np.log(30.0)), ke_c + 0.5)):
        x0 = np.r_[np.log(C.KM0), np.clip(ke0, lo[1] + 1e-6, hi[1] - 1e-6),
                   np.clip(th0, lo[2] + 1e-6, hi[2] - 1e-6), np.zeros(n)]
        if re:
            x0 = np.r_[x0, np.zeros(n)]
        s = least_squares(res, x0, jac=jac, bounds=(lo, hi), method="trf",
                          x_scale="jac", xtol=1e-10, ftol=1e-10, max_nfev=3000)
        if best is None or s.cost < best.cost:
            best = s
    _, lke, lth, _, eta = unpack(best.x)
    return lth, lke, lke + eta


def mc_re(omega, thetas, reps, seed0):
    """Replicated joint MAP; target is patient 0. Truth: ln ke_i = ln 5 + eta_i."""
    rng = np.random.default_rng(seed0)
    n = len(thetas)
    e_th, e_pop, e_own = [], [], []
    for k in range(reps):
        eta = rng.normal(0.0, omega, n) if omega > 0 else np.zeros(n)
        recs, ys, sgs = [], [], []
        for i in range(n):
            sd = 50_000 + 1000 * k + i
            r = rec_ke(thetas[i], C.KE * np.exp(eta[i]), seed=sd)
            y, sg = C.noisy(r, sd)
            recs.append(r); ys.append(y); sgs.append(sg)
        lth, lke, lke_i = joint_map_re(list(thetas), recs, ys, sgs, rng, omega)
        e_th.append(lth[0] - np.log(thetas[0]))
        e_pop.append(lke - np.log(C.KE))
        e_own.append(lke_i[0] - (np.log(C.KE) + eta[0]))
    e_th, e_pop, e_own = map(np.array, (e_th, e_pop, e_own))
    hw = lambda e: 0.5 * (np.percentile(e, 84) - np.percentile(e, 16))
    rmse = lambda e: float(np.sqrt(np.mean(e ** 2)))
    floor = np.log(1e-3) - np.log(thetas[0]) + 1e-3
    return dict(rmse=rmse(e_th), hw=float(hw(e_th)), bias=float(np.median(e_th)),
                off=float(np.mean(np.abs(e_th) > np.log(4.5))),
                at_floor=float(np.mean(e_th <= floor)),
                rmse_trim=rmse(e_th[e_th > floor]),
                pop_rmse=rmse(e_pop), own_rmse=rmse(e_own), own_hw=float(hw(e_own)))


def part_a():
    banner("A. TASK 2 READ AT N = 24 PATIENTS (FIM, t2_cohort.cohort), THEN A JOINT-MAP CHECK")
    print("  Each cell: target patient's sd(ln theta) | sd(ln k_e,pop) | sd of the target's OWN ln k_e.")
    print("  The target is patient 0 (a 5% patient, intake seed 0). The marginal prior on any one")
    print("  patient's ln k_e stays 0.8: sigma_pop^2 = 0.8^2 - omega^2 (Task 2's convention).")
    out = {}
    for name, fn in (("disease cohort, 24 at 5%", T2.disease_cohort),
                     ("q14 mixed cohort, 5/8/10/15/20% cycled", T2.mixed_cohort)):
        th, sd = fn(N24)
        print(f"\n  {name}, monthly (25 draws each)")
        print(f"     {'omega':>6} {'sd ln theta_0':>14} {'sd ln ke_pop':>13} {'sd own ln ke_0':>15}")
        for w in OMEGAS_A:
            a = T2.cohort(th, sd, w)
            out[(name[:7], w)] = a
            print(f"     {w:6.2f} {a[0]:14.3f} {a[1]:13.3f} {a[3]:15.3f}")
    one = T2.cohort(*T2.disease_cohort(1), 0.0)[0]
    print(f"\n  For reference, one patient alone: {one:.3f} (Task 2, N = 1).")

    print("\n  Monte Carlo joint MAP, 24 patients at 5% monthly, random effects eta_i ~ N(0, omega^2).")
    print("  Truth per replicate: eta_i drawn, ln k_e,i = ln 5 + eta_i, load 322.73 for everyone,")
    print("  fresh intake and noise seeds. Errors are in ln units; half-width = (P84 - P16)/2.")
    reps = 60
    print(f"  REPLICATES: {reps} per omega (the brief asked for about 40-60).")
    print(f"     {'omega':>6} {'FIM sd':>7} {'MC RMSE':>8} {'MC half-w':>10} {'median bias':>12} "
          f"{'>4.5x off':>10} {'at floor':>9} {'RMSE excl':>10} {'ke_pop RMSE':>12} {'own ke RMSE':>12} {'FIM own ke':>11}")
    th, _ = T2.disease_cohort(N24)
    mcr = {}
    for w, s0 in ((0.0, 9100), (0.33, 9133)):
        t0 = time.time()
        r = mc_re(w, th, reps, s0)
        mcr[w] = r
        f = out[("disease", w)]
        print(f"     {w:6.2f} {f[0]:7.3f} {r['rmse']:8.3f} {r['hw']:10.3f} {r['bias']:+12.3f} "
              f"{r['off']:10.1%} {r['at_floor']:9.1%} {r['rmse_trim']:10.3f} {r['pop_rmse']:12.3f} {r['own_rmse']:12.3f} {f[3]:11.3f}"
              f"   ({time.time() - t0:.0f} s)")
    print("  'at floor' = theta_hat_0 at the fit's lower bound 1e-3 (the one-sided likelihood of Task 1);")
    print("  'RMSE excl' = RMSE over the other fits.")
    print(f"  Reading: at omega = 0.33 the target's own ln k_e is known to {out[('disease', 0.33)][3]:.3f} "
          f"(FIM) / {mcr[0.33]['own_rmse']:.3f} (MC RMSE),")
    print(f"  against omega = 0.33 itself; its sd(ln theta) is {out[('disease', 0.33)][0]:.3f} (FIM), "
          f"against the 0.20 target.")
    return out, mcr


# ============================================================ B. plateau
def limit_cov(r, omega):
    """N -> infinity: k_e,pop and K_m known exactly, own ln k_e ~ N(pop, omega^2)."""
    F = C.fim(r, ke_prior=omega if omega > 0 else None)
    k = [C.I_TH, C.I_B, C.I_KE]
    return np.linalg.inv(F[np.ix_(k, k)])


def part_b():
    banner("B. THE PLATEAU APPROXIMATION  sqrt(0.39^2 + ((1-s)/s * omega)^2)  AGAINST COMPUTED VALUES")
    s = C.share(C.TH_5)
    amp = (1 - s) / s
    print(f"  s at 5% (untreated balance) = {s:.4f}, (1-s)/s = {amp:.3f}")
    recs = [C.record(C.TH_5, seed=sd) for sd in C.SEEDS]
    keep3 = [C.I_TH, C.I_KM, C.I_B]
    fl_free = np.array([np.sqrt(np.linalg.inv(C.fim(r)[np.ix_(keep3, keep3)])[0, 0]) for r in recs])
    keep2 = [C.I_TH, C.I_B]
    fl_known = np.array([np.sqrt(np.linalg.inv(C.fim(r)[np.ix_(keep2, keep2)])[0, 0]) for r in recs])
    print(f"  k_e-exact floor, K_m free (the 0.39): median {np.median(fl_free):.4f} "
          f"[{fl_free.min():.4f}, {fl_free.max():.4f}]; seed 0 {fl_free[0]:.4f}")
    print(f"  k_e-exact floor, K_m known (the N -> inf floor): median {np.median(fl_known):.4f} "
          f"[{fl_known.min():.4f}, {fl_known.max():.4f}]; seed 0 {fl_known[0]:.4f}")
    th1k, sd1k = T2.disease_cohort(1024)
    th4k, sd4k = T2.disease_cohort(4096)
    print("\n  Computed plateau (target seed 0 for the cohort columns; the limit column is the")
    print("  median over the 24 intake seeds, seed 0 in brackets):")
    print(f"     {'omega':>6} {'user approx':>12} {'N=1024':>8} {'N=4096':>8} {'N->inf':>16} "
          f"{'approx/limit':>13} {'A_eff':>7} {'|beta|':>7} {'own sd/omega':>13}")
    rows = {}
    for w in (0.2, 0.33, 0.5):
        approx = np.hypot(0.39, amp * w)
        c1 = T2.cohort(th1k, sd1k, w)[0]
        c4 = T2.cohort(th4k, sd4k, w)[0]
        lim, beta, own, aeff = [], [], [], []
        for r, fk in zip(recs, fl_known):
            P = limit_cov(r, w)
            lim.append(np.sqrt(P[0, 0]))
            beta.append(P[0, 2] / P[2, 2])
            own.append(np.sqrt(P[2, 2]))
            aeff.append(np.sqrt(P[0, 0] - fk ** 2) / w)
            # exact split: Var = Var(.|ke) + beta^2 Var(ke)
            assert abs(P[0, 0] - (fk ** 2 + (P[0, 2] / P[2, 2]) ** 2 * P[2, 2])) < 1e-9 * P[0, 0]
        lim, beta, own, aeff = map(np.array, (lim, beta, own, aeff))
        rows[w] = (approx, c1, c4, np.median(lim), np.median(aeff))
        print(f"     {w:6.2f} {approx:12.3f} {c1:8.3f} {c4:8.3f} {np.median(lim):8.3f} ({lim[0]:5.3f}) "
              f"{approx / np.median(lim):13.3f} {np.median(aeff):7.3f} {np.median(np.abs(beta)):7.3f} "
              f"{np.median(own) / w:13.3f}")
    print("  A_eff = sqrt(plateau^2 - floor_Km-known^2)/omega, the amplification the posterior actually")
    print("  applies; it factors exactly as |beta| * (own sd/omega), where beta is the posterior")
    print("  regression of ln theta on the patient's own ln k_e and own sd/omega < 1 is the shrinkage")
    print("  of the patient's own k_e by their own record (identity asserted per seed).")
    print("  The user's form overstates the plateau because it uses (1-s)/s instead of |beta|, and")
    print("  omega instead of the patient's posterior sd of ln k_e; see the two right-hand columns.")
    return rows, amp


# ============================================================ C. omega
def data_F(theta, seed, times):
    """One patient's DATA-ONLY information: no k_e prior; the 5% L load-bias prior kept."""
    return C.fim(C.record(theta, seed=seed, times=times), ke_prior=None)


def summary_V(thetas, seeds, times, km="shared", extra_ke_sd=None):
    """
    Covariance V of the per-patient data-only estimates of ln k_e_i, as
    fixed effects (no k_e prior of any kind).
      km = 'shared'     one ln K_m for the cohort (no prior): V is full
      km = 'perpatient' each patient profiles their own ln K_m: V diagonal
      km = 'known'      K_m known exactly: V diagonal
    extra_ke_sd: a hypothetical direct per-patient ln k_e reading of this sd.
    """
    Fs = [data_F(th, s, times) for th, s in zip(thetas, seeds)]
    if extra_ke_sd:
        for F in Fs:
            F[C.I_KE, C.I_KE] += extra_ke_sd ** -2
    n = len(Fs)
    if km == "perpatient":
        return np.diag([np.linalg.inv(F)[C.I_KE, C.I_KE] for F in Fs])
    loc = [C.I_TH, C.I_B, C.I_KE]
    if km == "known":
        return np.diag([np.linalg.inv(F[np.ix_(loc, loc)])[2, 2] for F in Fs])
    G, W, d = 0.0, [], []
    for F in Fs:
        A = F[np.ix_(loc, loc)]
        B = F[loc, C.I_KM]
        Ai = np.linalg.inv(A)
        G += F[C.I_KM, C.I_KM] - B @ Ai @ B
        W.append((Ai @ B)[2])
        d.append(Ai[2, 2])
    W = np.array(W)
    return np.diag(d) + np.outer(W, W) / G


def fisher_omega2(V, w2):
    """Expected information for omega^2 in z ~ N(mu 1, V + omega^2 I), ML and REML."""
    S = V + w2 * np.eye(V.shape[0])
    Si = np.linalg.inv(S)
    one = np.ones(V.shape[0])
    P = Si - np.outer(Si @ one, Si @ one) / (one @ Si @ one)
    return 0.5 * np.trace(Si @ Si), 0.5 * np.trace(P @ P)


W_GRID = np.r_[0.0, np.geomspace(1e-3, 30.0, 1500)]


def profile_dev(V, Z):
    """
    -2 log L(omega) profiled over mu (GLS), for each row of Z (replicates x N),
    on W_GRID. Returns (ML deviance, REML deviance), each (n_omega, n_rep).
    """
    lam, Q = np.linalg.eigh(V)
    zt = Z @ Q
    ut = Q.T @ np.ones(V.shape[0])
    Wt = 1.0 / (lam[None, :] + W_GRID[:, None] ** 2)
    logdet = np.log(lam[None, :] + W_GRID[:, None] ** 2).sum(1)
    a1 = Wt @ (zt ** 2).T
    a2 = Wt @ (zt * ut[None, :]).T
    a3 = Wt @ ut ** 2
    dml = logdet[:, None] + a1 - a2 ** 2 / a3[:, None]
    return dml, dml + np.log(a3)[:, None]


def interval(dev, crit=CHI95):
    """95% profile-likelihood interval [lo, hi] on W_GRID, column by column."""
    rel = dev - dev.min(0)
    inside = rel <= crit
    n = dev.shape[1]
    lo, hi = np.zeros(n), np.zeros(n)
    cens = np.zeros(n, bool)
    for j in range(n):
        idx = np.flatnonzero(inside[:, j])
        i0, i1 = idx[0], idx[-1]
        if i0 == 0:
            lo[j] = 0.0
        else:
            x0, x1, y0, y1 = W_GRID[i0 - 1], W_GRID[i0], rel[i0 - 1, j], rel[i0, j]
            lo[j] = x0 + (crit - y0) * (x1 - x0) / (y1 - y0)
        if i1 == W_GRID.size - 1:
            hi[j], cens[j] = W_GRID[-1], True
        else:
            x0, x1, y0, y1 = W_GRID[i1], W_GRID[i1 + 1], rel[i1, j], rel[i1 + 1, j]
            hi[j] = x0 + (crit - y0) * (x1 - x0) / (y1 - y0)
    what = W_GRID[np.argmin(dev, 0)]
    lr0 = rel[0]
    return lo, hi, cens, what, lr0


def simulate(V, omega, reps, seed):
    rng = np.random.default_rng(seed)
    n = V.shape[0]
    Lc = np.linalg.cholesky(V)
    Z = np.log(C.KE) + omega * rng.standard_normal((reps, n)) + rng.standard_normal((reps, n)) @ Lc.T
    res = []
    for k in range(0, reps, 500):
        dml, dre = profile_dev(V, Z[k:k + 500])
        res.append((interval(dml), interval(dre)))
    cat = lambda which, i: np.concatenate([r[which][i] for r in res])
    ml = [cat(0, i) for i in range(5)]
    reml = [cat(1, i) for i in range(5)]
    return ml, reml


def n_required(v, w, power=0.80, alpha=0.05):
    """
    Asymptotic cohort size for the boundary LR/Wald test of omega^2 = 0 at
    one-sided alpha with the given power, from the per-patient expected
    information 1/(2 (v_i + omega^2)^2) cycled over the composition.
    """
    v = np.asarray(v, float)
    h0 = np.mean(1.0 / v ** 2) ** -0.5
    h1 = np.mean(1.0 / (v + w ** 2) ** 2) ** -0.5
    za, zb = norm.ppf(1 - alpha), norm.ppf(power)
    return 2.0 * ((za * h0 + zb * h1) / w ** 2) ** 2


def v_required(n, w, power=0.80, alpha=0.05):
    """Per-patient data-only variance of ln k_e that makes n patients enough (equal v)."""
    za, zb = norm.ppf(1 - alpha), norm.ppf(power)
    c = np.sqrt(2.0 / n)
    return (w ** 2 * (1 - zb * c)) / (c * (za + zb))


def design_list():
    th5, sd5 = T2.disease_cohort(N24)
    thm, sdm = T2.mixed_cohort(N24)
    return [
        ("24 at 5%, monthly (25 draws)", th5, sd5, C.MONTHLY, None),
        ("q14 mix 5-20%, monthly", thm, sdm, C.MONTHLY, None),
        ("24 at 5%, every 2 d (365 draws)", th5, sd5, C.DENSE, None),
        ("q14 mix 5-20%, every 2 d", thm, sdm, C.DENSE, None),
        ("HYPOTHETICAL 24 at 5% monthly + own ln k_e read to sd 0.2", th5, sd5, C.MONTHLY, 0.2),
        ("HYPOTHETICAL 24 at 5% monthly + own ln k_e read to sd 0.1", th5, sd5, C.MONTHLY, 0.1),
    ]


def part_c():
    banner("C. COULD 24 RECORDS IDENTIFY OMEGA?  GAUSSIAN SUMMARY OF EACH PATIENT'S DATA-ONLY ln k_e")
    print("  NOTE. The repository's 24 records are NOT patients: they share one k_e = 5, so their")
    print("  between-record spread of k_e is 0 by construction and no estimate from them is evidence")
    print("  about patients. Everything below is about a HYPOTHETICAL cohort of 24 real patients")
    print("  simulated on the repository's footing (each patient's own intake seed 0..23).")
    print("\n  Model: z_i = data-only estimate of ln k_e_i (theta_i, b_i profiled, b prior 5% L kept,")
    print("  no k_e prior); z ~ N(mu 1, V + omega^2 I), mu profiled with no prior. V is evaluated at")
    print("  the typical patient (eta = 0) and held fixed across replicates.")

    designs = design_list()
    Vs = {}
    print("\n  (i) per-patient data-only sd of ln k_e_i, and the expected (Fisher) information for omega^2")
    print(f"     {'design':<60} {'sd z_i: Km shared':>18} {'Km per pt':>10} {'Km known':>9}")
    for lab, th, sd, tt, ex in designs:
        V = summary_V(th, sd, tt, "shared", ex)
        Vp = summary_V(th, sd, tt, "perpatient", ex)
        Vk = summary_V(th, sd, tt, "known", ex)
        Vs[lab] = (V, Vp, Vk)
        dsd = lambda M: np.sqrt(np.diag(M))
        print(f"     {lab:<60} {np.median(dsd(V)):8.3f} [{dsd(V).min():.2f},{dsd(V).max():.2f}]"
              f" {np.median(dsd(Vp)):10.3f} {np.median(dsd(Vk)):9.3f}")
    print("     (median over the 24 patients, [min, max] for the shared-K_m column)")
    print("     Shared K_m makes the z_i correlated; the largest correlation in V is printed below.")
    print(f"\n     {'design':<60} {'max corr':>8} {'sd(w^2)|w=0':>12} {'REML':>7}"
          + "".join(f" {'z at w=' + str(w):>10}" for w in OMEGAS_C[1:]))
    fish = {}
    for lab, *_ in designs:
        V = Vs[lab][0]
        dV = np.sqrt(np.diag(V))
        Cr = V / np.outer(dV, dV)
        mc = np.max(np.abs(Cr - np.eye(V.shape[0])))
        i0ml, i0re = fisher_omega2(V, 0.0)
        sd0 = i0ml ** -0.5
        zs = [w ** 2 / sd0 for w in OMEGAS_C[1:]]
        fish[lab] = sd0
        print(f"     {lab:<60} {mc:8.3f} {sd0:12.4f} {i0re ** -0.5:7.4f}"
              + "".join(f" {z:10.2f}" for z in zs))
    print("     sd(w^2)|w=0 = 1/sqrt(I(omega^2)) at omega = 0 (ML; REML beside it). 'z at w' is")
    print("     omega^2 / sd(omega^2 | omega = 0): the expected Wald signal of a true omega against the")
    print("     null's sd. omega = 0.33 is omega^2 = 0.109; omega = 0.5 is 0.25. I(omega) itself is")
    print("     4 omega^2 I(omega^2), which is ZERO at omega = 0: omega is a boundary parameter there.")
    print("\n     Fisher sd of omega at the true omega (delta method, ML, shared K_m):")
    print(f"     {'design':<60}" + "".join(f" {'w=' + str(w):>8}" for w in OMEGAS_C[1:]))
    for lab, *_ in designs:
        V = Vs[lab][0]
        cells = []
        for w in OMEGAS_C[1:]:
            iml, _ = fisher_omega2(V, w ** 2)
            cells.append(iml ** -0.5 / (2 * w))
        print(f"     {lab:<60}" + "".join(f" {c:8.3f}" for c in cells))

    # ------------------------------------------------------- (ii) simulation
    reps = 4000
    print(f"\n  (ii) {reps} simulated cohorts of 24 per (design, true omega); each cohort's 95% profile-")
    print("  likelihood interval for omega = {omega : D(omega) - D_min <= 3.84} (ML, mu profiled).")
    print("  'excl 0' = interval excludes omega = 0; 'test' = boundary-corrected 5% LR test of omega = 0")
    print("  (D(0) - D_min > 2.71, the 50:50 chi2_0/chi2_1 mixture); 'both' = interval holds 0 AND 0.33;")
    print("  'excl .33' = interval excludes 0.33. Upper limit: median [P10, P90]; censored = upper")
    print(f"  limit beyond the grid ({W_GRID[-1]:.0f}).")
    hdr = (f"     {'design':<60} {'w':>5} {'median w_hat':>12} {'upper: med [P10,P90]':>23} "
           f"{'med lower':>9} {'excl 0':>7} {'test':>6} {'both':>6} {'excl .33':>9} {'cens':>5}")
    sims = {}
    for km_i, km_lab in ((0, "K_m shared (primary)"), (1, "K_m per patient (diagonal V)")):
        print(f"\n   -- {km_lab}")
        print(hdr)
        for d_i, (lab, *_) in enumerate(designs):
            if km_i == 1 and "HYPOTHETICAL" in lab:
                continue
            V = Vs[lab][km_i]
            for w_i, w in enumerate(OMEGAS_C):
                ml, reml = simulate(V, w, reps, seed=700_000 + 1000 * d_i + 100 * km_i + w_i)
                lo, hi, cens, what, lr0 = ml
                both = (lo <= 0.0) & (hi >= 0.33)
                excl33 = (lo > 0.33) | (hi < 0.33)
                sims[(km_i, lab, w)] = dict(hi=hi, lo=lo, ex0=np.mean(lo > 0), test=np.mean(lr0 > CHI90),
                                            both=np.mean(both), reml=reml)
                print(f"     {lab:<60} {w:5.2f} {np.median(what):12.3f} "
                      f"{np.median(hi):7.3f} [{np.percentile(hi, 10):6.3f},{np.percentile(hi, 90):6.3f}] "
                      f"{np.median(lo):9.3f} {np.mean(lo > 0):7.1%} {np.mean(lr0 > CHI90):6.1%} "
                      f"{np.mean(both):6.1%} {np.mean(excl33):9.1%} {np.mean(cens):5.1%}")
    print("\n   REML check (shared K_m), same cohorts: excl 0 | median upper limit")
    for d_i, (lab, *_) in enumerate(designs):
        cells = []
        for w in OMEGAS_C:
            lo, hi = sims[(0, lab, w)]["reml"][0], sims[(0, lab, w)]["reml"][1]
            cells.append(f"w={w}: {np.mean(lo > 0):5.1%} | {np.median(hi):6.3f}")
        print(f"     {lab:<60} " + "   ".join(cells))

    # ------------------------------------------------ what would identify it
    print("\n  (iii) What would identify omega = 0.33 (boundary test, one-sided 5%, 80% power).")
    print("  ASYMPTOTIC formula from the per-patient expected information 1/(2 (v_i + omega^2)^2),")
    print("  diagonal V (K_m correlation ignored); NOT checked by simulation at the large N it returns.")
    print(f"     {'design':<60} {'N for w=0.33':>13} {'N for w=0.5':>12} {'N for w=0.2':>12}")
    for lab, *_ in designs:
        v = np.diag(Vs[lab][0])
        print(f"     {lab:<60} {n_required(v, 0.33):13.0f} {n_required(v, 0.5):12.0f} "
              f"{n_required(v, 0.2):12.0f}")
    vr = v_required(N24, 0.33)
    print(f"     With N = 24 fixed, the per-patient data-only sd of ln k_e needed for omega = 0.33: "
          f"{np.sqrt(vr):.3f} (equal v_i, same formula).")
    Vr = np.eye(N24) * vr
    for w in (0.0, 0.33):
        ml, _ = simulate(Vr, w, reps, seed=880_000 + int(100 * w))
        lo, hi, cens, what, lr0 = ml
        print(f"     simulated check at that sd, true omega = {w}: test rejects omega = 0 in "
              f"{np.mean(lr0 > CHI90):.1%}, interval excludes 0 in {np.mean(lo > 0):.1%}, "
              f"median upper limit {np.median(hi):.3f}")

    print("\n  SUMMARY of (ii), K_m shared, Gaussian summary: the typical 95% interval for omega")
    print(f"     {'design':<60} {'true w = 0':>20} {'true w = 0.33':>20}")
    for lab, *_ in designs:
        a, b = sims[(0, lab, 0.0)], sims[(0, lab, 0.33)]
        print(f"     {lab:<60} [{np.median(a['lo']):5.3f}, {np.median(a['hi']):6.3f}]     "
              f"  [{np.median(b['lo']):5.3f}, {np.median(b['hi']):6.3f}]")
    print("     (median lower limit, median upper limit over the simulated cohorts)")
    return Vs, sims, fish


# ----------------------------------------- C (iv): is the Gaussian summary right?
def _res_k(x, lke, r, y, sg, sb):
    lt, lk, b = x
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    m, inc = r["m"], r["inc"]
    pred = (th + inc) * m / (K + m) + ke * m + b
    return np.r_[(y - pred) / sg, b / sb]


def _jac_k(x, lke, r, y, sg, sb):
    lt, lk, b = x
    th, K = np.exp(lt), np.exp(lk)
    m, inc = r["m"], r["inc"]
    J = np.zeros((m.size + 1, 3))
    J[:m.size, 0] = -th * m / (K + m) / sg
    J[:m.size, 1] = (th + inc) * K * m / (K + m) ** 2 / sg
    J[:m.size, 2] = -1.0 / sg
    J[m.size, 2] = 1.0 / sb
    return J


LO3 = np.array([np.log(1e-6), np.log(1e-2), -np.inf])
HI3 = np.array([np.log(1e6), np.log(1e5), np.inf])
KE_OFF = np.array([-8, -6, -4, -3, -2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.1, 0.2, 0.3,
                   0.4, 0.5, 0.75, 1.0, 1.5, 2.0])


def ke_profile(r, y, sg, grid):
    """
    Data-only profile deviance in ln k_e for ONE patient: (ln theta, ln K_m, b)
    profiled, b prior 5% L, no k_e prior. Multi-start at the first point,
    warm start plus three starts afterwards.
    """
    sb = C.LOAD_FRAC * r["L"]
    out, th, prev = [], [], None
    for j, lke in enumerate(grid):
        if prev is None:
            cands = [np.array([np.log(t), np.log(k), 0.0])
                     for t in (1.0, 150.0, 1000.0) for k in (30.0, 100.0, 300.0)]
        else:
            cands = [prev] + [np.array([np.log(t), np.log(100.0), 0.0]) for t in (1.0, 150.0, 1000.0)]
        best = None
        for x0 in cands:
            s = least_squares(_res_k, x0, jac=_jac_k, args=(lke, r, y, sg, sb), bounds=(LO3, HI3),
                              method="trf", xtol=1e-12, ftol=1e-12, max_nfev=2000)
            if best is None or s.cost < best.cost:
                best = s
        prev = best.x
        out.append(2.0 * best.cost)
        th.append(np.exp(best.x[0]))
    return np.array(out), np.array(th)


def hier_from_profiles(cohorts, grid, mu_grid, w_grid, dx=0.01):
    """
    Profile-likelihood analogue of the Gaussian-summary model: each patient's
    data-only profile likelihood exp(-dev_i/2) in ln k_e stands in for their
    likelihood, integrated against N(mu, omega^2) on a fine x grid (linear
    interpolation in dev; held CONSTANT below the grid floor, an assumption),
    mu profiled on mu_grid. `cohorts` is a list of per-patient deviance lists.
    Returns D(omega) on w_grid, one row per cohort.
    """
    span = 5.0 * max(w_grid)
    x = np.arange(mu_grid[0] - span, mu_grid[-1] + span, dx)
    n_c, n_p = len(cohorts), len(cohorts[0])
    flat = [d for c in cohorts for d in c]
    E = np.array([np.exp(-0.5 * (np.interp(x, grid, d) - d.min())) for d in flat]).T
    E0 = np.array([np.exp(-0.5 * (np.interp(mu_grid, grid, d) - d.min())) for d in flat]).T
    shift = np.array([d.min() for d in flat]).reshape(n_c, n_p).sum(1)
    D = np.zeros((n_c, w_grid.size))
    for j, w in enumerate(w_grid):
        if w == 0:
            Lmu = E0
        else:
            K = norm.pdf((x[None, :] - mu_grid[:, None]) / w) / w * dx
            Lmu = K @ E
        with np.errstate(divide="ignore"):
            ll = np.log(Lmu).reshape(mu_grid.size, n_c, n_p).sum(2)
        D[:, j] = (-2.0 * ll).min(0) + shift
    return D


def interval_1d(rel, wg, crit=CHI95):
    inside = np.flatnonzero(rel <= crit)
    i0, i1 = inside[0], inside[-1]
    if i1 == wg.size - 1:
        hi, cens = wg[-1], True
    else:
        hi = wg[i1] + (crit - rel[i1]) * (wg[i1 + 1] - wg[i1]) / (rel[i1 + 1] - rel[i1])
        cens = False
    return i0 > 0, hi, cens


def shape_table(lab, times):
    """Per-patient data-only profile in ln k_e over the 24 intake seeds, against the Gaussian."""
    grid = np.log(C.KE) + KE_OFF
    show = [-8, -4, -2, -1, -0.5, 0, 0.2, 0.3, 0.5, 1.0]
    cols = [int(np.flatnonzero(np.isclose(KE_OFF, o))[0]) for o in show]
    t0 = time.time()
    devs, vs, ths = [], [], []
    for s in C.SEEDS:
        r = C.record(C.TH_5, seed=s, times=times)
        y, sg = C.noisy(r, s)
        d, th = ke_profile(r, y, sg, grid)
        devs.append(d); ths.append(th)
        vs.append(np.linalg.inv(data_F(C.TH_5, s, times))[C.I_KE, C.I_KE])
    devs, vs, ths = np.array(devs), np.array(vs), np.array(ths)
    rel = devs - devs.min(1, keepdims=True)
    r0 = C.record(C.TH_5, seed=0, times=times)
    y0 = (r0["theta"] + r0["inc"]) * r0["m"] / (r0["Km"] + r0["m"]) + r0["ke"] * r0["m"]
    d0, _ = ke_profile(r0, y0, C.sigmas(r0), grid)
    d0 = d0 - d0.min()
    print(f"\n  {lab}: deviance relative to the minimum, 24 intake seeds")
    print(f"     {'ln k_e offset':>17}" + "".join(f" {o:>7}" for o in show))
    print(f"     {'noiseless, seed 0':>17}" + "".join(f" {d0[c]:7.2f}" for c in cols))
    print(f"     {'noisy: median':>17}" + "".join(f" {np.median(rel[:, c]):7.2f}" for c in cols))
    print(f"     {'noisy: max':>17}" + "".join(f" {rel[:, c].max():7.2f}" for c in cols))
    print(f"     {'Gauss, median v':>17}" + "".join(f" {o ** 2 / np.median(vs):7.2f}" for o in show))
    print(f"     {'theta_hat median':>17}" + "".join(f" {np.median(ths[:, c]):7.1f}" for c in cols))
    inside_floor = np.mean(rel[:, 0] <= CHI95)
    up = []
    for d in rel:
        i = np.flatnonzero(d <= CHI95)[-1]
        up.append(KE_OFF[i] + (CHI95 - d[i]) * (KE_OFF[i + 1] - KE_OFF[i]) / (d[i + 1] - d[i]))
    up = np.array(up)
    print(f"     The 95% region (deviance <= 3.84) reaches the grid floor (k_e = 5 e^-8 = "
          f"{C.KE * np.exp(-8):.4f}) in {inside_floor:.0%} of the 24 records.")
    print(f"     Its upper end: ln k_e offset median {np.median(up):+.3f} [{up.min():+.3f}, {up.max():+.3f}], "
          f"i.e. k_e up to {C.KE * np.exp(np.median(up)):.2f} /d, a factor {np.exp(up.min()):.2f}-"
          f"{np.exp(up.max()):.2f} above the truth;")
    print(f"     the Gaussian summary puts it at +{1.96 * np.sqrt(np.median(vs)):.2f}. Plateau below the "
          f"truth (deviance at offset -8): noiseless {d0[0]:.2f}, noisy median {np.median(rel[:, 0]):.2f}.")
    print(f"     ({time.time() - t0:.0f} s)")


def hier_check(lab, times, reps, seed0):
    grid = np.log(C.KE) + KE_OFF
    mu_grid = np.log(C.KE) + np.arange(-6.0, 2.0001, 0.02)
    w_grid = np.r_[0.0, np.geomspace(0.02, 4.0, 80)]
    print(f"\n  {lab}: {reps} simulated cohorts of 24 per omega (REDUCED count: 24 profiles per cohort)")
    print(f"     {'omega':>6} {'excl 0':>7} {'test':>6} {'upper: med [min,max]':>24} {'both 0,.33':>11} "
          f"{'censored':>9} | {'Gauss: excl 0':>13} {'upper med':>10} {'z_i at floor':>13}")
    t0 = time.time()
    for w in (0.0, 0.33):
        rng = np.random.default_rng(seed0 + int(100 * w))
        cohorts, zs, vs = [], [], []
        for k in range(reps):
            eta = rng.normal(0.0, w, N24) if w > 0 else np.zeros(N24)
            dv, z, v = [], [], []
            for i in range(N24):
                sd = seed0 + 1000 * k + i
                r = rec_ke(C.TH_5, C.KE * np.exp(eta[i]), seed=sd, times=times)
                y, sg = C.noisy(r, sd)
                d, _ = ke_profile(r, y, sg, grid)
                dv.append(d)
                z.append(grid[np.argmin(d)])
                v.append(np.linalg.inv(C.fim(r, ke_prior=None))[C.I_KE, C.I_KE])
            cohorts.append(dv); zs.append(z); vs.append(v)
        D = hier_from_profiles(cohorts, grid, mu_grid, w_grid)
        his, ex0, tst, both, cen = [], [], [], [], []
        for row in D:
            rel_w = row - row.min()
            e0, hi, c = interval_1d(rel_w, w_grid)
            his.append(hi); ex0.append(e0); tst.append(rel_w[0] > CHI90); cen.append(c)
            both.append((not e0) and np.interp(0.33, w_grid, rel_w) <= CHI95)
        his = np.array(his)
        g_ex, g_hi = [], []
        for z, v in zip(zs, vs):
            dml, _ = profile_dev(np.diag(v), np.array(z)[None, :])
            lo, hi, _, _, _ = interval(dml)
            g_ex.append(lo[0] > 0); g_hi.append(hi[0])
        zs = np.array(zs)
        print(f"     {w:6.2f} {np.mean(ex0):7.0%} {np.mean(tst):6.0%} {np.median(his):8.3f} "
              f"[{his.min():5.3f},{his.max():6.3f}] {np.mean(both):11.0%} {np.mean(cen):9.0%} | "
              f"{np.mean(g_ex):13.0%} {np.median(g_hi):10.3f} {np.mean(zs <= grid[0] + 1e-9):13.1%}")
    print(f"     ({time.time() - t0:.0f} s)")


def part_c_check():
    banner("C (iv). CHECK OF THE GAUSSIAN SUMMARY: EACH PATIENT'S DATA-ONLY PROFILE LIKELIHOOD IN ln k_e")
    print("  One patient at 5%, n6-style noisy ordinate (core.noisy); (ln theta, ln K_m, b) profiled,")
    print("  b prior 5% L, no k_e prior. Columns are ln k_e offsets from the truth ln 5. 'Gauss' is the")
    print("  Gaussian summary's offset^2 / v_i (K_m per patient); 'theta_hat' is the profiled theta.")
    shape_table("24 at 5%, monthly", C.MONTHLY)
    shape_table("24 at 5%, every 2 d", C.DENSE)
    print("\n  Where theta_hat is ~1e-6 (the fit's floor) the wall is the theta >= 0 bound of Task 4.")
    print("\n  Profile-based hierarchy: each patient's profile likelihood exp(-dev_i/2) in ln k_e")
    print("  integrated against N(mu, omega^2), mu profiled on ln 5 + [-6, 2], K_m per patient. An")
    print("  approximation: a profile, not an integrated, likelihood per patient, and held constant")
    print("  below the profile grid's floor. Right of the bar, the SAME cohorts summarised the way a")
    print("  two-stage analysis would: z_i = the patient's profile ML of ln k_e (on the grid), v_i from")
    print("  the FIM, diagonal V. 'z_i at floor' = share of patients whose ML sits at the grid floor.")
    hier_check("24 at 5%, monthly", C.MONTHLY, 25, 300_000)
    hier_check("24 at 5%, every 2 d", C.DENSE, 15, 400_000)


# ============================================================ D. sources
def quote(path, a, b):
    lines = (ROOT / path).read_text().splitlines()
    print(f"\n  {path}, lines {a}-{b}:")
    for k in range(a, b + 1):
        print(f"    {k:4d}: {lines[k - 1]}")
    return "\n".join(lines[a - 1:b])


def part_d():
    banner("D. THE REPOSITORY'S ONLY INFORMATION ON OMEGA, QUOTED FROM THE FILES")
    q1 = quote("audit/q14_hierarchical.py", 19, 30)
    assert 'probably a third or' in q1 and "OAT1" in q1 and "eGFR" in q1
    q2 = quote("audit/PHASE10_SENSORS.md", 486, 494)
    assert "MARKER'S CLEARANCE" in q2 and "k_e/alpha" in q2
    print("\n  Two further places restate or border on it (found by grep; not asked for):")
    q3 = quote("audit/PHASE8_SIMULATIONS.md", 721, 723)
    assert "R1" in q3 and "no 3-HP figure" in q3
    q4 = quote("src/common/params.py", 363, 367)
    assert "51" in q4 and "mice" in q4
    print("\n  Reading: every between-patient figure traces to brief R1, which (PHASE8_SIMULATIONS.md")
    print("  721-723) is for 'the nearest comparable disposal rate' in healthy adults and has no 3-HP")
    print("  figure; PHASE10_SENSORS.md 487-494 adds that it is for the marker's clearance, not the")
    print("  model's k_e/alpha, so the model's omega could be larger. The params.py 51% is between")
    print("  MICE, on the 3-HP baseline level, not on k_e. No number in the repository measures omega.")
    readme = (ROOT / "README.md").read_text()
    assert "Parameter values are illustrative, not fitted" in readme
    print("  README.md: 'Parameter values are illustrative, not fitted' (asserted present).")


# ============================================================ figure
def figure(out_a, sims):
    from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.2))

    # (a) sd(ln theta_0) at N = 24 against omega
    a = ax[0]
    for key, col, lab in (("disease", C1, "24 at 5%"), ("q14 mix", C4, "q14 mix 5-20%")):
        y = [out_a[(key, w)][0] for w in OMEGAS_A]
        a.plot(OMEGAS_A, y, "-o", color=col, ms=4, lw=2, label=lab)
    a.axhline(0.20, color=INK, ls="--", lw=1)
    a.text(0.02, 0.23, "target 0.20", fontsize=8, color=INK)
    a.axhline(3.082, color=MUTED, ls=":", lw=1)
    a.text(0.02, 3.13, "one patient alone: 3.08", fontsize=8, color=MUTED)
    a.set_xlabel(r"between-patient sd $\omega$ of $\ln k_e$")
    a.set_ylabel(r"target patient's sd($\ln\hat\theta$)")
    a.set_title("(a) N = 24 cohort vs $\\omega$")
    a.set_ylim(0, 3.3)
    a.legend(fontsize=8, frameon=False, loc="lower right")
    finish(a)

    # (b) the per-patient profile likelihood in ln k_e: one-sided wall
    b = ax[1]
    grid = np.log(C.KE) + KE_OFF
    for times, col, lab in ((C.MONTHLY, C2, "monthly (25 draws)"), (C.DENSE, C3, "every 2 d (365)")):
        rels = []
        for s in C.SEEDS:
            r = C.record(C.TH_5, seed=s, times=times)
            y, sg = C.noisy(r, s)
            d, _ = ke_profile(r, y, sg, grid)
            rels.append(d - d.min())
        med = np.median(np.array(rels), 0)
        b.plot(KE_OFF, med, "-", color=col, lw=2, label=lab + " data")
        v = np.median([np.linalg.inv(data_F(C.TH_5, s, times))[C.I_KE, C.I_KE] for s in C.SEEDS])
        b.plot(KE_OFF, KE_OFF ** 2 / v, "--", color=col, lw=1.2, label=lab + " Gaussian")
    b.axhline(CHI95, color=INK, ls=":", lw=1)
    b.text(-3.8, CHI95 + 0.3, "95% (3.84)", fontsize=8, color=INK)
    b.set_xlim(-4, 1.0)
    b.set_ylim(0, 12)
    b.set_xlabel(r"$\ln k_e$ offset from the truth")
    b.set_ylabel("deviance (median over 24 records)")
    b.set_title("(b) one patient's data on $k_e$: a wall, not a Gaussian")
    b.legend(fontsize=7, frameon=False, loc="upper center")
    finish(b)

    # (c) distribution of the 95% upper limit for omega, true omega = 0 vs 0.33
    c = ax[2]
    designs = [("24 at 5%, monthly (25 draws)", "5% monthly"),
               ("q14 mix 5-20%, monthly", "mix monthly"),
               ("24 at 5%, every 2 d (365 draws)", "5% every 2 d"),
               ("q14 mix 5-20%, every 2 d", "mix every 2 d")]
    pos = np.arange(len(designs))
    for off, w, col in ((-0.17, 0.0, C1), (0.17, 0.33, C4)):
        meds = [np.median(sims[(0, lab, w)]["hi"]) for lab, _ in designs]
        p10 = [np.percentile(sims[(0, lab, w)]["hi"], 10) for lab, _ in designs]
        p90 = [np.percentile(sims[(0, lab, w)]["hi"], 90) for lab, _ in designs]
        c.errorbar(pos + off, meds, yerr=[np.array(meds) - p10, np.array(p90) - meds],
                   fmt="o", color=col, ms=5, lw=1.5, capsize=3, label=f"true $\\omega$ = {w}")
    c.axhline(0.33, color=MUTED, ls=":", lw=1)
    c.text(-0.4, 0.35, r"$\omega$ = 0.33", fontsize=8, color=MUTED)
    c.set_xticks(pos)
    c.set_xticklabels([s for _, s in designs], fontsize=8, rotation=20, ha="right")
    c.set_ylabel(r"95% upper limit for $\omega$ (median, P10-P90)")
    c.set_title("(c) what a cohort of 24 says about $\\omega$")
    c.set_ylim(0, 1.8)
    c.legend(fontsize=8, frameon=False, loc="upper right")
    finish(c)

    save(fig, "t9_omega")


# ============================================================ main
def main():
    banner("FOLLOW-UP ITEM 2.  N = 24 AND OMEGA")
    seeds = list(C.SEEDS)
    kes = {C.record(C.TH_5, seed=s)["ke"] for s in seeds}
    assert seeds == list(range(24)) and kes == {C.KE}
    print("  The repository's '24 records' are design.run.SEEDS = range(24): 24 intake sequences of ONE")
    print(f"  simulated patient. Distinct k_e values among them: {sorted(kes)}. Their between-record")
    print("  spread of k_e is therefore exactly 0 by construction; 'estimating omega from them' can")
    print("  only recover that construction (up to noise) and says nothing about patients. There is no")
    print("  patient data in the repository. Below, '24 patients' always means a SIMULATED cohort.")
    print(f"  Footing: theta_5% = {C.TH_5:.0f}, K_m = {C.KM0:.0f}, k_e = {C.KE:.0f}, L = {L_REC:.2f}, "
          f"m = {C.record(C.TH_5)['m_bar']:.2f}, s = {C.share(C.TH_5):.3f}.")

    out_a, _ = part_a()
    part_b()
    _, sims, _ = part_c()
    part_c_check()
    part_d()
    figure(out_a, sims)
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
