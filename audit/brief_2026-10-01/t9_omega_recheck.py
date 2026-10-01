"""
Adversarial recheck of t9_omega.py (follow-up item 2: N = 24 and omega).

Shares NO code with t9_omega.py or t2_cohort.py (neither is imported). core.py
is imported ONLY to assert that the footing rebuilt here (record generator and
sigma) reproduces core.record / core.sigmas bit for bit; nothing computed below
calls it.

Routes that differ from t9's:
  (a) the N = 24 cohort information matrix in the MARGINAL-PRIOR form: one
      parameter vector (ln K_m, {ln theta_i, b_i, ln k_e,i}), the population
      level integrated out analytically so the prior on the ln k_e vector is
      N(ln 5 1, omega^2 I + (0.8^2 - omega^2) 1 1'), the full 73-square matrix
      inverted directly (t2/t9 use (ln k_e,pop, eta_i) and a Schur complement).
      Sensitivities by complex-step differentiation of the prediction, not by
      hand-written rows.
  (b) per-patient data-only sd of ln k_e from the full fixed-effects cohort
      matrix (K_m shared) inverted directly.
  (c) omega from SIMULATED RAW RECORDS by the exact joint marginal likelihood.
      Given K_m the model is LINEAR in (theta_i, k_e,i, b_i), so b_i (Gaussian
      prior) is integrated analytically, theta_i >= 0 is either profiled (ML)
      or integrated under a flat prior on [0, inf) (REML-like), both in closed
      form, and ln k_e,i ~ N(mu, omega^2) is integrated numerically on a fine
      grid. mu and the shared K_m are profiled (likelihood) or integrated
      (Bayesian posterior for omega, half-normal(1) prior). No Gaussian
      per-patient summary is used anywhere in (c).
  Also: the Gaussian-summary model re-simulated with independent code; the
  asymptotic "N for 80% power" formula re-derived and checked by simulation at
  the N it returns; raw-likelihood power at larger N.

Assumptions shared with t9 (flagged, not tested unless said): every patient
has theta = 150 (5%), K_m = 100, the load 322.73 (the marker solves the balance
at the patient's own k_e), sigma evaluated at the patient's own untreated
level and treated as KNOWN in the fit. Different from t9's C(iv): the true
load-record bias b_i is DRAWN from its N(0, (5% L)^2) prior (t9's C(iv) used
b_i = 0 with the prior centred on it).

Run from this directory:
  python3 -W error t9_omega_recheck.py > out/t9_omega_recheck.txt
"""
import time

import numpy as np
from scipy.integrate import quad
from scipy.optimize import lsq_linear
from scipy.special import log_ndtr, logsumexp, ndtr
from scipy.stats import norm

T0 = time.time()


def el():
    return f"[{time.time() - T0:.0f} s]"


def banner(s):
    print("\n" + "=" * 100)
    print("  " + s)
    print("=" * 100)


# =====================================================================
#  0. Footing, rebuilt from the repository's definitions
# =====================================================================
TH0, KM0, KE0, MB0 = 3000.0, 100.0, 5.0, 10.0
L_REC = TH0 * MB0 / (KM0 + MB0) + KE0 * MB0       # convention A, held
TH5 = 0.05 * TH0
LN5 = np.log(KE0)
SB = 0.05 * L_REC                                  # load-bias prior sd
KE_SD = 0.8
D_AMP, LAM, TAU = 6000.0, np.log(2.0) / 3.5, 14.0
MONTHLY = np.arange(0.0, 600 * 30.4, 30.4)
MONTHLY = MONTHLY[MONTHLY < 730.0]
DENSE = np.arange(0.0, 730.0, 2.0)
CHI95, CHI90 = 3.841458820694124, 2.705543454095404


def dose(t):
    t = np.asarray(t, float)
    out = np.zeros_like(t)
    for k in range(200):
        on = t >= k * TAU
        out[on] += D_AMP * np.exp(-LAM * (t[on] - k * TAU))
    return out


INC = {"monthly": dose(MONTHLY), "2d": dose(DENSE)}
TIMES = {"monthly": MONTHLY, "2d": DENSE}


def level(th, K, ke, L, inc=0.0):
    B = ke * K + th + inc - L
    return (-B + np.sqrt(B * B + 4.0 * ke * L * K)) / (2.0 * ke)


def intake(seed, n, cv=0.10):
    rng = np.random.default_rng(seed)
    s2 = np.log1p(cv ** 2)
    return L_REC * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))


def sigma(th, K, ke):
    """hypot(5% of the balanced load, C'(m_bar) * 5% * m_bar) at the untreated level."""
    mb = level(th, K, ke, L_REC)
    return float(np.hypot(0.05 * (th * mb / (K + mb) + ke * mb),
                          (th * K / (K + mb) ** 2 + ke) * 0.05 * mb)), float(mb)


def check_footing():
    import core as C                          # assertion only
    worst = 0.0
    for des in ("monthly", "2d"):
        for s in (0, 5, 23):
            r = C.record(TH5, seed=s, times=TIMES[des])
            m = level(TH5, KM0, KE0, intake(s, TIMES[des].size), INC[des])
            worst = max(worst, np.max(np.abs(r["m"] / m - 1)), np.max(np.abs(r["inc"] - INC[des])))
            assert abs(C.sigmas(r)[0] - sigma(TH5, KM0, KE0)[0]) < 1e-12
    assert worst < 1e-13
    sg, mb = sigma(TH5, KM0, KE0)
    print(f"  rebuilt footing: L = {L_REC:.3f}, m_bar(5%) = {mb:.3f}, sigma = {sg:.4f}, b prior sd = {SB:.3f};")
    print(f"  {MONTHLY.size} monthly draws, {DENSE.size} every-2-d draws; max deviation from core.record "
          f"over 6 records = {worst:.1e} (asserted < 1e-13), sigma equal to core.sigmas (asserted).")


# =====================================================================
#  1. Information matrices (complex step), marginal-prior cohort form
# =====================================================================
def pred_c(p, m, inc):
    lt, lk, b, lke = p
    return (np.exp(lt) + inc) * m / (np.exp(lk) + m) + np.exp(lke) * m + b


def patient_F(seed, des, th=TH5, ke=KE0, h=1e-30):
    """Data-only information in (ln theta, ln K_m, b, ln k_e), b prior kept, no k_e prior."""
    inc = INC[des]
    m = level(th, KM0, ke, intake(seed, inc.size), inc)
    sg, _ = sigma(th, KM0, ke)
    p0 = np.array([np.log(th), np.log(KM0), 0.0, np.log(ke)], complex)
    J = np.zeros((m.size, 4))
    for j in range(4):
        p = p0.copy()
        p[j] += 1j * h
        J[:, j] = pred_c(p, m, inc).imag / h
    F = J.T @ J / sg ** 2
    F[2, 2] += SB ** -2
    return F


def cohort_marginal(Fs, omega, target=0):
    """sd ln theta_target and sd of its own ln k_e, population level integrated out."""
    n = len(Fs)
    if omega == 0.0:                          # one shared ln k_e, prior 0.8
        H = np.zeros((2 + 2 * n, 2 + 2 * n))
        for i, F in enumerate(Fs):
            idx = [2 + 2 * i, 0, 3 + 2 * i, 1]
            H[np.ix_(idx, idx)] += F
        H[1, 1] += KE_SD ** -2
        Ci = np.linalg.inv(H)
        return np.sqrt(Ci[2 + 2 * target, 2 + 2 * target]), np.sqrt(Ci[1, 1])
    H = np.zeros((1 + 3 * n, 1 + 3 * n))
    kei = []
    for i, F in enumerate(Fs):
        idx = [1 + 3 * i, 0, 2 + 3 * i, 3 + 3 * i]
        H[np.ix_(idx, idx)] += F
        kei.append(3 + 3 * i)
    S = omega ** 2 * np.eye(n) + (KE_SD ** 2 - omega ** 2) * np.ones((n, n))
    H[np.ix_(kei, kei)] += np.linalg.inv(S)
    Ci = np.linalg.inv(H)
    return np.sqrt(Ci[1 + 3 * target, 1 + 3 * target]), np.sqrt(Ci[kei[target], kei[target]])


def fe_cohort_V(Fs):
    """Covariance of the data-only ln k_e_i estimates, K_m shared, full matrix inverted."""
    n = len(Fs)
    H = np.zeros((1 + 3 * n, 1 + 3 * n))
    kei = []
    for i, F in enumerate(Fs):
        idx = [1 + 3 * i, 0, 2 + 3 * i, 3 + 3 * i]
        H[np.ix_(idx, idx)] += F
        kei.append(3 + 3 * i)
    Ci = np.linalg.inv(H)
    return Ci[np.ix_(kei, kei)], float(np.sqrt(Ci[0, 0]))


def part_a_b():
    banner("(a) N = 24 DISEASE COHORT, FIM IN THE MARGINAL-PRIOR FORM (independent algebra)")
    Fm = [patient_F(s, "monthly") for s in range(24)]
    print("  target = patient 0 (intake seed 0). t9 printed 1.172 (omega 0) and 1.817 (omega 0.33).")
    print(f"     {'omega':>7} {'sd ln theta_0':>14} {'sd own ln ke_0':>15}")
    A = {}
    for w in (0.0, 1e-3, 0.1, 0.2, 0.33, 0.5):
        a = cohort_marginal(Fm, w)
        A[w] = a
        print(f"     {w:7.3f} {a[0]:14.4f} {a[1]:15.4f}")
    print(f"  continuity omega 0 -> 1e-3: {abs(A[1e-3][0] - A[0.0][0]):.1e} (two different parameterisations)")
    one = cohort_marginal(Fm[:1], 0.0)[0]
    print(f"  one patient alone: {one:.4f}")

    # plateau: N -> infinity, K_m and mu known, own ln k_e ~ N(mu, omega^2)
    print("\n  Plateau check (t9 part B). N -> inf limit: K_m and the population mean known, own ln k_e")
    print("  carries prior sd omega. 'user approx' = sqrt(0.39^2 + ((1-s)/s omega)^2).")
    sg, mb = sigma(TH5, KM0, KE0)
    s = TH5 * mb / (KM0 + mb) / L_REC
    amp = (1 - s) / s
    print(f"  s = {s:.4f}, (1-s)/s = {amp:.3f}")
    keep = [0, 2, 3]
    print(f"     {'omega':>6} {'user approx':>12} {'limit median':>13} {'seed 0':>8} {'N=1024 seed 0':>14}")
    F1024 = [patient_F(si, "monthly") for si in range(1024)]
    for w in (0.2, 0.33, 0.5):
        lim = []
        for F in Fm:
            G = F[np.ix_(keep, keep)].copy()
            G[2, 2] += w ** -2
            lim.append(np.sqrt(np.linalg.inv(G)[0, 0]))
        lim = np.array(lim)
        big = cohort_marginal(F1024, w)[0]
        print(f"     {w:6.2f} {np.hypot(0.39, amp * w):12.3f} {np.median(lim):13.3f} {lim[0]:8.3f} {big:14.3f}")

    banner("(b) PER-PATIENT DATA-ONLY sd OF ln k_e (no k_e prior; b prior 5% L)")
    print("  t9 printed: monthly 1.064 [0.73, 1.33] (K_m shared), 1.349 (per patient), 1.052 (known);")
    print("              every 2 d 0.378 [0.35, 0.41], 0.379, 0.378.")
    Fd = [patient_F(si, "2d") for si in range(24)]
    Vs = {}
    for lab, FF in (("monthly", Fm), ("2d", Fd)):
        V, sdk = fe_cohort_V(FF)
        Vs[lab] = V
        sh = np.sqrt(np.diag(V))
        per = np.array([np.sqrt(np.linalg.inv(F)[3, 3]) for F in FF])
        kn = np.array([np.sqrt(np.linalg.inv(F[np.ix_(keep, keep)])[2, 2]) for F in FF])
        pk = np.median([np.sqrt(np.linalg.inv(F)[1, 1]) for F in FF])
        dv = np.sqrt(np.diag(V))
        mc = np.max(np.abs(V / np.outer(dv, dv) - np.eye(24)))
        print(f"  {lab:>8}: K_m shared {np.median(sh):.4f} [{sh.min():.3f}, {sh.max():.3f}]; per patient "
              f"{np.median(per):.4f}; K_m known {np.median(kn):.4f}; max corr {mc:.3f}; "
              f"sd ln K_m: cohort {sdk:.4f}, per patient {pk:.4f}")
    return Vs


# =====================================================================
#  2. Gaussian-summary model, re-simulated with independent code
# =====================================================================
W_G = np.r_[0.0, np.geomspace(1e-3, 30.0, 1200)]


def gauss_dev(V, Z):
    """ML deviance in omega (mu by GLS) for each row of Z, by Cholesky per omega."""
    n = V.shape[0]
    one = np.ones(n)
    D = np.empty((W_G.size, Z.shape[0]))
    for k, w in enumerate(W_G):
        Lc = np.linalg.cholesky(V + w * w * np.eye(n))
        ld = 2.0 * np.log(np.diag(Lc)).sum()
        A = np.linalg.solve(Lc, np.c_[one, Z.T])
        a1, az = A[:, 0], A[:, 1:]
        mu = (a1 @ az) / (a1 @ a1)
        R = az - np.outer(a1, mu)
        D[k] = ld + (R * R).sum(0)
    return D


def interval_from(D, wg):
    """95% profile interval (lo, hi), LR at omega = 0, omega_hat; D is (n_w, n_rep)."""
    rel = D - D.min(0)
    lo, hi = np.zeros(D.shape[1]), np.zeros(D.shape[1])
    for j in range(D.shape[1]):
        idx = np.flatnonzero(rel[:, j] <= CHI95)
        i0, i1 = idx[0], idx[-1]
        if i0 > 0:
            lo[j] = np.interp(CHI95, [rel[i0, j], rel[i0 - 1, j]], [wg[i0], wg[i0 - 1]])
        if i1 < wg.size - 1:
            hi[j] = np.interp(CHI95, [rel[i1, j], rel[i1 + 1, j]], [wg[i1], wg[i1 + 1]])
        else:
            hi[j] = wg[-1]
    return lo, hi, rel[0], wg[np.argmin(D, 0)]


def part_gauss(Vs):
    banner("GAUSSIAN-SUMMARY MODEL z ~ N(mu 1, V + omega^2 I), RE-SIMULATED (independent code)")
    print("  V = the K_m-shared data-only covariance from (b). 2000 cohorts per cell. t9 printed")
    print("  (4000 cohorts): monthly [0, 0.894] / [0, 1.000]; every 2 d [0, 0.325] / [0, 0.542].")
    print(f"     {'design':>8} {'w':>5} {'upper med [P10,P90]':>24} {'excl 0':>7} {'test':>6} {'covers w':>9}")
    G = {}
    for d_i, lab in enumerate(("monthly", "2d")):
        V = Vs[lab]
        Lc = np.linalg.cholesky(V)
        for w in (0.0, 0.33):
            rng = np.random.default_rng(31_000 + 10 * d_i + int(100 * w))
            Z = LN5 + w * rng.standard_normal((2000, 24)) + rng.standard_normal((2000, 24)) @ Lc.T
            lo, hi, lr0, _ = interval_from(gauss_dev(V, Z), W_G)
            G[(lab, w)] = (lo, hi)
            print(f"     {lab:>8} {w:5.2f} {np.median(hi):7.3f} [{np.percentile(hi, 10):.3f},{np.percentile(hi, 90):.3f}]"
                  f" {np.mean(lo > 0):7.1%} {np.mean(lr0 > CHI90):6.1%} {np.mean((lo <= w) & (hi >= w)):9.1%}")

    # asymptotic N for 80% power, re-derived: Wald on omega^2, var 2 (v + w^2)^2 per patient
    print("\n  Asymptotic N for the one-sided 5% test of omega = 0 with 80% power at omega = 0.33,")
    print("  re-derived: omega2_hat ~ N(omega^2, 2 h(omega)^2 / N), h^-2 = mean (v_i + omega^2)^-2,")
    print("  so N = 2 (z_.95 h(0) + z_.80 h(w))^2 / w^4 (diagonal V). t9 printed 1200 (monthly), 33 (2 d).")
    za, zb = norm.ppf(0.95), norm.ppf(0.80)
    Nreq = {}
    for lab in ("monthly", "2d"):
        v = np.diag(Vs[lab])
        h = lambda w: np.mean((v + w * w) ** -2) ** -0.5
        Nreq[lab] = 2.0 * ((za * h(0.0) + zb * h(0.33)) / 0.33 ** 2) ** 2
        print(f"     {lab:>8}: N = {Nreq[lab]:.0f}")
    # simulation check of that N inside the Gaussian model (diagonal V, v cycled)
    print("  Simulation check inside the Gaussian model at that N (diagonal V, the 24 v_i cycled),")
    print("  boundary LR test D(0) - D_min > 2.71; 1000 cohorts each:")
    wg = np.r_[0.0, np.geomspace(0.01, 3.0, 160)]
    for lab in ("monthly", "2d"):
        n = int(round(Nreq[lab]))
        v = np.resize(np.diag(Vs[lab]), n)
        for w in (0.0, 0.33):
            rng = np.random.default_rng(32_000 + n + int(100 * w))
            rej = []
            for blk in range(4):
                Z = LN5 + rng.standard_normal((250, n)) * np.sqrt(v + w * w)
                D = np.empty((wg.size, 250))
                for k, ww in enumerate(wg):
                    s2 = v + ww * ww
                    mu = (Z / s2).sum(1) / (1 / s2).sum()
                    D[k] = np.log(s2).sum() + (((Z - mu[:, None]) ** 2) / s2).sum(1)
                rej.append(D[0] - D.min(0) > CHI90)
            print(f"     {lab:>8}, N = {n}, true omega {w}: rejects omega = 0 in {np.mean(np.concatenate(rej)):.1%}")
    return G


# =====================================================================
#  3. Raw-record joint marginal likelihood for omega
# =====================================================================
#  u = ln k_e grid: coarse, medium, fine (trapezoid weights); analytic tail below.
U = np.r_[LN5 + np.arange(-12.0, -5.0, 0.1), LN5 + np.arange(-5.0, -2.0, 0.04),
          LN5 + np.arange(-2.0, 2.5 + 1e-9, 0.01)]
WU = np.zeros(U.size)
_d = np.diff(U)
WU[:-1] += _d / 2
WU[1:] += _d / 2
MU = LN5 + np.arange(-4.5, 2.0 + 1e-9, 0.04)
OM = np.r_[0.0, 0.01, 0.02, np.arange(0.03, 0.7, 0.03), np.arange(0.75, 1.5, 0.05),
           np.arange(1.5, 3.0 + 1e-9, 0.1), 3.25, 3.5, 3.75, 4.0]
OM_GH = OM[(OM > 0) & (OM < 0.025)]
OM_GR = OM[OM >= 0.025]
GHX, GHW = np.polynomial.hermite_e.hermegauss(10)
GHW = GHW / np.sqrt(2 * np.pi)


def kernels(mu, om):
    """Stacked (n_om * n_mu, n_u) kernel N(u; mu, w) * trapezoid weight, and lower-tail masses."""
    Kr = np.exp(-0.5 * ((U[None, None, :] - mu[None, :, None]) / om[:, None, None]) ** 2) \
        / (om[:, None, None] * np.sqrt(2 * np.pi)) * WU[None, None, :]
    tail = ndtr((U[0] - mu[None, :]) / om[:, None])
    return Kr.reshape(om.size * mu.size, U.size), tail.reshape(-1)


KERN, TAIL = kernels(MU, OM_GR)


def stats(ms, inc, ys, sgs, lnK):
    """GLS cross-products with b integrated: W = (sigma^2 I + sb^2 11')^-1. Arrays (n_pat, n_K)."""
    K = np.exp(lnK)[:, None]
    g = ms[:, None, :] / (K[None] + ms[:, None, :])
    mm = np.broadcast_to(ms[:, None, :], g.shape)
    r = ys[:, None, :] - inc * g
    n = ms.shape[1]
    c = (SB ** 2 / (sgs ** 2 + n * SB ** 2))[:, None]
    s2 = (sgs ** 2)[:, None]

    def dot(a, b):
        return ((a * b).sum(-1) - c * a.sum(-1) * b.sum(-1)) / s2
    return dict(gg=dot(g, g), gm=dot(g, mm), mm=dot(mm, mm), gr=dot(g, r), mr=dot(mm, r), rr=dot(r, r))


def loglik(S, u, mode, zread=None, zsd=None):
    """Per-patient log-likelihood in u = ln k_e, theta >= 0 profiled ('prof') or integrated ('int')."""
    ke = np.exp(u)
    gg, gm, mm, gr, mr, rr = (S[k][..., None] for k in ("gg", "gm", "mm", "gr", "mr", "rr"))
    th = (gr - ke * gm) / gg
    rss0 = rr - 2 * ke * mr + ke * ke * mm
    Q = rss0 - th * th * gg
    if mode == "prof":
        out = -0.5 * np.where(th >= 0, Q, rss0)
    else:
        out = -0.5 * Q + log_ndtr(th * np.sqrt(gg)) + 0.5 * np.log(2 * np.pi / gg)
    if zread is not None:
        out = out - 0.5 * ((u - zread[:, None, None]) / zsd) ** 2
    return out


def para_min(D):
    """Min along the last axis with parabolic refinement; also whether it sat on the edge."""
    j = np.argmin(D, -1)
    n = D.shape[-1]
    jj = np.clip(j, 1, n - 2)
    a = np.take_along_axis(D, (jj - 1)[..., None], -1)[..., 0]
    b = np.take_along_axis(D, jj[..., None], -1)[..., 0]
    c = np.take_along_axis(D, (jj + 1)[..., None], -1)[..., 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        den = a - 2 * b + c
        ref = np.where(den > 0, b - (c - a) ** 2 / (8 * den), b)
    ref = np.where(np.isfinite(ref), np.minimum(ref, D.min(-1)), D.min(-1))
    edge = (j == 0) | (j == n - 1)
    return np.where(edge, D.min(-1), ref), edge


def D_cube(S, mode, zread=None, zsd=None, mu=MU, kern=KERN, tail=TAIL, om_gr=OM_GR):
    """-2 log L(omega, mu, K) summed over patients; omega on OM (0, GH values, grid values)."""
    n, nK = S["gg"].shape
    out = np.empty((1 + OM_GH.size + om_gr.size, mu.size, nK))
    out[0] = -2 * loglik(S, mu, mode, zread, zsd).sum(0).T
    for k, w in enumerate(OM_GH):
        pts = (mu[:, None] + w * GHX[None, :]).ravel()
        lp = loglik(S, pts, mode, zread, zsd).reshape(n, nK, mu.size, GHX.size)
        out[1 + k] = -2 * logsumexp(lp, axis=-1, b=GHW).sum(0).T
    ll = loglik(S, U, mode, zread, zsd)
    M = ll.max(-1)
    E = np.exp(ll - M[..., None])
    I = kern @ E.reshape(n * nK, U.size).T + tail[:, None] * E[..., 0].reshape(1, -1)
    with np.errstate(divide="ignore"):
        lI = np.log(I).reshape(om_gr.size, mu.size, n, nK) + M[None, None]
    out[1 + OM_GH.size:] = -2 * lI.sum(2)
    return out, float(E[..., -1].max())


def make_cohort(rng, omega, des, n, seedbase, theta_sd=0.0, b_random=True):
    eta = rng.normal(0.0, omega, n) if omega > 0 else np.zeros(n)
    th = TH5 * np.exp(rng.normal(0.0, theta_sd, n)) if theta_sd else np.full(n, TH5)
    b = rng.normal(0.0, SB, n) if b_random else np.zeros(n)
    inc = INC[des]
    ms = np.empty((n, inc.size))
    ys = np.empty_like(ms)
    sgs = np.empty(n)
    for i in range(n):
        ke = KE0 * np.exp(eta[i])
        m = level(th[i], KM0, ke, intake(seedbase + i, inc.size), inc)
        sg, _ = sigma(th[i], KM0, ke)
        ys[i] = (th[i] + inc) * m / (KM0 + m) + ke * m + b[i] + rng.normal(0.0, sg, m.size)
        ms[i], sgs[i] = m, sg
    return ms, ys, sgs, eta


def k_center(ms, inc, ys, sgs):
    """The cohort's fixed-effects K_m estimate (each patient's own k_e), to centre the K_m window."""
    lk = np.log(KM0) + np.arange(-0.3, 0.3 + 1e-9, 0.01)
    fe = np.zeros(lk.size)
    for a in range(0, ms.shape[0], 48):
        S = stats(ms[a:a + 48], inc, ys[a:a + 48], sgs[a:a + 48], lk)
        fe += -2 * loglik(S, U, "prof").max(-1).sum(0)
    return lk[np.argmin(fe)], bool(np.argmin(fe) in (0, lk.size - 1))


LOG_PRIOR_MU = norm.logpdf(MU, LN5, 2.0)


def summarise_D(D):
    """Profile over K (if present) then mu: deviance curve in omega, plus edge flags."""
    if D.ndim == 3:
        DK, eK = para_min(D)
    else:
        DK, eK = D, np.zeros(D.shape, bool)
    Dw, eM = para_min(DK)
    rel = Dw - Dw.min()
    inside = rel <= CHI95
    idx = np.flatnonzero(inside)
    i0, i1 = idx[0], idx[-1]
    if i0 == 0:
        lo = 0.0
    elif np.isfinite(rel[i0 - 1]):
        lo = float(np.interp(CHI95, [rel[i0], rel[i0 - 1]], [OM[i0], OM[i0 - 1]]))
    else:
        lo = float(OM[i0])
    if i1 == OM.size - 1:
        hi = float(OM[-1])
    elif np.isfinite(rel[i1 + 1]):
        hi = float(np.interp(CHI95, [rel[i1], rel[i1 + 1]], [OM[i1], OM[i1 + 1]]))
    else:
        hi = float(OM[i1])
    jm = np.argmin(DK, -1)
    mu_edge = bool(np.any((jm[inside] == 0) | (jm[inside] == MU.size - 1)))
    k_edge = bool(np.any(eK[inside, :][np.arange(inside.sum()), jm[inside]]))
    return dict(lo=lo, hi=hi, lr0=float(rel[0]), what=float(OM[np.argmin(Dw)]),
                cens=i1 == OM.size - 1, mu_edge=mu_edge, k_edge=k_edge, rel=rel)


def posterior(D, dK):
    """Posterior for omega: half-normal(1) prior (and uniform(0,4)), mu ~ N(ln 5, 2^2), ln K flat."""
    lml = logsumexp(-0.5 * D + LOG_PRIOR_MU[None, :, None], axis=(1, 2)) + np.log(0.04 * dK)
    out = {}
    for name, lp in (("hn1", norm.logpdf(OM, 0, 1.0) + np.log(2)), ("u4", np.zeros(OM.size))):
        lw = lml + lp
        p = np.exp(lw - lw.max())
        cdf = np.r_[0.0, np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(OM))]
        cdf /= cdf[-1]
        q = lambda a: float(np.interp(a, cdf, OM))
        out[name] = dict(med=q(0.5), up95=q(0.95), p02=1 - float(np.interp(0.2, OM, cdf)),
                         p01=float(np.interp(0.1, OM, cdf)))
    i33 = np.argmin(np.abs(OM - 0.33))
    out["lbf"] = float(lml[i33] - lml[0])           # log Bayes factor omega = 0.33 vs 0
    return out


def analyse(ms, ys, sgs, des, zread=None, zsd=None, perpatient=True, bayes=True):
    inc = INC[des]
    kc, kc_edge = k_center(ms, inc, ys, sgs)
    lK = kc + np.arange(-0.06, 0.06 + 1e-9, 0.015)
    S = stats(ms, inc, ys, sgs, lK)
    res = dict(kc_edge=kc_edge)
    for mode in ("prof", "int"):
        D, etop = D_cube(S, mode, zread, zsd)
        res[mode] = summarise_D(D)
        res["etop"] = max(res.get("etop", 0.0), etop)
        if mode == "int" and bayes:
            res["post"] = posterior(D, 0.015)
    if perpatient:
        # K_m per patient, theta profiled: profile K inside each patient's likelihood
        lKp = np.log(KM0) + np.arange(-0.4, 0.4 + 1e-9, 0.02)
        Sp = stats(ms, inc, ys, sgs, lKp)

        llu = loglik(Sp, U, "prof", zread, zsd)              # (n, nK, nu)
        prof_u, _ = para_min(np.moveaxis(-llu, 1, -1))       # (n, nu) of -ll
        llm = loglik(Sp, MU, "prof", zread, zsd)
        prof_m, _ = para_min(np.moveaxis(-llm, 1, -1))
        n = ms.shape[0]
        M = (-prof_u).max(-1)
        E = np.exp(-prof_u - M[:, None])
        D = np.empty((OM.size, MU.size))
        D[0] = 2 * prof_m.sum(0)
        for k, w in enumerate(OM_GH):
            pts = (MU[:, None] + w * GHX[None, :]).ravel()
            lp = loglik(Sp, pts, "prof", zread, zsd)
            pp, _ = para_min(np.moveaxis(-lp, 1, -1))
            D[1 + k] = -2 * logsumexp(-pp.reshape(n, MU.size, GHX.size), axis=-1, b=GHW).sum(0)
        I = KERN @ E.T + TAIL[:, None] * E[:, 0][None, :]
        with np.errstate(divide="ignore"):
            lI = np.log(I).reshape(OM_GR.size, MU.size, n) + M[None, None]
        D[1 + OM_GH.size:] = -2 * lI.sum(-1)
        res["pp"] = summarise_D(D)
    return res


def validate():
    """Check the closed forms against brute force on one simulated patient."""
    banner("(c0) VALIDATION OF THE CLOSED FORMS AGAINST BRUTE FORCE (one simulated cohort)")
    rng = np.random.default_rng(424242)
    ms, ys, sgs, eta = make_cohort(rng, 0.33, "monthly", 3, 990_000)
    inc = INC["monthly"]
    lK = np.array([np.log(KM0)])
    S = stats(ms, inc, ys, sgs, lK)
    worst_p, worst_i, worst_q = 0.0, 0.0, 0.0
    for i in range(3):
        m, y, sg = ms[i], ys[i], sgs[i]
        g = m / (KM0 + m)
        for uu in (LN5 - 3, LN5 - 0.5, LN5, LN5 + 0.2, LN5 + 0.45):
            ke = np.exp(uu)
            r = y - inc * g - ke * m
            # brute profile: min over theta >= 0, b of |(r - th g - b)/sg|^2 + (b/sb)^2
            A = np.r_[np.c_[g / sg, np.ones_like(g) / sg], [[0.0, 1.0 / SB]]]
            bb = np.r_[r / sg, 0.0]
            sol = lsq_linear(A, bb, bounds=([0.0, -np.inf], [np.inf, np.inf]), tol=1e-14, lsmr_tol="auto")
            brute = -0.5 * 2 * sol.cost
            mine = loglik({k: v[i:i + 1] for k, v in S.items()}, np.array([uu]), "prof")[0, 0, 0]
            worst_p = max(worst_p, abs(brute - mine))
            # integrated theta: quad over theta of exp(-0.5 Q(theta)), b profiled (= integrated up to const)
            def negq(t):
                A2 = np.c_[np.ones_like(g) / sg]
                A2 = np.r_[A2, [[1.0 / SB]]]
                b2 = np.r_[(r - t * g) / sg, 0.0]
                bh = np.linalg.lstsq(A2, b2, rcond=None)[0]
                return 0.5 * np.sum((b2 - A2 @ bh) ** 2)
            ref = negq(max(0.0, (S["gr"][i, 0] - ke * S["gm"][i, 0]) / S["gg"][i, 0]))  # value at the mode
            tc = max(0.0, (S["gr"][i, 0] - ke * S["gm"][i, 0]) / S["gg"][i, 0])
            ts = 1.0 / np.sqrt(S["gg"][i, 0])
            val = quad(lambda t: np.exp(ref - negq(t)), 0.0, tc + 40 * ts, points=[tc], limit=400)[0]
            bi = np.log(val) - ref
            mi = loglik({k: v[i:i + 1] for k, v in S.items()}, np.array([uu]), "int")[0, 0, 0]
            worst_i = max(worst_i, abs(bi - mi))
        # u-integral: quad against the grid, two (omega, mu) points
        Si = {k: v[i:i + 1] for k, v in S.items()}
        ll = loglik(Si, U, "prof")[0, 0]
        Mx = ll.max()
        f = lambda uu: np.exp(loglik(Si, np.array([uu]), "prof")[0, 0, 0] - Mx)
        for w, mu in ((0.33, LN5), (0.9, LN5 - 0.7), (0.06, LN5 + 0.1), (2.0, LN5 - 2.0)):
            a, b = max(-60.0, mu - 12 * w), min(U[-1], mu + 12 * w)
            pts = [x for x in (mu, U[np.argmax(ll)]) if a < x < b]
            q = quad(lambda uu: f(uu) * norm.pdf(uu, mu, w), a, b, points=pts, limit=500,
                     epsabs=0.0, epsrel=1e-10)[0]
            gsum = np.sum(np.exp(ll - Mx) * norm.pdf(U, mu, w) * WU) + np.exp(ll[0] - Mx) * ndtr((U[0] - mu) / w)
            worst_q = max(worst_q, abs(np.log(gsum) - np.log(q)))
    print(f"  theta-profiled log-lik vs bounded least squares (lsq_linear): max |diff| = {worst_p:.2e}")
    print(f"  theta-integrated log-lik vs adaptive quadrature over theta:   max |diff| = {worst_i:.2e}")
    print(f"  grid integral over ln k_e vs adaptive quadrature (4 (omega, mu) x 3 patients): max |log diff| = {worst_q:.2e}")
    print("  (the brute force profiles b under its prior; the closed form uses r'(sigma^2 I + sb^2 11')^-1 r,")
    print("   which is the same minimum exactly; both drop the constant log-determinant.)")


def shape_table():
    banner("(c1) EACH PATIENT'S DATA-ONLY PROFILE IN ln k_e, CLOSED FORM (cross-check of t9 C(iv))")
    print("  t9's C(iv) data: intake seed s, noise seed 10000 + s, b = 0; K_m and theta profiled per")
    print("  patient (theta >= 0), b prior 5% L, no k_e prior. Deviance relative to the minimum over")
    print("  t9's 20 offsets. t9 printed, monthly median: 1.41 1.38 1.20 0.86 0.67 0.23 0.64 2.92 24.36 265.42")
    print("                            every 2 d median: 5.66 5.44 4.13 2.04 0.64 0.39 1.06 9.04 80.96 877.86")
    off20 = np.array([-8, -6, -4, -3, -2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.1, 0.2, 0.3,
                      0.4, 0.5, 0.75, 1.0, 1.5, 2.0])
    show = [-8, -4, -2, -1, -0.5, 0, 0.2, 0.3, 0.5, 1.0]
    cols = [int(np.flatnonzero(np.isclose(off20, o))[0]) for o in show]
    lKp = np.log(KM0) + np.arange(-3.0, 3.0 + 1e-9, 0.005)
    for des in ("monthly", "2d"):
        inc = INC[des]
        rel, walls, kedge = [], [], 0
        fine = LN5 + np.arange(-1.0, 1.0 + 1e-9, 0.002)
        for s in range(24):
            m = level(TH5, KM0, KE0, intake(s, inc.size), inc)
            sg, _ = sigma(TH5, KM0, KE0)
            y = (TH5 + inc) * m / (KM0 + m) + KE0 * m + np.random.default_rng(10_000 + s).normal(0.0, sg, m.size)
            S = stats(m[None], inc, y[None], np.array([sg]), lKp)
            p, ek = para_min(np.moveaxis(-2 * loglik(S, LN5 + off20, "prof"), 1, -1))
            p = p[0]
            kedge += int(ek.any())
            rel.append(p - p.min())
            pf, _ = para_min(np.moveaxis(-2 * loglik(S, fine, "prof"), 1, -1))
            pf = pf[0] - min(pf[0].min(), p.min())
            walls.append(fine[np.flatnonzero(pf <= CHI95)[-1]] - LN5)
        rel = np.array(rel)
        print(f"  {des:>8} median:" + "".join(f" {np.median(rel[:, c]):7.2f}" for c in cols))
        print(f"  {des:>8}    max:" + "".join(f" {rel[:, c].max():7.2f}" for c in cols))
        print(f"  {'':>8} 95% upper end in ln k_e offset: median {np.median(walls):+.3f} "
              f"[{np.min(walls):+.3f}, {np.max(walls):+.3f}]; region reaches offset -8 in "
              f"{np.mean(rel[:, 0] <= CHI95):.0%} of records (t9: +0.319 [+0.227, +0.467], 96% monthly;"
              f" +0.233 [-0.195, +0.312], 29% every 2 d); K_m on its +-3 grid edge in {kedge} records")


def run_cells(cells, label):
    rows = {}
    for (des, w, n_coh, zsd, seed0) in cells:
        t0 = time.time()
        R = []
        for k in range(n_coh):
            rng = np.random.default_rng([seed0, k])
            ms, ys, sgs, eta = make_cohort(rng, w, des, 24, seed0 * 10 ** 6 + 1000 * k)
            zr = (LN5 + eta + zsd * rng.standard_normal(24)) if zsd else None
            R.append(analyse(ms, ys, sgs, des, zr, zsd))
        rows[(des, w, zsd)] = R
        print(f"     ... {des} omega={w} read={zsd}: {n_coh} cohorts in {time.time() - t0:.0f} s {el()}")
    return rows


def print_cells(rows):
    print("\n  Profile-likelihood 95% intervals for omega, {D(omega) - D_min <= 3.84}, mu and K_m profiled.")
    print("  'excl 0' = lower end > 0; 'test' = D(0) - D_min > 2.71 (boundary 5% test); 'covers' = interval")
    print("  holds the true omega; 'both' = holds 0 and 0.33. Upper end median [P10, P90].")
    hdr = (f"     {'design':<22} {'w':>5} {'like.':>6} {'med w_hat':>9} {'upper med [P10,P90]':>23} "
           f"{'excl 0':>7} {'test':>6} {'covers':>7} {'both':>6} {'cens':>5}")
    print(hdr)
    for (des, w, zsd), R in rows.items():
        lab = des + ("" if not zsd else f" + ke read {zsd}")
        for key, nm in (("prof", "ML"), ("int", "REML"), ("pp", "ML-ppK")):
            if key not in R[0]:
                continue
            hi = np.array([r[key]["hi"] for r in R])
            lo = np.array([r[key]["lo"] for r in R])
            lr = np.array([r[key]["lr0"] for r in R])
            wh = np.array([r[key]["what"] for r in R])
            cen = np.array([r[key]["cens"] for r in R])
            relv = np.array([r[key]["rel"] for r in R])
            r33 = np.array([np.interp(0.33, OM, rv) for rv in relv])
            both = (lo == 0) & (r33 <= CHI95)
            cov = (lo <= w) & (np.array([np.interp(w, OM, rv) for rv in relv]) <= CHI95)
            print(f"     {lab:<22} {w:5.2f} {nm:>6} {np.median(wh):9.3f} {np.median(hi):7.3f} "
                  f"[{np.percentile(hi, 10):6.3f},{np.percentile(hi, 90):6.3f}] {np.mean(lo > 0):7.0%} "
                  f"{np.mean(lr > CHI90):6.0%} {np.mean(cov):7.0%} {np.mean(both):6.0%} {np.mean(cen):5.0%}")
    flags = [(r["kc_edge"], r["prof"]["k_edge"], r["prof"]["mu_edge"], r["int"]["k_edge"], r["int"]["mu_edge"],
              r["etop"]) for R in rows.values() for r in R]
    fl = np.array(flags, float)
    print(f"  grid checks over all cohorts: K centre at edge {fl[:, 0].mean():.0%}; K_m profile on the grid edge "
          f"inside the interval ML {fl[:, 1].mean():.0%} / REML {fl[:, 3].mean():.0%}; mu on the edge ML "
          f"{fl[:, 2].mean():.0%} / REML {fl[:, 4].mean():.0%}; largest scaled likelihood at the top of the "
          f"ln k_e grid {fl[:, 5].max():.1e}")
    print("\n  Bayesian posterior for omega (theta_i flat on [0, inf), b_i ~ N(0, (5% L)^2), mu ~ N(ln 5, 2^2),")
    print("  ln K_m flat on the cohort's +-0.06 window, omega ~ half-normal(1); 'u4' = omega ~ U(0, 4)).")
    print(f"  Prior alone, half-normal(1): median {norm.ppf(0.75):.3f}, 95% upper {norm.ppf(0.975):.3f}, "
          f"P(omega > 0.2) = {2 * norm.sf(0.2):.3f}.")
    print(f"     {'design':<22} {'w':>5} {'post med':>9} {'95% up: med [P10,P90]':>25} {'P(w>.2)':>8} "
          f"{'P(w>.2)>.95':>11} {'covers':>7} {'u4: 95% up':>11} {'log BF .33:0 med':>17}")
    for (des, w, zsd), R in rows.items():
        lab = des + ("" if not zsd else f" + ke read {zsd}")
        P = [r["post"] for r in R]
        up = np.array([p["hn1"]["up95"] for p in P])
        md = np.array([p["hn1"]["med"] for p in P])
        p02 = np.array([p["hn1"]["p02"] for p in P])
        u4 = np.array([p["u4"]["up95"] for p in P])
        lbf = np.array([p["lbf"] for p in P])
        print(f"     {lab:<22} {w:5.2f} {np.median(md):9.3f} {np.median(up):7.3f} [{np.percentile(up, 10):6.3f},"
              f"{np.percentile(up, 90):6.3f}] {np.mean(p02):8.2f} {np.mean(p02 > 0.95):11.0%} "
              f"{np.mean(up >= w):7.0%} {np.median(u4):11.3f} {np.median(lbf):+17.2f}")


def power_large_n(Ns, n_coh, seed0):
    """Raw-likelihood boundary test at larger N, monthly: reduced grids, theta profiled and integrated."""
    mu = LN5 + np.arange(-1.5, 0.8 + 1e-9, 0.03)
    om_gr = np.r_[np.arange(0.03, 0.7, 0.03), 0.8, 1.0, 1.3]
    kern, tail = kernels(mu, om_gr)
    out = {}
    for N in Ns:
        for w in (0.0, 0.33):
            t0 = time.time()
            rej = {"prof": [], "int": []}
            what = {"prof": [], "int": []}
            edge = 0
            for k in range(n_coh):
                rng = np.random.default_rng([seed0, N, int(100 * w), k])
                ms, ys, sgs, _ = make_cohort(rng, w, "monthly", N,
                                            (seed0 * 10 ** 6 + N * 10 ** 3 + int(100 * w)) * 10 ** 5 + 1000 * k)
                inc = INC["monthly"]
                kc, _ = k_center(ms, inc, ys, sgs)
                hw = max(0.03, 0.06 * np.sqrt(24.0 / N))
                lK = kc + np.linspace(-hw, hw, 7)
                S = stats(ms, inc, ys, sgs, lK)
                for mode in ("prof", "int"):
                    D, _ = D_cube(S, mode, mu=mu, kern=kern, tail=tail, om_gr=om_gr)
                    DK, eK = para_min(D)
                    Dw, eM = para_min(DK)
                    rej[mode].append(Dw[0] - Dw.min() > CHI90)
                    omg = np.r_[0.0, OM_GH, om_gr]
                    what[mode].append(omg[np.argmin(Dw)])
                    jm = np.argmin(DK, -1)
                    edge += int(np.any((jm == 0) | (jm == mu.size - 1)) or np.any(eK[np.arange(omg.size), jm]))
            out[(N, w)] = rej
            print(f"     N = {N:4d}, true omega {w:4.2f}: test rejects omega = 0 in ML {np.mean(rej['prof']):4.0%}, "
                  f"REML {np.mean(rej['int']):4.0%} of {n_coh} cohorts; median omega_hat ML "
                  f"{np.median(what['prof']):.3f}, REML {np.median(what['int']):.3f}; grid-edge hits {edge} "
                  f"({time.time() - t0:.0f} s) {el()}")
    return out


def main():
    banner("RECHECK OF t9_omega.py: N = 24 AND OMEGA")
    check_footing()
    Vs = part_a_b()
    print(f"  {el()}")
    part_gauss(Vs)
    print(f"  {el()}")
    validate()
    shape_table()
    print(f"  {el()}")

    banner("(c) OMEGA FROM SIMULATED RAW RECORDS: JOINT MARGINAL LIKELIHOOD AND POSTERIOR, 24 PATIENTS")
    print("  Each cohort: 24 patients at theta = 150, eta_i ~ N(0, omega^2), k_e,i = 5 e^eta_i, b_i ~ N(0, (5% L)^2),")
    print("  fresh intake seeds, Gaussian noise sigma_i. 'ML' = theta_i profiled (>= 0), K_m shared; 'REML' =")
    print("  theta_i integrated flat on [0, inf), K_m shared; 'ML-ppK' = theta_i and K_m profiled per patient.")
    print(f"  Grids: ln k_e {U.size} points (step 0.01 on ln 5 + [-2, 2.5]), mu {MU.size} points (step 0.04,"
          f" parabolic refinement), omega {OM.size} points to {OM[-1]:.0f}, K_m 9 points of 0.015 about the")
    print("  cohort's own fixed-effects estimate (parabolic refinement).")
    cells = [("monthly", 0.0, 80, None, 101), ("monthly", 0.33, 80, None, 102),
             ("2d", 0.0, 60, None, 103), ("2d", 0.33, 60, None, 104),
             ("monthly", 0.0, 40, 0.2, 105), ("monthly", 0.33, 40, 0.2, 106)]
    rows = run_cells(cells, "main")
    print_cells(rows)

    banner("(c2) SENSITIVITY: CAPACITY VARYING BETWEEN PATIENTS (ln theta_i ~ N(ln 150, 0.5^2)), monthly")
    print("  Tests whether the raw likelihood's information on omega depends on every patient sharing")
    print("  theta = 150 (so that each patient's theta >= 0 wall sits the same distance above its k_e).")
    rows2 = {}
    for w, seed0 in ((0.0, 201), (0.33, 202)):
        R = []
        for k in range(40):
            rng = np.random.default_rng([seed0, k])
            ms, ys, sgs, eta = make_cohort(rng, w, "monthly", 24, seed0 * 10 ** 6 + 1000 * k, theta_sd=0.5)
            R.append(analyse(ms, ys, sgs, "monthly", perpatient=False))
        rows2[("monthly", w, None)] = R
    print_cells(rows2)
    print(f"  {el()}")

    banner("(c3) RAW-LIKELIHOOD POWER AT LARGER N (monthly): the claim 'about 1200 patients for 80% power'")
    print("  Same model and estimators as (c), coarser grids (mu step 0.03 on ln 5 + [-1.5, 0.8], omega to 1.3).")
    power_large_n((24, 96, 384), 20, 301)
    print(f"\n  total runtime {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
