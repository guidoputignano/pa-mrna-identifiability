"""
Adversarial recheck of t10_factorial.py (follow-up item 3), by a route that
shares NO code with it: nothing is imported from core, dyn, t6_estimators,
t7_design, n4_theta_scope or the repository's src/. Every model quantity is
rebuilt here from the stated footing.

  integrator   a fixed-step classical Runge-Kutta 4 written out below (not
               LSODA), with the step grid aligned to every load step, the known
               load evaluated per step and the drug kernel at the stage times
  sensitivity  COMPLEX-STEP differentiation of that integrator (perturb one
               parameter by i*1e-30, read d m/dp off the imaginary part), not
               forward-sensitivity ODEs. Exact for the discrete solution; the
               only error left is the RK4 truncation, checked by halving h.
  record       the 25-draw monthly record rebuilt from its definition (dose
               train summed by hand, lognormal 10% intake with the same
               default_rng(seed) stream, sigma from the tangent slope at the
               untreated level), its information written in the transient's
               coordinates directly
  effects      the 2^4 effects from a least-squares fit of the full saturated
               +-1 model, the Shapley shares by brute force over all 24 orders
  estimator    a batched Levenberg-Marquardt written here (own damping, own
               box handling, complex-step Jacobian), fitted from starts that
               include and exclude the truth, run to a Gauss-Newton step check,
               on new seeds and with more replicates than t10

Footing (from the brief, re-derived, not imported): normal theta 3000, K_m 100,
k_e 5 /d, m = 10 at normal so L = 3000*10/110 + 50 = 322.73 (convention A, held
at 5% and 20%); unrecorded L0 = 0.5 L constant; ln y = ln m + N(0, cv^2); the
drug kernel rises to its peak at 13 h and decays with a 3.5 d half-life; the
record's sigma = hypot(0.05 L, C'(m_bar) 0.05 m_bar) with m_bar the untreated
level; priors k_e ~ 0.8 in ln, record bias b ~ N(0, (0.05 L)^2).

Run from this directory:
    python3 -W error t10_factorial_recheck.py > out/t10_factorial_recheck.txt 2>&1
"""
import itertools
import time

import numpy as np

T_START = time.time()

# ------------------------------------------------------------------ footing
TH0, KM, KE, M_NORMAL = 3000.0, 100.0, 5.0, 10.0
L = TH0 * M_NORMAL / (KM + M_NORMAL) + KE * M_NORMAL
L0 = 0.5 * L
TH5, TH20 = 0.05 * TH0, 0.20 * TH0
GATE = 0.20
KD = np.log(2) / 3.5
T_PEAK = 13 / 24


def _bisect(f, a, b, n=200):
    fa = f(a)
    for _ in range(n):
        c = 0.5 * (a + b)
        fc = f(c)
        if fa * fc <= 0:
            b = c
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


KR = _bisect(lambda k: np.log(k / KD) / (k - KD) - T_PEAK, KD * 1.01, 200.0)
KNORM = np.exp(-KD * T_PEAK) - np.exp(-KR * T_PEAK)


def steady(thp, ke, K, Lt):
    """Positive root of thp m/(K+m) + ke m = Lt (complex-safe)."""
    B = ke * K + thp - Lt
    return (-B + np.sqrt(B * B + 4 * ke * Lt * K)) / (2 * ke)


# ------------------------------------------------------------------ designs
class Design:
    """Known load u(t) as (t_start_days, u) steps; drug kernel peak Dpk (0 = none)."""

    def __init__(self, steps, Dpk=0.0, T_h=48.0, label=""):
        self.steps = sorted(steps)
        self.Dpk = Dpk
        self.T_h = T_h
        self.label = label

    def u(self, t):
        v = self.steps[0][1]
        for ts, us in self.steps:
            if t >= ts:
                v = us
        return v

    def delta(self, t):
        t = np.asarray(t, float)
        return np.where(t < 0, 0.0, self.Dpk * (np.exp(-KD * t) - np.exp(-KR * t)) / KNORM)


def d_withdraw(r, T_h):
    return Design([(-1, 0.5 * L), (0, (0.5 - r) * L)], T_h=T_h, label=f"withdraw {r:.0%}")


def d_up(x, T_h, hold_h=None):
    hold_h = T_h / 2 if hold_h is None else hold_h
    return Design([(-1, 0.5 * L), (0, (x - 0.5) * L), (hold_h / 24, 0.0)], T_h=T_h,
                  label=f"load x{x:g} for {hold_h:g} h")


def d_drug(Dpk, T_h):
    return Design([(-1, 0.5 * L)], Dpk=Dpk, T_h=T_h, label=f"drug {Dpk}")


# ------------------------------------------------------------------ RK4
def path(th, ke, K, l0, kap, des, h_min, rec_min):
    """m at every rec_min over [0, T] for arrays of parameters (any dtype)."""
    h = h_min / 1440.0
    n = int(round(des.T_h * 60 / h_min))
    every = int(round(rec_min / h_min))
    assert abs(n * h_min - des.T_h * 60) < 1e-9 and abs(every * h_min - rec_min) < 1e-9
    for ts, _ in des.steps:
        if ts > 0:
            assert abs(ts * 1440 / h_min - round(ts * 1440 / h_min)) < 1e-9
    us = np.array([des.u((k + 0.5) * h) for k in range(n)]) + 0.0
    dn = des.delta(np.arange(2 * n + 1) * 0.5 * h)
    m = steady(th + kap * float(des.delta(0.0)), ke, K, des.u(-1e-9) + l0)
    out = [m]
    for k in range(n):
        U = us[k] + l0
        a0, a1, a2 = th + kap * dn[2 * k], th + kap * dn[2 * k + 1], th + kap * dn[2 * k + 2]
        k1 = U - a0 * m / (K + m) - ke * m
        mm = m + 0.5 * h * k1
        k2 = U - a1 * mm / (K + mm) - ke * mm
        mm = m + 0.5 * h * k2
        k3 = U - a1 * mm / (K + mm) - ke * mm
        mm = m + h * k3
        k4 = U - a2 * mm / (K + mm) - ke * mm
        m = m + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if (k + 1) % every == 0:
            out.append(m)
    return np.array(out)


EPS = 1e-30


def sens(theta, des, h_min=1.0, rec_min=5.0, kfree=False):
    """d ln m / d(ln th, ln ke, ln Km, L0[, ln kap]) at every rec_min, complex step."""
    npar = 5 if kfree else 4
    x0 = np.array([np.log(theta), np.log(KE), np.log(KM), L0, 0.0])[:npar]
    X = np.tile(x0.astype(complex)[:, None], (1, npar))
    X[np.arange(npar), np.arange(npar)] += 1j * EPS
    kap = np.exp(X[4]) if kfree else np.ones(npar, complex)
    m = path(np.exp(X[0]), np.exp(X[1]), np.exp(X[2]), X[3], kap, des, h_min, rec_min)
    return (m.imag / EPS) / m.real, m.real[:, 0]


def sd_th(F):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


# ------------------------------------------------------------------ the monthly record
REC_T = 30.4 * np.arange(25)
LAM_DOSE = np.log(2) / 3.5
REC_INC = np.array([sum(6000.0 * np.exp(-LAM_DOSE * (t - 14.0 * k)) for k in range(200) if t >= 14.0 * k)
                    for t in REC_T])


def record(theta, seed):
    rng = np.random.default_rng(seed)
    s2 = np.log1p(0.10 ** 2)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), REC_T.size))
    m = steady(theta + REC_INC, KE, KM, loads)
    mbar = steady(theta, KE, KM, L)
    slope = theta * KM / (KM + mbar) ** 2 + KE
    sig = np.hypot(0.05 * L, slope * 0.05 * mbar)
    return m, sig


def record_info(theta, seed, kfree=False):
    """Record information in z = (ln th, ln ke, ln Km, L0, b[, ln kap]); b prior included."""
    m, sig = record(theta, seed)
    rows = np.zeros((m.size, 6 if kfree else 5))
    rows[:, 0] = theta * m / (KM + m)
    rows[:, 1] = KE * m
    rows[:, 2] = -(theta + REC_INC) * KM * m / (KM + m) ** 2
    rows[:, 4] = 1.0
    rows /= sig
    Z = rows.T @ rows
    Z[4, 4] += (0.05 * L) ** -2
    return Z


# knowledge rungs: indices in (ln th, ln ke, ln Km, L0, [b], [ln kap])
RUNGS = (("none", {}), ("k_e 0.8", {1: 0.8}), ("+K_m 0.062", {1: 0.8, 2: 0.062}),
         ("+L0 5%", {1: 0.8, 2: 0.062, 3: 0.05 * L}), ("+L0 1%", {1: 0.8, 2: 0.062, 3: 0.01 * L}))
KCOST = (0, 0, 0, 1, 2)


def alone(F, pri):
    F = F.copy()
    for j, s in pri.items():
        F[j, j] += s ** -2
    return F


def joined(F, pri, Zrec):
    """Transient info F in (th, ke, Km, L0[, kap]) joined to the record; K_m prior dropped."""
    n = F.shape[0]
    Z = Zrec.copy()
    idx = [0, 1, 2, 3] + ([5] if n == 5 else [])
    Z[np.ix_(idx, idx)] += F
    for j, s in pri.items():
        if j == 2:
            continue
        jj = 5 if j == 4 else j
        Z[jj, jj] += s ** -2
    return Z


def banner(s):
    print("\n" + "=" * 100 + "\n  " + s + "\n" + "=" * 100)


# t10's printed values, copied from out/t10_factorial.txt, for side-by-side comparison only
T10 = dict(
    cells_np=[383.3133, 3.1140, 76.6627, 0.6228, 199.6758, 1.6088, 39.9352, 0.3218,
              280.2237, 0.8666, 56.0447, 0.1733, 148.3284, 0.4464, 29.6657, 0.0893],
    cells_r=[4.1260, 0.5651, 1.3931, 0.1181, 2.9875, 0.2929, 0.7397, 0.0674,
             0.9753, 0.1316, 0.5815, 0.0520, 0.8854, 0.0785, 0.3531, 0.0435],
    main_np=[-0.3444, -0.2836, -0.6990, -2.3038], main_r=[-0.4507, -0.1882, -0.4599, -0.9831],
    capxexc_np=-0.2118, capxcv_r=0.1391,
    shap_np=[-0.3450, -0.2843, -0.6990, -2.3045], shap_r=[-0.4366, -0.1618, -0.4279, -0.9509],
    var_main_np=0.993, kmsd5=0.0617, kmsd20=0.0606,
    n_pass=1422, cls_cv=[[0, 0, 0, 36], [0, 9, 33, 71], [5, 30, 59, 84], [20, 52, 77, 84],
                         [27, 64, 78, 89], [37, 67, 83, 102], [37, 75, 91, 112]],
    pareto=(27, 452))


# =========================================================== (a) factorial
def part_a():
    banner("(a) THE 2^4 FACTORIAL, RECOMPUTED (RK4 + complex step; own effects and Shapley)")
    print(f"  footing: L = {L:.4f}, L0 = {L0:.4f}; m at 5% untreated {steady(TH5, KE, KM, L):.4f}, "
          f"at 20% {steady(TH20, KE, KM, L):.4f}; drug kernel rise {KR:.4f} /d, peak {24*T_PEAK:.1f} h")
    cells = list(itertools.product((0, 1), repeat=4))
    y0, y2 = [], []
    t0 = time.time()
    for x in cells:
        th = (TH5, TH20)[x[0]]
        tau = (60, 15)[x[1]]
        cv = (0.05, 0.01)[x[2]]
        des = (Design([(-1, 0.5 * L), (0, 0.0)], T_h=48) if x[3] == 0 else
               Design([(-1, 0.5 * L), (0, 5.5 * L), (1.0, 0.0)], T_h=48))
        R, _ = sens(th, des, h_min=1.0, rec_min=tau)
        F = R.T @ R / cv ** 2
        y0.append(sd_th(F))
        y2.append(sd_th(alone(F, {1: 0.8, 2: 0.062})))
    # step-size check on the two corners
    chk = []
    for x, th, tau, cv, des in (((0,) * 4, TH5, 60, 0.05, Design([(-1, 0.5 * L), (0, 0.0)], T_h=48)),
                                ((1,) * 4, TH20, 15, 0.01, Design([(-1, 0.5 * L), (0, 5.5 * L), (1.0, 0.0)], T_h=48))):
        v = []
        for h in (0.5, 1.0, 5.0):
            R, _ = sens(th, des, h_min=h, rec_min=tau)
            v.append(sd_th(R.T @ R / cv ** 2))
        chk.append(v)
    print(f"  RK4 step check, sd no priors at h = 0.5 / 1 / 5 min: TB-A {chk[0][0]:.4f} / {chk[0][1]:.4f} / "
          f"{chk[0][2]:.4f};  TB-B {chk[1][0]:.5f} / {chk[1][1]:.5f} / {chk[1][2]:.5f}")
    print(f"\n  {'cell (cap,samp,CV,exc)':24} {'mine, none':>11} {'t10':>10} {'rel':>8}   {'mine, rung':>11} "
          f"{'t10':>8} {'rel':>8}")
    for i, x in enumerate(cells):
        a, b = y0[i], T10["cells_np"][i]
        c, d = y2[i], T10["cells_r"][i]
        print(f"  {str(x):24} {a:11.4f} {b:10.4f} {a/b-1:+8.1e}   {c:11.4f} {d:8.4f} {c/d-1:+8.1e}")
    print(f"  ({time.time()-t0:.1f} s)")
    # effects by least squares on the saturated +-1 model
    Xc = 2 * np.array(cells) - 1
    terms = [S for k in range(1, 5) for S in itertools.combinations(range(4), k)]
    M = np.column_stack([np.ones(16)] + [np.prod(Xc[:, list(S)], axis=1) for S in terms])
    names = ("capacity", "sampling", "assay CV", "excursion")
    out = {}
    for lab, y, mk, sk in (("no priors", np.log10(y0), "main_np", "shap_np"),
                           ("k_e 0.8 + K_m 0.062", np.log10(y2), "main_r", "shap_r")):
        beta = np.linalg.lstsq(M, y, rcond=None)[0]
        eff = {S: 2 * b for S, b in zip(terms, beta[1:])}
        sst = float(np.sum((y - y.mean()) ** 2))
        main = sum(16 * (e / 2) ** 2 for S, e in eff.items() if len(S) == 1) / sst
        two = sum(16 * (e / 2) ** 2 for S, e in eff.items() if len(S) == 2) / sst
        # Shapley by brute force over the 24 orders
        pos = {c: i for i, c in enumerate(cells)}
        phi = np.zeros(4)
        for order in itertools.permutations(range(4)):
            cur = [0, 0, 0, 0]
            for f in order:
                before = y[pos[tuple(cur)]]
                cur[f] = 1
                phi[f] += (y[pos[tuple(cur)]] - before) / 24
        dAB = y[-1] - y[0]
        print(f"\n  [{lab}] log10 sd(B) - log10 sd(A) = {dAB:+.4f}; main-effect share of the variance {main:.1%},"
              f" two-way {two:.1%}")
        print("     main effects (mine | t10): " + ", ".join(
            f"{names[i]} {eff[(i,)]:+.4f} | {T10[mk][i]:+.4f}" for i in range(4)))
        two_way = sorted(((abs(e), S, e) for S, e in eff.items() if len(S) == 2), reverse=True)
        print("     largest two-way terms: " + ", ".join(
            f"{names[S[0]]} x {names[S[1]]} {e:+.4f}" for _, S, e in two_way[:3]))
        print("     Shapley (mine | t10): " + ", ".join(
            f"{names[i]} {phi[i]:+.4f} ({phi[i]/dAB:.1%}) | {T10[sk][i]:+.4f}" for i in range(4)))
        out[lab] = dict(eff=eff, phi=phi, dAB=dAB, main=main)
    # K_m from the record, 24 seeds, with and without the k_e prior
    print("\n  sd(ln K_m) the monthly record alone gives, over seeds 0..23 (median [min, max]):")
    km = {}
    for th in (TH5, TH20):
        for kp in (0.8, None):
            v = []
            for s in range(24):
                Z = record_info(th, s)[np.ix_([0, 1, 2, 4], [0, 1, 2, 4])]
                if kp:
                    Z[1, 1] += kp ** -2
                v.append(np.sqrt(np.linalg.inv(Z)[2, 2]))
            v = np.array(v)
            km[(th, kp)] = v
            print(f"     theta {th:5.0f}, k_e prior {str(kp):4}: {np.median(v):.4f} [{v.min():.4f}, {v.max():.4f}]")
    print(f"     (t10 printed 0.0617 at 5% and 0.0606 at 20%, both with the k_e prior)")
    k20 = float(np.median(km[(TH20, 0.8)]))
    dmax = 0.0
    for i, x in enumerate(cells):
        if x[0] != 1:
            continue
        des = (Design([(-1, 0.5 * L), (0, 0.0)], T_h=48) if x[3] == 0 else
               Design([(-1, 0.5 * L), (0, 5.5 * L), (1.0, 0.0)], T_h=48))
        R, _ = sens(TH20, des, h_min=1.0, rec_min=(60, 15)[x[1]])
        F = R.T @ R / (0.05, 0.01)[x[2]] ** 2
        dmax = max(dmax, abs(np.log10(sd_th(alone(F, {1: 0.8, 2: k20}))) - np.log10(y2[i])))
    print(f"     using the 20% record's own {k20:.4f} instead of 0.062 in the 8 cells at 20% moves log10 sd by at "
          f"most {dmax:.4f} [t10 0.0074]")
    return y0, y2, out, km


# =========================================================== (b) the 5% grid
EXC = (("drug 450", 0, 0, "drug", 450), ("drug 1500", 0, 1, "drug", 1500), ("drug 2850", 0, 2, "drug", 2850),
       ("withdraw 25%", 0, 3, "out", 0.25), ("withdraw 50%", 0, 4, "out", 0.5),
       ("load x1.5", 1, 0, "up", 1.5), ("load x2", 2, 0, "up", 2), ("load x3", 3, 0, "up", 3),
       ("load x4", 4, 0, "up", 4), ("load x6", 5, 0, "up", 6), ("load x8", 6, 0, "up", 8))
TAUS = (5, 15, 30, 60)
TS = (24, 48, 72)
CVS = (0.05, 0.02, 0.01, 0.005)


def make(kind, arg, T, hold=None):
    if kind == "drug":
        return d_drug(arg, T)
    if kind == "out":
        return d_withdraw(arg, T)
    return d_up(arg, T, hold)


def unit_info(kind, arg, T, tau, kfree=False, cache={}):
    """Unit-CV information of one transient at 5%."""
    Tkey = 72 if kind in ("drug", "out") else T          # no T dependence: subsample a 72 h run
    key = (kind, arg, Tkey, kfree)
    if key not in cache:
        cache[key] = sens(TH5, make(kind, arg, Tkey), h_min=1.0, rec_min=5.0, kfree=kfree)[0]
    R = cache[key]
    step = tau // 5
    n = T * 60 // tau + 1
    Rs = R[::step][:n]
    assert Rs.shape[0] == n
    return Rs.T @ Rs


def part_b():
    banner("(b) THE 5% GRID, RECOMPUTED (same axes as t10; own FIM, own record, own ordering and Pareto)")
    t0 = time.time()
    Z0 = record_info(TH5, 0)
    Z0k = record_info(TH5, 0, kfree=True)
    SD = np.full((len(EXC), len(TS), len(TAUS), len(CVS), 5, 2), np.nan)
    for e, (nm, cls, rk, kind, arg) in enumerate(EXC):
        for a, T in enumerate(TS):
            for b, tau in enumerate(TAUS):
                F1 = unit_info(kind, arg, T, tau)
                for c, cv in enumerate(CVS):
                    F = F1 / cv ** 2
                    for r, (_, pri) in enumerate(RUNGS):
                        SD[e, a, b, c, r, 0] = sd_th(alone(F, pri))
                        SD[e, a, b, c, r, 1] = sd_th(joined(F, pri, Z0))
    print(f"  {SD.size} grid points in {time.time()-t0:.1f} s")
    passing = []
    for idx in zip(*np.nonzero(SD <= GATE)):
        e, a, b, c, r, j = (int(v) for v in idx)
        passing.append(dict(e=e, cls=EXC[e][1], rk=EXC[e][2], T=TS[a], tau=TAUS[b], c=c, r=r, j=j,
                            kc=KCOST[r], N=TS[a] * 60 // TAUS[b] + 1, sd=float(SD[idx])))
    print(f"  passing points: {len(passing)} of {SD.size} (t10 printed {T10['n_pass']})")
    print(f"  passing by safety class x CV (mine | t10):")
    for cls in range(7):
        mine = [sum(1 for d in passing if d['cls'] == cls and d['c'] == c) for c in range(4)]
        print(f"     class {cls}: " + "  ".join(f"{m:4d} | {t:4d}" for m, t in zip(mine, T10['cls_cv'][cls])))
    # key single values
    E = {nm: e for e, (nm, *_) in enumerate(EXC)}
    ia = {T: a for a, T in enumerate(TS)}
    ib = {t: b for b, t in enumerate(TAUS)}
    ic = {cv: c for c, cv in enumerate(CVS)}
    def g(nm, T, tau, cv, r, j):
        return SD[E[nm], ia[T], ib[tau], ic[cv], r, j]
    print("\n  key values (mine; t10's printed value in brackets):")
    keyv = [("x4 hourly 48 h, L0 1% alone, CV 5%", g("load x4", 48, 60, 0.05, 4, 0), 0.2054),
            ("x6 hourly 48 h, L0 1% alone, CV 5%", g("load x6", 48, 60, 0.05, 4, 0), 0.1757),
            ("x4 hourly 72 h, L0 1% alone, CV 5%", g("load x4", 72, 60, 0.05, 4, 0), 0.182),
            ("x4 30 min 48 h, L0 1% alone, CV 5%", g("load x4", 48, 30, 0.05, 4, 0), 0.174),
            ("drug 2850 15m/72h CV 0.5% none alone", g("drug 2850", 72, 15, 0.005, 0, 0), 0.1942),
            ("drug 1500 5m/72h CV 0.5% none alone", g("drug 1500", 72, 5, 0.005, 0, 0), 0.1235),
            ("withdraw 50% 5m/48h CV 0.5% L0 1% alone", g("withdraw 50%", 48, 5, 0.005, 4, 0), 0.1950),
            ("x2 15m/72h CV 5% L0 1% joined", g("load x2", 72, 15, 0.05, 4, 1), 0.1989),
            ("x3 5m/72h CV 5% +K_m alone", g("load x3", 72, 5, 0.05, 2, 0), 0.1933),
            ("x4 5m/48h CV 5% none joined", g("load x4", 48, 5, 0.05, 0, 1), 0.1987),
            ("x6 15m/72h CV 5% +K_m alone", g("load x6", 72, 15, 0.05, 2, 0), 0.1977),
            ("x8 15m/72h CV 5% +K_m alone", g("load x8", 72, 15, 0.05, 2, 0), 0.1725),
            ("x1.5 5m/72h CV 2% none joined", g("load x1.5", 72, 5, 0.02, 0, 1), 0.1998),
            ("drug 2850 15m/72h CV 5% L0 1% alone", g("drug 2850", 72, 15, 0.05, 4, 0), 1.481),
            ("drug 2850 15m/72h CV 5% L0 1% joined", g("drug 2850", 72, 15, 0.05, 4, 1), 1.260)]
    for lab, v, t in keyv:
        print(f"     {lab:44} {v:8.4f}  [{t}]")
    # class 0 at CV 5/2/1%: the best value anywhere, and at L0 1% in particular
    print("\n  class 0 (no added substrate): the best sd anywhere in the grid, per CV, over all rungs and modes")
    for e in range(5):
        print(f"     {EXC[e][0]:14}" + "".join(f"  CV {100*cv:g}%: {np.nanmin(SD[e, :, :, c]):.3f}"
                                         for c, cv in enumerate(CVS)))
    best0 = min(np.nanmin(SD[e, :, :, c]) for e in range(5) for c in range(3))
    print(f"     the best class-0 sd at CV 5, 2 or 1% is {best0:.4f}")
    # ordering (i) and (ii), and the Pareto set, by my own code
    k1 = lambda d: (d["cls"], d["kc"], d["c"], d["N"], d["T"], 5 * d["j"] + d["r"], d["rk"])
    k2 = lambda d: (d["c"], d["cls"], d["kc"], d["N"], d["T"], 5 * d["j"] + d["r"], d["rk"])
    def desc(d):
        return (f"{EXC[d['e']][0]:13} {d['tau']:2d} min/{d['T']} h (N {d['N']:3d}), CV {100*CVS[d['c']]:g}%, "
                f"{RUNGS[d['r']][0]} {('alone', 'joined')[d['j']]}: {d['sd']:.4f}")
    for lab, key in (("(i) safety > knowledge > CV > N", k1), ("(ii) CV > safety > knowledge > N", k2)):
        srt = sorted(passing, key=key)
        print(f"\n  ordering {lab}: first three distinct excursions")
        seen = []
        for d in srt:
            if d["e"] not in seen:
                seen.append(d["e"])
                print("     " + desc(d))
            if len(seen) == 3:
                break
    vecs = {}
    for d in passing:
        v = (d["cls"], d["kc"], d["c"], d["N"])
        vecs[v] = vecs.get(v, d)
    V = list(vecs)
    nd = [v for v in V if not any(w != v and all(wi <= vi for wi, vi in zip(w, v)) for w in V)]
    print(f"\n  Pareto set over (class, knowledge, CV, N): {len(nd)} non-dominated of {len(V)} distinct vectors "
          f"(t10 printed {T10['pareto'][0]} of {T10['pareto'][1]}); class-0 members: "
          + ", ".join(f"{v}" for v in sorted(nd) if v[0] == 0))
    print("\n  load x1.5, fewest samples reaching 0.20 with free knowledge (rungs none/k_e/K_m, alone or joined):")
    for c, cv in enumerate(CVS):
        cand = [(d["N"], d["T"], d["sd"], d["tau"]) for d in passing
                if d["e"] == E["load x1.5"] and d["c"] == c and d["kc"] == 0]
        if cand:
            n, T, v, tau = min(cand)
            print(f"     CV {100*cv:g}%: {n} ({tau} min/{T} h) sd {v:.4f}")
        else:
            print(f"     CV {100*cv:g}%: none; best {np.nanmin(SD[E['load x1.5'], :, :, c, :3, :]):.4f}")
    print("     [t10: none (0.442); 865 5m/72h 0.200; 289 15m/72h 0.183; 73 60m/72h 0.180]")
    # kappa
    print("\n  drug 2850 with kappa uncertain (15 min / 72 h, CV 0.5%, no priors, alone):")
    F = unit_info("drug", 2850, 72, 15, kfree=True) / 0.005 ** 2
    print(f"     kappa to 12%: {sd_th(alone(F, {4: 0.12})):.4f};  kappa free: {sd_th(F):.4f}  [t10 0.3450, 0.3701]")
    best = None
    for T in TS:
        for tau in TAUS:
            Fk = unit_info("drug", 2850, T, tau, kfree=True) / 0.005 ** 2
            for r, (_, pri) in enumerate(RUNGS[:3]):
                for j in range(2):
                    q = {**pri, 4: 0.12}
                    v = sd_th(alone(Fk, q)) if j == 0 else sd_th(joined(Fk, q, Z0k))
                    if v <= GATE:
                        cand = (T * 60 // tau + 1, T, v, tau)
                        best = cand if best is None or cand < best else best
    print(f"     kappa 12%, free knowledge, CV 0.5%: fewest samples {best[0]} ({best[3]} min/{best[1]} h) sd "
          f"{best[2]:.4f}  [t10 865 5m/72h 0.182]")
    return SD, passing, Z0


def part_b_probes(SD):
    banner("(b') PROBES t10 DID NOT RUN: the load hold confound, record-seed dependence, the K_m prior")
    print("  In t10's grid (as in t7_design.washout) the load is held for the FIRST HALF of the experiment,")
    print("  so T = 24/48/72 h holds it 12/24/36 h. Separating sampling duration from load exposure:")
    Z0 = record_info(TH5, 0)
    rows = (("x4 hourly, L0 1% alone, CV 5%", 4, 60, 0.05, 4, 0),
            ("x2 15 min, L0 1% joined, CV 5%", 2, 15, 0.05, 4, 1),
            ("x3 5 min, +K_m alone, CV 5%", 3, 5, 0.05, 2, 0),
            ("x6 15 min, +K_m alone, CV 5%", 6, 15, 0.05, 2, 0),
            ("x8 15 min, +K_m alone, CV 5%", 8, 15, 0.05, 2, 0))
    print(f"     {'design':34} {'T 48 h hold 24 h':>17} {'T 72 h hold 36 h':>17} {'T 72 h hold 24 h':>17}")
    res = {}
    for lab, x, tau, cv, r, j in rows:
        v = []
        for T, hold in ((48, 24), (72, 36), (72, 24)):
            R, _ = sens(TH5, d_up(x, T, hold), h_min=1.0, rec_min=tau)
            F = R.T @ R / cv ** 2
            pri = RUNGS[r][1]
            v.append(sd_th(alone(F, pri)) if j == 0 else sd_th(joined(F, pri, Z0)))
        res[lab] = v
        print(f"     {lab:34} {v[0]:17.4f} {v[1]:17.4f} {v[2]:17.4f}")
    print("\n  The grid's load excursions with the hold CAPPED at 24 h (T = 24 h keeps 12 h), the same exposure")
    print("  at 48 and 72 h: passing points at CV 5% by class, and ordering (ii)'s first three at CV 5%")
    passing = []
    for e, (nm, cls, rk, kind, arg) in enumerate(EXC):
        if kind != "up":
            continue
        for T in TS:
            hold = min(T / 2, 24)
            Rf, _ = sens(TH5, d_up(arg, T, hold), h_min=1.0, rec_min=5.0)
            for tau in TAUS:
                Rs = Rf[::tau // 5]
                F = Rs.T @ Rs / 0.05 ** 2
                for r, (_, pri) in enumerate(RUNGS):
                    for j in range(2):
                        v = sd_th(alone(F, pri)) if j == 0 else sd_th(joined(F, pri, Z0))
                        if v <= GATE:
                            passing.append((cls, KCOST[r], T * 60 // tau + 1, T, 5 * j + r, nm, tau, v,
                                            RUNGS[r][0], ("alone", "joined")[j]))
    cnt = [sum(1 for d in passing if d[0] == c) for c in range(1, 7)]
    print("     passing at CV 5%, classes 1..6 (hold capped): " + ", ".join(f"{c}" for c in cnt)
          + "   [t10, hold T/2: 0, 5, 20, 27, 37, 37]")
    seen = []
    for d in sorted(passing):
        if d[5] not in seen:
            seen.append(d[5])
            print(f"     {len(seen)}. {d[5]:10} {d[6]:2d} min/{d[3]} h (N {d[2]}), {d[8]} {d[9]}: {d[7]:.4f}")
        if len(seen) == 3:
            break

    print("\n  The joined designs use ONE record (intake seed 0). The borderline joined passes over seeds 0..23:")
    Zs = [record_info(TH5, s) for s in range(24)]
    probes = (("x2 15m/72h CV 5% L0 1% joined", "up", 2, 72, 15, 0.05, 4),
              ("x4 5m/48h CV 5% none joined", "up", 4, 48, 5, 0.05, 0),
              ("x1.5 5m/72h CV 2% none joined", "up", 1.5, 72, 5, 0.02, 0),
              ("x3 60m/72h CV 5% L0 1% joined", "up", 3, 72, 60, 0.05, 4),
              ("x4 60m/48h CV 5% L0 1% joined", "up", 4, 48, 60, 0.05, 4))
    for lab, kind, arg, T, tau, cv, r in probes:
        F = unit_info(kind, arg, T, tau) / cv ** 2
        v = np.array([sd_th(joined(F, RUNGS[r][1], Z)) for Z in Zs])
        print(f"     {lab:34} seed 0 {v[0]:.4f}; median {np.median(v):.4f} [{v.min():.4f}, {v.max():.4f}];"
              f" pass in {int(np.sum(v <= GATE))}/24 seeds")
    print("\n  The '+K_m from the record (0.062)' rung, used alone, adds the record's MARGINAL sd(ln K_m),")
    print("  computed with the k_e prior in force, and then applies the k_e prior again. Against the proper")
    print("  joint analysis (transient joined to the record, k_e prior once), for the alone passes at 5% CV:")
    for lab, kind, arg, T, tau in (("x3 5m/72h", "up", 3, 72, 5), ("x6 15m/72h", "up", 6, 72, 15),
                                   ("x8 15m/72h", "up", 8, 72, 15)):
        F = unit_info(kind, arg, T, tau) / 0.05 ** 2
        a = sd_th(alone(F, RUNGS[2][1]))
        b = sd_th(joined(F, RUNGS[1][1], Z0))
        print(f"     {lab:12} alone with K_m 0.062 + k_e 0.8: {a:.4f};  joined, k_e 0.8 once: {b:.4f}")
    return res


# =========================================================== (c) Monte Carlo, own estimator
LO = np.array([np.log(1e-2), np.log(1e-3), np.log(1.0), -1e4, -1e4])
HI = np.array([np.log(1e5), np.log(1e3), np.log(1e5), 1e5, 1e4])


class Problem:
    """Batched residuals and complex-step Jacobian for R replicates at once."""

    def __init__(self, des, theta, tau, cv, pri, rec_seed, h_min):
        self.des, self.tau, self.cv, self.pri, self.h = des, tau, cv, dict(pri), h_min
        self.joined = rec_seed is not None
        self.rec_seed = rec_seed
        self.np_ = 5 if self.joined else 4
        self.truth = np.array([np.log(theta), np.log(KE), np.log(KM), L0, 0.0])[:self.np_]
        if self.joined:
            self.pri.pop(2, None)
            self.rm, sig = record(theta, rec_seed)
            self.rsig = np.full(self.rm.size, float(sig))

    def model(self, X):
        """ln m at the samples for parameter rows X (R, >=4), any dtype."""
        m = path(np.exp(X[:, 0]), np.exp(X[:, 1]), np.exp(X[:, 2]), X[:, 3], np.ones(X.shape[0]),
                 self.des, self.h, self.tau)
        return m

    def rec_pred(self, X):
        th, ke, K = np.exp(X[:, 0]), np.exp(X[:, 1]), np.exp(X[:, 2])
        rm = self.rm[None, :]
        return (th[:, None] + REC_INC[None, :]) * rm / (K[:, None] + rm) + ke[:, None] * rm + X[:, 4:5]

    def resid(self, X, data):
        with np.errstate(all="ignore"):
            m = self.model(X).T                                   # (R, N)
            ok = np.all(np.isfinite(m) & (m > 0), axis=1)
            r = [(np.log(np.where(m > 0, m, 1.0)) - data["ly"]) / self.cv]
        for j, s in self.pri.items():
            r.append(((X[:, j] - data["cen"][j]) / s)[:, None])
        if self.joined:
            r.append((self.rec_pred(X) - data["yrec"]) / self.rsig[None, :])
            r.append(((X[:, 4] - data["cb"]) / (0.05 * L))[:, None])
        r = np.concatenate(r, axis=1)
        r[~ok] = np.nan
        return r

    def jac(self, X):
        R = X.shape[0]
        Xc = np.repeat(X[:, None, :4].astype(complex), 4, axis=1)   # (R, 4, 4)
        Xc[:, np.arange(4), np.arange(4)] += 1j * EPS
        Xc = Xc.reshape(R * 4, 4)
        with np.errstate(all="ignore"):
            m = self.model(Xc)                                    # (N, 4R)
            N = m.shape[0]
            Jd = ((m.imag / EPS) / m.real).T.reshape(R, 4, N).transpose(0, 2, 1) / self.cv   # (R, N, 4)
        blocks = [np.concatenate([Jd, np.zeros((R, N, self.np_ - 4))], axis=2)]
        for j, s in self.pri.items():
            P = np.zeros((R, 1, self.np_))
            P[:, 0, j] = 1 / s
            blocks.append(P)
        if self.joined:
            th, ke, K = np.exp(X[:, 0]), np.exp(X[:, 1]), np.exp(X[:, 2])
            rm = self.rm[None, :]
            Jr = np.zeros((R, rm.shape[1], 5))
            Jr[:, :, 0] = th[:, None] * rm / (K[:, None] + rm)
            Jr[:, :, 1] = ke[:, None] * rm
            Jr[:, :, 2] = -(th[:, None] + REC_INC[None, :]) * K[:, None] * rm / (K[:, None] + rm) ** 2
            Jr[:, :, 4] = 1.0
            blocks.append(Jr / self.rsig[None, :, None])
            B = np.zeros((R, 1, 5))
            B[:, 0, 4] = 1 / (0.05 * L)
            blocks.append(B)
        return np.concatenate(blocks, axis=1)

    def fim(self):
        """The information at the truth, at this integrator step."""
        Xc = np.repeat(self.truth[None, :4].astype(complex), 4, axis=0)
        Xc[np.arange(4), np.arange(4)] += 1j * EPS
        m = self.model(Xc)
        Rd = (m.imag / EPS) / m.real
        F = Rd.T @ Rd / self.cv ** 2
        if not self.joined:
            return alone(F, self.pri)
        return joined(F, self.pri, record_info(np.exp(self.truth[0]), self.rec_seed))


def draw(prob, R, rng):
    m = prob.model(prob.truth[None, :]).real[:, 0]
    data = dict(ly=np.log(m)[None, :] + rng.normal(0, prob.cv, (R, m.size)))
    data["cen"] = {j: prob.truth[j] + s * rng.normal(size=R) for j, s in prob.pri.items()}
    if prob.joined:
        yr = prob.rec_pred(prob.truth[None, :])[0]
        data["yrec"] = yr[None, :] + rng.normal(size=(R, yr.size)) * prob.rsig[None, :]
        data["cb"] = 0.05 * L * rng.normal(size=R)
    return data


def draw_t10(prob, seeds):
    """t10's replicate draws, regenerated from its seed convention (one default_rng per replicate):
    transient noise, then each prior centre in rung order, then (joined) the record noise and the
    load-bias centre. Only the convention is shared; the data are made by this file's integrator."""
    m = prob.model(prob.truth[None, :]).real[:, 0]
    ly, cen, yrec, cb = [], {j: [] for j in prob.pri}, [], []
    yr = prob.rec_pred(prob.truth[None, :])[0] if prob.joined else None
    for s in seeds:
        rng = np.random.default_rng(s)
        ly.append(np.log(m) + rng.normal(0, prob.cv, m.size))
        for j, sd in prob.pri.items():
            cen[j].append(prob.truth[j] + sd * rng.normal())
        if prob.joined:
            yrec.append(yr + rng.normal(0, prob.rsig))
            cb.append(0.05 * L * rng.normal())
    data = dict(ly=np.array(ly), cen={j: np.array(v) for j, v in cen.items()})
    if prob.joined:
        data["yrec"], data["cb"] = np.array(yrec), np.array(cb)
    return data


def lm(prob, data, X0, iters=80):
    """Batched Levenberg-Marquardt with box clipping; returns x, cost, final GN step in ln theta."""
    X = np.clip(X0, LO[:prob.np_] + 1e-9, HI[:prob.np_] - 1e-9)
    r = prob.resid(X, data)
    cost = 0.5 * np.sum(np.nan_to_num(r) ** 2, axis=1)
    cost[~np.all(np.isfinite(r), axis=1)] = np.inf
    J = np.nan_to_num(prob.jac(X))
    mu = np.full(X.shape[0], 1e-3)
    I = np.eye(prob.np_)[None]
    idle = 0
    for _ in range(iters):
        live = np.isfinite(cost)
        JTJ = np.einsum("rnp,rnq->rpq", J, J)
        gr = np.einsum("rnp,rn->rp", J, np.nan_to_num(r))
        Dg = np.einsum("rpp->rp", JTJ)
        A = JTJ + mu[:, None, None] * (Dg[:, :, None] * I + 1e-12 * Dg.max(axis=1)[:, None, None] * I)
        with np.errstate(all="ignore"):
            dX = -(np.linalg.pinv(A) @ gr[:, :, None])[:, :, 0]
        Xt = np.clip(X + np.nan_to_num(dX), LO[:prob.np_] + 1e-9, HI[:prob.np_] - 1e-9)
        rt = prob.resid(Xt, data)
        ct = 0.5 * np.sum(np.nan_to_num(rt) ** 2, axis=1)
        ct[~np.all(np.isfinite(rt), axis=1)] = np.inf
        acc = (ct < cost) & live
        moved = acc & (np.abs(Xt - X).max(axis=1) > 1e-10)
        mu = np.where(acc, np.maximum(mu / 3, 1e-12), np.minimum(mu * 4, 1e12))
        if acc.any():
            X[acc], cost[acc], r[acc] = Xt[acc], ct[acc], rt[acc]
            J[acc] = np.nan_to_num(prob.jac(X[acc]))
        idle = 0 if moved.any() else idle + 1
        if idle >= 8:
            break
    JTJ = np.einsum("rnp,rnq->rpq", J, J)
    gr = np.einsum("rnp,rn->rp", J, np.nan_to_num(r))
    with np.errstate(all="ignore"):
        gn = -(np.linalg.pinv(JTJ) @ gr[:, :, None])[:, :, 0]
    return X, cost, np.abs(gn[:, 0])


def chi_band(R):
    """5-95% band of RMSE/sd for R draws of an exactly efficient, unbiased Gaussian estimator."""
    from scipy.stats import chi2
    return float(np.sqrt(chi2.ppf(0.05, R) / R)), float(np.sqrt(chi2.ppf(0.95, R) / R))


def mc(lab, des, theta, tau, cv, pri, rec_seed, R, seed, h_min=5.0, t10=None, t10_seeds=None):
    prob = Problem(des, theta, tau, cv, pri, rec_seed, h_min)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    # the FIM at this integrator step, against the 1-min step
    Ff = prob.fim()
    prob1 = Problem(des, theta, tau, cv, pri, rec_seed, 1.0)
    f_sd, f_sd1 = sd_th(Ff), sd_th(prob1.fim())
    if t10_seeds is not None:
        data = draw_t10(prob, t10_seeds)
        R = len(t10_seeds)
    else:
        data = draw(prob, R, rng)
    tr = prob.truth
    dp = np.array([0.7, -0.3, 0.5, -0.2 * L0, 0.0])[:prob.np_]
    starts = {"truth": np.tile(tr, (R, 1)), "truth+dp": np.tile(tr + dp, (R, 1)),
              "truth-dp": np.tile(tr - dp, (R, 1))}
    rs = tr + rng.normal(size=(R, prob.np_)) * np.array([1.0, 0.5, 0.7, 0.3 * L0, 5.0])[:prob.np_]
    starts["random"] = rs
    if prob.pri:
        pc = np.tile(tr, (R, 1))
        for j in prob.pri:
            pc[:, j] = data["cen"][j]
        starts["prior centres"] = pc
    fits = {k: lm(prob, data, X0) for k, X0 in starts.items()}
    keys = list(fits)
    C_ = np.array([fits[k][1] for k in keys])                      # (S, R)
    def pick(names):
        ii = [keys.index(k) for k in names]
        b = np.argmin(C_[ii], axis=0)
        x = np.array([fits[keys[ii[bb]]][0][r, 0] for r, bb in enumerate(b)])
        gn = np.array([fits[keys[ii[bb]]][2][r] for r, bb in enumerate(b)])
        return x - tr[0], gn
    e_all, gn_all = pick(keys)
    e_not, _ = pick([k for k in keys if k != "truth"])
    e_tru, _ = pick(["truth"])
    def stats(e):
        e = e[np.isfinite(e)]
        rm = float(np.sqrt(np.mean(e ** 2)))
        bs = np.random.default_rng(seed + 1).choice(e, (4000, e.size))
        lo, hi = np.percentile(np.sqrt(np.mean(bs ** 2, axis=1)), [5, 95])
        rob = 0.5 * (np.percentile(e, 84) - np.percentile(e, 16))
        return rm, lo, hi, rob, float(np.median(e)), e.size
    sa, sn, st = stats(e_all), stats(e_not), stats(e_tru)
    conv = float(np.mean(gn_all < 1e-3))
    agree = float(np.mean(np.abs(e_all - e_not) < 1e-4))
    print(f"\n  {lab}")
    print(f"     FIM sd at h = {h_min:g} min {f_sd:.4f} (h = 1 min {f_sd1:.4f}); {R} replicates, seed {seed}, "
          f"{time.time()-t0:.0f} s")
    print(f"     best of all starts:       RMSE {sa[0]:.4f} [90% boot {sa[1]:.4f}, {sa[2]:.4f}]  robust {sa[3]:.4f}"
          f"  median bias {sa[4]:+.4f}  MC/FIM {sa[0]/f_sd:.3f}")
    print(f"     best excluding the truth: RMSE {sn[0]:.4f} [{sn[1]:.4f}, {sn[2]:.4f}]; same optimum as "
          f"all-starts in {agree:.0%} of replicates")
    print(f"     from the truth only:      RMSE {st[0]:.4f}")
    print(f"     remaining Gauss-Newton step in ln theta < 1e-3 at the chosen optimum: {conv:.0%} of replicates")
    cb_ = chi_band(sa[5])
    print(f"     an efficient Gaussian estimator would give MC/FIM in [{cb_[0]:.2f}, {cb_[1]:.2f}] (5-95%) at "
          f"{sa[5]} replicates")
    if t10:
        print(f"     t10 printed: {t10}")
    return dict(fim=f_sd, rmse=sa[0], lo=sa[1], hi=sa[2], robust=sa[3], rn=sn[0], conv=conv, n=sa[5])


def part_c():
    banner("(c) MONTE CARLO WITH AN INDEPENDENT ESTIMATOR (batched LM, complex-step Jacobian, new seeds)")
    print("  Priors are pseudo-observations with centres drawn at truth + sd N(0,1) per replicate; joined designs")
    print("  redraw the record's ordinate noise and its load-bias centre (record intake seed 0, as t10/t7).")
    print("  Starts: the truth, truth +- t6's offsets, one random start, and the prior centres when there are")
    print("  priors; the lowest cost wins. Integrator step 5 min (data and fit), checked against 1 min via the FIM.")
    D0 = (("[0] drug 2850, 15 min / 72 h, CV 0.5%, no priors, alone", d_drug(2850, 72), TH5, 15, 0.005, {}, None,
           "FIM 0.1942, MC RMSE 0.1686 [0.144, 0.193], MC/FIM 0.87 (40 reps)"),
          ("[1] drug 1500, 5 min / 72 h, CV 0.5%, no priors, alone", d_drug(1500, 72), TH5, 5, 0.005, {}, None,
           "FIM 0.1235, MC RMSE 0.1127 [0.094, 0.130], MC/FIM 0.91 (40 reps)"),
          ("[2] withdraw 50%, 5 min / 48 h, CV 0.5%, k_e 0.8 + K_m 0.062 + L0 1%, alone", d_withdraw(0.5, 48),
           TH5, 5, 0.005, RUNGS[4][1], None, "FIM 0.1950, MC RMSE 0.2370 [0.198, 0.273], MC/FIM 1.22 (40 reps)"),
          ("[3] load x2, 15 min / 72 h, CV 5%, k_e 0.8 + L0 1%, joined to the record", d_up(2, 72), TH5, 15,
           0.05, RUNGS[4][1], 0, "FIM 0.1989, MC RMSE 0.2369 [0.181, 0.291], MC/FIM 1.19 (40 reps)"),
          ("[4] load x3, 5 min / 72 h, CV 5%, k_e 0.8 + K_m 0.062, alone", d_up(3, 72), TH5, 5, 0.05,
           RUNGS[2][1], None, "FIM 0.1933, MC RMSE 0.2033 [0.167, 0.238], MC/FIM 1.05 (40 reps)"))
    reps = (300, 150, 150, 300, 150)
    out = {}
    print("\n  (c1) new draws, more replicates than t10's 40")
    for i, ((lab, des, th, tau, cv, pri, rs, t10), R) in enumerate(zip(D0, reps)):
        out[i] = mc(lab, des, th, tau, cv, pri, rs, R, 7001 + i, t10=t10)
    print("\n  (c2) t10's own 40 replicate draws (seeds 52000 + 1000 id + rep, its draw order) refitted with")
    print("       this file's estimator, to separate estimator behaviour from the luck of 40 draws")
    out2 = {}
    for i, (lab, des, th, tau, cv, pri, rs, t10) in enumerate(D0):
        out2[i] = mc(lab + "  [t10's draws]", des, th, tau, cv, pri, rs, 40, 0, t10=t10,
                     t10_seeds=[52_000 + 1000 * i + r for r in range(40)])
    print("\n  (c3) TB-B and its reverse cell 'sampling at A', no priors (t10's part (a) Monte Carlo)")
    dB = Design([(-1, 0.5 * L), (0, 5.5 * L), (1.0, 0.0)], T_h=48)
    out["B"] = mc("TB-B (20%, 15 min, CV 1%, x6), new draws", dB, TH20, 15, 0.01, {}, None, 200, 7006,
                  t10="FIM 0.089, MC RMSE 0.100 [0.080, 0.119] (30 reps); t6 printed NLS RMSE 0.079 (80 reps)")
    out["Bt"] = mc("TB-B, t10's own 30 draws (seeds 36000 + rep)", dB, TH20, 15, 0.01, {}, None, 30, 0,
                   t10="MC RMSE 0.100", t10_seeds=[31_000 + 1000 * 5 + r for r in range(30)])
    out["Bs"] = mc("TB-B with hourly sampling (sampling at A), new draws", dB, TH20, 60, 0.01, {}, None, 200, 7007,
                   t10="FIM 0.173, MC RMSE 0.205 [0.154, 0.252] (20 reps)")
    print("\n  (c4) 'B, excursion at A' (20%, 15 min, CV 1%, diet out; FIM sd 29.7): t10's MC RMSE 1.303 (20 reps).")
    print("       The same t10 draws (seeds 40000 + rep) fitted here, converged, to show whether 1.303 measures")
    print("       the data or the estimator:")
    dBe = Design([(-1, 0.5 * L), (0, 0.0)], T_h=48)
    out["Be"] = mc("B, excursion at A, t10's 20 draws", dBe, TH20, 15, 0.01, {}, None, 20, 0,
                   t10="FIM 29.666, MC RMSE 1.303 [1.194, 1.405], off >2x 100%",
                   t10_seeds=[31_000 + 1000 * 9 + r for r in range(20)])
    print("\n  SUMMARY (c): FIM | t10 MC (40) | this file, new draws (n) | this file on t10's draws")
    t10mc = (0.1686, 0.1127, 0.2370, 0.2369, 0.2033)
    for i in range(5):
        a, b = out[i], out2[i]
        print(f"     [{i}] {a['fim']:.4f} | {t10mc[i]:.4f} | {a['rmse']:.4f} [{a['lo']:.4f}, {a['hi']:.4f}] "
              f"MC/FIM {a['rmse']/a['fim']:.3f} (n {a['n']}) | {b['rmse']:.4f}")
    return out, out2


def main():
    print("=" * 100)
    print("  RECHECK OF t10_factorial.py BY AN INDEPENDENT ROUTE (no shared code)")
    print("=" * 100)
    part_a()
    SD, passing, Z0 = part_b()
    part_b_probes(SD)
    part_c()
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
