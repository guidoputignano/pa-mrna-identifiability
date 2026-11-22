"""
Adversarial recheck of t8_floor.py (follow-up item 1, the 0.39 "floor").

Shares no code with t8_floor.py, core.py, dyn.py, t*.py, n4_theta_scope.py or
design.run: every piece of the footing is rebuilt here from the stated model.

    steady state   (theta + Delta) m/(K_m + m) + k_e m = load, solved by
                   vectorised BISECTION (the repository uses the quadratic root)
    dose train     direct sum of decaying doses D exp(-lambda (t - 14 k)) for
                   14 k <= t, summed over the last 80 doses (the repository
                   sums all 200 from k = 0; the investigator uses a geometric
                   closed form); optionally capped at 200 doses
    sigma          hypot(scatter * L, C'(m_bar) * assay * m_bar), with C'(m_bar)
                   by COMPLEX STEP, not written down
    sensitivities  complex-step derivatives of the predicted recorded load
                   (theta + Delta) g(m) + k_e m + b in (ln theta, ln K_m, b,
                   ln k_e), not written down
    Var(ln theta)  Frisch-Waugh: 1 / ||M_Z x_0||^2 on the whitened design with
                   the priors appended as pseudo-rows (no matrix inverse, no
                   coordinate deletion: an exact nuisance is simply a column
                   left out)
    zero noise     computed at sigma scaled by 1e-8, not asserted
    Shapley        by explicit enumeration of the 24 orders AND from the Mobius
                   dividends (sum_{T contains i} d(T)/|T|); both are printed
    asymptotics    stationary expected information per draw by quadrature over
                   the dose phase lattice and Gauss-Hermite over the load
    Monte Carlo    MAP fits (multi-start, bounded trust region) on simulated
                   records, b drawn from its prior, k_e fixed at the truth

The 24 intake records are rebuilt with the repository's convention
(default_rng(seed), lognormal CV 10%) so that seed-for-seed comparison is
possible; a second, disjoint seed set (1000-1023) tests whether the medians
depend on which 24 records were drawn.

The investigator's printed values are carried in CLAIM below only so that the
comparison table can print them beside the recomputed values. They are the
claims under test, not expected results.
"""
import itertools
import time

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import brentq, least_squares

T_START = time.time()

# ------------------------------------------------------------------ footing
TH_N, KM, KE, M_N = 3000.0, 100.0, 5.0, 10.0
L_REC = TH_N * M_N / (KM + M_N) + KE * M_N           # convention A
DOSE, TAU, HALF = 6000.0, 14.0, 3.5
LAM = np.log(2.0) / HALF
LOAD_CV, LF, ACV, KE_SD = 0.10, 0.05, 0.05, 0.8
SEEDS = list(range(24))
ALT_SEEDS = list(range(1000, 1024))
GATE = 0.20
TH5 = 0.05 * TH_N
EDGE = 2472.0
MONTH = 30.4

CLAIM = {  # the investigator's printed values (out/t8_floor.txt), under test
    "sd|ke 5% 25m": 0.3879, "Var|ke 5%": 0.15047, "sd 5% 25m": 3.0673,
    "scatter term": 16.136, "assay term": 15.214, "sigma": 22.177,
    "S{b}": 0.2979, "S{K_m}": 0.3786, "S{scatter}": 0.3139, "S{assay}": 0.3237,
    "S{b,K_m}": 0.1943, "S{b,scatter}": 0.2044, "S{b,K_m,scatter}": 0.1333,
    "S{b,K_m,assay}": 0.1413,
    "Sh b": 0.292, "Sh K_m": 0.093, "Sh scatter": 0.325, "Sh assay": 0.291,
    "first b": 0.410, "first K_m": 0.047, "max b": 0.702, "max K_m": 0.339,
    "D{b,scatter}": 0.425, "D{b,assay}": 0.385, "D{all}": -0.243,
    "Sh100 b": 0.109, "Sh100 K_m": 0.266,
    "cont 100": 0.2767, "cont 600": 0.1624, "cont 6000": 0.0587,
    "grid 100": 0.2633, "grid 600": 0.1762,
    "dense 365": 0.1947, "dense 730": 0.1507, "dense 100": 0.2728,
    "bex 25": 0.2979, "bex 6000": 0.0189, "loadex 25": 0.2044,
    "kmex 6000": 0.0282,
    "floor": 0.1889, "c": 1.7323, "pred 6000": 0.1902,
    "slope 25-50": -0.222, "slope 400-600": -0.330, "slope 2400-6000": -0.475,
    "nVar bex 6000": 2.153, "c_inf ke": 4.5431, "c_inf bex": 1.4674,
    "c_inf loadex": 1.0066,
    "bratio 25": 0.9588, "bratio 365": 0.6379, "bratio 6000": 0.2034,
    "n monthly ke": 328, "n monthly bex": 58, "n monthly loadex": 29,
    "n monthly kmex": 114, "n monthly kmbex": 24,
    "first 2y ke": 320, "stays 2y ke": 324,
    "c 365 ke1": 0.1998, "c 730 ke2": 0.1762, "c 730 ke5 best": 0.2095,
    "c 365 ke2 b2": 0.1561, "c 100 ke2 b2 sc2": 0.1746,
    "c 25 kmex ke2 b1 sc1": 0.1964,
    "p 365 now": 0.0102, "p 730 now": 0.0291, "b need 100 keex": 0.0262,
    "share 5%": 0.9838, "share edge": 0.8425, "share 100%": 0.7724,
    "share 1000%": 0.0022, "share 1%": 0.9840, "share 50.12%": 0.9352,
}

COMPARE = []


def cmp(key, val, fmt="{:.4f}"):
    """Record a recomputed value beside the investigator's claim."""
    COMPARE.append((key, CLAIM.get(key), val, fmt))


# ------------------------------------------------------------ the record
def dose_inc(t, cap=None, keep=80):
    """Direct sum over the last `keep` doses with 14 k <= t (k < cap if capped)."""
    t = np.asarray(t, float)
    kmax = np.floor(t / TAU).astype(int) + 1
    out = np.zeros_like(t)
    for j in range(keep + 2):
        k = kmax - j
        on = (k >= 0) & (k * TAU <= t)
        if cap is not None:
            on &= k < cap
        out[on] += DOSE * np.exp(-LAM * (t[on] - k[on] * TAU))
    return out


def clearance(m, theta, K=KM, ke=KE, inc=0.0):
    return (theta + inc) * m / (K + m) + ke * m


def m_solve(theta, load, inc, K=KM, ke=KE):
    lo = np.zeros_like(load)
    hi = load / ke
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f = clearance(mid, theta, K, ke, inc) - load
        hi = np.where(f > 0, mid, hi)
        lo = np.where(f > 0, lo, mid)
    return 0.5 * (lo + hi)


def m_bar_of(theta):
    return brentq(lambda x: clearance(x, theta) - L_REC, 0.0, L_REC / KE, xtol=1e-14, rtol=1e-15)


def cstep(f, x, h=1e-30):
    return np.imag(f(x + 1j * h)) / h


class Rec:
    def __init__(self, theta, times, seed, cap=None):
        times = np.asarray(times, float)
        rng = np.random.default_rng(seed)
        s2 = np.log1p(LOAD_CV ** 2)
        self.loads = L_REC * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
        self.inc = dose_inc(times, cap)
        self.theta, self.times = theta, times
        self.m = m_solve(theta, self.loads, self.inc)
        self.mbar = m_bar_of(theta)
        self._J = None

    def J(self):
        """d(predicted recorded load)/d(ln th, ln K_m, b, ln k_e), complex step."""
        if self._J is None:
            p0 = np.array([np.log(self.theta), np.log(KM), 0.0, np.log(KE)], complex)
            J = np.empty((self.m.size, 4))
            for j in range(4):
                p = p0.copy()
                p[j] += 1e-30j
                th, K, b, ke = np.exp(p[0]), np.exp(p[1]), p[2], np.exp(p[3])
                pred = (th + self.inc) * self.m / (K + self.m) + ke * self.m + b
                J[:, j] = pred.imag / 1e-30
            self._J = J
        return self._J

    def sigma(self, sc=LF, acv=ACV, noise="repo"):
        if noise == "repo":
            cp = cstep(lambda x: clearance(x, self.theta), self.mbar)
            return np.full(self.m.size, np.hypot(sc * L_REC, cp * acv * self.mbar))
        cp = cstep(lambda x: clearance(x, self.theta, inc=self.inc), self.m)
        return np.hypot(sc * L_REC, cp * acv * self.m)


def var_th(r, ke_sd=0.0, km_exact=False, b_frac=LF, sc=LF, acv=ACV, noise="repo",
           want_b=False):
    """
    Var(ln theta_hat) by Frisch-Waugh on the whitened design.
    ke_sd = 0: k_e exact (column left out); b_frac = 0: b exact; b_frac = inf: no b prior.
    """
    sig = r.sigma(sc, acv, noise)
    A = r.J() / sig[:, None]
    cols = [0] + ([] if km_exact else [1]) + ([] if b_frac == 0 else [2]) + ([] if ke_sd == 0 else [3])
    X = A[:, cols]
    extra = []
    if 2 in cols and np.isfinite(b_frac):
        row = np.zeros(len(cols)); row[cols.index(2)] = 1.0 / (b_frac * L_REC); extra.append(row)
    if 3 in cols:
        row = np.zeros(len(cols)); row[cols.index(3)] = 1.0 / ke_sd; extra.append(row)
    if extra:
        X = np.vstack([X, np.array(extra)])

    def fw(k):
        x0 = X[:, k]
        Z = np.delete(X, k, axis=1)
        if Z.shape[1]:
            coef, *_ = np.linalg.lstsq(Z, x0, rcond=None)
            res = x0 - Z @ coef
        else:
            res = x0
        return 1.0 / float(res @ res)
    if want_b:
        return fw(0), fw(cols.index(2))
    return fw(0)


def med(v):
    return float(np.median(np.asarray(v, float)))


def sd_med(recs, **kw):
    return med([np.sqrt(var_th(r, **kw)) for r in recs])


def recs_for(theta, times, seeds=SEEDS, cap=None):
    return [Rec(theta, times, s, cap) for s in seeds]


def monthly(n):
    return np.arange(n) * MONTH


def hr(title):
    print("\n" + "-" * 96)
    print("  " + title)
    print("-" * 96)


# =================================================================== (0)
print("=" * 96)
print("  RECHECK OF t8_floor.py: independent footing, FIM by complex step + Frisch-Waugh, Monte Carlo")
print("=" * 96)
hr("(0) FOOTING REBUILT FROM THE STATED MODEL")
mb5 = m_bar_of(TH5)
r0 = Rec(TH5, monthly(25), 0)
cp5 = cstep(lambda x: clearance(x, TH5), mb5)
print(f"  L = {L_REC:.4f}; m at normal = {m_bar_of(TH_N):.4f}; m_bar at 5% = {mb5:.4f}; "
      f"s = theta g / L at 5% = {TH5 * mb5 / (KM + mb5) / L_REC:.4f}")
print(f"  sigma terms at 5%: scatter 0.05 L = {LF * L_REC:.3f}; assay C'(m_bar) 0.05 m_bar = "
      f"{cp5 * ACV * mb5:.3f}; total {r0.sigma()[0]:.3f}")
cmp("scatter term", LF * L_REC, "{:.3f}")
cmp("assay term", cp5 * ACV * mb5, "{:.3f}")
cmp("sigma", r0.sigma()[0], "{:.3f}")
g = np.arange(0.0, 600 * MONTH, MONTH)
inc_cap = dose_inc(g, cap=200)
inc_cont = dose_inc(g)
k_last = int(np.sum(g <= 199 * TAU))
print(f"  dose train capped at 200 doses: last dose day {199 * TAU:.0f}; monthly draws on or before it: "
      f"{k_last}; increment at draw {k_last + 1} (t = {g[k_last]:.1f}) = {inc_cap[k_last]:.1f}; "
      f"max increment from draw {k_last + 2} on = {inc_cap[k_last + 1:].max():.3f}")
print(f"  capped and continued trains agree before day 2786 to "
      f"{np.max(np.abs(inc_cap[:k_last] - inc_cont[:k_last]) / inc_cont[:k_last]):.1e} (relative)")
tt = np.arange(0.0, 730.0, 2.0)
print(f"  dose increment range over 2 years every 2 d: {dose_inc(tt).min():.1f} to {dose_inc(tt).max():.1f}")

# =================================================================== (a)
hr("(a) Var(ln theta | k_e exact), 25 monthly draws: subsets, Shapley, orders, Mobius")
FAC = ("b", "K_m", "scatter", "assay")
SUBS = [frozenset(c) for k in range(5) for c in itertools.combinations(FAC, k)]
TINY = 1e-8


def kw_removed(S):
    sc = LF * (TINY if "scatter" in S else 1.0)
    ac = ACV * (TINY if "assay" in S else 1.0)
    return dict(ke_sd=0.0, km_exact="K_m" in S, b_frac=0.0 if "b" in S else LF, sc=sc, acv=ac)


def game(recs):
    per = {S: np.array([var_th(r, **kw_removed(S)) for r in recs]) for S in SUBS}
    return {S: float(np.median(per[S])) for S in SUBS}, per


def shapley_orders(V):
    phi = {f: 0.0 for f in FAC}
    contrib = {f: [] for f in FAC}
    for order in itertools.permutations(FAC):
        S = frozenset()
        for f in order:
            d = V[S] - V[S | {f}]
            phi[f] += d / 24.0
            contrib[f].append(d)
            S = S | {f}
    return phi, contrib


def mobius(V):
    full = frozenset(FAC)
    w = {T: V[full - T] for T in SUBS}
    return {T: sum((-1) ** (len(T) - len(S)) * w[S] for S in SUBS if S <= T) for T in SUBS}


def lab(S):
    return "{" + ", ".join(f for f in FAC if f in S) + "}"


games = {}
for pct in (5, 10, 20, 50, 100):
    recs = recs_for(TH_N * pct / 100, monthly(25))
    games[pct] = game(recs)
V5, per5 = games[5]
tot5 = V5[frozenset()]
print(f"  Var(ln theta | k_e exact) at 5% (median over 24 records) = {tot5:.5f}, sd = {np.sqrt(tot5):.4f}")
print(f"  with both noise terms scaled by {TINY:.0e}: largest Var over subsets = "
      f"{max(V5[S] for S in SUBS if {'scatter', 'assay'} <= S):.2e}")
cmp("Var|ke 5%", tot5, "{:.5f}")
cmp("sd|ke 5% 25m", np.sqrt(tot5))
print("  subset sds (k_e exact + the listed factors removed):")
print("    " + "".join(f"{p:>9d}%" for p in games) + "   removed")
for S in SUBS:
    print("    " + "".join(f"{np.sqrt(games[p][0][S]):10.4f}" for p in games) + f"   {lab(S)}")
for S, key in ((frozenset({"b"}), "S{b}"), (frozenset({"K_m"}), "S{K_m}"),
               (frozenset({"scatter"}), "S{scatter}"), (frozenset({"assay"}), "S{assay}"),
               (frozenset({"b", "K_m"}), "S{b,K_m}"), (frozenset({"b", "scatter"}), "S{b,scatter}"),
               (frozenset({"b", "K_m", "scatter"}), "S{b,K_m,scatter}"),
               (frozenset({"b", "K_m", "assay"}), "S{b,K_m,assay}")):
    cmp(key, np.sqrt(V5[S]))

print("\n  Shapley shares (median game), by enumeration of orders and from Mobius dividends:")
print(f"    {'pct':>4} " + "".join(f"{f:>18}" for f in FAC) + "   sigma total   |orders - Mobius|")
for p in games:
    V = games[p][0]
    phi, contrib = shapley_orders(V)
    D = mobius(V)
    phi2 = {f: sum(D[T] / len(T) for T in SUBS if f in T) for f in FAC}
    t = V[frozenset()] - V[frozenset(FAC)]
    gap = max(abs(phi[f] - phi2[f]) for f in FAC) / t
    print(f"    {p:4d} " + "".join(f"{phi[f] / t:18.1%}" for f in FAC)
          + f"   {(phi['scatter'] + phi['assay']) / t:10.1%}   {gap:.1e}")
    if p == 5:
        for f in FAC:
            cmp(f"Sh {f}", phi[f] / t, "{:.3f}")
        cmp("first b", (V[frozenset()] - V[frozenset({'b'})]) / t, "{:.3f}")
        cmp("first K_m", (V[frozenset()] - V[frozenset({'K_m'})]) / t, "{:.3f}")
        cmp("max b", max(contrib["b"]) / t, "{:.3f}")
        cmp("max K_m", max(contrib["K_m"]) / t, "{:.3f}")
        cmp("D{b,scatter}", D[frozenset({"b", "scatter"})] / t, "{:.3f}")
        cmp("D{b,assay}", D[frozenset({"b", "assay"})] / t, "{:.3f}")
        cmp("D{all}", D[frozenset(FAC)] / t, "{:.3f}")
        print("      5% Mobius dividends (share of Var):  " + ", ".join(
            f"{lab(T)} {D[T] / t:.1%}" for T in SUBS[1:] if abs(D[T] / t) > 5e-4))
    if p == 100:
        cmp("Sh100 b", phi["b"] / t, "{:.3f}")
        cmp("Sh100 K_m", phi["K_m"] / t, "{:.3f}")

# A different but equally defensible game: the noise is SCALED, not removed. Var is
# homogeneous of degree 2 in (sigma, sd_b) jointly, so in the k_e-exact problem with b
# at its prior the noise share depends on how the b prior is tied to the noise.
print("\n  Sensitivity of the attribution to the game: b prior tied to the scatter term")
print("  (core.fim's default coupling, load_frac moves both): removing 'scatter' then also")
print("  makes b exact. Var with both together removed vs separately:")
recs5 = recs_for(TH5, monthly(25))
v_tied = med([var_th(r, ke_sd=0, b_frac=0.0, sc=0.0) for r in recs5])
print(f"    b exact + scatter 0 (load measured exactly): sd {np.sqrt(v_tied):.4f}; share of Var removed "
      f"{1 - v_tied / tot5:.1%}; assay alone then carries {v_tied / tot5:.1%}")
print("  Elasticity of Var to each source at the footing (d ln Var / d ln source):")
for nm, kw_hi in (("b prior sd", dict(b_frac=LF * 1.01)), ("scatter", dict(sc=LF * 1.01)),
                  ("assay CV", dict(acv=ACV * 1.01))):
    v1 = med([var_th(r, ke_sd=0, **kw_hi) for r in recs5])
    print(f"    {nm:12} {np.log(v1 / tot5) / np.log(1.01):.3f}")

# =================================================================== (b)
hr("(b) n-DEPENDENCE AT 5% WITH k_e EXACT")
SC = (("k_e exact", dict(ke_sd=0.0)),
      ("+ b exact", dict(ke_sd=0.0, b_frac=0.0)),
      ("+ load exact", dict(ke_sd=0.0, b_frac=0.0, sc=0.0)),
      ("+ K_m exact", dict(ke_sd=0.0, km_exact=True)),
      ("+ K_m + b exact", dict(ke_sd=0.0, km_exact=True, b_frac=0.0)))
NS = [6, 8, 12, 16, 25, 36, 50, 75, 100, 150, 200, 300, 400, 600, 1200, 2400, 6000]
cont = {}
print(f"    {'n':>5}" + "".join(f"{nm:>17}" for nm, _ in SC) + "   (monthly, therapy continued)")
for n in NS:
    recs = recs_for(TH5, monthly(n))
    cont[n] = [sd_med(recs, **kw) for _, kw in SC]
    print(f"    {n:5d}" + "".join(f"{v:17.4f}" for v in cont[n]))
cmp("cont 100", cont[100][0]); cmp("cont 600", cont[600][0]); cmp("cont 6000", cont[6000][0])
cmp("bex 25", cont[25][1]); cmp("bex 6000", cont[6000][1]); cmp("loadex 25", cont[25][2])
cmp("kmex 6000", cont[6000][3])
for n in (100, 600):
    v = sd_med(recs_for(TH5, monthly(n), cap=200), ke_sd=0.0)
    print(f"    repository's capped train (untreated after day 2786), n = {n}: k_e exact {v:.4f}")
    cmp(f"grid {n}", v)

print("\n  n * Var (median sd squared times n):")
for j, (nm, _) in enumerate(SC):
    print(f"    {nm:17}" + "".join(f"  n={n}: {n * cont[n][j] ** 2:7.3f}" for n in (25, 100, 600, 6000)))
cmp("nVar bex 6000", 6000 * cont[6000][1] ** 2, "{:.3f}")
print("  local slopes d ln sd / d ln n:")
for a, b in ((25, 50), (400, 600), (2400, 6000)):
    sl = [np.log(cont[b][j] / cont[a][j]) / np.log(b / a) for j in range(len(SC))]
    print(f"    {a}-{b}: " + "  ".join(f"{nm} {v:.3f}" for (nm, _), v in zip(SC, sl)))
    if a in (25, 400, 2400):
        cmp(f"slope {a}-{b}", sl[0], "{:.3f}")

# floor fits, two different criteria
nfit = np.array([n for n in NS if n <= 600], float)
sfit = np.array([cont[int(n)][0] for n in nfit])
res = least_squares(lambda p: np.sqrt(p[0] ** 2 + p[1] ** 2 / nfit) - sfit, x0=[0.3, 2.0],
                    bounds=([0, 0], [10, 100]))
fl, c = res.x
A = np.stack([np.ones_like(nfit), 1 / nfit], 1)
(f2, c2), *_ = np.linalg.lstsq(A, sfit ** 2, rcond=None)
print(f"\n  floor fit on sd (nonlinear LS), n 6-600: floor {fl:.4f}, c {c:.4f}; predicts at 6000 "
      f"{np.sqrt(fl ** 2 + c ** 2 / 6000):.4f} vs FIM {cont[6000][0]:.4f}")
print(f"  floor fit on sd^2 (linear LS), n 6-600: floor^2 {f2:.5f} (floor {np.sqrt(max(f2, 0)):.4f}), "
      f"c {np.sqrt(max(c2, 0)):.4f}; predicts at 6000 {np.sqrt(max(f2 + c2 / 6000, 0)):.4f}")
cmp("floor", fl); cmp("c", c); cmp("pred 6000", np.sqrt(fl ** 2 + c ** 2 / 6000))

# closed-form asymptotics: stationary expected information per monthly draw
print("\n  closed-form asymptotics: expected information per draw, monthly sampling, therapy in")
print("  steady state, dose phase on its 35-point lattice (30.4 mod 14 = 2.4, gcd 0.4 d),")
print("  load by 60-node Gauss-Hermite. c_inf = sqrt(lim n Var).")
z, wz = hermegauss(60)
wz = wz / wz.sum()
s2 = np.log1p(LOAD_CV ** 2)
loads_q = L_REC * np.exp(-0.5 * s2 + np.sqrt(s2) * z)
q = np.exp(-LAM * TAU)
sig5 = r0.sigma()[0]


def info_per_draw(phase_zero_as_peak=True):
    phases = np.arange(35) * 0.4
    if not phase_zero_as_peak:
        phases[0] = TAU - 1e-12
    E = np.zeros((4, 4))
    for ph in phases:
        dlt = DOSE * np.exp(-LAM * ph) / (1 - q)
        m = m_solve(TH5, loads_q, np.full_like(loads_q, dlt))
        p0 = np.array([np.log(TH5), np.log(KM), 0.0, np.log(KE)], complex)
        J = np.empty((m.size, 4))
        for j in range(4):
            p = p0.copy(); p[j] += 1e-30j
            th, K, b, ke = np.exp(p[0]), np.exp(p[1]), p[2], np.exp(p[3])
            J[:, j] = ((th + dlt) * m / (K + m) + ke * m + b).imag / 1e-30
        E += (J * wz[:, None]).T @ J / 35.0
    return E / sig5 ** 2


for peak in (True, False):
    E = info_per_draw(peak)
    c_bfree = np.sqrt(np.linalg.inv(E[:3, :3])[0, 0])
    c_bex = np.sqrt(np.linalg.inv(E[np.ix_([0, 1], [0, 1])])[0, 0])
    sb2 = 1 / (LF * L_REC) ** 2
    pr = np.diag([0.0, 0.0, sb2])
    curve = {n: np.sqrt(np.linalg.inv(n * E[:3, :3] + pr)[0, 0]) for n in (25, 100, 600, 6000, 60000)}
    print(f"    phase 0 {'at the dose peak' if peak else 'just before the dose'}: c_inf, b free (no prior) "
          f"{c_bfree:.4f}; b exact {c_bex:.4f}")
    print("      expected-information sd with the b prior: " + ", ".join(
        f"n={n}: {v:.4f}" for n, v in curve.items()) + f"; sqrt(6000) x sd(6000) = {np.sqrt(6000) * curve[6000]:.4f}")
    if peak:
        cmp("c_inf bex", c_bex)
        c_free_peak = c_bfree
print(f"    record-based sqrt(n Var) at n = 6000: k_e exact {np.sqrt(6000) * cont[6000][0]:.4f}, "
      f"b exact {np.sqrt(6000) * cont[6000][1]:.4f}, load exact {np.sqrt(6000) * cont[6000][2]:.4f}")
cmp("c_inf ke", np.sqrt(6000) * cont[6000][0])
cmp("c_inf loadex", np.sqrt(6000) * cont[6000][2])

# what the record teaches about b
print("\n  posterior sd(b) / prior sd(b), k_e exact:")
for lbl, times, key in (("25 monthly", monthly(25), "bratio 25"),
                        ("365 every 2 d", np.arange(0.0, 730.0, 2.0), "bratio 365"),
                        ("6000 monthly", monthly(6000), "bratio 6000")):
    v = med([np.sqrt(var_th(r, ke_sd=0.0, want_b=True)[1]) / (LF * L_REC)
             for r in recs_for(TH5, times)])
    print(f"    {lbl:16} {v:.4f}")
    cmp(key, v)

# dense within two years
print("\n  denser sampling within two years (therapy on):")
dense = {}
for dt in (30.4, 14.5, 7.3, 2.0, 1.0):
    times = np.arange(0.0, 730.0, dt)
    recs = recs_for(TH5, times)
    dense[dt] = [sd_med(recs, **kw) for _, kw in SC]
    print(f"    dt {dt:5.1f} n {times.size:4d}" + "".join(f"  {nm} {v:.4f}" for (nm, _), v in zip(SC, dense[dt])))
cmp("dense 365", dense[2.0][0]); cmp("dense 730", dense[1.0][0]); cmp("dense 100", dense[7.3][0])

# draws needed, monthly therapy continued (own bisection)
print("\n  monthly draws (therapy continued) needed for 0.20, own bisection:")
for (nm, kw), key in zip(SC, ("n monthly ke", "n monthly bex", "n monthly loadex",
                              "n monthly kmex", "n monthly kmbex")):
    f = lambda n: sd_med(recs_for(TH5, monthly(n)), **kw)
    lo, hi = 2, 1000
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if f(mid) <= GATE else (mid, hi)
    print(f"    {nm:17} n = {hi} (sd {f(hi):.4f}; at n-1 {f(hi - 1):.4f}); {hi * MONTH / 365.25:.1f} years")
    cmp(key, hi, "{:d}")

print("\n  within two years at spacing 730/n, k_e exact, scan n = 300-800:")
scan = {}
for n in range(300, 801):
    scan[n] = sd_med(recs_for(TH5, np.arange(n) * (730.0 / n)), ke_sd=0.0)
ok = [n for n in scan if scan[n] <= GATE]
bad = [n for n in scan if scan[n] > GATE]
print(f"    first passing n = {min(ok)} (sd {scan[min(ok)]:.4f}); last failing n = {max(bad)} "
      f"(sd {scan[max(bad)]:.4f}); max sd over n 330-800 = {max(scan[n] for n in range(330, 801)):.4f}")
cmp("first 2y ke", min(ok), "{:d}")
cmp("stays 2y ke", max(bad) + 1, "{:d}")

# =================================================================== (c)
hr("(c) REQUIREMENT CELLS AT 5% (two-year records)")
D730, D365, D100, D25 = (np.arange(0.0, 730.0, 1.0), np.arange(0.0, 730.0, 2.0),
                         np.arange(0.0, 730.0, 7.3), np.arange(0.0, 730.0, 30.4))
R = {k: recs_for(TH5, t) for k, t in (("730", D730), ("365", D365), ("100", D100), ("25", D25))}
cells = (("c 365 ke1", "365", dict(ke_sd=np.log(1.01))),
         ("c 730 ke2", "730", dict(ke_sd=np.log(1.02))),
         ("c 730 ke5 best", "730", dict(ke_sd=np.log(1.05), b_frac=0.0, sc=0.0)),
         ("c 365 ke2 b2", "365", dict(ke_sd=np.log(1.02), b_frac=0.02)),
         ("c 100 ke2 b2 sc2", "100", dict(ke_sd=np.log(1.02), b_frac=0.02, sc=0.02)),
         ("c 25 kmex ke2 b1 sc1", "25", dict(ke_sd=np.log(1.02), b_frac=0.01, sc=0.01, km_exact=True)))
for key, dk, kw in cells:
    v = sd_med(R[dk], **kw)
    print(f"    {key:22} {v:.4f}")
    cmp(key, v)
print("  k_e 5% across the 64-cell (b, scatter) grid at each design: minimum sd")
for dk in R:
    vals = [sd_med(R[dk], ke_sd=np.log(1.05), b_frac=b, sc=s) for b in (0.0, 0.01, 0.02, 0.05)
            for s in (0.0, 0.01, 0.02, 0.05)]
    print(f"    {dk:>4} draws: min {min(vals):.4f}")


def loosest_ke(dk, b, s):
    f = lambda lp: sd_med(R[dk], ke_sd=np.exp(lp), b_frac=b, sc=s) - GATE
    if f(np.log(1e-5)) > 0:
        return None
    return float(np.exp(brentq(f, np.log(1e-5), np.log(0.8), xtol=1e-6)))


for dk, key in (("365", "p 365 now"), ("730", "p 730 now")):
    sdk = loosest_ke(dk, LF, LF)
    print(f"    loosest k_e precision, {dk} draws, load as recorded: sd(ln k_e) = {sdk:.5f}, "
          f"p = {np.exp(sdk) - 1:.4f}")
    cmp(key, np.exp(sdk) - 1)
fb = lambda lb: sd_med(R["100"], ke_sd=0.0, b_frac=np.exp(lb)) - GATE
bneed = float(np.exp(brentq(fb, np.log(1e-4), np.log(0.05), xtol=1e-7)))
print(f"    100 every 7.3 d, k_e exact, scatter 5%: b prior needed {bneed:.4f} L")
cmp("b need 100 keex", bneed)

# =================================================================== (d)
hr("(d) k_e VARIANCE SHARE 1 - Var(ln th | k_e)/Var(ln th), 25 monthly draws")
pcts = np.unique(np.r_[np.geomspace(1.0, 1000.0, 31), 5.0, 100.0 * EDGE / TH_N, 100.0])
share = {}
print(f"    {'pct':>8} {'sd':>9} {'sd|ke':>8} {'share':>8}")
for pct in pcts:
    recs = recs_for(TH_N * pct / 100, monthly(25))
    tot = np.array([var_th(r, ke_sd=KE_SD) for r in recs])
    cond = np.array([var_th(r, ke_sd=0.0) for r in recs])
    sh = med(1 - cond / tot)
    share[round(pct, 2)] = sh
    print(f"    {pct:8.2f} {med(np.sqrt(tot)):9.4f} {med(np.sqrt(cond)):8.4f} {sh:8.2%}")
    if abs(pct - 5) < 1e-9:
        cmp("sd 5% 25m", med(np.sqrt(tot)))
        cmp("share 5%", sh)
    if abs(pct - 100 * EDGE / TH_N) < 1e-9:
        cmp("share edge", sh)
    if abs(pct - 100) < 1e-9:
        cmp("share 100%", sh)
    if abs(pct - 1000) < 1e-6:
        cmp("share 1000%", sh)
    if abs(pct - 1) < 1e-9:
        cmp("share 1%", sh)
    if abs(pct - 50.12) < 0.01:
        cmp("share 50.12%", sh)
ks = sorted(share)
for target in (0.9, 0.5, 0.1):
    k = next(i for i, p in enumerate(ks) if share[p] < target)
    print(f"    share first below {target:.0%} between {ks[k - 1]}% and {ks[k]}%")

# =================================================================== (e)
hr("(e) ROBUSTNESS: a disjoint seed set (1000-1023) and per-draw noise")
rows = (("25 monthly, k_e exact", monthly(25), dict(ke_sd=0.0)),
        ("25 monthly, k_e + load exact", monthly(25), dict(ke_sd=0.0, b_frac=0.0, sc=0.0)),
        ("25 monthly, k_e + K_m + b exact", monthly(25), dict(ke_sd=0.0, b_frac=0.0, km_exact=True)),
        ("365 every 2 d, k_e exact", D365, dict(ke_sd=0.0)),
        ("365 every 2 d, k_e 1%", D365, dict(ke_sd=np.log(1.01))),
        ("730 daily, k_e 2%", D730, dict(ke_sd=np.log(1.02))),
        ("100 every 7.3 d, k_e 2%, b 2%, sc 2%", D100, dict(ke_sd=np.log(1.02), b_frac=0.02, sc=0.02)))
print(f"    {'design':40} {'seeds 0-23':>11} {'[min, max]':>17} {'seeds 1000-23':>14} {'per-draw sigma':>15}")
for nm, times, kw in rows:
    a = [np.sqrt(var_th(r, **kw)) for r in recs_for(TH5, times)]
    b = sd_med(recs_for(TH5, times, ALT_SEEDS), **kw)
    c_ = sd_med(recs_for(TH5, times), noise="perdraw", **kw)
    print(f"    {nm:40} {med(a):11.4f}  [{min(a):.4f}, {max(a):.4f}] {b:14.4f} {c_:15.4f}")
for nm, times in (("25 monthly", monthly(25)), ("365 every 2 d", D365)):
    print(f"    fraction of the 24 records at or below 0.20, {nm}, k_e exact: "
          f"{np.mean([np.sqrt(var_th(r, ke_sd=0.0)) <= GATE for r in recs_for(TH5, times)]):.2f}")
sh_alt = med([1 - var_th(r, ke_sd=0.0) / var_th(r, ke_sd=KE_SD) for r in recs_for(TH5, monthly(25), ALT_SEEDS)])
sh_alt_e = med([1 - var_th(r, ke_sd=0.0) / var_th(r, ke_sd=KE_SD) for r in recs_for(EDGE, monthly(25), ALT_SEEDS)])
print(f"    k_e share, seeds 1000-1023: 5% {sh_alt:.4f}; band edge {sh_alt_e:.4f}")
ph = [sd_med(recs_for(TH5, np.arange(25) * dt), ke_sd=0.0, b_frac=0.0, sc=0.0)
      for dt in (28.0, 29.2, 30.0, 30.4, 31.0, 32.0)]
print("    25 draws, k_e + load exact, by spacing 28.0/29.2/30.0/30.4/31.0/32.0 d: "
      + " / ".join(f"{v:.4f}" for v in ph))

# =================================================================== (f)
hr("(f) MONTE CARLO: MAP fits with k_e fixed at the truth, b drawn from its prior")
REPS = 25


STARTS6 = ((0, 0), (-1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-2.5, 0.0))
STARTS40 = tuple((a, b) for a in (-5.0, -3.5, -2.5, -1.5, -0.75, 0.0, 0.75, 1.5, 2.5, 3.5)
                 for b in (-2.0, -0.7, 0.7, 2.0))


def mc(times, b_known=False, sc=LF, seeds=SEEDS, reps=REPS, seed0=777, starts=STARTS6,
       return_raw=False):
    rng = np.random.default_rng(seed0)
    est, crb, at_bound, costs = [], [], 0, []
    lt0, lk0 = np.log(TH5), np.log(KM)
    sb = LF * L_REC
    for s in seeds:
        r = Rec(TH5, times, s)
        sig = r.sigma(sc=sc)
        crb.append(var_th(r, ke_sd=0.0, b_frac=0.0 if b_known else LF, sc=sc))
        truth = clearance(r.m, TH5, inc=r.inc)
        for _ in range(reps):
            bt = 0.0 if b_known else rng.normal(0.0, sb)
            y = truth + bt + rng.normal(0.0, sig)

            def resid(x):
                th, K = np.exp(x[0]), np.exp(x[1])
                b = 0.0 if b_known else x[2]
                pr = (th + r.inc) * r.m / (K + r.m) + KE * r.m + b
                out = (y - pr) / sig
                return out if b_known else np.r_[out, b / sb]
            lo = [lt0 - 7, lk0 - 5] + ([] if b_known else [-15 * sb])
            hi = [lt0 + 5, lk0 + 5] + ([] if b_known else [15 * sb])
            best = None
            for dlt, dlk in starts:
                x0 = [lt0 + dlt, lk0 + dlk] + ([] if b_known else [0.0])
                sol = least_squares(resid, x0, bounds=(lo, hi), method="trf",
                                    xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=400)
                if best is None or sol.cost < best.cost:
                    best = sol
            est.append(best.x[0] - lt0)
            costs.append(best.cost)
            at_bound += int(np.any(np.isclose(best.x, lo, atol=1e-6)) or np.any(np.isclose(best.x, hi, atol=1e-6)))
    e = np.array(est)
    if return_raw:
        return e, np.array(costs)
    q16, q50, q84 = np.percentile(e, [16, 50, 84])
    bs = np.random.default_rng(seed0 + 1).integers(0, e.size, (400, e.size))
    eb = e[bs]
    se_rmse = float(np.std(np.sqrt(np.mean(eb ** 2, axis=1))))
    qb = np.percentile(eb, [16, 84], axis=1)
    se_half = float(np.std((qb[1] - qb[0]) / 2))
    return dict(rmse=float(np.sqrt(np.mean(e ** 2))), sd=float(e.std()), bias=float(np.mean(e)),
                med=q50, half=(q84 - q16) / 2, crb_rms=float(np.sqrt(np.mean(crb))),
                crb_med=float(np.median(np.sqrt(crb))), n=e.size, at_bound=at_bound / e.size,
                frac_off2=float(np.mean(np.abs(e) > np.log(2))), se_rmse=se_rmse, se_half=se_half)


print("  Replicates per record are listed; 6 starts per fit unless stated. Bootstrap standard")
print("  errors (400 resamples of the fits) are given for the RMSE and the 16-84% half-width.")
print(f"    {'design':40} {'fits':>5} {'CRB med':>8} {'CRB rms':>8} {'MC rmse (se)':>15} {'bias':>7} "
      f"{'16-84 half (se)':>16} {'|err|>ln2':>9} {'bound':>6} {'time':>6}")
MC = {}
for nm, times, bk, sc, reps in (("25 monthly, k_e exact", monthly(25), False, LF, 50),
                                ("25 monthly, k_e + load exact", monthly(25), True, 0.0, 25),
                                ("328 monthly (cont.), k_e exact", monthly(328), False, LF, 25),
                                ("365 every 2 d, k_e exact", D365, False, LF, 100),
                                ("730 daily, k_e exact", D730, False, LF, 25)):
    t0 = time.time()
    o = mc(times, bk, sc, reps=reps)
    MC[nm] = o
    print(f"    {nm:40} {o['n']:5d} {o['crb_med']:8.4f} {o['crb_rms']:8.4f} "
          f"{o['rmse']:8.4f} ({o['se_rmse']:.4f}) {o['bias']:+7.4f} {o['half']:8.4f} ({o['se_half']:.4f}) "
          f"{o['frac_off2']:9.3f} {o['at_bound']:6.3f} {time.time() - t0:5.0f}s")

print("\n  start-set check, 365 every 2 d, k_e exact, records 0-5 x 20 replicates: 6 starts vs a")
print("  40-start grid (ln theta offsets -5..+3.5, ln K_m offsets -2..+2), same simulated data")
e6, c6 = mc(D365, seeds=SEEDS[:6], reps=20, seed0=4242, return_raw=True)
e40, c40 = mc(D365, seeds=SEEDS[:6], reps=20, seed0=4242, starts=STARTS40, return_raw=True)
print(f"    fits where 40 starts found a lower cost by > 1e-6: {int(np.sum(c6 - c40 > 1e-6))} of {e6.size}; "
      f"RMSE 6 starts {np.sqrt(np.mean(e6 ** 2)):.4f}, 40 starts {np.sqrt(np.mean(e40 ** 2)):.4f}")

# =================================================================== (g)
hr("(g) DOES THE 200-DOSE CAP (F8) TOUCH REPORT.md TASK 3? band edge, k_e prior 0.8")
print("  REPORT Task 3 fits sd(n) = sqrt(floor^2 + c^2/n) over n = 8..400 on GRID prefixes at")
print("  theta = 2472. Draws beyond ~93 there are untreated. Same fit, capped vs continued train:")
ns_t3 = np.array([8, 12, 16, 25, 36, 50, 75, 100, 150, 200, 300, 400], float)
for cap, tag in ((200, "capped (repository GRID)"), (None, "therapy continued")):
    sds, nv = [], {}
    for n in ns_t3.astype(int):
        recs = recs_for(EDGE, monthly(n), cap=cap)
        sds.append(sd_med(recs, ke_sd=KE_SD))
        if n in (25, 100, 400):
            nv[n] = n * med([var_th(r, ke_sd=0.0) for r in recs])
    sds = np.array(sds)
    rs = least_squares(lambda p: np.sqrt(p[0] ** 2 + p[1] ** 2 / ns_t3) - sds, x0=[0.18, 0.4],
                       bounds=([0, 0], [10, 100]))
    fl3, c3 = rs.x
    print(f"    {tag:26} sd(400) {sds[-1]:.4f}; floor {fl3:.4f}, c {c3:.4f}; floor share at n = 25 "
          f"{fl3 ** 2 / (fl3 ** 2 + c3 ** 2 / 25):.1%}; n Var|k_e at 25/100/400: "
          + "/".join(f"{nv[k]:.4f}" for k in (25, 100, 400)))

# =================================================================== summary
hr("COMPARISON WITH THE INVESTIGATOR'S PRINTED VALUES")
print(f"    {'quantity':24} {'claimed':>10} {'recheck':>10} {'rel diff':>10}")
worst = 0.0
for key, cl, val, fmt in COMPARE:
    rd = (val - cl) / abs(cl) if cl not in (None, 0) else float("nan")
    worst = max(worst, abs(rd)) if np.isfinite(rd) else worst
    print(f"    {key:24} {fmt.format(cl) if cl is not None else '-':>10} {fmt.format(val):>10} {rd:10.2e}")
print(f"  largest relative difference: {worst:.2e}")
print(f"\n  runtime {time.time() - T_START:.0f} s")
