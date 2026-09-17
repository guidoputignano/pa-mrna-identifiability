"""
Phase 1, item 1: the estimator.

A variable-projection fit, solved GLOBALLY over a grid in K_m, checked against
the repository's own nonlinear fit, and used to measure the one thing the
repository flags and does not measure: the errors-in-variables gap.

    The balance is  (theta + Delta_k) g(m_k) + k_e m_k = L_k + b,
                    g(m) = m/(K_m + m).

    Hold K_m FIXED and the balance is LINEAR in (theta, k_e, b):

        s_k  =  theta g_k  +  k_e m_k  +  b ,        g_k = m_k/(K_m + m_k),

    with s_k the observed ordinate less Delta_k g_k. So there is exactly ONE
    nonlinear parameter. Scanning it on a grid and solving the 3x3 normal
    equations at each node gives the global minimiser of the same objective
    the repository's least_squares minimises locally, with no start, no
    convergence test and no local minimum to fall into.

    The `beta` form of P3 is the other route: multiplying through by
    (K_m + m) makes the problem linear in
    beta = (K_m, K_m b, b - theta - k_e K_m, -k_e) against (L_k, 1, m_k, m_k^2).
    That is linear in FOUR parameters with no grid at all, at the price of
    reweighting every residual by (K_m + m_k). Both are fitted here and the
    price is reported rather than assumed.

Parts
    A  the three estimators agree on noiseless data, at a varying load
    B  and on noisy data, where the repository's own noise model is used
    C  the errors-in-variables bias at the paper's assay CV
    D  at a CONSTANT load the profile over K_m is FLAT: the estimator returns
       the whole fibre, not a point

Footing, stated once. Local model constants as src/common/params.py carries
them: theta = 3000, K_m = 100, k_e = 5, L = load(0) = 322.7273, D = 6000,
T = 14 d. Design: the paper's own monthly record, np.arange(0, 730, 30.4),
n = 25. Noise: the repository's, ordinate noise of sd sigma_total(load_frac,
assay_cv), load record 5%, assay CV 5%, which is ASSAY_CV/100 in the
repository's units. Parameters are reported in the paper's coordinates,
ln theta and ln K_m. Every Monte Carlo figure is the median over 8 independent
seeds of a batch statistic, with the min-max range beside it.

Nothing here records an expected result in advance.

Run:  python audit/i1_estimator.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from design.run import (P, delta, load, LOAD_FRAC, ASSAY_CV, M_BAR,
                        clearance_A, DOSE_INTERVAL, D_AMP, LAMBDA)
import model_A.run as LOC

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
#: design.run.ASSAY_CV is a FRACTION (0.05) while model_A.run.ASSAY_CV is a
#: PERCENTAGE (5.0). Both mean 5%, and sigma_total's own default divides the
#: second by 100, so sigma_total(LOAD_FRAC, ASSAY_CV) and sigma_total(LOAD_FRAC)
#: agree at 21.9612. Dividing design.run's value by 100 again does not, and an
#: earlier draft of this script did exactly that.
CV = float(ASSAY_CV)
SIG = LOC.sigma_total(LOAD_FRAC, CV)
#: the paper's own priors, which are what make its 0.1697 a number at all.
#: Without them the four-parameter fit on this design is not identified: the
#: Fisher bound with no priors is sd(ln theta) = 5.54. The k_e prior is
#: Gaussian on ln k_e with sd 0.8 in the filter; the linear inner problem
#: carries its delta-method equivalent, Gaussian on k_e with sd 0.8 k_e, and
#: the nonlinear fit carries the same one so that the two are the SAME
#: estimator and any difference between them is optimisation, not specification.
S_B = LOAD_FRAC * LB
S_KE = 0.8 * KE
#: the local slope of the balance in m, which is what turns an assay error on
#: the abscissa into an ordinate error. design/run.py:623's b_loc.
BLOC = TH * KM / (KM + M_BAR) ** 2 + KE
MONTHLY = np.arange(0.0, 730.0, 30.4)
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
SEEDS = range(8)
#: the K_m grid. Three decades around the truth, log-spaced, so the scan is
#: global in any sense a reader would mean by it.
GRID = np.geomspace(1.0, 1000.0, 601)


def banner(t):
    print("\n" + "=" * 96); print("  " + t); print("=" * 96)


def spread(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if not v.size:
        return "all non-finite"
    return f"{np.median(v):.6g} (8 seeds {v.min():.6g} to {v.max():.6g})"


# --------------------------------------------------------------- the truth
def levels(times, loads, b=0.0, ke=None, theta=None, Km=None):
    """The levels that solve the balance exactly, at the given nuisances."""
    ke = KE if ke is None else ke
    theta = TH if theta is None else theta
    Km = KM if Km is None else Km
    inc = np.asarray(delta(times), float)
    tot = loads + b
    B = ke * Km + theta + inc - tot
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * tot * Km)) / (2.0 * ke)


def loads_for(times, load_cv, seed):
    if not load_cv:
        return np.full(len(times), LB)
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    return LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), len(times)))


# ------------------------------------------------------------ the estimators
def _normal(Km, m, inc, s, pen):
    """
    The 3x3 penalised normal equations at a fixed K_m, batched over Km.

    Model     s_k = theta g_k + k_e m_k + b + Delta_k g_k,  g = m/(K_m+m)
    Objective ||X p - y||^2/sigma^2 + (k_e - k_e0)^2/s_ke^2 + b^2/s_b^2
    """
    Km = np.atleast_1d(np.asarray(Km, float))
    g = m[None, :] / (Km[:, None] + m[None, :])            # (G, n)
    X = np.stack([g, np.broadcast_to(m, g.shape), np.ones_like(g)], axis=-1)
    y = s[None, :] - inc[None, :] * g                      # (G, n)
    A = X.transpose(0, 2, 1) @ X / SIG ** 2
    rhs = np.einsum("gnj,gn->gj", X, y) / SIG ** 2
    if pen:
        D = np.diag([0.0, 1.0 / S_KE ** 2, 1.0 / S_B ** 2])
        A = A + D
        rhs = rhs + D @ np.array([0.0, KE, 0.0])
    return X, y, A, rhs


def varpro_profile(m, inc, s, grid=GRID, pen=True):
    """The penalised objective and the linear solution at every grid node."""
    X, y, A, rhs = _normal(grid, m, inc, s, pen)
    p = np.linalg.solve(A, rhs[..., None])[..., 0]          # (G, 3)
    r = np.einsum("gnj,gj->gn", X, p) - y
    J = np.sum(r ** 2, axis=1) / SIG ** 2
    if pen:
        J = J + (p[:, 1] - KE) ** 2 / S_KE ** 2 + p[:, 2] ** 2 / S_B ** 2
    return J, p


def varpro(m, inc, s, grid=GRID, pen=True, refine=True):
    """
    The global variable-projection fit: scan K_m, solve the rest exactly.

    There is one nonlinear parameter, so a grid scan IS a global search. No
    start, no convergence test, no local minimum.
    """
    J, p = varpro_profile(m, inc, s, grid, pen)
    j = int(np.argmin(J))
    Km_hat = float(grid[j])
    if refine and 0 < j < len(grid) - 1:
        f = lambda lk: float(varpro_profile(m, inc, s,
                                            np.array([np.exp(lk)]), pen)[0][0])
        r = minimize_scalar(f, bounds=(np.log(grid[j - 1]), np.log(grid[j + 1])),
                            method="bounded", options={"xatol": 1e-12})
        Km_hat = float(np.exp(r.x))
    Jh, ph = varpro_profile(m, inc, s, np.array([Km_hat]), pen)
    return dict(theta=float(ph[0, 0]), Km=Km_hat, ke=float(ph[0, 1]),
                b=float(ph[0, 2]), rss=float(Jh[0])), J


def beta_fit(m, inc, s):
    """
    The P3 linear form, for comparison. Multiplying the balance by (K_m + m)
    makes it linear in beta = (K_m, K_m c, c - theta - k_e K_m, -k_e) against
    the regressors (s_k, 1, m_k, m_k^2), response (Delta_k - s_k) m_k.

    Two things it cannot do, and both are reported rather than worked around.
    It reweights every residual by (K_m + m_k), so it is not minimising the
    paper's objective. And the paper's priors on b and k_e are not linear
    functions of beta, so it cannot carry them: this fit is UNPENALISED
    whatever the caller wants.
    """
    X = np.stack([s, np.ones_like(m), m, m ** 2], axis=1)
    y = (inc - s) * m
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    Km_hat = beta[0]
    if not np.isfinite(Km_hat) or abs(Km_hat) < 1e-12:
        return dict(theta=np.nan, Km=np.nan, ke=np.nan, b=np.nan, rss=np.nan)
    c, ke = beta[1] / Km_hat, -beta[3]
    r = X @ beta - y
    return dict(theta=c - beta[2] - ke * Km_hat, Km=float(Km_hat),
                ke=float(ke), b=float(c), rss=float(r @ r))


def nonlinear_fit(m, inc, s, seed=0, jitter=0.5, pen=True):
    """
    The repository's own fit, reproduced from design/run.py:628-641, with the
    same priors as varpro appended as residuals so the two are one estimator.
    """
    rng = np.random.default_rng(seed)

    def res(q):
        t, K, k = np.exp(np.clip(q[:3], -25.0, 25.0))
        r = ((t + inc) * m / (K + m) + k * m + q[3] - s) / SIG
        if pen:
            r = np.concatenate([r, [(k - KE) / S_KE, q[3] / S_B]])
        return r

    start = np.concatenate([np.log([TH, KM, KE]) + rng.normal(0.0, jitter, 3),
                            [0.0]])
    try:
        r = least_squares(res, start, method="lm", max_nfev=5000)
    except Exception:
        return dict(theta=np.nan, Km=np.nan, ke=np.nan, b=np.nan, rss=np.nan)
    t, K, k = np.exp(np.clip(r.x[:3], -25.0, 25.0))
    return dict(theta=t, Km=K, ke=k, b=float(r.x[3]), rss=float(2.0 * r.cost))


def filter_sd(m, inc, pen=True):
    """
    The information filter on the SAME record and the SAME priors, so the
    Monte Carlo has a bound to be compared against rather than a folklore one.
    """
    g = m / (KM + m)
    Q = np.stack([TH * g, -(TH + inc) * KM * m / (KM + m) ** 2,
                  np.ones_like(m), KE * m], axis=1)
    F = Q.T @ Q / SIG ** 2
    if pen:
        F = F + np.diag([0.0, 0.0, 1.0 / S_B ** 2, 1.0 / 0.8 ** 2])
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


# ------------------------------------------------------------------ part A
def part_a():
    banner("(A) NOISELESS DATA, VARYING LOAD: DO THE THREE AGREE WITH THE TRUTH?")
    print("  Monthly n=25, load CV 10%, b_true = 0, no observation noise at all.")
    print("  Unpenalised, so nothing but the data is pulling. Any departure from")
    print("  the truth is the estimator's and not the record's.\n")
    L = loads_for(MONTHLY, 0.10, 3)
    m = levels(MONTHLY, L)
    inc = np.asarray(delta(MONTHLY), float)
    s = L.copy()
    print(f"    {'estimator':>22} {'theta':>14} {'K_m':>12} {'k_e':>10} "
          f"{'b':>12} {'objective':>11}")
    print("  " + "-" * 92)
    print(f"    {'the truth':>22} {TH:14.6f} {KM:12.6f} {KE:10.6f} "
          f"{0.0:12.6f} {'':>11}")
    out = {}
    for name, fn in (
            ("varpro (global grid)", lambda: varpro(m, inc, s, pen=False)[0]),
            ("beta (P3 linear)", lambda: beta_fit(m, inc, s)),
            ("nonlinear (the repo's)",
             lambda: nonlinear_fit(m, inc, s, seed=0, pen=False))):
        r = fn()
        out[name] = r
        print(f"    {name:>22} {r['theta']:14.6f} {r['Km']:12.6f} "
              f"{r['ke']:10.6f} {r['b']:12.6f} {r['rss']:11.3e}")
    print(f"\n    worst relative error in theta "
          f"{max(abs(r['theta']/TH - 1) for r in out.values()):.3e}, "
          f"in K_m {max(abs(r['Km']/KM - 1) for r in out.values()):.3e}")
    print("\n    All three are exact on noiseless data at a varying load, which")
    print("    is the only claim this part makes. It rules out a coding error")
    print("    in any of the three before any is asked a harder question.")
    return out


# ------------------------------------------------------------------ part B
def one_record(batch_seed, rep, load_cv=0.10, eiv=False, cv=CV,
               load_frac=LOAD_FRAC, draw_nuisance=False):
    """
    One synthetic record.

    The LOAD RECORD is fixed within a batch and varies between batches: it is
    part of the design, not of the observation noise, and mixing the two would
    report design variation as estimator variance. The 8-seed spread therefore
    covers load records as well as noise draws.

    eiv=False is the REPOSITORY'S model (design/run.py:598-604): the levels are
    exact and the assay error is placed on the ORDINATE as b_loc cv m_bar.
    eiv=True is the physical model its docstring says it is not doing: the
    assay perturbs the LEVEL, the record perturbs the LOAD, and the fit sees
    only the perturbed pair.
    """
    L = loads_for(MONTHLY, load_cv, batch_seed)
    rng = np.random.default_rng(1_000_003 * batch_seed + rep)
    # the SAME standard normals feed both noise models, drawn in the same
    # order, so that `ordinate` and `EIV` differ only in WHERE the assay error
    # is placed. Without that pairing the second-order gap being measured is
    # far smaller than the Monte Carlo noise it is measured through.
    z_load = rng.standard_normal(len(MONTHLY))
    z_assay = rng.standard_normal(len(MONTHLY))
    if draw_nuisance:
        # the nuisances drawn FROM the priors the estimator carries, which is
        # the footing on which the filter's number is a posterior sd. theta is
        # held at the truth; only b and k_e are nuisances.
        b_t = rng.normal(0.0, S_B)
        ke_t = KE * np.exp(rng.normal(0.0, 0.8))
    else:
        b_t, ke_t = 0.0, KE
    m = levels(MONTHLY, L, b=b_t, ke=ke_t)
    inc = np.asarray(delta(MONTHLY), float)
    if not eiv:
        s = ((L + b_t) * (1.0 + load_frac * z_load)
             + BLOC * cv * M_BAR * z_assay)
        return m, inc, s
    return (m * (1.0 + cv * z_assay), inc,
            (L + b_t) * (1.0 + load_frac * z_load))


def batch(fit, reps, batch_seed, **kw):
    v = []
    for i in range(reps):
        m, inc, s = one_record(batch_seed, i, **kw)
        r = fit(m, inc, s, i)
        t = r["theta"]
        v.append(np.log(t) if np.isfinite(t) and t > 0 else np.nan)
    return np.array(v, float)


def fits(pen):
    return {
        "varpro (global grid)": lambda m, i, s, k: varpro(m, i, s, pen=pen)[0],
        "beta (P3 linear)": lambda m, i, s, k: beta_fit(m, i, s),
        "nonlinear (the repo's)":
            lambda m, i, s, k: nonlinear_fit(m, i, s, seed=k, pen=pen),
    }


def part_b(reps=60):
    banner("(B) NOISY DATA, THE REPOSITORY'S OWN NOISE MODEL")
    print(f"  Monthly n=25, load CV 10%, ordinate noise sd "
          f"sigma_total({LOAD_FRAC:g}, {CV:g}) = {SIG:.4f}.")
    print(f"  {reps} replicates per seed, 8 seeds, load record fixed within a")
    print("  seed. Reported in ln theta; `bias` is mean(ln theta_hat) - ln theta,")
    print("  and `rmse` is sqrt(bias^2 + sd^2) so a biased-but-tight estimator")
    print("  cannot look good by hiding in the sd column.\n")
    store = {}
    for pen, dn, lab in (
            (False, False, "UNPENALISED, four free parameters, nuisances at "
                           "their prior means"),
            (True, False, "PENALISED with the paper's own priors, nuisances AT "
                          "their prior means"),
            (True, True, "PENALISED, nuisances DRAWN FROM those priors, which "
                         "is the footing on which the filter's number is a "
                         "posterior sd")):
        fb = [filter_sd(*one_record(s0, 0, draw_nuisance=dn)[:2], pen=pen)
              for s0 in SEEDS]
        print(f"  {lab}")
        print(f"    the filter's bound on the same records: {spread(fb)}")
        print(f"    {'estimator':>22} {'bias':>28} {'sd':>28} {'rmse':>10} "
              f"{'sd/bound':>9}")
        print("  " + "-" * 102)
        for name, fit in fits(pen).items():
            if name.startswith("beta") and pen:
                print(f"    {name:>22}   cannot carry the priors: they are not "
                      f"linear in beta")
                continue
            bs, sds = [], []
            for s0 in SEEDS:
                v = batch(fit, reps, s0, draw_nuisance=dn)
                v = v[np.isfinite(v)]
                bs.append(np.mean(v) - np.log(TH))
                sds.append(np.std(v, ddof=1))
            store[(pen, dn, name)] = (bs, sds)
            rmse = np.median([np.hypot(a, b) for a, b in zip(bs, sds)])
            print(f"    {name:>22} {spread(bs):>28} {spread(sds):>28} "
                  f"{rmse:10.4f} {np.median(sds)/np.median(fb):9.2f}")
        print()
    return store


def part_b_starts(reps=40):
    banner("(B, continued) WHAT THE GRID BUYS: START DEPENDENCE")
    print("  The variable-projection fit has one nonlinear parameter and scans")
    print("  it, so it has no start to depend on. The repository's fit has four")
    print("  and is started at the truth plus N(0, jitter) in log space. Below,")
    print("  the same records fitted at several jitters, penalised, against the")
    print("  grid fit's own objective on the same record.\n")
    print(f"    {'jitter (log sd)':>16} {'worse than the grid':>22} "
          f"{'median excess objective':>26} {'bias in ln theta':>28}")
    print("  " + "-" * 98)
    for jit in (0.0, 0.25, 0.5, 1.0, 2.0):
        frac, exc, bias = [], [], []
        for s0 in SEEDS:
            w, e, v = 0, [], []
            for i in range(reps):
                m, inc, sobs = one_record(s0, i)
                gr, _ = varpro(m, inc, sobs, pen=True)
                nl = nonlinear_fit(m, inc, sobs, seed=i, jitter=jit, pen=True)
                if np.isfinite(nl["rss"]):
                    d = nl["rss"] - gr["rss"]
                    e.append(d)
                    if d > 1e-6 * max(gr["rss"], 1.0):
                        w += 1
                if np.isfinite(nl["theta"]) and nl["theta"] > 0:
                    v.append(np.log(nl["theta"]) - np.log(TH))
            frac.append(w / reps)
            exc.append(np.median(e) if e else np.nan)
            bias.append(np.mean(v) if v else np.nan)
        print(f"    {jit:16g} {np.median(frac):21.1%} "
              f"{np.median(exc):26.3e} {spread(bias):>28}")
    print("\n    `worse than the grid` is the fraction of replicates whose local")
    print("    fit ends at a strictly larger objective than the grid's global")
    print("    minimum on the same data, so it is a lower bound on how often the")
    print("    local search stopped somewhere it should not have.")


# ------------------------------------------------------------------ part C
def part_c(reps=60, pair_reps=1500):
    banner("(C) THE ERRORS-IN-VARIABLES BIAS AT THE PAPER'S ASSAY CV")
    print("  design/run.py:598-604 states the gap and declines to measure it.")
    print("  The assay error is placed on the ORDINATE, where the Fisher matrix")
    print("  puts it, rather than on the abscissa where it physically sits, and")
    print("  the docstring says fitting against a perturbed m `would make the")
    print("  fit look better than it is`. Both halves of that are testable.\n")
    print("  To FIRST order the two models agree by construction: a level error")
    print("  dm enters the balance as b_loc dm, and sigma_total is exactly")
    print(f"  hypot(load_frac L, b_loc cv m_bar) with b_loc = {BLOC:.4f}. So any")
    print("  gap is second order and should scale as cv^2. Penalised varpro,")
    print("  monthly n=25, load CV 10%.\n")
    print(f"    {'assay CV':>9} {'model':>10} {'bias in ln theta':>30} "
          f"{'sd(ln theta_hat)':>30}")
    print("  " + "-" * 88)
    rows = {}
    for cv in (0.0, 0.02, CV, 0.10, 0.20):
        for eiv, lab in ((False, "ordinate"), (True, "EIV")):
            bs, sds = [], []
            for s0 in SEEDS:
                v = batch(fits(True)["varpro (global grid)"], reps, s0,
                          eiv=eiv, cv=cv)
                v = v[np.isfinite(v)]
                bs.append(np.mean(v) - np.log(TH))
                sds.append(np.std(v, ddof=1))
            rows[(cv, lab)] = (np.array(bs), np.array(sds))
            star = "  <- the paper's" if abs(cv - CV) < 1e-12 else ""
            print(f"    {cv:9.0%} {lab:>10} {spread(bs):>30} "
                  f"{spread(sds):>30}{star}")
        print()
    print("    The two models above share their random draws replicate by")
    print("    replicate, so the difference can be taken PAIRED, which is the")
    print("    only way a second-order effect is visible through first-order")
    print(f"    noise. {pair_reps} pairs per seed, 8 seeds = {8*pair_reps} pairs")
    print("    per row; `+/-` is 1 standard error of the paired mean.\n")
    print(f"    {'assay CV':>9} {'paired gap in ln theta':>28} {'gap/cv^2':>10} "
          f"{'resolved?':>10} {'sd ratio EIV/ord':>18}")
    print("  " + "-" * 82)
    for cv in (0.02, CV, 0.10, 0.20):
        d = []
        for s0 in SEEDS:
            for i in range(pair_reps):
                m_o, inc, s_o = one_record(s0, i, eiv=False, cv=cv)
                m_e, _, s_e = one_record(s0, i, eiv=True, cv=cv)
                a = varpro(m_o, inc, s_o, pen=True)[0]["theta"]
                b = varpro(m_e, inc, s_e, pen=True)[0]["theta"]
                if a > 0 and b > 0:
                    d.append(np.log(b) - np.log(a))
        d = np.array(d)
        se = d.std(ddof=1) / np.sqrt(d.size)
        g = d.mean()
        sr = (np.median(rows[(cv, "EIV")][1])
              / np.median(rows[(cv, "ordinate")][1]))
        print(f"    {cv:9.0%} {g:15.6f} +/- {se:8.6f} {g/cv**2:10.4f} "
              f"{('yes' if abs(g) > 2*se else 'NO'):>10} {sr:18.4f}")
    return rows


# ------------------------------------------------------------------ part D
def part_d():
    banner("(D) AT A CONSTANT LOAD THE PROFILE IS FLAT: THE WHOLE FIBRE")
    print("  Noiseless, load CV 0%, so every recorded load is L exactly. The")
    print("  variable-projection profile over K_m is the exact profile of the")
    print("  objective, because the other three parameters are solved and not")
    print("  searched. If P5 is right this profile is not shallow, it is")
    print("  CONSTANT, and the (theta, k_e, b) at each node lies on the fibre")
    print("      K_m(b) = K_0 L/(L+b),  theta(b) = theta_0 + b + k_e(K_0 - K_m),")
    print("  with k_e unchanged.\n")
    L = np.full(len(MONTHLY), LB)
    m = levels(MONTHLY, L)
    inc = np.asarray(delta(MONTHLY), float)
    s = L.copy()
    J0, _ = varpro_profile(m, inc, s, GRID, pen=False)
    print(f"    UNPENALISED profile over {len(GRID)} nodes in "
          f"[{GRID[0]:g}, {GRID[-1]:g}]:")
    print(f"      min {J0.min():.6e}   max {J0.max():.6e}   "
          f"range {J0.max()-J0.min():.6e}")
    print(f"      for scale, one unit of this objective is one sigma of one "
          f"residual.\n")
    print("    b is reported in the estimator's own sign convention, an offset")
    print("    on the model flux, so it is the NEGATIVE of the fibre's b.\n")
    print(f"    {'K_m node':>9} {'theta_hat':>13} {'k_e_hat':>9} {'b_hat':>13} "
          f"{'objective':>11} {'fibre b':>13} {'fibre theta':>13} {'match':>7}"
          f" {'admissible':>11}")
    print("  " + "-" * 98)
    worst = 0.0
    for k in (10.0, 30.0, 100.0, 300.0, 1000.0):
        J, ph = varpro_profile(m, inc, s, np.array([k]), pen=False)
        th, ke, b = ph[0]
        # the fibre is stated with b added to the LOAD; this estimator's b is
        # added to the MODEL FLUX, which is design/run.py:630's convention, so
        # the two differ by a sign. Getting that wrong reads a perfect match as
        # a total failure, which an earlier draft of this script did.
        b_bal = LB * (KM / k - 1.0)
        th_fib = TH + b_bal + KE * (KM - k)
        b_fib = -b_bal
        d = max(abs(th - th_fib) / TH, abs(b - b_fib) / max(abs(b_fib), 1.0),
                abs(ke - KE) / KE)
        worst = max(worst, d)
        adm = "yes" if (th_fib > 0 and LB + b_bal > 0) else "NO"
        print(f"    {k:9g} {th:13.5f} {ke:9.5f} {b:13.5f} {J[0]:11.3e} "
              f"{b_fib:13.5f} {th_fib:13.5f} {'exact' if d < 1e-7 else 'NO':>7}"
              f" {adm:>11}")
    print(f"\n    worst relative departure from the fibre: {worst:.3e}")
    # where on this grid the fibre is PHYSICALLY admissible, from Phase 2's
    # own bounds on b: theta(b) > 0 below and the capacity ceiling above.
    # solved, not scanned: reading these off the grid would report the grid's
    # own spacing as the bound, which an earlier draft of this script did.
    from scipy.optimize import brentq
    th_b = lambda bb: TH + bb + KE * (KM - KM * LB / (LB + bb))
    b_lo = brentq(th_b, -LB + 1e-9, 0.0, xtol=1e-12)
    for frac, lab in ((0.10, "10%"), (0.05, "5%")):
        b_hi = brentq(lambda bb: th_b(bb) + PEAK - TH / frac, 0.0, 1e7,
                      xtol=1e-9)
        k_hi, k_lo = KM * LB / (LB + b_lo), KM * LB / (LB + b_hi)
        print(f"    admissible b in [{b_lo:.4f}, {b_hi:.4f}], so admissible "
              f"K_m at {lab}: [{k_lo:.4f}, {k_hi:.4f}]")
    print("    The profile is flat over the WHOLE grid. Only that sub-interval")
    print("    describes a patient, and it is h1 check 6's bound on b read in")
    print("    K_m instead: the two are the same constraint in two coordinates.")
    print("\n    So the estimator does not fail at a constant load and does not")
    print("    return a wrong answer. It returns the CORRECT answer, which is a")
    print("    CURVE. Software that reports a point here is reporting its start,")
    print("    its grid's first node, or its stopping rule.")

    J1, p1 = varpro_profile(m, inc, s, GRID, pen=True)
    j = int(np.argmin(J1))
    print(f"\n    PENALISED, the same data. The profile is no longer flat:")
    print(f"      min {J1.min():.6e} at K_m = {GRID[j]:.4f}, "
          f"max {J1.max():.6e}")
    print(f"      the minimiser's b is {p1[j,2]:.6e} and the truth's is 0")
    ok = GRID[(J1 - J1.min()) <= 1.0]
    print(f"      the objective rises by 1 over K_m in "
          f"[{ok.min():.4f}, {ok.max():.4f}]")
    pen_only = J1 - J0
    print(f"      the penalty's own share of the objective, J_pen - J_data,")
    print(f"      differs from J_pen by at most {np.max(np.abs(pen_only - J1)):.3e}")
    print(f"      across the whole grid, because J_data is {J0.max():.3e} at worst.")
    print("\n    The curvature that picks K_m out here is the PRIOR'S, not the")
    print("    data's: the data contributed a flat profile and the b penalty")
    print("    supplied everything else. That is the estimator-side statement of")
    print("    what PHASE3 section 4(a) found on the design side, and it is why")
    print("    a constant-load record must not be reported as an estimate with")
    print("    an interval, whatever a fitting routine returns.")
    return J0, J1


if __name__ == "__main__":
    print(__doc__)
    part_a()
    part_b()
    part_b_starts()
    part_c()
    part_d()
