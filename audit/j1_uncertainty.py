"""
Phase 1, step 2: items 2 and 2b, uncertainty.

Four ways of putting an interval on ln theta, compared on a deliberately small
grid of conditions, and one recommendation at the end.

    filter          Wald, theta_hat +/- 1.96 sd from the inverse information
                    matrix at the fit
    profile         {theta : J(theta) - J_min <= 3.84}, the chi^2_1 threshold
    calibrated      the same set with the threshold replaced by the 95th
                    percentile of a PARAMETRIC BOOTSTRAP LR distribution
                    (test inversion)
    percentile      the 2.5 and 97.5 percentiles of ln theta_hat* over the
                    same bootstrap

The parametric bootstrap is run two ways, which the brief asks to be compared:
resimulating with the nuisances FIXED at their estimates, and resimulating with
them REDRAWN from their priors.

Footing, stated once and repeated beside anything that departs from it. Local
model constants as src/common/params.py carries them. Design: the paper's own
monthly record np.arange(0, 730, 30.4), n = 25, and nothing else. Ordinate
noise sd sigma_total(0.05, 0.05) = 21.9612. Parameters in the paper's
coordinates. `b` is the BALANCE's b, an unrecorded addition to the load, so
b >= 0 means the record under-states intake.

The main setting draws the nuisances FROM their priors, which is the footing on
which the filter's number is a posterior sd (i1 part B). `at prior means` is
kept only as a reference row.

`load prior off` means the load-bias prior only. The k_e log-prior always
stays: i1's ladder shows that flattening it alone costs a factor of 23 and the
fully prior-free design is not identified at all, so a study of interval
quality there would be a study of nothing.

The priors on k_e and on the dose scale are carried in their DELTA-METHOD form,
Gaussian on the parameter rather than on its log, because the inner problem has
to stay linear. i1's `part_b_prior_form` measures what that costs on this
design: a median 5.9e-5 and a worst 2.0e-3 in ln theta, against interval widths
here of order 0.7. The same objective is used for J_min, for J(theta) and for
every bootstrap replicate, so the LR statistics are internally consistent
whatever the prior's form.

Every figure is over 8 seeds with the load record fixed within a seed, reported
as the median with the min-max range. Nothing here records an expected result
in advance.

Run:  python audit/j1_uncertainty.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from design.run import P, delta, load, LOAD_FRAC, ASSAY_CV, M_BAR
import model_A.run as LOC

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
CV = float(ASSAY_CV)
SIG = LOC.sigma_total(LOAD_FRAC, CV)
BLOC = TH * KM / (KM + M_BAR) ** 2 + KE
S_B = LOAD_FRAC * LB
KE_PRIOR_LOG = 0.8
S_KE = KE_PRIOR_LOG * KE
MONTHLY = np.arange(0.0, 730.0, 30.4)
SEEDS = range(8)
GRID = np.geomspace(3.0, 400.0, 121)
#: the ln theta grid the profile is inverted on, +/- 1.5 in log around the
#: truth, which is about nine times the reported sd.
LT = np.log(TH) + np.linspace(-1.5, 1.5, 161)
CHI2_1_95 = 3.841458820694124
NREP, NBOOT = 30, 199


def banner(t):
    print("\n" + "=" * 110); print("  " + t); print("=" * 110)


def spread(v, fmt="{:.4f}"):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if not v.size:
        return "all non-finite"
    return (fmt.format(np.median(v)) + " (" + fmt.format(v.min()) + " to "
            + fmt.format(v.max()) + ")")


# --------------------------------------------------------------- the model
def loads_for(times, load_cv, seed):
    if not load_cv:
        return np.full(len(times), LB)
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    return LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), len(times)))


def levels(times, loads, b=0.0, ke=None, theta=None, Km=None, kappa=1.0):
    ke = KE if ke is None else ke
    theta = TH if theta is None else theta
    Km = KM if Km is None else Km
    inc = kappa * np.asarray(delta(times), float)
    tot = loads + b
    B = ke * Km + theta + inc - tot
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * tot * Km)) / (2.0 * ke)


# --------------------------------------------- the penalised bounded solver
def build(Km, m, inc, s, cfg, theta_fixed=None):
    """
    The design matrix, response, prior precisions and prior means at one or
    many K_m, with theta optionally held fixed.

    Parameter order: [theta] [kappa] k_e b, where b is the BALANCE's b and so
    enters with a column of -1. The residual is

        theta g + kappa Delta g + k_e m - b - s .
    """
    Km = np.atleast_1d(np.asarray(Km, float))
    s = np.asarray(s, float)
    g = m / (Km[..., None] + m)
    # every array is broadcast to one common (..., n) shape first, so that a
    # grid of K_m, a grid of theta and a stack of bootstrap resamples can all
    # be in flight at once without the caller reshaping anything.
    shape = np.broadcast_shapes(g.shape, s.shape)
    if theta_fixed is not None:
        shape = np.broadcast_shapes(shape, np.shape(theta_fixed) + (1,) * 0)
        shape = np.broadcast_shapes(shape, (np.asarray(theta_fixed) * g).shape)
    g = np.broadcast_to(g, shape)
    y = np.broadcast_to(s, shape).astype(float).copy()
    cols, sds, p0 = [], [], []
    if theta_fixed is None:
        cols.append(g); sds.append(np.inf); p0.append(0.0)
    else:
        y = y - np.broadcast_to(np.asarray(theta_fixed) * (m / (Km[..., None] + m)),
                                shape)
    if cfg["kappa"] is None:
        y = y - inc * g
    else:
        cols.append(inc * g); sds.append(cfg["kappa"]); p0.append(1.0)
    kp = cfg.get("ke_prior", KE_PRIOR_LOG)
    cols.append(np.broadcast_to(m, shape))
    sds.append(np.inf if kp is None else kp * KE); p0.append(KE)
    cols.append(np.broadcast_to(-1.0, shape))
    sds.append(S_B if cfg["load_prior"] else np.inf); p0.append(0.0)
    return (np.stack(cols, axis=-1), y, np.array(sds, float),
            np.array(p0, float))


def solve(X, y, sds, p0, bound_b=False):
    """
    Penalised least squares, batched over every leading axis, with an optional
    b >= 0 bound.

    b is the LAST parameter. One inequality has an exact active-set solution:
    solve free, and where b came out negative re-solve with b held at zero.
    Returns (p, J, on_bound).
    """
    D = np.diag(np.where(np.isinf(sds), 0.0, 1.0 / np.maximum(sds, 1e-300) ** 2))
    A = np.einsum("...ni,...nj->...ij", X, X) / SIG ** 2 + D
    rhs = np.einsum("...ni,...n->...i", X, y) / SIG ** 2 + D @ p0
    p = np.linalg.solve(A, rhs[..., None])[..., 0]
    on = np.zeros(p.shape[:-1], bool)
    if bound_b:
        bad = p[..., -1] < 0.0
        if bad.any():
            Xr, yr = X[..., :-1], y
            Dr, p0r = D[:-1, :-1], p0[:-1]
            Ar = (np.einsum("...ni,...nj->...ij", Xr, Xr) / SIG ** 2 + Dr)
            rr = np.einsum("...ni,...n->...i", Xr, yr) / SIG ** 2 + Dr @ p0r
            pr = np.linalg.solve(Ar, rr[..., None])[..., 0]
            p = np.where(bad[..., None],
                         np.concatenate([pr, np.zeros(pr.shape[:-1] + (1,))],
                                        axis=-1), p)
            on = bad
    r = np.einsum("...ni,...i->...n", X, p) - y
    d = p - p0
    J = (np.sum(r ** 2, axis=-1) / SIG ** 2
         + np.sum(np.where(np.isinf(sds), 0.0,
                           d ** 2 / np.maximum(sds, 1e-300) ** 2), axis=-1))
    return p, J, on


def fit(m, inc, s, cfg, grid=GRID):
    """The global fit: scan K_m, solve the rest exactly at each node."""
    X, y, sds, p0 = build(grid, m, inc, s, cfg)
    p, J, on = solve(X, y, sds, p0, cfg["bound_b"])
    j = int(np.argmin(J))
    out = dict(theta=float(p[j, 0]), Km=float(grid[j]), J=float(J[j]),
               on_bound=bool(on[j]),
               on_edge=bool(j == 0 or j == len(grid) - 1))
    k = 1
    if cfg["kappa"] is not None:
        out["kappa"] = float(p[j, k]); k += 1
    out["ke"] = float(p[j, k]); out["b"] = float(p[j, k + 1])
    return out


def profile(m, inc, s, cfg, lt=LT, grid=GRID):
    """J profiled over K_m and the linear parameters, at each ln theta."""
    th = np.exp(lt)
    X, y, sds, p0 = build(grid, m, inc, s, cfg,
                          theta_fixed=th[:, None, None])
    _, J, _ = solve(X, y, sds, p0, cfg["bound_b"])
    return J.min(axis=1)


def interval_from(lt, lr, thr):
    """The connected set {lt : lr <= thr} around its minimum, by interpolation."""
    ok = lr <= thr
    if not ok.any():
        return np.nan, np.nan
    j = int(np.argmin(lr))
    i = j
    while i > 0 and ok[i - 1]:
        i -= 1
    k = j
    while k < len(lt) - 1 and ok[k + 1]:
        k += 1
    def cross(a, b):
        if lr[b] == lr[a]:
            return lt[a]
        return lt[a] + (thr - lr[a]) * (lt[b] - lt[a]) / (lr[b] - lr[a])
    lo = lt[0] if i == 0 else cross(i, i - 1)
    hi = lt[-1] if k == len(lt) - 1 else cross(k, k + 1)
    # an interval that runs off the end of the ln theta grid is NOT an
    # interval of that width, it is an unbounded one truncated by the grid.
    # Reporting its width as the grid's would read a failure as a result.
    censored = bool(i == 0 or k == len(lt) - 1)
    return float(lo), float(hi), censored


# ------------------------------------------------------------ the bootstrap
def bootstrap(m, inc, s, cfg, f, nboot=NBOOT, redraw=False, rng=None):
    """
    A parametric bootstrap at the fit.

    redraw=False resimulates with the nuisances FIXED at their estimates,
    which is the textbook parametric bootstrap.
    redraw=True redraws b and k_e from their PRIORS for each resample, which
    asks a different question: not `how would this estimate vary if the
    patient were re-measured` but `how would it vary across the population the
    prior describes`.

    Returns (ln theta_hat* array, LR* array). LR*_j = J*(theta_hat) - J*_min,
    the statistic whose 95th percentile calibrates the profile threshold.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    kap = f.get("kappa", 1.0)
    if redraw:
        bs = rng.normal(0.0, S_B, nboot)
        if cfg["bound_b"]:
            bs = np.abs(bs)
        kes = KE * np.exp(rng.normal(0.0, KE_PRIOR_LOG, nboot))
    else:
        bs = np.full(nboot, f["b"])
        kes = np.full(nboot, f["ke"])
    g = m / (f["Km"] + m)
    mu = (f["theta"] * g + kap * inc * g
          + kes[:, None] * m[None, :] - bs[:, None])
    S = mu + rng.normal(0.0, SIG, (nboot, len(m)))

    X, y, sds, p0 = build(GRID, m, inc, S[:, None, :], cfg)
    p, J, _ = solve(X, y, sds, p0, cfg["bound_b"])          # (nboot, G)
    jmin = J.min(axis=1)
    Xr, yr, sdr, p0r = build(GRID, m, inc, S[:, None, :], cfg,
                             theta_fixed=f["theta"])
    _, Jr, _ = solve(Xr, yr, sdr, p0r, cfg["bound_b"])
    jres = Jr.min(axis=1)
    k = np.argmin(J, axis=1)
    th = p[np.arange(nboot), k, 0]
    lt = np.where(th > 0, np.log(np.maximum(th, 1e-300)), np.nan)
    return lt, jres - jmin


# ------------------------------------------------------------- one scenario
def record(seed, rep, cfg, load_cv, b_mode, kappa_true=1.0):
    """
    One synthetic record and the truth behind it.

    b_mode `draw`  b and k_e drawn from their priors (the main setting)
           `mean`  b = 0, k_e = KE (the reference row)
           a float b fixed at that value, k_e drawn
    """
    L = loads_for(MONTHLY, load_cv, seed)
    rng = np.random.default_rng(7_000_003 * seed + rep)
    if b_mode == "draw":
        b_t = rng.normal(0.0, S_B)
        if cfg["bound_b"]:
            b_t = abs(b_t)
        ke_t = KE * np.exp(rng.normal(0.0, KE_PRIOR_LOG))
    elif b_mode == "mean":
        b_t, ke_t = 0.0, KE
    else:
        b_t = float(b_mode)
        ke_t = KE * np.exp(rng.normal(0.0, KE_PRIOR_LOG))
    m = levels(MONTHLY, L, b=b_t, ke=ke_t, kappa=kappa_true)
    inc = np.asarray(delta(MONTHLY), float)
    # at the truth the model flux less b IS the recorded load, so the
    # observation is L + noise, with the noise sd the filter's own combined
    # sigma_total rather than its two components separately.
    g = m / (KM + m)
    s = (TH * g + kappa_true * inc * g + ke_t * m - b_t
         + rng.normal(0.0, SIG, len(m)))
    return m, inc, s, b_t, ke_t


def scenario(cfg, load_cv, b_mode, redraw=False, nrep=NREP, nboot=NBOOT,
             seeds=SEEDS):
    """One configuration: coverage, width and diagnostics for all four methods."""
    METH = ("filter", "profile", "calibrated", "percentile")
    per_seed = {k: {"cov": [], "wid": [], "cen": []} for k in METH}
    bias, rmse, lrq, lrt, bnd, edge = [], [], [], [], [], []
    for s0 in seeds:
        hit = {k: [] for k in METH}
        wid = {k: [] for k in METH}
        cns = {k: [] for k in METH}
        est, lr95, lrtrue, nb, ne = [], [], [], 0, 0
        for i in range(nrep):
            m, inc, s, b_t, ke_t = record(s0, i, cfg, load_cv, b_mode)
            f = fit(m, inc, s, cfg)
            if not (f["theta"] > 0):
                continue
            nb += f["on_bound"]; ne += f["on_edge"]
            lth = np.log(f["theta"])
            est.append(lth)
            # filter: Wald from the curvature of the same objective
            lr = profile(m, inc, s, cfg)
            # the statistic whose 95th percentile IS the calibrated threshold:
            # the profile LR evaluated at the TRUE ln theta, over replicates.
            lrtrue.append(float(np.interp(np.log(TH), LT, lr - lr.min())))
            sd_w = _wald_sd(m, inc, s, cfg, f)
            iv = {"filter": (lth - 1.96 * sd_w, lth + 1.96 * sd_w),
                  "profile": interval_from(LT, lr - lr.min(), CHI2_1_95)}
            iv["filter"] = iv["filter"] + (False,)
            rng = np.random.default_rng(90_000 + 1000 * s0 + i)
            bl, blr = bootstrap(m, inc, s, cfg, f, nboot, redraw, rng)
            q = float(np.nanquantile(blr, 0.95)) if np.isfinite(blr).any() else np.nan
            lr95.append(q)
            iv["calibrated"] = (interval_from(LT, lr - lr.min(), q)
                                if np.isfinite(q) else (np.nan, np.nan, True))
            bl = bl[np.isfinite(bl)]
            iv["percentile"] = ((float(np.quantile(bl, 0.025)),
                                 float(np.quantile(bl, 0.975)), False)
                                if bl.size > 10 else (np.nan, np.nan, True))
            for k in METH:
                lo, hi, cen = iv[k]
                if np.isfinite(lo) and np.isfinite(hi):
                    hit[k].append(lo <= np.log(TH) <= hi)
                    wid[k].append(hi - lo)
                    cns[k].append(cen)
        e = np.array(est)
        bias.append(e.mean() - np.log(TH))
        rmse.append(np.sqrt(np.mean((e - np.log(TH)) ** 2)))
        lrq.append(np.median(lr95))
        lrt.append(float(np.quantile(lrtrue, 0.95)) if lrtrue else np.nan)
        bnd.append(nb / max(len(est), 1)); edge.append(ne / max(len(est), 1))
        for k in METH:
            per_seed[k]["cov"].append(np.mean(hit[k]) if hit[k] else np.nan)
            per_seed[k]["wid"].append(np.median(wid[k]) if wid[k] else np.nan)
            per_seed[k]["cen"].append(np.mean(cns[k]) if cns[k] else np.nan)
    return dict(per_seed=per_seed, bias=bias, rmse=rmse, lrq=lrq,
                lrtrue=lrt, bound=bnd, edge=edge)


def _wald_sd(m, inc, s, cfg, f):
    """sd(ln theta) from the information matrix at the FIT, in log coordinates."""
    Km, th = f["Km"], f["theta"]
    ke, b = f["ke"], f["b"]
    kap = f.get("kappa", 1.0)
    g = m / (Km + m)
    # EVERY column here is a derivative in the coordinate its prior is stated
    # in: log for theta, K_m, k_e and kappa, natural for b. Mixing the two,
    # by pairing the d/d ln k_e column with a prior sd of 0.8 k_e, inflates
    # sd(ln theta) by about a factor of four. An earlier draft did exactly
    # that, and the bootstrap caught it.
    cols = [th * g, -(th + kap * inc) * Km * m / (Km + m) ** 2,
            np.ones_like(m), ke * m]
    kp = cfg.get("ke_prior", KE_PRIOR_LOG)
    sds = [np.inf, np.inf, S_B if cfg["load_prior"] else np.inf,
           np.inf if kp is None else kp]
    if cfg["kappa"] is not None:
        cols.append(kap * inc * g); sds.append(cfg["kappa"])
    Q = np.stack(cols, axis=1)
    D = np.diag([0.0 if np.isinf(x) else 1.0 / x ** 2 for x in sds])
    F = Q.T @ Q / SIG ** 2 + D
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


# ------------------------------------------------------------------ report
METH = ("filter", "profile", "calibrated", "percentile")


def show(tag, r, extra=""):
    ps = r["per_seed"]
    print(f"    {tag:<44}{extra}")
    for k in METH:
        c = np.array(ps[k]["cov"], float)
        w = np.array(ps[k]["wid"], float)
        cen = np.array(ps[k]["cen"], float)
        flag = "" if 0.90 <= np.nanmedian(c) <= 0.99 else "   <-"
        if np.nanmedian(cen) > 0.02:
            flag += f"  CENSORED {np.nanmedian(cen):.0%} of intervals run off the grid"
        print(f"      {k:>12}  coverage {spread(c, '{:.3f}'):>26}"
              f"   width {spread(w):>26}{flag}")
    print(f"      {'':>12}  bias {spread(r['bias']):>21}   "
          f"rmse {spread(r['rmse']):>21}")
    print(f"      {'':>12}  EMPIRICAL LR 95% at the true theta "
          f"{spread(r['lrtrue'], '{:.3f}'):>24}  against chi^2_1 = 3.841")
    print(f"      {'':>12}  bootstrap LR 95% (calibration threshold) "
          f"{spread(r['lrq'], '{:.3f}'):>17}")
    print(f"      {'':>12}  fits on the b bound {spread(r['bound'], '{:.3f}'):>20}"
          f"   on a K_m grid edge {spread(r['edge'], '{:.3f}'):>20}")
    print()


def main():
    print(__doc__)
    KAP = 0.1168      # i1 part B: the repository's own required dose-scale CV,
                      # recomputed over 8 seeds. Its printed value is 0.1167.
    base = dict(load_prior=False, kappa=None, bound_b=False)

    banner("(A) THE TWO PARAMETRIC BOOTSTRAPS. WHICH ONE RESAMPLES THE RIGHT THING")
    print("  Monthly n=25, load CV 10%, nuisances drawn from their priors, load")
    print("  prior OFF, dose scale known, no bound. The only difference between")
    print("  the two blocks is what the bootstrap holds fixed when it")
    print(f"  resimulates. {NREP} replicates x 8 seeds, {NBOOT} resamples each.\n")
    for rd, lab in ((False, "nuisances FIXED at their estimates"),
                    (True, "nuisances REDRAWN from their priors")):
        show(lab, scenario(base, 0.10, "draw", redraw=rd))
    print("  `filter` and `profile` do not use the bootstrap, so their rows are")
    print("  the same in both blocks up to the Monte Carlo; only `calibrated`")
    print("  and `percentile` move. Everything below uses the REDRAWN bootstrap,")
    print("  for the reason the table gives.")

    banner("(B) THE MAIN GRID. NUISANCES DRAWN FROM THEIR PRIORS")
    print("  Monthly n=25 throughout. `load prior` is the prior on the load")
    print("  BIAS b only; the k_e log-prior of 0.8 is always on. `dose` is the")
    print(f"  log-prior on the dose scale, either exact or the repository's own")
    print(f"  required {KAP:.4f}.\n")
    for cv in (0.02, 0.10):
        for lp in (False, True):
            for kap in (None, KAP):
                cfg = dict(load_prior=lp, kappa=kap, bound_b=False)
                show(f"load CV {cv:.0%}, load prior {'ON ' if lp else 'off'}, "
                     f"dose {'exact' if kap is None else f'CV {kap:.1%}'}",
                     scenario(cfg, cv, "draw", redraw=True))

    banner("(C) THE REFERENCE ROW: NUISANCES AT THEIR PRIOR MEANS")
    print("  The same estimator on a patient who sits exactly at the prior mean,")
    print("  b = 0 and k_e = k_e0. Kept only as a reference: i1 part B shows the")
    print("  filter's number is a posterior sd and this is not the footing it is")
    print("  stated on.\n")
    show("load CV 10%, load prior off, dose exact, b = 0, k_e = k_e0",
         scenario(base, 0.10, "mean", redraw=True))

    banner("(D) THE b >= 0 BOUND, AT TWO TRUE VALUES OF b")
    print("  b is the BALANCE's b, an unrecorded addition to the load, so b >= 0")
    print("  says the diet record cannot over-state intake. Imposing a bound")
    print("  that the truth sits ON is the classic boundary problem, and b = 0")
    print("  is exactly that case. k_e is still drawn from its prior.\n")
    for b_true in (0.0, S_B):
        for bd in (False, True):
            cfg = dict(load_prior=False, kappa=None, bound_b=bd)
            show(f"b_true = {b_true:.4f}, bound {'ON ' if bd else 'off'}, "
                 f"load CV 10%, load prior off",
                 scenario(cfg, 0.10, b_true, redraw=True))

    banner("(E) THE CONSTANT-LOAD CASE, WHERE THE BOOTSTRAP HAS NOTHING TO STAND ON")
    print("  Load CV 0%. i1 part D showed the profile over K_m is exactly flat")
    print("  here without the load prior, so the estimator returns a curve. What")
    print("  an interval method does with that is the question.\n")
    for lp in (True, False):
        cfg = dict(load_prior=lp, kappa=None, bound_b=False)
        show(f"load CV 0%, load prior {'ON ' if lp else 'off'}, dose exact",
             scenario(cfg, 0.0, "draw", redraw=True))


if __name__ == "__main__":
    main()
