"""
Phase 1, step 2, addendum: proper test inversion, the boundary cutoff, the
Monte Carlo standard errors, and what a seed actually varies.

Four things the first pass got wrong or left out.

(1) WHAT `calibrated` WAS, AND WHAT TEST INVERSION IS. j1's `calibrated`
    bootstraps at theta_hat, forms LR* = J*(theta_hat) - J*_min once, and uses
    that single 95th percentile as the cutoff at EVERY theta_0. That is not
    test inversion. Test inversion simulates under the CONSTRAINED fit at each
    theta_0 separately, so the cutoff q(theta_0) varies along the grid:

        accept theta_0  iff  J(theta_0) - J_min  <=  q(theta_0),
        q(theta_0) = 95th percentile of  J*(theta_0) - J*_min  over datasets
                     simulated at the constrained fit (theta_0, K_m, k_e, b).

    It is also orthogonal to `fixed vs redrawn`, which is about what the
    bootstrap does with the NUISANCES. Test inversion's textbook form fixes
    them at the CONSTRAINED estimates, which is a third thing again.

(2) THE BOUNDARY CUTOFF. With b >= 0 imposed and b_true = 0, a fit lands on
    the bound about 38% of the time. The usual correction for a boundary is a
    50:50 mixture of chi^2_0 and chi^2_1, whose 95% point is 2.706 rather than
    3.841. Whether that is the right correction HERE is a separate question,
    since the parameter on the boundary is a nuisance and not the parameter of
    interest, so it is measured rather than assumed.

(3) WHAT A SEED VARIES. In j1 and here, a seed is a DIET RECORD: the load
    realisation is drawn once per seed and held across that seed's replicates,
    which vary the observation noise and the drawn nuisances. So the seed
    spread is a DESIGN quantity, not Monte Carlo error, and the right summary
    of it is the fraction of diet patterns that reach the target.

(4) MONTE CARLO ERROR. Every coverage figure gets its replicate count and a
    binomial standard error.

Footing as j1's and repeated where it changes: monthly n=25, ordinate noise sd
sigma_total(0.05, 0.05) = 21.9612, nuisances drawn from their priors, k_e
log-prior 0.8 always on, dose scale exact unless stated. Nothing here records
an expected result in advance.

Run:  python audit/j3_testinversion.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import j1_uncertainty as J

#: the theta_0 grid the test is inverted on, in ln theta around the fit.
NT = 21
HALF = 1.1
NBOOT_TI = 99
CHI2_MIX_95 = 2.705543454095404      # 95% point of 0.5 chi^2_0 + 0.5 chi^2_1


def banner(t):
    print("\n" + "=" * 104); print("  " + t); print("=" * 104)


def spread(v, fmt="{:.4f}"):
    return J.spread(v, fmt)


def mc_se(p, n):
    return float(np.sqrt(max(p * (1.0 - p), 0.0) / max(n, 1)))


# ------------------------------------------------------- proper test inversion
def invert(m, inc, s, cfg, f, nboot=NBOOT_TI, rng=None, redraw=False):
    """
    The interval {theta_0 : LR(theta_0) <= q(theta_0)}, with q simulated
    separately at each theta_0 under the CONSTRAINED fit there.

    Returns (lo, hi, censored, the q grid, the LR grid).
    """
    rng = np.random.default_rng(0) if rng is None else rng
    lt0 = np.log(f["theta"]) + np.linspace(-HALF, HALF, NT)
    lr_obs = J.profile(m, inc, s, cfg, lt=lt0) - f["J"]
    q = np.empty(NT)
    for i, lt in enumerate(lt0):
        th0 = float(np.exp(lt))
        # the constrained fit at theta_0 on the real data
        X, y, sds, p0 = J.build(J.GRID, m, inc, s, cfg, theta_fixed=th0)
        p, Jr, _ = J.solve(X, y, sds, p0, cfg["bound_b"])
        k = int(np.argmin(Jr))
        pc = p[k]
        Km_c = J.GRID[k]
        kk = 0
        kap_c = 1.0
        if cfg["kappa"] is not None:
            kap_c = float(pc[kk]); kk += 1
        ke_c, b_c = float(pc[kk]), float(pc[kk + 1])
        # datasets simulated AT that constrained fit
        if redraw:
            bs = rng.normal(0.0, J.S_B, nboot)
            if cfg["bound_b"]:
                bs = np.abs(bs)
            kes = J.KE * np.exp(rng.normal(0.0, J.KE_PRIOR_LOG, nboot))
        else:
            bs = np.full(nboot, b_c)
            kes = np.full(nboot, ke_c)
        g = m / (Km_c + m)
        mu = (th0 * g + kap_c * inc * g + kes[:, None] * m[None, :]
              - bs[:, None])
        S = mu + rng.normal(0.0, J.SIG, (nboot, len(m)))
        Xu, yu, su, pu = J.build(J.GRID, m, inc, S[:, None, :], cfg)
        _, Ju, _ = J.solve(Xu, yu, su, pu, cfg["bound_b"])
        Xr, yr, sr, pr = J.build(J.GRID, m, inc, S[:, None, :], cfg,
                                 theta_fixed=th0)
        _, Jr2, _ = J.solve(Xr, yr, sr, pr, cfg["bound_b"])
        q[i] = float(np.quantile(np.maximum(Jr2.min(axis=1) - Ju.min(axis=1),
                                            0.0), 0.95))
    ok = lr_obs <= q
    if not ok.any():
        return np.nan, np.nan, True, q, lr_obs
    j = int(np.argmin(lr_obs - q))
    i0 = j
    while i0 > 0 and ok[i0 - 1]:
        i0 -= 1
    i1 = j
    while i1 < NT - 1 and ok[i1 + 1]:
        i1 += 1
    d = lr_obs - q

    def cross(a, b):
        if d[b] == d[a]:
            return lt0[a]
        return lt0[a] - d[a] * (lt0[b] - lt0[a]) / (d[b] - d[a])

    lo = lt0[0] if i0 == 0 else cross(i0, i0 - 1)
    hi = lt0[-1] if i1 == NT - 1 else cross(i1, i1 + 1)
    return float(lo), float(hi), bool(i0 == 0 or i1 == NT - 1), q, lr_obs


# ------------------------------------------------------------------ the study
def run(cfg, load_cv, b_mode="draw", nrep=20, seeds=range(8), thr_extra=()):
    """
    Coverage and width for: the chi^2_1 profile, j1's theta_hat-bootstrap
    `calibrated`, proper test inversion, and any extra fixed thresholds.
    """
    names = ["profile (chi^2 3.841)", "j1 calibrated (one threshold)",
             "TEST INVERSION (per theta_0)"] + [n for n, _ in thr_extra]
    hit = {n: [] for n in names}
    wid = {n: [] for n in names}
    cen = {n: [] for n in names}
    per_seed_cov = {n: [] for n in names}
    rmse, bnd, ntot = [], [], 0
    for s0 in seeds:
        h = {n: [] for n in names}
        est = []
        nb = 0
        for i in range(nrep):
            m, inc, s, b_t, ke_t = J.record(s0, i, cfg, load_cv, b_mode)
            f = J.fit(m, inc, s, cfg)
            if not (f["theta"] > 0):
                continue
            ntot += 1
            nb += f["on_bound"]
            est.append(np.log(f["theta"]))
            lr = J.profile(m, inc, s, cfg)
            lrm = lr - lr.min()
            iv = {"profile (chi^2 3.841)":
                  J.interval_from(J.LT, lrm, J.CHI2_1_95)}
            for n, t in thr_extra:
                iv[n] = J.interval_from(J.LT, lrm, t)
            rng = np.random.default_rng(300_000 + 1000 * s0 + i)
            _, blr = J.bootstrap(m, inc, s, cfg, f, 199, True, rng)
            qq = float(np.nanquantile(blr, 0.95))
            iv["j1 calibrated (one threshold)"] = J.interval_from(J.LT, lrm, qq)
            lo, hi, cz, _, _ = invert(m, inc, s, cfg, f,
                                      rng=np.random.default_rng(
                                          400_000 + 1000 * s0 + i))
            iv["TEST INVERSION (per theta_0)"] = (lo, hi, cz)
            for n in names:
                a, b, c = iv[n]
                if np.isfinite(a) and np.isfinite(b):
                    h[n].append(a <= np.log(J.TH) <= b)
                    hit[n].append(a <= np.log(J.TH) <= b)
                    wid[n].append(b - a); cen[n].append(c)
        e = np.array(est)
        rmse.append(np.sqrt(np.mean((e - np.log(J.TH)) ** 2)))
        bnd.append(nb / max(len(est), 1))
        for n in names:
            per_seed_cov[n].append(np.mean(h[n]) if h[n] else np.nan)
    return dict(names=names, hit=hit, wid=wid, cen=cen, rmse=rmse,
                bound=bnd, n=ntot, nseed=len(list(seeds)), nrep=nrep,
                per_seed=per_seed_cov)


def show(tag, r):
    print(f"    {tag}")
    print(f"      {'method':>30} {'coverage':>10} {'+/- MC se':>10} "
          f"{'per-seed range':>20} {'median width':>13} {'censored':>9}")
    print("    " + "-" * 98)
    for n in r["names"]:
        h = np.array(r["hit"][n], float)
        c = float(h.mean()) if h.size else np.nan
        ps = np.array(r["per_seed"][n], float)
        print(f"      {n:>30} {c:10.3f} {mc_se(c, h.size):10.3f} "
              f"{f'{np.nanmin(ps):.3f} to {np.nanmax(ps):.3f}':>20} "
              f"{np.median(r['wid'][n]):13.4f} "
              f"{np.mean(r['cen'][n]):9.1%}")
    print(f"      {'':>30} {r['n']} replicates = {r['nseed']} diet "
          f"records x {r['nrep']}; rmse {spread(r['rmse'])}; "
          f"on the b bound {spread(r['bound'], '{:.3f}')}")
    print()


# ------------------------------------------------- (3) what a seed varies
def diet_patterns(cfg, load_cv, nseed=32, nrep=60, target=0.20):
    """
    RMSE per DIET RECORD. A seed draws one load realisation and holds it; the
    replicates inside it vary the observation noise and the drawn nuisances.
    So this is a design question, `for how many diet patterns does the design
    reach the target`, and not Monte Carlo error.
    """
    out = []
    for s0 in range(nseed):
        e = []
        for i in range(nrep):
            m, inc, s, _, _ = J.record(s0, i, cfg, load_cv, "draw")
            f = J.fit(m, inc, s, cfg)
            if f["theta"] > 0:
                e.append(np.log(f["theta"]) - np.log(J.TH))
        e = np.array(e)
        out.append(np.sqrt(np.mean(e ** 2)))
    out = np.array(out)
    return out, float(np.mean(out <= target))


def main():
    print(__doc__)
    base_off = dict(load_prior=False, kappa=None, bound_b=False)
    base_on = dict(load_prior=True, kappa=None, bound_b=False)

    banner("(1) TEST INVERSION AGAINST THE chi^2_1 PROFILE, ON THE MAIN ROWS")
    print("  Monthly n=25, nuisances drawn from their priors, dose exact. Test")
    print(f"  inversion uses {NT} theta_0 nodes and {NBOOT_TI} datasets per node,")
    print("  simulated at the CONSTRAINED fit there, so the cutoff varies along")
    print("  the grid. `j1 calibrated` is the single-threshold version for")
    print("  comparison. Nominal coverage 0.95.\n")
    for cv in (0.10, 0.02):
        for cfg, lab in ((base_off, "load prior off"), (base_on, "load prior ON ")):
            show(f"load CV {cv:.0%}, {lab}, dose exact",
                 run(cfg, cv))

    banner("(2) THE BOUNDARY CUTOFF, b >= 0 WITH b_true ON THE BOUND")
    print("  The usual correction for a parameter on a boundary is a 50:50")
    print(f"  mixture of chi^2_0 and chi^2_1, 95% point {CHI2_MIX_95:.4f} instead of")
    print("  3.8415. Here the parameter ON the boundary is a NUISANCE and the")
    print("  one being tested is not, so whether the correction applies is the")
    print("  question rather than the assumption.\n")
    cfgb = dict(load_prior=False, kappa=None, bound_b=True)
    for b_true, lab in ((0.0, "b_true = 0, ON the bound"),
                        (J.S_B, f"b_true = {J.S_B:.4f}, one prior sd away")):
        show(f"{lab}, load CV 10%, load prior off",
             run(cfgb, 0.10, b_mode=b_true,
                 thr_extra=(("profile (mixture 2.706)", CHI2_MIX_95),)))

    banner("(3) WHAT A SEED VARIES, AND HOW MANY DIET PATTERNS REACH THE TARGET")
    print("  A seed is a DIET RECORD. The load realisation is drawn once per")
    print("  seed and held; the replicates inside it vary the observation noise")
    print("  and the nuisances drawn from their priors. The spread across seeds")
    print("  is therefore a property of the DESIGN and not Monte Carlo error,")
    print("  and the summary that matters is the fraction of diet patterns for")
    print("  which the design reaches the 0.20 target on RMSE.\n")
    print(f"    {'setting':>40} {'median rmse':>12} {'range':>22} "
          f"{'fraction <= 0.20':>17}")
    print("    " + "-" * 96)
    for cv in (0.10, 0.02):
        for cfg, lab in ((base_off, "load prior off"),
                         (base_on, "load prior ON ")):
            v, frac = diet_patterns(cfg, cv)
            print(f"    {f'load CV {cv:.0%}, {lab}, dose exact':>40} "
                  f"{np.median(v):12.4f} "
                  f"{f'{v.min():.4f} to {v.max():.4f}':>22} "
                  f"{f'{frac:.0%} of 32':>17}")
    for kap, lab in ((0.1168, "dose CV 11.7%"),):
        cfg = dict(load_prior=False, kappa=kap, bound_b=False)
        v, frac = diet_patterns(cfg, 0.10)
        print(f"    {f'load CV 10%, load prior off, {lab}':>40} "
              f"{np.median(v):12.4f} {f'{v.min():.4f} to {v.max():.4f}':>22} "
              f"{f'{frac:.0%} of 32':>17}")
    print("\n    32 diet records, 60 replicates each. The 8-seed ranges quoted")
    print("    throughout j1 are the same quantity at a quarter of the")
    print("    resolution, so they understate how often a design misses.")


if __name__ == "__main__":
    main()
