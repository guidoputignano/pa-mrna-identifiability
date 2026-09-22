"""
The headline precision is a curvature. What interval does a record support?

Every precision this paper quotes, 0.1699 monthly and the rest, is
sqrt([F^-1]_11): the Cramer-Rao bound, which is the CURVATURE of the
log-likelihood at the optimum. A reader converts it to an interval by
multiplying by 1.96, and that step is only valid if the likelihood is
quadratic over that range.

`audit/e9_rank_vs_precision.py` and Remark rem:second_order already establish
that this surface is NOT quadratic: there is an exact fibre at a constant
recorded load along which theta and K_m trade off with no change to the
marker, and the specification's day-to-day load variation only bends it rather
than closing it. The paper says, correctly, that the diagnostic that works is
the profile and not a ray. It has never run that profile at the design it
recommends.

So this script profiles. For each simulated record it fixes ln theta on a grid
and re-optimises (ln K_m, b, ln k_e) with the paper's own priors, then reads
the 95% interval off the profile at the usual 1.92 drop in the negative log
likelihood.

WHAT THIS IS NOT. It is not a claim that the published number is wrong. The
CRB correctly describes the sampling spread of theta-hat over records, and
this script reports that spread too, so the two can be compared directly. What
it bounds is the different quantity a reader infers: the interval ONE record
supports.

RE-CHECK. On a NOISELESS record the profile's curvature at the truth must
equal the CRB exactly, since both are the same second derivative. Running that
is how this file demonstrates its optimiser is not the thing being measured;
if the two disagreed, nothing below would mean anything.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import minimize

import n4_theta_scope as N4
from n4_theta_scope import KM0, KE, D0, MONTHLY, footing, delta_at, TH0, TH_DISEASE
from design.run import LOAD_FRAC, ASSAY_CV
from common.models import ParamsA
import model_A.run as LOC

DROP = 1.92          # chi2(1) at 95%, halved, the usual profile cut
STARTS = ([np.log(KM0), 0.0, np.log(KE)],
          [np.log(KM0) + 0.6, 8.0, np.log(KE) - 0.4],
          [np.log(KM0) - 0.6, -8.0, np.log(KE) + 0.4])


def record(theta, seed, D=D0, times=MONTHLY, load_cv=0.10, noisy=True):
    L, m_bar, _, _ = footing(theta, convention="A")
    p = ParamsA(theta=theta, Km=KM0, ke=KE)
    sigma = float(LOC.sigma_total(LOAD_FRAC, ASSAY_CV, p=p, m_bar=m_bar))
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    inc = delta_at(times, D)
    B = KE * KM0 + theta + inc - loads
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * loads * KM0)) / (2.0 * KE)
    y = (theta + inc) * m / (KM0 + m) + KE * m
    if noisy:
        y = y + rng.normal(0.0, sigma, times.size)
    return m, inc, y, sigma, L


def nll(par, m, inc, y, sigma, L, lt_fixed=None):
    """Negative log likelihood, with the paper's own two priors."""
    lt, lk, b, lke = par if lt_fixed is None else (lt_fixed, *par)
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    pred = (th + inc) * m / (K + m) + ke * m + b
    return (0.5 * np.sum(((y - pred) / sigma) ** 2)
            + 0.5 * (b / (LOAD_FRAC * L)) ** 2
            + 0.5 * ((lke - np.log(KE)) / 0.8) ** 2)


def profile_at(lt, args):
    """Re-optimise the other three at this fixed ln theta, from three starts."""
    return min(minimize(nll, np.array(st), args=(*args, lt),
                        method="Nelder-Mead",
                        options=dict(maxiter=40000, maxfev=40000,
                                     xatol=1e-11, fatol=1e-13)).fun
               for st in STARTS)


def interval(theta, seed, span=(-2.5, 1.5), n=161, **kw):
    args = record(theta, seed, **kw)
    grid = np.log(theta) + np.linspace(*span, n)
    v = np.array([profile_at(lt, args) for lt in grid])
    v -= v.min()
    inside = grid[v <= DROP]
    edge = inside.min() <= grid[0] + 1e-9 or inside.max() >= grid[-1] - 1e-9
    return (float(np.exp(inside.min())), float(np.exp(inside.max())),
            float(np.exp(grid[int(np.argmin(v))])), edge)


def recheck():
    print("\n  RE-CHECK. On a noiseless record the profile's curvature at the")
    print("  truth is the same second derivative as the CRB, so they must agree.")
    worst = 0.0
    for theta in (TH0, TH_DISEASE):
        crb = N4.over_seeds(theta=theta, D=D0, convention="A")[0]
        vals = []
        for seed in (500, 503, 507):
            args = record(theta, seed, noisy=False)
            lt0, h = np.log(theta), 0.02
            d2 = (profile_at(lt0 + h, args) - 2 * profile_at(lt0, args)
                  + profile_at(lt0 - h, args)) / h ** 2
            vals.append(1.0 / np.sqrt(d2))
        got = float(np.median(vals))
        worst = max(worst, abs(got - crb) / crb)
        print(f"    theta {theta:6.0f}: profile curvature sd {got:.4f} "
              f"against CRB {crb:.4f}, ratio {got/crb:.4f}")
    print(f"    worst relative disagreement {worst:.3%} against 1%")
    ok = worst < 0.01
    print(f"    RE-CHECK {'PASSES' if ok else 'FAILS'}")
    return ok


def main():
    print("=" * 88)
    print("  WHAT INTERVAL DOES ONE RECORD SUPPORT, AGAINST WHAT THE CRB IMPLIES")
    print("=" * 88)
    for theta, span, n in ((TH0, (-2.5, 1.5), 161), (TH_DISEASE, (-4.0, 4.0), 81)):
        crb = N4.over_seeds(theta=theta, D=D0, convention="A")[0]
        implied = np.exp(2 * 1.96 * crb)
        print(f"\n  theta = {theta:.0f} ({100*theta/TH0:.0f}% of nominal). "
              f"CRB sd(ln theta) = {crb:.4f},")
        print(f"  so a reader infers a 95% interval a factor of "
              f"{implied:.2f} wide.")
        seeds = range(600, 616) if theta == TH0 else range(600, 606)
        w, mles, his, edges = [], [], [], 0
        for s in seeds:
            lo, hi, mle, edge = interval(theta, s, span=span, n=n)
            w.append(hi / lo)
            mles.append(mle)
            his.append(hi)
            edges += bool(edge)
            print(f"    seed {s}: MLE {mle:9.1f}  95% profile "
                  f"[{lo:9.1f}, {hi:9.1f}]  factor {hi/lo:8.2f}"
                  f"{'  GRID EDGE' if edge else ''}")
        w = np.array(w)
        print(f"    over {len(w)} records: median width a factor of "
              f"{np.median(w):.2f} (range {w.min():.2f} to {w.max():.2f}),")
        print(f"    against the {implied:.2f} the CRB implies.")
        if theta == TH0:
            #  Remark 2 and the introduction both quote these, and both
            #  quoted them as literals
            from common.numbers import write_numbers
            write_numbers("profile",
                          {"profile_factor": float(np.median(w)),
                           "profile_factor_lo": float(w.min()),
                           "profile_factor_hi": float(w.max()),
                           "crb_factor": float(implied)},
                          {"profile_factor": "{:.2f}",
                           "profile_factor_lo": "{:.2f}",
                           "profile_factor_hi": "{:.2f}",
                           "crb_factor": "{:.2f}"})
        if edges:
            print(f"    {edges} of {len(w)} intervals reach the grid edge, so "
                  f"those widths are LOWER")
            print(f"    bounds and the ratio to the CRB is not meaningful here: "
                  f"the record does")
            print(f"    not bound theta from below AT ALL. The MLE itself runs "
                  f"{min(mles):.0f} to {max(mles):.0f}")
            print(f"    against a truth of {theta:.0f}, so the point estimate "
                  f"is not usable either.")
            print(f"    What the record DOES supply is an upper bound: it "
                  f"rules out capacity above")
            print(f"    about {np.median(his):.0f}, which is "
                  f"{100*np.median(his)/TH0:.0f}% of nominal. That is the half "
                  f"of the question a")
            print(f"    clinician can still be given an answer to.")
        else:
            print(f"    ratio of the honest interval to the CRB interval: "
                  f"{np.median(w)/implied:.2f}.")
            print(f"    MLE spread over records {min(mles):.0f} to "
                  f"{max(mles):.0f}, truth {theta:.0f}: the CRB describes")
            print(f"    THAT correctly, and understates the interval. The "
                  f"asymmetry is the fibre:")
            print(f"    lowering theta is compensated by raising K_m, so the "
                  f"profile has a long")
            print(f"    tail downward while its upper end sits near where the "
                  f"CRB puts it.")
    ok = recheck()
    print("\n" + "=" * 88)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
