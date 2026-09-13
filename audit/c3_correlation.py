"""
C3 -- is corr(ln theta, ln K_m) flipping +0.9988 -> -0.3163 a real reversal, or
a coincidence of one parameter set?

The mechanism, if it is structural: the ratio of the two corrected columns is

    q_Km / q_theta = -(1 + Delta/theta) * K_m/(K_m+m) ,

so the therapy moves the ratio directly, while the marker band moves it only
through phi(m) = K_m/(K_m+m).  The published column drops the first factor, so
the ratio varies only as phi does, which over a managed band is a few per cent
-- hence near-perfect collinearity.  If that is the mechanism, the sign of the
correlation should track the RANGE of (1+Delta/theta) against the range of
phi(m), and not the particular value of theta, K_m or k_e.

Swept below: operating point m_bar/K_m, theta, k_e, dose amplitude, dose
interval, sampling interval, and record length.  Also Monte Carlo across seeds,
since the covariance itself is deterministic and "seeds" only enter the fit.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import least_squares
import design.run as D, model_A.run as LOC
from common.models import ParamsA, clearance_A

LOAD_FRAC = D.LOAD_FRAC


def record(theta, Km, ke, m_bar, d_amp, dose_int, dt, T, lam=None):
    """Levels, increments and the recorded load for one design, from scratch."""
    lam = np.log(2.0) / 3.5 if lam is None else lam
    p = ParamsA(theta=theta, Km=Km, ke=ke)
    Lb = clearance_A(m_bar, p)
    t = np.arange(0.0, T, dt)
    # closed-form geometric dose train, as design/run.py:delta
    n = np.floor(t / dose_int)
    frac = t - n * dose_int
    r = np.exp(-lam * dose_int)
    dd = d_amp * np.exp(-lam * frac) * (1.0 - r ** (n + 1)) / (1.0 - r)
    B = ke * Km + theta + dd - Lb
    lv = (-B + np.sqrt(B ** 2 + 4.0 * ke * Lb * Km)) / (2.0 * ke)
    return lv, dd, Lb, p


def cov(lv, dd, theta, Km, ke, Lb, sigma, use_delta, ke_free, b_frac=LOAD_FRAC):
    pri = [0.8 ** 2, 0.8 ** 2, (b_frac * Lb) ** 2] + ([0.8 ** 2] if ke_free else [])
    Finv = np.linalg.inv(np.diag(pri))
    for m, dl in zip(lv, dd):
        eff = theta + (dl if use_delta else 0.0)
        q = [theta * m / (Km + m), -eff * Km * m / (Km + m) ** 2, 1.0]
        if ke_free:
            q.append(ke * m)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / sigma ** 2
    C = np.linalg.inv(Finv)
    return float(np.sqrt(C[0, 0])), float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]))


def row(label, theta, Km, ke, m_bar, d_amp=6000.0, dose_int=14.0, dt=2.0, T=730.0,
        ke_free=False):
    lv, dd, Lb, p = record(theta, Km, ke, m_bar, d_amp, dose_int, dt, T)
    sigma = LOC.sigma_total(LOAD_FRAC, D.ASSAY_CV, p=p, m_bar=m_bar)
    a, ra = cov(lv, dd, theta, Km, ke, Lb, sigma, False, ke_free)
    b, rb = cov(lv, dd, theta, Km, ke, Lb, sigma, True, ke_free)
    rng_phi = Km / (Km + lv.min()) - Km / (Km + lv.max())
    rng_dlt = (dd.max() - dd.min()) / theta
    print(f"{label:>28} | {m_bar/Km:7.3f} {rng_dlt:8.3f} {rng_phi:8.4f} | "
          f"{a:8.4f} {ra:+8.4f} | {b:8.4f} {rb:+8.4f} | "
          f"{'FLIPS' if rb < 0 < ra else ('drops' if rb < ra - 0.05 else 'no')}")
    return ra, rb


def main():
    print(__doc__)
    hdr = (f"{'case':>28} | {'mbar/Km':>7} {'rng D/th':>8} {'rng phi':>8} | "
           f"{'sd pub':>8} {'corr pub':>8} | {'sd corr':>8} {'corr corr':>8} | verdict")
    for ke_free in (False, True):
        print("=" * len(hdr))
        print(f"k_e {'FREE (0.8 log-prior) -- the footing A1 recommends' if ke_free else 'KNOWN -- the published footing'}")
        print("=" * len(hdr)); print(hdr); print("-" * len(hdr))
        row("published base case", 3000., 100., 5., 10., ke_free=ke_free)
        print(f"{'--- operating point':>28} |")
        for f in (0.02, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0):
            row(f"m_bar/K_m = {f}", 3000., 100., 5., f * 100., ke_free=ke_free)
        print(f"{'--- theta':>28} |")
        for th in (300., 1000., 3000., 10000., 30000.):
            row(f"theta = {th:.0f}", th, 100., 5., 10., ke_free=ke_free)
        print(f"{'--- k_e':>28} |")
        for k in (0.5, 2., 5., 20., 50.):
            row(f"k_e = {k}", 3000., 100., k, 10., ke_free=ke_free)
        print(f"{'--- dose amplitude':>28} |")
        for a in (600., 1500., 3000., 6000., 12000.):
            row(f"D_AMP = {a:.0f} (D/th {a/3000:.1f})", 3000., 100., 5., 10.,
                d_amp=a, ke_free=ke_free)
        print(f"{'--- schedule':>28} |")
        for di, dt_, T in ((7., 2., 730.), (14., 2., 730.), (28., 2., 730.),
                           (14., 1., 730.), (14., 7., 730.), (14., 2., 28.),
                           (14., 2., 14.)):
            row(f"dose {di:.0f}d, sample {dt_:.0f}d, {T:.0f}d",
                3000., 100., 5., 10., dose_int=di, dt=dt_, T=T, ke_free=ke_free)
        print()

    print("=" * len(hdr))
    print("ACROSS SEEDS -- the covariance is deterministic, so seeds can only enter")
    print("the fit.  Monte Carlo of the corrected estimator, 5 seeds, n=365:")
    print("=" * len(hdr))
    lv, dd, Lb, p = record(3000., 100., 5., 10., 6000., 14., 2., 730.)
    sigma = LOC.sigma_total(LOAD_FRAC, D.ASSAY_CV)
    s_true = (3000. + dd) * lv / (100. + lv) + 5. * lv
    fisher_sd, fisher_r = cov(lv, dd, 3000., 100., 5., Lb, sigma, True, False)
    for seed in range(5):
        rng = np.random.default_rng(seed)
        ests = []
        for _ in range(200):
            s = s_true + rng.normal(0., sigma, lv.size)
            # MAP: data residual stacked with the module's own priors
            def res(x):
                th, km = 3000. * np.exp(x[0]), 100. * np.exp(x[1])
                pred = (th + dd) * lv / (km + lv) + 5. * lv + x[2]
                return np.concatenate([(pred - s) / sigma,
                                       [x[0] / 0.8, x[1] / 0.8,
                                        x[2] / (LOAD_FRAC * Lb)]])
            r = least_squares(res, np.zeros(3), method="lm", xtol=1e-14, ftol=1e-14)
            ests.append(r.x[:2])
        e = np.array(ests)
        print(f"  seed {seed}: MC sd(ln theta) = {e[:,0].std(ddof=1):.4f}  "
              f"corr = {np.corrcoef(e.T)[0,1]:+.4f}   "
              f"(Fisher: {fisher_sd:.4f}, {fisher_r:+.4f})")


if __name__ == "__main__":
    main()
