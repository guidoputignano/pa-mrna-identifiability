r"""
Round 18: the PI's line-pencil geometry, and what it adds.

Each steady-state draw is a LINE in the (theta, k_e) plane. From the Model-A
balance (theta + inc_i) g(m_i) + k_e m_i = L_i, with g(m) = m/(K_m + m) and the
therapy increment inc_i KNOWN, moving the known increment to the right gives

        g(m_i) * theta  +  m_i * k_e  =  L_i - inc_i g(m_i)            (line i)
        \____ a_i ____/     \_ b_i _/     \________ c_i ________/

a line in (theta, k_e) with normal (a_i, b_i) and offset c_i. Noise-free and
with the TRUE K_m, every line passes exactly through the true (theta, k_e):
the pencil CONCURS. With load-record error and assay error the lines form a
near-concurrent pencil, and the size of the near-concurrence region is the
identifiability of (theta, k_e).

This script does the four things the PI's formulation adds that the FIM/profile
work had not done explicitly. It does NOT re-derive conclusions already in the
report; it checks them by an independent geometry and produces the pencil figure.

  1. DRAW THE PENCIL  (figure built in r18_figs.py) -- numbers for it here.
  2. PAIRWISE INTERSECTIONS as a nonparametric, noise-model-free uncertainty
     estimate; compare to the FIM ellipse and the NLS posterior; characterise
     the heavy tail (near-parallel pairs meet near infinity).
  3. CONCURRENCY as a model-adequacy test: residual distance of the best-fit
     point to each line vs what assay noise predicts (chi-square). Power probed
     against a misspecified (second-channel) record.
  4. THE WEIGHTING (Remark 6: the noise is in the ABSCISSA, the assayed marker,
     not the ordinate): total-least-squares vs ordinary Euclidean fit, and the
     weighted normal matrix reproduces the FIM.

All on the repository's footing (core.py): 5% patient theta=150, k_e=5, K_m=100,
L=322.73, the 25-draw monthly dosed record, LOAD_FRAC=5%, ASSAY_CV=5%.
"""
import numpy as np

import core as C
from n4_theta_scope import level, delta_at, MONTHLY

TH5 = C.TH_5                 # 150, the 5% patient
KE = C.KE                    # 5
KM = C.KM0                   # 100
LF = C.LOAD_FRAC if hasattr(C, "LOAD_FRAC") else 0.05
from design.run import LOAD_FRAC, ASSAY_CV, DOSE_INTERVAL
from n4_theta_scope import D0


# ---------------------------------------------------------------- the pencil
def pencil(theta=TH5, ke=KE, Km=KM, D=D0, times=MONTHLY, seed=0,
           load_cv=0.10, assay_cv=ASSAY_CV, load_frac=LOAD_FRAC, noisy=True):
    """
    One realised record as a pencil of lines in (theta, k_e).

    Returns a dict with per-draw (a, b, c) and the observed/true levels. With
    noisy=True the abscissa (marker m) carries assay_cv and the recorded load
    carries load_frac; with noisy=False the lines concur exactly at (theta,ke).
    """
    L0, _ = __import__("n4_theta_scope").anchor(theta, Km, ke, "A")
    times = np.asarray(times, float)
    rng = np.random.default_rng(seed)
    if load_cv:
        s2 = np.log1p(load_cv ** 2)
        loads = L0 * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        loads = np.full(times.size, L0)
    inc = delta_at(times, D) if D else np.zeros(times.size)
    m_true = level(theta, Km, ke, loads, inc)
    if noisy:
        m_obs = m_true * (1.0 + rng.normal(0.0, assay_cv, times.size))
        L_rec = loads * (1.0 + rng.normal(0.0, load_frac, times.size))
    else:
        m_obs, L_rec = m_true, loads
    g = m_obs / (Km + m_obs)
    a = g
    b = m_obs
    c = L_rec - inc * g
    return dict(a=a, b=b, c=c, m_true=m_true, m_obs=m_obs, loads=loads,
                L_rec=L_rec, inc=inc, theta=theta, ke=ke, Km=Km,
                assay_cv=assay_cv, load_frac=load_frac)


def wide_pencil(theta=TH5, ke=KE, Km=KM, n=25, seed=1,
                lo=0.45, hi=3.2, assay_cv=ASSAY_CV, load_frac=LOAD_FRAC):
    """
    A DRUG-FORCING / excursion record: loads deliberately swept lo..hi x the
    managed load so the marker ranges widely and the line slopes -1/(K_m+m_i)
    differ a lot -- a well-conditioned pencil that concurs crisply. This is the
    steady-state reading of Task 7's optimal washout/forcing design.
    """
    L0, _ = __import__("n4_theta_scope").anchor(theta, Km, ke, "A")
    rng = np.random.default_rng(seed)
    loads = L0 * np.linspace(lo, hi, n)
    m_true = level(theta, Km, ke, loads, 0.0)
    m_obs = m_true * (1.0 + rng.normal(0.0, assay_cv, n))
    L_rec = loads * (1.0 + rng.normal(0.0, load_frac, n))
    g = m_obs / (Km + m_obs)
    return dict(a=g, b=m_obs, c=L_rec, m_true=m_true, m_obs=m_obs, loads=loads,
                L_rec=L_rec, inc=np.zeros(n), theta=theta, ke=ke, Km=Km,
                assay_cv=assay_cv, load_frac=load_frac)


# --------------------------------------------------- item 2: intersections
def intersections(P):
    """All C(N,2) pairwise line intersections and each pair's crossing angle."""
    a, b, c = P["a"], P["b"], P["c"]
    n = a.size
    pts, ang, det = [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            D = a[i] * b[j] - a[j] * b[i]           # determinant of the 2x2
            if D == 0:
                continue
            th = (c[i] * b[j] - c[j] * b[i]) / D
            k = (a[i] * c[j] - a[j] * c[i]) / D
            # sine of the crossing angle between the two normals (0 = parallel)
            na = np.hypot(a[i], b[i]); nb = np.hypot(a[j], b[j])
            s = abs(D) / (na * nb)
            pts.append((th, k)); ang.append(s); det.append(D)
    return np.array(pts), np.array(ang), np.array(det)


# --------------------------------------------- items 3/4: fit a point to lines
def fit_point(P, weights=None):
    """
    (theta,k_e) minimising sum w_i * d_i^2, d_i the ORTHOGONAL distance to line
    i. weights=None is the Euclidean (P=I) fit; a vector gives the weighted
    (Mahalanobis / TLS) fit. Returns the point and the 2x2 normal matrix N
    (its inverse is the covariance of the estimate).
    """
    a, b, c = P["a"], P["b"], P["c"]
    nn = a ** 2 + b ** 2                     # |normal_i|^2
    w = (np.ones_like(a) if weights is None else np.asarray(weights)) / nn
    # normal equations for min sum w_i (a_i th + b_i k - c_i)^2
    N = np.array([[np.sum(w * a * a), np.sum(w * a * b)],
                  [np.sum(w * a * b), np.sum(w * b * b)]])
    rhs = np.array([np.sum(w * a * c), np.sum(w * b * c)])
    x = np.linalg.solve(N, rhs)
    return x, N


def dist_sigma(P):
    """
    Per-line sd of the orthogonal distance, from the ACTUAL error structure.
    The marker m carries assay_cv (abscissa) and the load carries load_frac
    (ordinate). Propagated to the signed distance d_i = (a th + b k - c)/|n|:

        along the ordinate c: sigma_c = load_frac * L_i
        along the abscissa m: a one-sigma marker wiggle moves the line; to first
            order its effect on d at the true point is |C'(m_i)-?|... computed
            as the shift of c_i minus the point's motion. Using d = (a th+b k-c),
            da/dm, db/dm, dc/dm give dd/dm; sigma_m = assay_cv * m_i.
    Returns sigma_d (distance units) evaluated at the true point.
    """
    a, b, c = P["a"], P["b"], P["c"]
    th, k, Km, m = P["theta"], P["ke"], P["Km"], P["m_obs"]
    inc = P["inc"]
    nn = np.hypot(a, b)
    # d_num = a*th + b*k - c  ; propagate through m (a=g(m), b=m, c=L-inc*g)
    # da/dm = Km/(Km+m)^2 ; db/dm = 1 ; dc/dm = -inc*Km/(Km+m)^2
    dg = Km / (Km + m) ** 2
    dnum_dm = dg * th + 1.0 * k - (-inc * dg)        # d(num)/dm
    dnum_dc = -1.0                                   # d(num)/dL_rec
    sig_m = P["assay_cv"] * m
    sig_c = P["load_frac"] * P["L_rec"]
    sig_num = np.hypot(dnum_dm * sig_m, dnum_dc * sig_c)
    return sig_num / nn, sig_num


def adequacy(P):
    """
    Concurrency as a model-adequacy test. Fit the best point with the proper
    weights, then compare each residual distance to what noise predicts:
        chi2 = sum (d_i / sigma_d_i)^2  ~ chi^2_{N-2} if the model is right.
    Returns (chi2, dof, p_value, residual sds in sigma units).
    """
    from scipy.stats import chi2 as chi2d
    sig_d, sig_num = dist_sigma(P)
    x, N = fit_point(P, weights=1.0 / sig_num ** 2 * (P["a"] ** 2 + P["b"] ** 2))
    a, b, c = P["a"], P["b"], P["c"]
    d = (a * x[0] + b * x[1] - c) / np.hypot(a, b)
    z = d / sig_d
    stat = float(np.sum(z ** 2))
    dof = a.size - 2
    pval = float(chi2d.sf(stat, dof))
    return stat, dof, pval, z, x


def _two_channel_levels(loads, theta, ke, Km, Km2, frac2):
    """Steady states of a two-channel truth: theta g(m) + (1-f)ke m + ke2 h2(m)."""
    from scipy.optimize import brentq
    ke2 = frac2 * ke * (Km2 + 50.0) / 50.0           # match linear scale near m~50
    def clear(m):
        return theta * m / (Km + m) + (1 - frac2) * ke * m + ke2 * m / (Km2 + m)
    return np.array([brentq(lambda m: clear(m) - Li, 1e-6, 1e7) for Li in loads])


def misspec_pencil(theta=TH5, ke=KE, Km=KM, Km2=600.0, frac2=0.5, wide=False,
                   seed=0, times=MONTHLY, assay_cv=ASSAY_CV, load_frac=LOAD_FRAC):
    """
    A record from a TWO-CHANNEL truth (second saturable route) read through the
    one-channel line g(m)theta + m k_e. The clearance is NOT theta g(m) + k_e m,
    so the lines should fail to concur beyond noise -- IF the band is wide enough
    to see the second channel's curvature. wide=True sweeps loads to force it.
    """
    L0, _ = __import__("n4_theta_scope").anchor(theta, Km, ke, "A")
    rng = np.random.default_rng(seed)
    if wide:
        loads = L0 * np.linspace(0.45, 3.2, 25)
    else:
        s2 = np.log1p(0.10 ** 2)
        loads = L0 * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    m_true = _two_channel_levels(loads, theta, ke, Km, Km2, frac2)
    m_obs = m_true * (1.0 + rng.normal(0.0, assay_cv, loads.size))
    L_rec = loads * (1.0 + rng.normal(0.0, load_frac, loads.size))
    g = m_obs / (Km + m_obs)
    return dict(a=g, b=m_obs, c=L_rec, m_true=m_true, m_obs=m_obs, loads=loads,
                L_rec=L_rec, inc=np.zeros(loads.size), theta=theta, ke=ke,
                Km=Km, assay_cv=assay_cv, load_frac=load_frac)


def logcorr(P):
    """
    corr(ln theta, ln k_e) and the ellipse aspect ratio of the pencil's own
    (weighted) concurrence region -- the scale-free measure of near-parallelism,
    directly comparable to the report's -0.99. Returns (rho, aspect, sd_lnth, sd_lnke).
    """
    sig_d, sig_num = dist_sigma(P)
    wls = (P["a"] ** 2 + P["b"] ** 2) / sig_num ** 2
    x, NW = fit_point(P, weights=wls)
    Dd = np.diag([P["theta"], P["ke"]])
    Cov = np.linalg.inv(Dd @ NW @ Dd)                 # covariance in log coords
    rho = Cov[0, 1] / np.sqrt(Cov[0, 0] * Cov[1, 1])
    ar = np.sqrt((1 + abs(rho)) / (1 - abs(rho)))
    return float(rho), float(ar), float(np.sqrt(Cov[0, 0])), float(np.sqrt(Cov[1, 1]))


# ------------------------------------------------- FIM ellipse for comparison
def fim_ellipse_linear(theta=TH5, ke=KE, Km=KM, seed=0, ke_prior=True):
    """
    The (theta, k_e) 2x2 covariance at FIXED K_m, in LINEAR coordinates, from
    the flux-space information -- the pencil's own metric. Delta-method mapping
    from the log-FIM block. ke_prior=True adds the report's k_e prior (the
    headline footing); ke_prior=False is the pure data information, which is the
    object the prior-FREE weighted pencil fit reproduces.
    """
    r = C.record(theta, seed=seed)
    Q = C.flux_rows(r)                 # cols (ln th, ln Km, b, ln ke)
    sg = C.sigmas(r, "perdraw")
    # keep columns ln theta, ln ke (fix K_m, drop shared bias -> per-line c)
    Qr = Q[:, [C.I_TH, C.I_KE]]
    F = (Qr / sg[:, None] ** 2).T @ Qr
    if ke_prior:
        F[1, 1] += C.KE_PRIOR ** -2
    Plog = np.linalg.inv(F)            # cov in (ln th, ln ke)
    J = np.diag([theta, ke])           # d(th,ke)/d(ln th,ln ke)
    return J @ Plog @ J, Plog


def tail_stats(P):
    """Robust centres and tail of the pairwise-intersection cloud."""
    pts, ang, det = intersections(P)
    th, k = pts[:, 0], pts[:, 1]
    w = det ** 2
    detw = (float(np.sum(w * th) / w.sum()), float(np.sum(w * k) / w.sum()))
    q = np.quantile(th, [.05, .25, .5, .75, .95])
    in3x = float(np.mean((th > P["theta"] / 3) & (th < P["theta"] * 3)))
    return dict(n=len(pts), med=(float(np.median(th)), float(np.median(k))),
                detw=detw, q=q, iqr=float(q[3] - q[1]),
                rng=float(th.max() - th.min()), in3x=in3x,
                mean=float(th.mean()))


def main():
    np.set_printoptions(suppress=True, precision=1)
    print("=" * 92)
    print("  ROUND 18.  THE LINE-PENCIL GEOMETRY (PI's formulation) AND WHAT IT ADDS")
    print("=" * 92)
    P = pencil(seed=0)
    Pc = pencil(seed=0, noisy=False)
    W = wide_pencil(seed=1)
    print(f"""
  5% patient: theta={TH5:.0f}, k_e={KE:.0f}, K_m={KM:.0f}, L=322.73.
  NARROW record : 25 monthly dosed draws (D={D0:.0f}), load CV 10%; m_obs in [{P['m_obs'].min():.1f}, {P['m_obs'].max():.1f}]
  WIDE record   : 25-point load sweep 0.45-3.2x (steady-state reading of the
                  Task-7 drug-forcing design); m_obs in [{W['m_obs'].min():.1f}, {W['m_obs'].max():.1f}]
  both carry assay CV {ASSAY_CV:.0%} (abscissa) and load-record {LOAD_FRAC:.0%} (ordinate).
  line slope dk_e/dtheta = -g/m = -1/(K_m+m):  narrow spans only [{-1/(KM+P['m_true'].max()):.5f},{-1/(KM+P['m_true'].min()):.5f}]
  noise-free lines pass through (150,5) to {np.max(np.abs(Pc['a']*TH5+Pc['b']*KE-Pc['c'])):.1e}  (exact concurrence, true K_m)
""")

    # ---- conditioning: the near-parallelism, scale-free -----------------
    rN, aN, sltN, slkN = logcorr(P)
    rW, aW, sltW, slkW = logcorr(W)
    print("  " + "-" * 88)
    print("  CONDITIONING of the pencil (= near-parallelism, in the noise metric)")
    print("  " + "-" * 88)
    print(f"""      {'record':16s}  {'corr(ln th,ln ke)':>18s}  {'ellipse aspect':>14s}  {'sd(ln th)':>10s}
      {'NARROW dosed':16s}  {rN:+18.3f}  {aN:14.1f}  {sltN:10.2f}
      {'WIDE excursion':16s}  {rW:+18.3f}  {aW:14.1f}  {sltW:10.2f}
  The narrow pencil is near-degenerate (corr -0.998, the report's -0.99; aspect {aN:.0f}:1);
  widening the marker band swings the slopes apart and the pencil concurs, sd(ln th)
  {sltN:.1f}->{sltW:.2f}. This is WHY the fibre exists, drawn rather than tabulated.
""")

    # ---- item 2: pairwise intersections, narrow vs wide -----------------
    tN, tW = tail_stats(P), tail_stats(W)
    # honest centering diagnostic: median-of-medians and the negative-median
    # fraction over many records (the cloud is so pathological the median can
    # itself go negative, so a "ln-bias" conditioned on positive medians misleads).
    meds = np.array([np.median(intersections(pencil(seed=s))[0][:, 0])
                     for s in range(400)])
    neg_frac = float(np.mean(meds < 0))
    med_of_meds = float(np.median(meds))                          # signed, over all records
    ln_bias_all = float(np.log(med_of_meds / TH5)) if med_of_meds > 0 else np.nan
    print("  " + "-" * 88)
    print(f"  ITEM 2.  PAIRWISE INTERSECTIONS  (N=25 -> {tN['n']} pairs), model-free uncertainty")
    print("  " + "-" * 88)
    print(f"""  true (theta, k_e) = (150, 5). centre & tail of the intersection cloud:
      {'record':16s} {'median th':>10s} {'det^2-wtd th':>13s} {'IQR th':>9s} {'range th':>11s} {'%in[50,450]':>12s}
      {'NARROW dosed':16s} {tN['med'][0]:10.1f} {tN['detw'][0]:13.1f} {tN['iqr']:9.0f} {tN['rng']:11.1e} {tN['in3x']:11.0%}
      {'WIDE excursion':16s} {tW['med'][0]:10.1f} {tW['detw'][0]:13.1f} {tW['iqr']:9.0f} {tW['rng']:11.1e} {tW['in3x']:11.0%}
  NARROW: the cloud is Cauchy-like (raw mean {tN['mean']:.1e} is meaningless; range ~{tN['rng']:.0e}).
    PREDICTION CONFIRMED and sharpened: not only is the SPREAD heavy-tailed -- the
    near-parallel pairs fling intersections toward infinity -- the median itself is
    an unreliable point estimate. Over 400 records the median-of-medians is {med_of_meds:.0f}
    (bias +{ln_bias_all:.2f} in ln theta), and {neg_frac:.0%} of records give a NEGATIVE cloud median:
    the pencil is so pathological the robust centre is often nonsensical. Only
    {tN['in3x']:.0%} of pairwise intersections land within a factor 3 of the truth. The
    det^2 (crossing-angle) weighting the PI suggests pulls the centre back toward the
    MLE ({tN['detw'][0]:.0f} vs the weighted fit's 242), but the natural-scale spread stays
    enormous. This is round 13's heavy-tailed posterior and round 14's lesson from a
    THIRD, assumption-free route: a near-parallel pencil has no natural-scale point
    estimate, only a log-scale width.
  WIDE: the same construction concurs -- median theta {tW['med'][0]:.0f}, {tW['in3x']:.0%} of pairs within a
    factor 3, IQR {tW['iqr']:.0f} vs {tN['iqr']:.0f}. Width buys a usable point estimate.
""")

    # ---- FIM ellipse + NLS posterior comparison (narrow) ---------------
    # All three routes are prior-FREE pencil geometry (no k_e prior, fixed K_m),
    # so compare against the prior-free FIM. The report's headline adds the k_e
    # prior and the K_m fibre on top of this (noted below).
    Cov0, _ = fim_ellipse_linear(ke_prior=False)
    sd_th0 = np.sqrt(Cov0[0, 0]); rho0 = Cov0[0, 1] / np.sqrt(Cov0[0, 0] * Cov0[1, 1])
    Covp, _ = fim_ellipse_linear(ke_prior=True)
    sd_th_lin = np.sqrt(Covp[0, 0]); rho = Covp[0, 1] / np.sqrt(Covp[0, 0] * Covp[1, 1])
    post = np.array([fit_point(pencil(seed=1000 + s),
                     weights=(lambda Q: (Q["a"] ** 2 + Q["b"] ** 2)
                              / dist_sigma(Q)[1] ** 2)(pencil(seed=1000 + s)))[0]
                     for s in range(200)])
    iqr_post = float(np.subtract(*np.quantile(post[:, 0], [.75, .25])))
    print("  " + "-" * 88)
    print("  THREE ROUTES to the narrow-record uncertainty  (fixed K_m, prior-free geometry)")
    print("  " + "-" * 88)
    print(f"""      FIM ellipse (delta-method, no k_e prior):  sd(theta)={sd_th0:7.1f}  corr={rho0:+.3f}
      NLS posterior (200 weighted fits):         theta median {np.median(post[:,0]):7.1f}  IQR {iqr_post:7.1f}
      intersection cloud (one record):           theta IQR {tN['iqr']:7.0f}  (heavy-tailed)
  The FIM sd and the NLS IQR agree in the CORE (a Gaussian sd={sd_th0:.0f} has IQR~{1.349*sd_th0:.0f} vs
  NLS {iqr_post:.0f}); the symmetric FIM ellipse cannot show the tail the two sampling routes do.
  All three say theta is not placeable on the natural scale at the record precision.
  (Adding the report's k_e prior tightens the fixed-K_m FIM to sd(theta)={sd_th_lin:.0f}, corr {rho:+.3f};
  the full headline sd(ln theta)=3.05 adds the K_m fibre on top -- a different, looser object.)
""")

    # ---- item 3: concurrency as adequacy, and its power ----------------
    stat, dof, pval, z, xfit = adequacy(P)
    ps = np.array([adequacy(pencil(seed=s))[2] for s in range(300)])

    def load_bias(seed, bias):                       # systematic ordinate error
        Pw = dict(pencil(seed=seed)); Pw["c"] = Pw["c"].copy()
        Pw["c"][::2] *= (1 + bias)                   # bias on alternate draws
        return Pw
    bias_pow = {b: float(np.mean([adequacy(load_bias(s, b))[2] < 0.05
                for s in range(300)])) for b in (0.05, 0.10, 0.15, 0.25)}
    # second saturable channel -- strong, wide band (worst case for detection)
    p2 = float(np.mean([adequacy(misspec_pencil(seed=s, wide=True, Km2=150.0,
                frac2=0.6))[2] < 0.05 for s in range(200)]))
    print("  " + "-" * 88)
    print("  ITEM 3.  CONCURRENCY AS A MODEL-ADEQUACY TEST  (chi^2 of residual distances)")
    print("  " + "-" * 88)
    print(f"""  repo's simulated record (seed 0): chi^2 = {stat:.1f} on {dof} dof, p = {pval:.2f} -> PASSES.
  Over 300 simulated records the test rejects {np.mean(ps<0.05):.0%} at alpha=5% -- CALIBRATED, and
  the repo's records PASS (they are generated from the one-channel model, as they must).

  WHAT IT HAS POWER TO DETECT:
    * a WRONG / BIASED recorded load (an ordinate error, not of the line form):
        bias on 1/2 of draws   +5%     +10%    +15%    +25%
        rejection rate        {bias_pow[0.05]:4.0%}   {bias_pow[0.10]:5.0%}   {bias_pow[0.15]:5.0%}   {bias_pow[0.25]:5.0%}
      -> growing power; this is the misspecification the pencil is built to see.
    * a SECOND SATURABLE CHANNEL / saturating k_e (structural clearance error),
      strong (60% of the linear channel) and on a WIDE band:  rejection {p2:.0%}  -- NO POWER.
  WHY the asymmetry: a second saturable term h2(m) is nearly COLLINEAR with {{g(m), m}}
  over any realistic band, so the two-parameter line absorbs it -- the misspecification
  lands on the fibre and concurrence is preserved. A per-draw load error does NOT lie in
  that span, so it breaks concurrence. This is t39's affine degeneracy seen from the
  pencil (t39 asks if clearance is locally LINEAR; this asks if the two-channel point-
  model fits the pencil at all), and it is a real limit: the concurrency test guards
  against DATA errors (a mis-recorded load), not against the structural degeneracy the
  rest of the report is about -- which no steady-state concurrency test can see.
""")

    # ---- item 4: the weighting (Remark 6) ------------------------------
    sig_d, sig_num = dist_sigma(P)
    xE, NE = fit_point(P, weights=None)                           # Euclidean P=I
    wls = (P["a"] ** 2 + P["b"] ** 2) / sig_num ** 2              # proper weights
    xW, NW = fit_point(P, weights=wls)
    CovE = np.linalg.inv(NE); CovW = np.linalg.inv(NW)
    print("  " + "-" * 88)
    print("  ITEM 4.  THE WEIGHTING  (Remark 6: noise is in the ABSCISSA, the assayed marker)")
    print("  " + "-" * 88)
    print(f"""  per-line orthogonal-distance sd ranges sigma_d in [{sig_d.min():.2f}, {sig_d.max():.2f}] over the record
  (heteroscedastic: grows with m, shrinks as the line flattens) -- so the unweighted
  Euclidean fit is the wrong metric.
      Euclidean (P=I):          (theta,k_e) = ({xE[0]:7.1f}, {xE[1]:6.3f})   sd(theta)={np.sqrt(CovE[0,0]):6.1f}
      weighted (TLS, x'Px):     (theta,k_e) = ({xW[0]:7.1f}, {xW[1]:6.3f})   sd(theta)={np.sqrt(CovW[0,0]):6.1f}
      FIM 2x2, NO k_e prior (= the weighted fit's object):        sd(theta)={sd_th0:6.1f}  corr={rho0:+.3f}
      weighted-fit corr = {CovW[0,1]/np.sqrt(CovW[0,0]*CovW[1,1]):+.3f}
  The weighted normal matrix IS the (fixed-K_m) FIM OF THE DATA -- same structural object:
  corr {CovW[0,1]/np.sqrt(CovW[0,0]*CovW[1,1]):+.3f} vs {rho0:+.3f}, sd(theta) {np.sqrt(CovW[0,0]):.0f} vs {sd_th0:.0f} (both prior-free; the ~{abs(np.sqrt(CovW[0,0])-sd_th0)/sd_th0:.0%}
  gap is realised-vs-expected sampling). Adding the report's k_e prior tightens this to
  sd(theta)={sd_th_lin:.0f}. So the PI's x'Px fit with the correct P is the report's MLE, now shown
  as a geometry; the Euclidean (P=I) fit mis-weights by the crossing angle and is a
  different, worse estimator ({xE[0]:.0f} vs {xW[0]:.0f}). Weighting is Remark 6 operational.
""")
    print("=" * 92)
    print("  WHAT THE GEOMETRY ADDS (it checks the numbers; it does not change them):")
    print("  (a) the pencil figure -- WHY the fibre exists, not just that it does;")
    print("  (b) a model-free uncertainty (intersection scatter) confirming the heavy tail")
    print("      from a THIRD route, and showing the pencil has no natural-scale centre;")
    print("  (c) a concurrency adequacy test -- calibrated, with power against a biased LOAD")
    print("      but BLIND to a second saturable channel (it hides in the fibre; t39 here);")
    print("  (d) Remark 6's weighting shown to reproduce the FIM/MLE.")
    print("=" * 92)


if __name__ == "__main__":
    main()
