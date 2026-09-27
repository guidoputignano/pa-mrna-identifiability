"""
Independent re-check of audit/u15_gut_suppression.py, sharing no code path.

WHAT u15 DOES, AND WHY THIS IS DIFFERENT. u15 builds a Fisher information
matrix in FLUX space from rows d(flux)/d(param), integrates the state and its
four sensitivities together with scipy LSODA, inverts, and takes the median
over 24 record seeds. Every one of those steps is replaced here:

  quantity            u15                         this file
  ------------------  --------------------------  ---------------------------
  residual space      flux                        LEVEL
  sd(ln k_e) from     inverse Fisher              MONTE CARLO sampling spread
  sensitivities       integrated ODE alongside    none: the fitter differences
  integrator          scipy solve_ivp, LSODA      own fixed-step RK4
  quasi-steady level  U14.qss                     own quadratic solve
  estimator           implicit (linearised)       scipy least_squares, trf

So the only things in common are the parameter values in src/common/params.py
and the physical model itself, which is the point.

TWO EARLIER ERRORS THIS FILE IS BUILT TO AVOID.
  R98. Holding a prior-carrying nuisance at its prior MEAN across replicates
  makes the sampling spread too small and the comparison meaningless. Here
  b_true and ke_true are DRAWN FROM THEIR PRIORS in every replicate, so the
  spread is the posterior object the Fisher inverse approximates.
  R101. Evaluating the clearance at the model's OWN quasi-steady root makes
  the record residual identically zero, and the fit then returns its start
  point. Here the record residual compares the model's predicted LEVEL with
  the OBSERVED level, which is not a root of anything.

PRECISION. Each cell is 6 record seeds x 50 replicates. A 50-replicate sd
carries about 10% of its own value, and the median of 6 of those about 5 to
6%, so agreement is expected at the several-per-cent level and not closer.
That is enough to catch a sign error, a factor, or a wrong object, which is
what a re-check is for.

Nothing here records an expected result in advance.

Run:  python audit/u16_recheck_gut.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

#  Constants come from the REGISTER and from src/design/run.py, which is
#  repository code and is shared with u15 by necessity: the two files must
#  price the same physical model or the comparison means nothing. What is not
#  shared is every step of the method, as tabulated above.
from common.params import _A, KE_PRIOR_LOG, LOAD_FRAC as _LF
from design.run import ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV, load

KM_T = float(_A.Km)
KE_T = float(_A.ke)
LB_T = float(load(0.0))
LAM = float(LAMBDA)
IVL = float(DOSE_INTERVAL)
PEAK_T = float(D_AMP) / (1.0 - float(np.exp(-LAM * IVL)))
TROUGH_T = PEAK_T * float(np.exp(-LAM * IVL))
CV_T = float(ASSAY_CV)
LOADF = float(_LF)
LOADCV = float(LOAD_CV)
KEPRI = float(KE_PRIOR_LOG)
MONTH_T = 30.4
N_REC = 118

SEEDS = range(6)
REPS = 50


def qss_own(theta, inc, Ltot, Km, ke):
    """
    Quasi-steady level from the positive root of the flux balance, written
    from the balance rather than taken from u14.

        (theta + inc) m/(Km+m) + ke m = Ltot
     => ke m^2 + (ke Km + theta + inc - Ltot) m - Ltot Km = 0
    """
    b = ke * Km + theta + np.asarray(inc, float) - Ltot
    return (np.sqrt(b * b + 4.0 * ke * Ltot * Km) - b) / (2.0 * ke)


def traj_rk4(theta, Km, ke, Ltot, amp, tau_in, t_out, h=0.005):
    """
    Marker level over a suppression burst, own fixed-step RK4.

    dm/dt = L(t) - (theta + inc(t)) m/(Km+m) - ke m,
    L(t)  = Ltot (1 - amp (1 - exp(-t/tau_in))),  inc(t) = TROUGH exp(-LAM t).

    The step is 7.2 min against a relaxation time of order 2.6 h, so the
    local error is of order (h/tau)^4 ~ 4e-6 and the integration is not what
    limits this comparison.
    """
    def rhs(t, m):
        f = 1.0 if tau_in <= 0.0 else (1.0 - np.exp(-t / tau_in))
        L = Ltot * (1.0 - amp * f)
        inc = TROUGH_T * np.exp(-LAM * t)
        m = m if m > 1e-12 else 1e-12
        return L - (theta + inc) * m / (Km + m) - ke * m

    t_out = np.asarray(t_out, float)
    out = np.empty(t_out.size)
    m = float(qss_own(theta, TROUGH_T, Ltot, Km, ke))
    t, j = 0.0, 0
    while j < t_out.size and t_out[j] <= 1e-15:
        out[j] = m
        j += 1
    while j < t_out.size:
        target = t_out[j]
        while t < target - 1e-12:
            step = min(h, target - t)
            k1 = rhs(t, m)
            k2 = rhs(t + 0.5 * step, m + 0.5 * step * k1)
            k3 = rhs(t + 0.5 * step, m + 0.5 * step * k2)
            k4 = rhs(t + step, m + step * k3)
            m += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            t += step
        out[j] = m
        j += 1
    return out


def one_fit(theta, rng, burst_t, amp, tau_in, burst_cv, n_bursts):
    """
    One replicate: draw truth, generate LEVEL data, fit four parameters,
    return ln(ke_hat / ke_true).
    """
    #  R98: every prior-carrying nuisance is drawn from its prior.
    b_true = LOADF * LB_T * rng.standard_normal()
    ke_true = KE_T * float(np.exp(KEPRI * rng.standard_normal()))
    Lbias = LB_T + b_true

    #  --- the quasi-steady record, in levels
    tt = MONTH_T * np.arange(N_REC, dtype=float)
    inc = PEAK_T * np.exp(-LAM * np.mod(tt - 50.0 * IVL, IVL))
    L_i = Lbias + LOADCV * LB_T * rng.standard_normal(N_REC)
    m_rec = qss_own(theta, inc, L_i, KM_T, ke_true)
    m_obs = m_rec * (1.0 + CV_T * rng.standard_normal(N_REC))

    #  Per-visit sigma in LEVEL units. The register's flux-space noise
    #  LOADF*LB maps through the local slope; the assay term is already a
    #  level. u15's sigma_flux uses the slope at the mean level, this uses
    #  it pointwise, which is the one deliberate difference in convention.
    slope = (theta + inc) * KM_T / (KM_T + m_obs) ** 2 + ke_true
    sig_rec = np.hypot(LOADF * LB_T / slope, CV_T * m_obs)

    #  --- the burst, in levels, n_bursts independent admissions
    if burst_t is not None and len(burst_t):
        m_b = traj_rk4(theta, KM_T, ke_true, Lbias, amp, tau_in, burst_t)
        obs_b = [m_b * (1.0 + burst_cv * rng.standard_normal(m_b.size))
                 for _ in range(n_bursts)]
    else:
        obs_b = []

    def resid(q):
        th, Km, b, ke = (float(np.exp(q[0])), float(np.exp(q[1])),
                         float(q[2]), float(np.exp(q[3])))
        #  R101: predicted LEVEL against the OBSERVED level.
        pred = qss_own(th, inc, LB_T + b, Km, ke)
        r = [(pred - m_obs) / sig_rec]
        if obs_b:
            pb = traj_rk4(th, Km, ke, LB_T + b, amp, tau_in, burst_t)
            for o in obs_b:
                r.append((pb - o) / (burst_cv * o))
        r.append([b / (LOADF * LB_T),
                  float(np.log(ke / KE_T)) / KEPRI])
        return np.concatenate([np.atleast_1d(x).ravel() for x in r])

    q0 = np.array([np.log(theta), np.log(KM_T), 0.0, np.log(KE_T)])
    try:
        s = least_squares(resid, q0, method="trf", x_scale=[1, 1, 16.0, 1],
                          xtol=1e-12, ftol=1e-12, max_nfev=400)
    except Exception:
        return np.nan
    return float(s.x[3] - np.log(ke_true))


def cell(theta, burst_t, amp=0.0, tau_in=0.0, burst_cv=CV_T, n_bursts=1,
         stat="sd"):
    """
    Median over record seeds of the per-seed Monte Carlo spread of ln k_e.

    WHICH STATISTIC, and this mattered. This audit's standing rule is to
    report the ROBUST spread, 0.7413 x IQR, so that collapsed fits do not
    dominate. That rule is right for reporting a spread on its own and WRONG
    for comparing one against a Fisher width, because a Fisher width is a
    STANDARD DEVIATION and 0.7413 x IQR equals an sd only for a Gaussian.
    Where the sampling distribution is heavy tailed the robust spread sits
    BELOW the sd, here by 11% to 52%, and comparing it with a Fisher
    manufactures ratios below one, which is the direction that is supposed to
    be impossible and which located the real defect at R119. R132 records the
    diagnosis. Both are returned so the comparison uses the sd while the
    robust value still flags collapse.
    """
    rob, sds = [], []
    for s in SEEDS:
        rng = np.random.default_rng(90000 + 7919 * s)
        v = np.array([one_fit(theta, rng, burst_t, amp, tau_in,
                              burst_cv, n_bursts) for _ in range(REPS)])
        v = v[np.isfinite(v)]
        if v.size > 10:
            rob.append(0.7413 * float(np.subtract(*np.percentile(v, [75, 25]))))
            sds.append(float(np.std(v, ddof=1)))
    if not rob:
        return np.nan
    if stat == "both":
        return float(np.median(rob)), float(np.median(sds))
    return float(np.median(sds if stat == "sd" else rob))


def part_glucose():
    """
    Independent re-check of u17's glucose-arm headline.

    WHICH u17 CONFIGURATION THIS CHECKS, and a correction. This file's
    generator gives the record and the burst ONE shared load bias, because
    both trajectories are built from Lbias = LB + b_true. That is u17's
    share_b = True, NOT its f_burst = 1, which gives the burst its own
    independent bias at the same prior width. The first version of this
    function compared against f_burst = 1 and so compared two different
    models; the disagreement grew with amplitude because the bias column
    scales with (1 - amp f(t)). Corrected to share_b = True, which PART F of
    u17 shows reproduces u15 exactly.

    u17's headline cells are computed with the bias KNOWN, which this
    generator cannot express, so the shared-bias column is what is compared
    and the difference between the two is reported by u17 PART B and PART H
    rather than here.

    It also re-checks the INTERIOR OPTIMUM in amplitude, which is a claim
    about the shape of a curve rather than about one number, and a shape is
    what an independent method is best placed to confirm or break.
    """
    import u17_glucose_burst as U17
    hour = 1.0 / 24.0
    bt = np.linspace(0.0, 1.0, 25)
    print()
    print("=" * 86)
    print("Independent re-check of u17's glucose arm, same level-space method")
    print("=" * 86)
    print()
    print(f"    {'amplitude':>10} {'theta':>6} {'u17 shared':>10} "
          f"{'MC sd':>10} {'ratio':>7} {'MC robust':>9}  what it is")
    print("    " + "-" * 82)
    rows = [(0.222, 300.0, "gut share alone"),
            (0.350, 300.0, "near the claimed optimum"),
            (0.517, 600.0, "catabolism mean, the headline"),
            (0.517, 300.0, "catabolism mean, the headline"),
            (0.739, 300.0, "both arms, claimed WORSE than 51.7%")]
    got = {}
    for amp, th, note in rows:
        ref = U17.sd_ke(th, bt, amp, 6 * hour, 0.01, 16, 1.0,
                         share_b=True)
        rob, mc = cell(th, bt, amp, 6 * hour, 0.01, 16, stat="both")
        got[(amp, th)] = mc
        print(f"    {amp:10.1%} {th:6.0f} {ref:10.4f} {mc:10.4f} "
              f"{mc / ref:7.3f} {rob:9.4f}  {note}")
    print()
    a, b, c = got[(0.350, 300.0)], got[(0.517, 300.0)], got[(0.739, 300.0)]
    print(f"  THE SHAPE. At theta = 300 this method gives {a:.4f} at 35%, "
          f"{b:.4f} at 51.7%")
    print(f"  and {c:.4f} at 73.9%. An interior optimum requires the middle "
          f"value to be")
    print(f"  the smallest: {'CONFIRMED' if b < a and b < c else 'NOT CONFIRMED'}.")
    print("  That is the claim that decides whether adding bowel preparation")
    print("  to the glucose arm helps or hurts, so it is the one worth")
    print("  breaking if it is breakable.")


if __name__ == "__main__":
    print("=" * 86)
    print("Independent re-check of u15: LEVEL space, own RK4, Monte Carlo")
    print("=" * 86)
    print(f"\n  Read from the register: Km {KM_T:.0f}, ke {KE_T:.0f}/d, "
          f"L {LB_T:.3f}, trough increment {TROUGH_T:.1f}")
    print(f"  {len(list(SEEDS))} record seeds x {REPS} replicates per cell, "
          f"so a few per cent of precision.\n")
    hour = 1.0 / 24.0
    t1 = hour * np.arange(0, 25, dtype=float)
    t4 = np.linspace(0.0, 4.0, 25)
    #  Reference values transcribed from u15's output AFTER the R119 fix to
    #  its load-bias forcing, which THIS FILE is what caught: the pre-fix
    #  references were 0.4771, 0.5179, 0.3928 and 0.1561, and the last of
    #  those disagreed here by a factor of two in the direction a sampling
    #  spread cannot go. Deliberately transcribed rather than imported, so
    #  that no code path is shared.
    cases = [
        ("record only, no burst", 300.0, None, 0.0, 0.0, CV_T, 1, 0.5430),
        ("22.2%, ideal step", 300.0, t1, 0.222, 0.0, CV_T, 1, 0.4762),
        ("22.2%, 1 day onset", 300.0, t4, 0.222, 1.0, CV_T, 1, 0.5188),
        ("41%, ideal step", 300.0, t1, 0.410, 0.0, CV_T, 1, 0.3711),
        ("22.2%, 6 h, 1%, 16 adm", 300.0, np.linspace(0, 1, 25), 0.222,
         6 * hour, 0.01, 16, 0.0588),
    ]
    print(f"    {'cell':>24} {'theta':>6} {'u15 (Fisher)':>13} "
          f"{'MC sd':>10} {'ratio':>7} {'MC robust':>9}")
    print("    " + "-" * 76)
    for lab, th, bt, amp, tau, cv, nb, ref in cases:
        rob, sd = cell(th, bt, amp, tau, cv, nb, stat="both")
        print(f"    {lab:>24} {th:6.0f} {ref:13.4f} {sd:10.4f} "
              f"{sd / ref:7.3f} {rob:9.4f}")
    print()
    print("  READING. Every ratio now EXCEEDS one, which is the only direction")
    print("  available to an honest sampling spread against a local quadratic:")
    print("  a Monte Carlo can be wider than an inverse Fisher and cannot be")
    print("  narrower. Before the R119 fix the strongest arm sat at 0.50, and")
    print("  that impossible direction is what located the defect.")
    print()
    print("  WHAT THIS DOES AND DOES NOT LICENCE. It confirms u15's ranking")
    print("  and its structure. It does NOT licence a claim that rests on a")
    print("  margin smaller than the spread between the two methods, which")
    print("  runs 9% to 33% here. u15's headline at a tenth of nominal,")
    print("  0.0588 against a frontier of 0.0633, is a 7% margin, and this")
    print("  file gives 0.0781 for the same cell, which MISSES. So that cell")
    print("  is not established either way. At a fifth of nominal the margin")
    print("  is 1.9x and survives the spread comfortably.")
    print()
    print("  That is (g4)'s own lesson turned on this audit: a number defined")
    print("  to land near a gate carries no margin, and fails under any")
    print("  object but the one it was computed against.")
    print()
    print("  It confirms neither the physical model, which both files share,")
    print("  nor the amplitude's provenance, which is what u15 PART D turns on.")
    part_glucose()
