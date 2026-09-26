"""
S1. The increment time course: analyst against truth.

Phase 8 item 1. The robustness review of 25 September 2026 left a gap it named
B: if the true therapeutic increment follows the published two-compartment
protein kinetics, or a 5 to 7 day half-life, while the analysis assumes the
committed 3.5 days, `theta_hat/theta` is 1.3 to 2.0 at a structural residual of
0.3 to 0.7 sigma. That is a bias with no warning in the fit. This script
measures the whole grid of truth-against-analyst pairs with the repository's
own estimators, adds the trial's real irregular dose train, and ends with the
design question the gap actually poses: what would have to be measured, and
when, to pin the shape.

WHAT IS BEING VARIED, and what is not. Only the SHAPE of `Delta(t)` and the
DOSE TRAIN that shape is convolved with. The marker model, the priors, the
noise model, the schedule and every constant in `src/common/params.py` are the
committed ones. Nothing in `src/` is edited. `design.run.LAMBDA` and
`design.run.D_AMP` are never read after Part 0, because every increment below
is assembled here from an explicit kernel and an explicit train rather than by
calling `design.run.delta`, and the committed shape is reproduced as one member
of that family and checked against `design.run.delta` to machine precision.

    THE SHAPES (truths). Each is a unit-impulse response `phi(u)`, u >= 0.

      exp(t_half)   single exponential, t_half in {2.5, 3.5, 5.0, 7.0} d.
                    3.5 d is the repository's committed LAMBDA and so is the
                    well-specified case. The others are the range brief R2
                    (the increment time course) leaves open: R2 finds NO human
                    hepatic PCC protein or activity time course at all, so
                    every shape constant below is a mouse value or an
                    assumption, and the half-life is carried as a sweep.
      attarwala     the two-compartment protein model of Attarwala et al. 2023
                    (doi:10.1089/nat.2022.0036), Suppl Table S1(B), driven
                    through its effect compartment:
                        d(Ce)/dt  = -ke0 Ce
                        d(PCC)/dt = Slope*Ce + kq (PCP - PCC) - kdeg PCC
                        d(PCP)/dt = kq (PCC - PCP)
                    with ke0, Slope, kq, kdeg from its Table 1. MOUSE
                    parameters throughout, as `audit/k1_sources.py` finding 2
                    established and as brief R2 restates.

    THE TRAINS. `nominal` is the committed uniform 14 d train. The irregular
    trains are built from what the Nature 2024 interim analysis publishes and
    `audit/PHASE6_REVIEWS_SEP2026.md` section 8 records: cohort 1 is Q3W and
    cohorts 2 to 5 Q2W, amplitudes 0.30 to 0.90 mg/kg, one participant had dose
    2 delayed, dose 3 skipped and a de-escalation from 0.90 to 0.60, and the one
    slippage the paper quantifies is a dose 3 at day 52 against a 42 d nominal.
    Every construction choice is printed with its source or marked ASSUMED.

    THE ANALYSTS.
      fixed 3.5 d       the committed shape. 4 free parameters.
      half-life free    a single exponential whose half-life is estimated.  5.
      two-compartment   the Attarwala family with kq and kdeg estimated.    6.
      amplitude free    the committed shape with a free dose scale kappa and
                        a 12% log-prior on it, which the review reports does
                        not help. Carried so that report is checked and not
                        inherited.

    THE AMPLITUDE QUESTION, which has no single right answer and is therefore
    carried as a footing sweep rather than settled. Comparing two shapes
    requires saying what is held equal. Three conventions are computed:
      peak    the steady-state peak increment of the nominal train matches the
              committed train's. Holds the maximum enzyme level equal, so the
              capacity ceiling of Section 7 is comparable across shapes.
      dose    the SINGLE-dose peak matches, so each administration deposits the
              same increment and a longer half-life accumulates more.
      mean    the steady-state cycle mean matches, i.e. equal exposure.
    `peak` is the headline because it is the one under which the admissibility
    check does not move. Every bias figure is reported under all three.

REPORTED, per (truth, analyst) cell: median `theta_hat/theta` with its seed
range, RMSE of `ln theta_hat`, the RMS standardised residual against the
`sqrt((n-p)/n)` a well-specified fit would give, and the coverage of the 95%
profile interval over seeds.

Nothing here records an expected result in advance. The independent re-check is
`audit/q9_recheck.py`, which shares no code path with this file.

Run:  python audit/q1_increment_shape.py
"""
import sys
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize

import design.run as DR
import i1_estimator as I1
import k1_sources as K1
import model_A.run as LOC
from common.models import ParamsA
from design.run import ASSAY_CV, LOAD_FRAC

# --------------------------------------------------------------------------
#  The committed footing, taken from the modules that own it. Nothing here
#  assigns a value that `src/` already carries.
# --------------------------------------------------------------------------
TH, KM, KE = I1.TH, I1.KM, I1.KE          # 3000, 100, 5
LB = I1.LB                                 # the recorded load at the operating point
SIG = I1.SIG                               # sigma_total(0.05, 0.05)
S_B, S_KE = I1.S_B, I1.S_KE                # the paper's two priors
M_BAR, BLOC, CV = I1.M_BAR, I1.BLOC, I1.CV
MONTHLY = I1.MONTHLY                       # np.arange(0, 730, 30.4), n = 25
T_DOSE = DR.DOSE_INTERVAL                  # 14 d
D_AMP = DR.D_AMP                           # 6000
LAMBDA0 = DR.LAMBDA                        # ln 2 / 3.5 d
TH_HALF0 = float(np.log(2.0) / LAMBDA0)    # 3.5 d

#: the log-prior width on a free dose scale. The review used 12%; it is an
#: ASSUMPTION about how well a delivered increment is known, not a measurement,
#: and `audit/g_design_common.py:sd_ln_theta` already carries kappa the same way.
KAPPA_PRIOR = 0.12

#: the seed set every design number in this repository is reported over.
SEEDS = tuple(range(24))

#: the K_m grid the variable-projection scan runs on, i1_estimator's own.
KM_GRID = I1.GRID
#: a coarser grid for the profile, whose only job is to find the basin before
#: the exact polish. Part 1 checks the profile against the fine grid.
KM_GRID_PROF = np.geomspace(1.0, 1000.0, 241)

#: chi2(1) at 95%, halved. n6_profile_vs_crb.py's DROP, and it is a cut on the
#: NEGATIVE LOG LIKELIHOOD. Everything in this file carries the objective in
#: `i1_estimator.obj_exact`'s scale, which is TWICE that, so the cut applied
#: here is 2 * 1.92. An earlier draft of this file cut at 1.92 on the doubled
#: scale and reported intervals about a third too narrow; `audit/q9_recheck.py`
#: check 5 found it against n6's own objective.
DROP_NLL = 1.92
DROP = 2.0 * DROP_NLL
#: how far the profile thins the outer shape grid. 1 is no thinning.
PROFILE_STRIDE = 5

#: WHICH SIDE OF AN INFUSION A DRAW AT THE SAME INSTANT IS READ ON, for the
#: truth and the analyst alike. The paper's monthly record puts one draw at a
#: dose instant, t = 0, and `design.run.delta` reads it on the POST-dose side
#: (`t >= t_k`). That draw is inside the marker's own settling window, 0.054 d
#: by `g_design_common.settling_days`, so the quasi-steady model does not hold
#: there; and what an increment kernel says at that instant is exactly where
#: kernels differ most, a jump against a rise. The first run of this script
#: read it post-dose and it alone carried the Attarwala cell from about 1.1 to
#: 5.6 (`q1b_reconcile.py`, R2 and R4). "pre-dose" reads it before the
#: infusion, as a clinic draws on a dose day. "post-dose" reproduces the
#: first run.
DOSE_DAY_DRAW = "pre-dose"
_PRE = 1e-9


def at_draw(times):
    """The instants the increment is evaluated at, per DOSE_DAY_DRAW."""
    times = np.asarray(times, float)
    return times - _PRE if DOSE_DAY_DRAW == "pre-dose" else times

#: burn-in before the increment train is read as steady. The slowest kernel
#: here has a terminal half-life near 12 d, so 200 d is 12 interval-lengths;
#: Part 0 prints the convergence actually achieved rather than assuming it.
BURN_IN = 200.0
FINE = np.arange(0.0, 900.0, 0.02)


def banner(t):
    print("\n" + "=" * 96)
    print("  " + t)
    print("=" * 96)


def spread(v, fmt="{:.4f}"):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if not v.size:
        return "all non-finite"
    return (fmt.format(np.median(v)) + " [" + fmt.format(v.min()) + ", "
            + fmt.format(v.max()) + "]")


@contextmanager
def rebind(mod, **kw):
    """Temporarily rebind module constants and restore them. k1_sources' `_set`."""
    old = {k: getattr(mod, k) for k in kw}
    for k, v in kw.items():
        setattr(mod, k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            setattr(mod, k, v)


# ==========================================================================
#  Part 0a. The kernels.
# ==========================================================================
class ExpSum:
    """
    A unit-impulse response written as `phi(u) = sum_i c_i exp(w_i u)`, u in days.

    Both shape families here are exactly of this form: a single exponential
    trivially, and the two-compartment protein model because a linear
    time-invariant system's impulse response is a sum of exponentials in its
    eigenvalues. Writing them this way makes every normalisation statistic
    CLOSED FORM rather than a quadrature, which is what lets the outer shape
    search be a grid scan rather than a local descent. `_attarwala_profile`
    below is the independent slow route the constructor is checked against.

    The cycle mean deserves a line, because it is the one that looks like an
    integral and is not: the steady-state increment averaged over one dose
    interval is the AREA of one impulse response divided by the interval, for
    any kernel, since every dose contributes its whole area exactly once per
    interval in steady state.
    """

    def __init__(self, c, w):
        self.c = np.asarray(c, float)
        self.w = np.asarray(w, float)
        if np.any(self.w > 0):
            raise ValueError("a growing increment is not a kernel")

    def __call__(self, u):
        u = np.asarray(u, float)
        e = np.exp(np.clip(self.w[None, :] * u[..., None], -700.0, 0.0))
        return (e * self.c[None, :]).sum(axis=-1)

    def area(self):
        return float(np.sum(self.c / -self.w))

    def steady(self, tau, period):
        """sum_j phi(tau + j*period), the steady increment at cycle phase tau."""
        tau = np.asarray(tau, float)
        d = self.c / (1.0 - np.exp(self.w * period))
        e = np.exp(np.clip(self.w[None, :] * tau[..., None], -700.0, 0.0))
        return (e * d[None, :]).sum(axis=-1)

    @staticmethod
    def _argmax(f, lo, hi, n=257, rounds=4):
        """A bracketed maximum of a smooth 1-D function, refined by zooming."""
        for _ in range(rounds):
            x = np.linspace(lo, hi, n)
            y = f(x)
            j = int(np.argmax(y))
            lo = x[max(j - 1, 0)]
            hi = x[min(j + 1, n - 1)]
        return float(0.5 * (lo + hi)), float(np.max(f(np.array([lo, hi]))))

    def single_peak(self, horizon=400.0):
        return self._argmax(self.__call__, 0.0, horizon)[1]

    def steady_peak(self, period):
        return self._argmax(lambda t: self.steady(t, period), 0.0, period)[1]

    def steady_mean(self, period):
        return self.area() / period


def kernel_exp(t_half):
    """Single-exponential unit-impulse response, in days."""
    return ExpSum([1.0], [-np.log(2.0) / float(t_half)])


def kernel_twocomp(kq, kdeg, ke0=None, slope=None):
    """
    The two-compartment protein model of Attarwala et al. 2023 at FREE rates,
    which is analyst option 3, as an exponential sum.

        d(Ce)/dt  = -ke0 Ce
        d(PCC)/dt = Slope*Ce + kq (PCP - PCC) - kdeg PCC
        d(PCP)/dt = kq (PCC - PCP)

    Rates are per HOUR in the source and the kernel is in DAYS, so the
    eigenvalues are multiplied by 24. `slope` scales the whole response and is
    normalised away, so it is held at its published value and never freed.
    """
    ke0 = K1.KE0 if ke0 is None else ke0
    slope = K1.SLOPE if slope is None else slope
    A = np.array([[-ke0, 0.0, 0.0],
                  [slope, -(kq + kdeg), kq],
                  [0.0, kq, -kq]], float)
    w, V = np.linalg.eig(A)
    if np.linalg.cond(V) > 1e10:
        #  a repeated eigenvalue makes the decomposition defective. Nudging
        #  ke0 by one part in 1e-7 moves the kernel by less than any figure
        #  below is quoted to and restores a basis. Reported, not silent.
        _DEFECTIVE.append((kq, kdeg, ke0))
        return kernel_twocomp(kq, kdeg, ke0 * (1.0 + 1e-7), slope)
    c = np.linalg.solve(V, np.array([1.0, 0.0, 0.0]))
    coef = np.real(V[1, :] * c)
    rate = np.real(w) * 24.0
    #  the zero mode of the kq -> 0 limit carries no PCC weight; drop the
    #  numerically zero rates rather than dividing by them in `area`.
    keep = rate < -1e-12
    if not keep.all():
        if np.any(np.abs(coef[~keep]) > 1e-9 * np.abs(coef).max()):
            _NONDECAY.append((kq, kdeg, float(np.abs(coef[~keep]).max())))
        coef, rate = coef[keep], rate[keep]
    return ExpSum(coef, rate)


#: kernels whose eigen-decomposition needed a nudge, and kernels with a
#: non-decaying mode of non-negligible weight. Printed by Part 0; a non-empty
#: second list would mean the family reached outside what a kernel can be.
_DEFECTIVE, _NONDECAY = [], []


def _attarwala_profile(t_end_days=900.0, n=180001):
    """
    PCC(t) after a unit bolus into the effect compartment, on a dense grid.

    The three states are Attarwala's Suppl Table S1(B) as `audit/k1_sources.py`
    reads them: the printed `Ksyn*Ce` with `Ksyn = Ce*Slope` is quadratic in Ce
    and dimensionally wrong, the text says synthesis rises linearly in mRNA, so
    the synthesis term is `Slope*Ce`. The eigenvalues do not depend on that
    reading; the amplitude does, and the amplitude is normalised away below.
    """
    hours = np.linspace(0.0, t_end_days * 24.0, n)

    def rhs(t, y):
        ce, pcc, pcp = y
        return [-K1.KE0 * ce,
                K1.SLOPE * ce + K1.KQ * (pcp - pcc) - K1.KDEG * pcc,
                K1.KQ * (pcc - pcp)]

    sol = solve_ivp(rhs, (0.0, hours[-1]), [1.0, 0.0, 0.0], t_eval=hours,
                    rtol=1e-11, atol=1e-14, method="LSODA")
    return hours / 24.0, np.asarray(sol.y[1], float)


_ATT_U, _ATT_Y = _attarwala_profile()


def kernel_attarwala():
    """The published two-compartment protein response, at Table 1's rates."""
    return kernel_twocomp(K1.KQ, K1.KDEG)


KERNELS = {
    "exp 2.5 d": kernel_exp(2.5),
    "exp 3.5 d": kernel_exp(3.5),
    "exp 5.0 d": kernel_exp(5.0),
    "exp 7.0 d": kernel_exp(7.0),
    "Attarwala 2-comp": kernel_attarwala(),
}


# ==========================================================================
#  Part 0b. The dose trains.
# ==========================================================================
def nominal_train(t_end=900.0, period=T_DOSE):
    """The committed uniform train: unit doses every 14 days from t = 0."""
    return np.arange(0.0, t_end, period), np.ones(int(np.ceil(t_end / period)))


#: The published per-regimen dose counts, Nature 2024 interim analysis, as
#: `audit/PHASE6_REVIEWS_SEP2026.md` section 8 records them. Used only to fix
#: the mg/kg that the committed D_AMP stands in for.
PUBLISHED_DOSES = (
    (0.30, "Q3W", 73), (0.30, "Q2W", 136), (0.45, "Q2W", 68),
    (0.60, "Q2W", 45), (0.90, "Q2W", 24),
)


def published_mean_mgkg():
    n = sum(c for _, _, c in PUBLISHED_DOSES)
    return sum(d * c for d, _, c in PUBLISHED_DOSES) / n, n


def irregular_train(kind, t_end=900.0, delay=7.0):
    """
    A single patient's dose train, in (times, mg/kg).

    PUBLISHED, from the Nature 2024 interim analysis:
      - cohort 1 is 0.30 mg/kg Q3W (n = 4); cohorts 2 to 5 are 0.30, 0.45,
        0.60 and 0.90 mg/kg Q2W (n = 3 each).
      - one participant had a dose decrease from 0.90 to 0.60 after a grade 3
        pancreatitis, SKIPPED DOSE 3 entirely, and had dose 2 delayed by
        hospitalisation.
      - two cohort-1 participants moved to Q2W.
      - for a cohort-1 participant the predose sample for dose 3 was 52 days
        after the first dose, against a 42 d nominal: a 10 d slip.

    ASSUMED, because the paper does not say:
      - the SIZE of the dose-2 delay. `delay` days, default 7, which is the
        published 10 d slip on a 21 d interval scaled to a 14 d one. Part 2
        sweeps it.
      - WHEN the cohort-1 switch to Q2W happened. Taken at dose 6.
      - that the record runs to `t_end` at the post-change regimen.
    """
    if kind == "C5 de-escalated":
        #  0.90 Q2W, dose 2 delayed, dose 3 skipped, 0.60 Q2W from dose 4
        t = [0.0, T_DOSE + delay]
        a = [0.90, 0.90]
        #  dose 3 would have been at 2T; it is skipped. Dosing resumes on the
        #  nominal grid at dose 4 and continues fortnightly at 0.60.
        k = 3
        while k * T_DOSE < t_end:
            t.append(k * T_DOSE)
            a.append(0.60)
            k += 1
        return np.array(t), np.array(a)
    if kind == "C1 Q3W then Q2W":
        #  0.30 Q3W with the published 10 d slip on dose 3, then Q2W from
        #  dose 6
        t = [0.0, 21.0, 52.0, 73.0, 94.0]
        a = [0.30] * 5
        u = 94.0 + T_DOSE
        while u < t_end:
            t.append(u)
            a.append(0.30)
            u += T_DOSE
        return np.array(t), np.array(a)
    if kind == "composite":
        #  every published feature in one record, which is a CONSTRUCTION and
        #  not a transcription of any one participant: Q3W start with the
        #  published slip, a switch to Q2W, an amplitude ladder across the five
        #  published levels, a delayed dose and a skipped one.
        t = [0.0, 21.0, 52.0]
        a = [0.30, 0.30, 0.30]
        levels = [0.45, 0.45, 0.60, 0.60, 0.90, 0.90, 0.60]
        u, i = 52.0 + T_DOSE, 0
        skipped, delayed = False, False
        while u < t_end:
            if not delayed and i == 2:
                u += delay
                delayed = True
            if not skipped and i == 5:
                skipped = True
                u += T_DOSE
                i += 1
                continue
            t.append(u)
            a.append(levels[i] if i < len(levels) else levels[-1])
            u += T_DOSE
            i += 1
        return np.array(t), np.array(a)
    raise ValueError(kind)


def assemble(kernel, times, weights, t):
    """Delta(t) = sum_k w_k phi(t - t_k) over the doses already given."""
    t = np.atleast_1d(np.asarray(t, float))
    out = np.zeros(t.shape, float)
    for tk, wk in zip(np.asarray(times, float), np.asarray(weights, float)):
        on = t >= tk
        if on.any():
            out[on] += wk * kernel(t[on] - tk)
    return out


def norm_stats(kernel, times, weights, burn=BURN_IN, grid=FINE):
    """
    Single-dose peak, steady peak and steady mean of one kernel-train pair, by
    direct assembly on a fine grid. This is the SLOW and obvious route; it is
    the reference `norm_stats_uniform` is checked against in Part 0.
    """
    single = float(np.max(kernel(grid[grid <= 400.0])))
    y = assemble(kernel, times, weights, grid)
    #  exactly whole dose cycles after the burn-in, so the mean is a cycle mean
    #  and not a window mean that happens to end mid-cycle.
    ncyc = int((grid[-1] - burn) // T_DOSE)
    win = (grid >= burn) & (grid <= burn + ncyc * T_DOSE)
    ys, xs = y[win], grid[win]
    mean = float(np.sum(0.5 * (ys[1:] + ys[:-1]) * np.diff(xs))
                 / (ncyc * T_DOSE))
    return single, float(np.max(ys)), mean


def norm_stats_uniform(kernel, period=T_DOSE):
    """
    The same three statistics for the NOMINAL uniform unit train, in closed
    form rather than by assembly. Part 0 checks it against `norm_stats`.
    """
    return (kernel.single_peak(), kernel.steady_peak(period),
            kernel.steady_mean(period))


_AMP_CACHE = {}


def amplitude(kernel, convention, ref, key=None):
    """
    The scale A that puts this kernel on the committed footing.

    `ref` is the (single, peak, mean) triple of the committed 3.5 d kernel on
    the nominal train, and the scale is always computed on the NOMINAL train so
    that it is a property of the kernel and not of the train. An irregular
    train then carries its own mg/kg weights on top.

    `key` caches the result for a kernel the outer search revisits, which the
    shape grid does at every profile node. Cached on the caller's own key, so a
    caller that supplies none pays full price and nothing is ever returned for
    a kernel that was not the one asked about.
    """
    if key is not None and (key, convention) in _AMP_CACHE:
        return _AMP_CACHE[(key, convention)]
    s, p, m = norm_stats_uniform(kernel)
    got = {"dose": s, "peak": p, "mean": m}[convention]
    want = {"dose": ref[0], "peak": ref[1], "mean": ref[2]}[convention]
    A = want / got
    if key is not None:
        _AMP_CACHE[(key, convention)] = A
    return A


# ==========================================================================
#  Part 0c. The record, and the estimator.
# ==========================================================================

# ==========================================================================
#  The footing: every constant a fit needs at one capacity, in one object.
# ==========================================================================
class Footing:
    """
    The operating point, the noise and the two priors at one capacity.

    This exists because Part 4 asks the same question at a disease-range
    capacity, where the operating level, sigma, the local slope and the
    load-bias prior all move together. Rebinding module constants would not
    reach a default argument, which binds at definition; passing a footing
    does. At theta = 3000 it must reproduce i1_estimator's committed constants
    exactly, and Part 0 checks that rather than asserting it.

    `convention` is n4_theta_scope's: "A" holds the recorded load at its
    theta = 3000 value and lets the marker rise, which is the physical reading
    and the paper's own.
    """

    def __init__(self, theta=TH, Km=KM, ke=KE, convention="A"):
        import n4_theta_scope as N4
        L, m_bar, sigma, s_b = N4.footing(theta, Km, ke, convention)
        p = ParamsA(theta=theta, Km=Km, ke=ke)
        self.theta, self.Km, self.ke = float(theta), float(Km), float(ke)
        self.L, self.m_bar = float(L), float(m_bar)
        self.sigma, self.s_b = float(sigma), float(s_b)
        self.s_ke = 0.8 * float(ke)
        self.bloc = float(LOC.fit_ab(p, m_bar)[1])
        self.convention = convention

    def __repr__(self):
        return (f"Footing(theta={self.theta:.0f}, L={self.L:.4f}, "
                f"m_bar={self.m_bar:.4f}, sigma={self.sigma:.4f}, "
                f"s_b={self.s_b:.4f}, b_loc={self.bloc:.4f})")


#: the committed footing, which every part but Part 4 runs on.
FOOT = Footing()


def level(theta, inc, L, Km=KM, ke=KE):
    """Positive root of (theta + inc) m/(Km + m) + ke m = L. The paper's balance."""
    B = ke * Km + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * L * Km)) / (2.0 * ke)


def make_record(inc_true, seed, F=None, load_cv=0.10, times=MONTHLY):
    """
    One synthetic record, on `i1_estimator.one_record`'s own noise convention.

    The load record is fixed within a seed and the assay error is placed on the
    ORDINATE as `b_loc cv m_bar`, which is `design/run.py:598-604`. The levels
    are exact given the TRUE increment; the analyst's increment never enters
    here, which is the whole point.
    """
    F = FOOT if F is None else F
    rng_load = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    L = F.L * np.exp(rng_load.normal(-0.5 * s2, np.sqrt(s2), times.size))
    rng = np.random.default_rng(1_000_003 * seed + 7)
    z_load = rng.standard_normal(times.size)
    z_assay = rng.standard_normal(times.size)
    m = level(F.theta, inc_true, L, Km=F.Km, ke=F.ke)
    s = L * (1.0 + LOAD_FRAC * z_load) + F.bloc * CV * F.m_bar * z_assay
    return m, s, L


def _linear_solve(cols, y, sigma, prior_sd, prior_mu):
    """
    A penalised linear least squares, batched over the leading axis.

    cols  (G, n, p)   design columns
    y     (G, n)      response
    Returns the solution (G, p) and the penalised objective (G,).
    """
    A = cols.transpose(0, 2, 1) @ cols / sigma ** 2
    rhs = np.einsum("gnj,gn->gj", cols, y) / sigma ** 2
    D = np.diag([0.0 if s is None else 1.0 / s ** 2 for s in prior_sd])
    mu = np.array([0.0 if s is None else m for s, m in zip(prior_sd, prior_mu)])
    A = A + D
    rhs = rhs + D @ mu
    p = np.linalg.solve(A, rhs[..., None])[..., 0]
    r = np.einsum("gnj,gj->gn", cols, p) - y
    J = np.sum(r ** 2, axis=1) / sigma ** 2
    for j, s in enumerate(prior_sd):
        if s is not None:
            J = J + (p[:, j] - mu[j]) ** 2 / s ** 2
    return p, J


def varpro_km(m, s, phi, free_kappa, grid=KM_GRID, F=None,
              theta_fixed=None, kappa_prior=None, refine=True):
    """
    Scan K_m, solve the rest exactly. One nonlinear parameter, so this is global.

    `phi` is the ANALYST'S assembled increment on the visit times, at unit
    scale. With `free_kappa` the dose scale joins the linear block; without it
    the increment is subtracted off at scale 1, which is `i1_estimator._normal`.
    `theta_fixed` holds ln theta for the profile and drops it from the block.
    """
    F = FOOT if F is None else F
    sigma, sb, ske, ke0 = F.sigma, F.s_b, F.s_ke, F.ke
    Km = np.atleast_1d(np.asarray(grid, float))
    g = m[None, :] / (Km[:, None] + m[None, :])
    one = np.ones_like(g)
    mm = np.broadcast_to(m, g.shape)
    cols, sd, mu = [], [], []
    if theta_fixed is None:
        cols.append(g)
        sd.append(None)
        mu.append(0.0)
    if free_kappa:
        cols.append(phi[None, :] * g)
        sd.append(kappa_prior)
        mu.append(1.0)
    cols += [mm, one]
    sd += [ske, sb]
    mu += [ke0, 0.0]
    X = np.stack(cols, axis=-1)
    y = np.broadcast_to(s, g.shape).copy()
    if theta_fixed is not None:
        y = y - theta_fixed * g
    if not free_kappa:
        y = y - phi[None, :] * g
    p, J = _linear_solve(X, y, sigma, sd, mu)
    j = int(np.argmin(J))
    if refine and 0 < j < Km.size - 1:
        lo, hi = np.log(Km[j - 1]), np.log(Km[j + 1])
        fine = np.exp(np.linspace(lo, hi, 61))
        p2, J2 = varpro_km(m, s, phi, free_kappa, fine, F,
                           theta_fixed, kappa_prior, refine=False)[:2]
        if J2.min() < J[j]:
            return p2, J2, float(fine[int(np.argmin(J2))])
    return p, J, float(Km[j])


#: a fit is called RUNAWAY when it returns a capacity this many times the
#: truth. It is not a physical bound on theta and nothing is clipped to it: it
#: is a label for the cells where the fit has left the range a capacity can be
#: read in, so that a median is not quoted as though it were an estimate.
RUNAWAY = 100.0
#: and DEGENERATE when it drives capacity to the floor instead.
FLOOR = 1e-4


def unpack(p_row, free_kappa, theta_fixed):
    """(theta, kappa, ke, b) from one solution row."""
    i = 0
    if theta_fixed is None:
        theta = float(p_row[i])
        i += 1
    else:
        theta = float(theta_fixed)
    kappa = 1.0
    if free_kappa:
        kappa = float(p_row[i])
        i += 1
    return theta, kappa, float(p_row[i]), float(p_row[i + 1])


def exact_obj(theta, Km, ke, b, kappa, m, s, phi, F=None,
              kappa_prior=None):
    """
    The PAPER'S penalised objective, priors in the coordinates the filter states
    them in: Gaussian on ln k_e, sd 0.8. `i1_estimator.obj_exact`.
    """
    F = FOOT if F is None else F
    sigma, sb, ke_ref = F.sigma, F.s_b, F.ke
    if theta <= 0 or Km <= 0 or ke <= 0:
        return np.inf
    pred = (theta + kappa * phi) * m / (Km + m) + ke * m + b
    J = float(np.sum(((pred - s) / sigma) ** 2))
    J += b ** 2 / sb ** 2 + (np.log(ke) - np.log(ke_ref)) ** 2 / 0.8 ** 2
    if kappa_prior is not None:
        #  on LN kappa, which is the repository's convention: the dose-scale
        #  column in g_design_common.sd_ln_theta is `inc g`, the derivative in
        #  ln kappa, and its prior is stated on that coordinate. An earlier
        #  draft carried (kappa - 1)^2 here, the delta-method form, which
        #  charges 68 for kappa = 0.01 where the log prior charges 1470, and
        #  let the fit switch the increment off. q1b_reconcile.py R2 found it
        #  against the probe's fitter on identical records.
        if kappa <= 0:
            return np.inf
        J += np.log(kappa) ** 2 / kappa_prior ** 2
    return J


def polish(theta, Km, ke, b, kappa, m, s, phi, free_kappa, theta_fixed=None,
           F=None, kappa_prior=None):
    """
    Re-minimise the EXACT objective from the grid solution.

    The linear inner problem must carry the delta-method form of the k_e
    log-prior; this step starts at its answer and minimises the exact one, so
    the grid supplies the global basin and the polish the exact stationary
    point. `i1_estimator.polish`.

    K_m IS BOUNDED TO THE SCAN'S OWN BOX, which `i1_estimator.polish` does not
    do because it is never run far enough from a well-specified fit to need it.
    Here it is: on the badly misspecified cells the polish walks K_m out past
    the top of the grid, where `theta g(m)` degenerates into `(theta/K_m) m`
    and is no longer distinguishable from the `k_e m` column at all, and
    returns capacities of order 1e7. That is not an estimate, it is the
    estimator leaving its own domain. Bounding it makes the scan and the polish
    one estimator, and `cell` counts how often the bound binds rather than
    hiding it.
    """
    F = FOOT if F is None else F
    sigma, sb, ke_ref = F.sigma, F.s_b, F.ke
    lo_km, hi_km = np.log(KM_GRID[0]), np.log(KM_GRID[-1])
    free = []
    if theta_fixed is None:
        free.append(("theta", np.log(max(theta, 1e-12))))
    free.append(("Km", np.log(max(Km, 1e-12))))
    #  The linear inner solve carries the k_e prior in its delta-method form,
    #  which does not keep k_e positive, and with theta held high on a profile
    #  it can return k_e <= 0. Starting the polish from ln(1e-12) there leaves
    #  it stranded where exp(ln k_e) has no gradient, and it was accepted
    #  because the start's own objective is infinite: the profile jumped from
    #  1.5 to 553 between two adjacent theta, and upper interval edges came
    #  out short. q9_recheck.py check 5 found it against n6. So a non-positive
    #  or negligible k_e restarts from the prior's own centre instead.
    ke_start = ke if ke > 1e-3 * ke_ref else ke_ref
    free.append(("ke", np.log(ke_start)))
    free.append(("b", b))
    if free_kappa:
        #  ln kappa, so the exact polish can carry the log prior and kappa
        #  stays positive whether or not it has a prior
        free.append(("kappa", np.log(max(kappa, 1e-6))))
    names = [n for n, _ in free]
    q0 = np.array([v for _, v in free])

    def unpack_q(q):
        d = dict(zip(names, q))
        th = np.exp(np.clip(d["theta"], -25, 25)) if "theta" in d else theta_fixed
        return (th, np.exp(np.clip(d["Km"], -25, 25)),
                np.exp(np.clip(d["ke"], -25, 25)), d["b"],
                float(np.exp(np.clip(d["kappa"], -25, 25)))
                if "kappa" in d else 1.0)

    def res(q):
        th, K, k, bb, kap = unpack_q(q)
        r = ((th + kap * phi) * m / (K + m) + k * m + bb - s) / sigma
        extra = [bb / sb, (np.log(k) - np.log(ke_ref)) / 0.8]
        if kappa_prior is not None:
            extra.append(np.log(kap) / kappa_prior)
        return np.concatenate([r, extra])

    lb = np.array([lo_km if n == "Km" else -np.inf for n in names])
    ub = np.array([hi_km if n == "Km" else np.inf for n in names])
    try:
        out = least_squares(res, np.clip(q0, lb + 1e-9, ub - 1e-9),
                            bounds=(lb, ub), method="trf", max_nfev=8000,
                            xtol=1e-14, ftol=1e-14, gtol=1e-14)
    except Exception:
        return (theta, Km, ke, b, kappa,
                exact_obj(theta, Km, ke, b, kappa, m, s, phi, F, kappa_prior))
    th, K, k, bb, kap = unpack_q(out.x)
    J_new = float(2.0 * out.cost)
    J_old = exact_obj(theta, Km, ke, b, kappa, m, s, phi, F, kappa_prior)
    if not np.isfinite(J_new) or J_new > J_old:
        return theta, Km, ke, b, kappa, J_old
    return th, K, k, bb, kap, J_new


def rms_residual(theta, Km, ke, b, kappa, m, s, phi, F=None):
    """RMS standardised residual at the fit, in sigma units."""
    F = FOOT if F is None else F
    pred = (theta + kappa * phi) * m / (Km + m) + ke * m + b
    return float(np.sqrt(np.mean(((pred - s) / F.sigma) ** 2)))


# ==========================================================================
#  Part 0d. The four analysts.
# ==========================================================================
class Analyst:
    """
    One analyst: an increment family, the parameters it frees, and a fit.

    `phi_of(q)` returns the analyst's increment on the visit times at unit
    scale, for shape parameters `q` (a vector on the log scale, possibly
    empty). `n_par` is the total free parameter count, which is what the
    expected residual sqrt((n-p)/n) is formed with.
    """

    def __init__(self, name, train, shape_kind, free_kappa=False,
                 kappa_prior=None, t_end=900.0):
        self.name = name
        self.train = train                       # (times, weights)
        self.shape_kind = shape_kind             # 'fixed' | 'exp' | '2comp'
        self.free_kappa = free_kappa
        self.kappa_prior = kappa_prior
        self.t_end = t_end
        self.n_shape = {"fixed": 0, "exp": 1, "2comp": 2}[shape_kind]
        #  theta, K_m, k_e, b, plus the shape and the scale
        self.n_par = 4 + self.n_shape + int(free_kappa)
        self._ref = None

    def set_reference(self, ref, convention):
        """The committed (single, peak, mean) triple and which one is held."""
        self._ref = (ref, convention)

    def kernel(self, q):
        if self.shape_kind == "fixed":
            return kernel_exp(TH_HALF0)
        if self.shape_kind == "exp":
            return kernel_exp(float(np.exp(q[0])))
        return kernel_twocomp(float(np.exp(q[0])), float(np.exp(q[1])))

    def phi(self, q, times):
        k = self.kernel(q)
        key = (self.shape_kind, tuple(np.round(np.asarray(q, float), 12)))
        A = amplitude(k, self._ref[1], self._ref[0], key=key)
        return A * assemble(k, self.train[0], self.train[1], at_draw(times))

    def q0(self):
        if self.shape_kind == "exp":
            return np.array([np.log(TH_HALF0)])
        if self.shape_kind == "2comp":
            return np.log(np.array([K1.KQ, K1.KDEG]))
        return np.zeros(0)

    #: the coarse SHAPE grid the outer search is made global on. A grid scan
    #: followed by a local polish, rather than a multistart local search, so
    #: that the outer problem has the same standing as the K_m scan: no start
    #: and no basin to miss. Half-lives from 1 to 30 d; the two-compartment
    #: rates over two decades either side of the published pair.
    SHAPE_GRID = {
        "fixed": [np.zeros(0)],
        "exp": [np.array([v]) for v in np.log(np.geomspace(0.5, 40.0, 25))],
        "2comp": [np.array([a, b]) for a in np.log(np.geomspace(
            K1.KQ / 100.0, K1.KQ * 100.0, 13))
            for b in np.log(np.geomspace(K1.KDEG / 100.0, K1.KDEG * 100.0,
                                         13))],
    }

    # ---------------------------------------------------------------- fit
    def _inner(self, q, m, s, times, theta_fixed, grid, do_polish=True,
               phi=None, F=None):
        """
        The inner problem at a fixed shape: K_m scanned globally, the rest
        solved exactly. `do_polish` re-minimises the EXACT objective, which
        differs from the grid's only in carrying the k_e prior on ln k_e rather
        than its delta-method form. It is off during the outer shape search,
        where it costs a nonlinear solve per evaluation and moves the shape
        optimum by less than the tolerance; `fit` polishes once at the end and
        Part 1 reports the objective the polish recovers.
        """
        phi = self.phi(q, times) if phi is None else phi
        p, J, Km = varpro_km(m, s, phi, self.free_kappa, grid=grid, F=F,
                             theta_fixed=theta_fixed,
                             kappa_prior=self.kappa_prior,
                             refine=do_polish)
        j = int(np.argmin(J))
        th, kap, ke, b = unpack(p[j], self.free_kappa, theta_fixed)
        Km = Km if np.isfinite(Km) else float(grid[j])
        if not do_polish:
            return (th, Km, ke, b, kap, float(J[j])), phi
        out = polish(th, Km, ke, b, kap, m, s, phi, self.free_kappa,
                     theta_fixed=theta_fixed, F=F,
                     kappa_prior=self.kappa_prior)
        return out, phi

    def _search(self, m, s, times, theta_fixed, grid, warm=None, stride=1,
                F=None):
        """
        Scan the shape grid, then polish the best node locally.

        `stride` thins the grid for the profile, where the shape is already
        warm-started from the full fit and the scan is a guard against a
        distant basin rather than the search itself. Part 1 reports what
        thinning costs by rerunning one column at stride 1.
        """
        nodes = list(self.SHAPE_GRID[self.shape_kind])[::stride]
        if warm is not None and self.n_shape:
            nodes = [np.asarray(warm, float)] + nodes
        best_q, best_J = None, np.inf
        for q in nodes:
            J = self._inner(q, m, s, times, theta_fixed, grid,
                            do_polish=False, F=F)[0][5]
            if J < best_J:
                best_q, best_J = q, J
        if self.n_shape:
            #  bounded, so the local step cannot leave the box the scan was
            #  global over and land in a corner the grid never saw.
            lo = np.min(np.stack(self.SHAPE_GRID[self.shape_kind]), axis=0)
            hi = np.max(np.stack(self.SHAPE_GRID[self.shape_kind]), axis=0)
            r = minimize(lambda q: self._inner(q, m, s, times, theta_fixed,
                                               grid, do_polish=False,
                                               F=F)[0][5],
                         np.clip(best_q, lo, hi), method="Nelder-Mead",
                         bounds=list(zip(lo, hi)),
                         options=dict(maxiter=60 if stride > 1 else 250,
                                      xatol=1e-4, fatol=1e-7))
            if r.fun < best_J:
                best_q, best_J = r.x, float(r.fun)
        return best_q, best_J

    def fit(self, m, s, times=MONTHLY, grid=KM_GRID, F=None):
        q, _ = self._search(m, s, times, None, grid, F=F)
        (th, Km, ke, b, kap, J), phi = self._inner(q, m, s, times, None, grid,
                                                   F=F)
        best = type("R", (), dict(x=q, fun=J))
        return dict(theta=th, Km=Km, ke=ke, b=b, kappa=kap, J=J, q=best.x,
                    phi=phi)

    # ------------------------------------------------------------ profile
    def profile_at(self, lt, m, s, times, q_warm, grid=KM_GRID_PROF, F=None):
        """
        The profile objective at a fixed ln theta: K_m scanned globally, the
        shape re-optimised, and then the EXACT objective minimised from there.

        The polish is not optional here, as it nearly is at the optimum. The
        linear inner problem must carry the delta-method form of the k_e prior,
        Gaussian on k_e with sd 0.8 k_e, and the paper's prior is Gaussian on
        ln k_e. Down the profile's low tail the compensating direction is k_e
        RISING, to five times its reference and beyond, and there the two
        priors are not the same prior at all: the delta-method form charges
        about seven times what the log prior does, so an unpolished profile
        climbs a wall that is not in the paper's objective. Dropping the polish
        here is the error q9_recheck.py check 5 found.
        """
        theta = float(np.exp(lt))
        q = (np.zeros(0) if self.n_shape == 0
             else self._search(m, s, times, theta, grid, warm=q_warm,
                               stride=PROFILE_STRIDE, F=F)[0])
        return self._inner(q, m, s, times, theta, grid, do_polish=True,
                           F=F)[0][5]

    def interval(self, m, s, fit, times=MONTHLY, span=6.0, bisect=9, F=None):
        """
        The 95% profile interval on theta, found by bracketing rather than on a
        fixed grid: expand from the MLE until the objective has risen by DROP,
        then bisect. `span` bounds the expansion in ln theta; an edge reached
        is reported as an OPEN side, which is the honest answer when the record
        does not bound theta at all.
        """
        lt0 = np.log(max(fit["theta"], 1e-12))
        #  every objective below is the UNPOLISHED one, so the 1.92 drop is
        #  measured on one scale. The polish shifts the whole profile by a
        #  constant to the precision it matters at, and Part 1 prints what it
        #  recovers at the optimum.
        J0 = self.profile_at(lt0, m, s, times, fit["q"], F=F)
        out = []
        for sgn in (-1.0, +1.0):
            lo, hi, step = 0.0, None, 0.25
            while step <= span:
                if self.profile_at(lt0 + sgn * step, m, s, times, fit["q"],
                                   F=F) - J0 >= DROP:
                    hi = step
                    break
                lo, step = step, step * 1.7
            if hi is None:
                out.append((sgn * span, True))
                continue
            for _ in range(bisect):
                mid = 0.5 * (lo + hi)
                if self.profile_at(lt0 + sgn * mid, m, s, times,
                                   fit["q"], F=F) - J0 >= DROP:
                    hi = mid
                else:
                    lo = mid
            out.append((sgn * 0.5 * (lo + hi), False))
        return (float(np.exp(lt0 + out[0][0])), float(np.exp(lt0 + out[1][0])),
                out[0][1], out[1][1])


# ==========================================================================
#  The cells: one (truth, analyst) pair over seeds.
# ==========================================================================
REF_KEYS = ("dose", "peak", "mean")


def committed_ref():
    """The (single-dose peak, steady peak, steady cycle mean) of the committed
    3.5 d kernel on the nominal train, at the committed D_AMP."""
    k = kernel_exp(TH_HALF0)
    return tuple(np.array(norm_stats_uniform(k)) * D_AMP)


def truth_increment(kernel_name, train, convention, ref, times):
    k = KERNELS[kernel_name]
    A = amplitude(k, convention, ref, key=("truth", kernel_name))
    return A * assemble(k, train[0], train[1], at_draw(times))


def make_analysts(train, ref, convention, which=None):
    out = {}
    spec = (("fixed 3.5 d", "fixed", False),
            ("fixed + free amp", "fixed", True),
            ("half-life free", "exp", False),
            ("two-comp free", "2comp", False))
    for name, kind, fk in spec:
        if which is not None and name not in which:
            continue
        a = Analyst(name, train, kind, free_kappa=fk,
                    kappa_prior=KAPPA_PRIOR if fk else None)
        a.set_reference(ref, convention)
        out[name] = a
    return out


def cell(analyst, inc_true, seeds=SEEDS, F=None, times=MONTHLY,
         coverage=True, **rec_kw):
    """
    One (truth, analyst) pair over seeds.

    Returns the ratio, the residual, the log error and the interval outcome per
    seed. Nothing is averaged here; the caller decides what summary to print,
    so a median and its range always travel together.
    """
    F = FOOT if F is None else F
    theta = F.theta
    ratio, resid, lerr, cov, open_lo, hl = [], [], [], [], [], []
    for sd in seeds:
        m, s, _ = make_record(inc_true, sd, F=F, times=times, **rec_kw)
        f = analyst.fit(m, s, times=times, F=F)
        ratio.append(f["theta"] / theta)
        lerr.append(np.log(max(f["theta"], 1e-300)) - np.log(theta))
        resid.append(rms_residual(f["theta"], f["Km"], f["ke"], f["b"],
                                  f["kappa"], m, s, f["phi"], F=F))
        hl.append(np.exp(f["q"][0]) if analyst.shape_kind == "exp" else np.nan)
        if coverage:
            lo, hiv, o_lo, o_hi = analyst.interval(m, s, f, times=times, F=F)
            cov.append(bool(lo <= theta <= hiv))
            open_lo.append(o_lo)
        else:
            cov.append(None)
            open_lo.append(None)
    n_par = analyst.n_par
    ratio = np.array(ratio)
    n_run = int(np.sum(ratio > RUNAWAY))
    n_floor = int(np.sum(ratio < FLOOR))
    #  an RMSE formed over seeds where the fit left its domain is a reading of
    #  the floor and the bound, not of the estimator, so it is withheld rather
    #  than printed with a caveat nobody carries forward.
    rmse = (float(np.sqrt(np.mean(np.array(lerr) ** 2)))
            if n_run == 0 and n_floor == 0 else np.nan)
    return dict(ratio=ratio, resid=np.array(resid),
                lerr=np.array(lerr), cov=cov, open_lo=open_lo,
                hl=np.array(hl), n_run=n_run, n_floor=n_floor, rmse=rmse,
                expected_resid=float(np.sqrt(max(times.size - n_par, 0)
                                             / times.size)))


def _g(x, w, d=3):
    """A number wide enough to read, in fixed point while it fits."""
    return (f"{x:{w}.{d}f}" if abs(x) < 10 ** (w - d - 2)
            else f"{x:{w}.1e}")


def cell_row(name, c, width=18):
    rat, res = c["ratio"], c["resid"]
    covs = [x for x in c["cov"] if x is not None]
    cv = (f"{100 * sum(covs) / len(covs):5.1f}%" if covs else "    .")
    flags = []
    if c["n_run"]:
        flags.append(f"{c['n_run']} runaway")
    if c["n_floor"]:
        flags.append(f"{c['n_floor']} at floor")
    n_open = sum(1 for x in c["open_lo"] if x)
    if n_open:
        flags.append(f"{n_open} open below")
    rmse = "    n/a" if not np.isfinite(c["rmse"]) else f"{c['rmse']:7.3f}"
    return (f"    {name:{width}s} {_g(float(np.median(rat)), 7)} "
            f"[{_g(float(rat.min()), 6)},{_g(float(rat.max()), 7)}] {rmse} "
            f"{np.median(res):7.3f} [{res.min():5.3f},{res.max():5.3f}] "
            f"{cv}" + ("  " + ", ".join(flags) if flags else ""))


def cell_header(width=18, expected=None):
    e = "" if expected is None else f"   (a well-specified fit: {expected:.3f})"
    return (f"    {'analyst':{width}s} {'theta_hat/theta, median [min,max]':>33s} "
            f"{'RMSE':>7s} {'RMS resid / sigma':>23s} {'cover':>6s}" + e)


# ==========================================================================
#  Part 0. The footing, and the checks that the shapes are what they say.
# ==========================================================================
def part_0():
    banner("(0) FOOTING, AND THE CHECKS ON THE SHAPES THEMSELVES")
    print(f"  Local model constants, from src/common/params.py through "
          f"i1_estimator:")
    print(f"    theta = {TH:.0f}, K_m = {KM:.0f}, k_e = {KE:.0f}, "
          f"L = {LB:.4f}, sigma = {SIG:.4f}")
    print(f"    priors: b ~ N(0, {S_B:.4f}); ln k_e ~ N(ln {KE:.0f}, 0.8). "
          f"Both the paper's own.")
    print(f"    design: monthly, n = {MONTHLY.size}, over "
          f"{MONTHLY[-1] - MONTHLY[0]:.1f} d; load CV 10%; dose interval "
          f"{T_DOSE:.0f} d.")
    print(f"    committed increment: D_AMP = {D_AMP:.0f}, half-life "
          f"{TH_HALF0:.2f} d, steady peak {D_AMP / (1 - np.exp(-LAMBDA0 * T_DOSE)):.0f}.")
    print(f"    seeds: {len(SEEDS)}. Every figure below is a median with its "
          f"seed range beside it.")

    print("\n  CHECK 1. The committed kernel and train reproduce "
          "design.run.delta, which is")
    print("  what makes 'the committed shape' one member of this family "
          "rather than a")
    print("  second spelling of it.")
    tn, wn = nominal_train()
    t = np.linspace(0.0, 600.0, 6001)
    mine = (D_AMP / norm_stats_uniform(kernel_exp(TH_HALF0))[0]) * \
        assemble(kernel_exp(TH_HALF0), tn, wn, t)
    theirs = np.asarray(DR.delta(t), float)
    rel = float(np.max(np.abs(mine - theirs)) / np.max(np.abs(theirs)))
    print(f"    worst relative difference over 0 to 600 d: {rel:.3e}")

    print("\n  CHECK 2. The two-compartment kernel by eigen-decomposition "
          "against the same")
    print("  system integrated by solve_ivp, which is the route "
          "audit/k1_sources.py takes.")
    u = _ATT_U[_ATT_U <= 300.0]
    rel2 = float(np.max(np.abs(kernel_attarwala()(u) - _ATT_Y[:u.size]))
                 / np.max(np.abs(_ATT_Y[:u.size])))
    print(f"    worst relative difference over 0 to 300 d: {rel2:.3e}")

    print("\n  CHECK 3. The closed-form normalisation against direct assembly "
          "on a 0.02 d")
    print("  grid. They are NOT expected to agree exactly on the cycle mean: "
          "a trapezoid")
    print("  across the jump a bolus makes in the increment is wrong by about "
          "half the")
    print("  jump times the step, per dose, so the gap should halve when the "
          "grid does.")
    print(f"    {'kernel':18s} {'single peak':>24s} {'steady peak':>24s} "
          f"{'cycle mean':>24s}")
    for nm, k in KERNELS.items():
        sl, fa = norm_stats(k, tn, wn), norm_stats_uniform(k)
        cells = "".join(f" {a:11.5f}/{b:11.5f}" for a, b in zip(sl, fa))
        print(f"    {nm:18s}{cells}")
    fine2 = np.arange(0.0, 900.0, 0.01)
    for nm in ("exp 3.5 d", "Attarwala 2-comp"):
        k = KERNELS[nm]
        g1 = abs(norm_stats(k, tn, wn)[2] / norm_stats_uniform(k)[2] - 1.0)
        g2 = abs(norm_stats(k, tn, wn, grid=fine2)[2]
                 / norm_stats_uniform(k)[2] - 1.0)
        print(f"    {nm:18s} cycle-mean gap {g1:.3e} at h=0.02 d, "
              f"{g2:.3e} at h=0.01 d, ratio {g1 / g2 if g2 else np.nan:.2f}")

    print("\n  CHECK 4. Are any kernels in the analyst's two-compartment box "
          "degenerate or")
    print("  non-decaying? Both would mean the family reaches outside what an "
          "increment")
    print("  can be, and both are recorded rather than silently repaired.")
    for q in Analyst.SHAPE_GRID["2comp"]:
        kernel_twocomp(*np.exp(q))
    print(f"    defective eigen-decompositions on the grid: {len(_DEFECTIVE)}")
    print(f"    non-decaying modes of non-negligible weight: {len(_NONDECAY)}")
    hl = [float(np.log(2.0) / -kernel_twocomp(*np.exp(q)).w.max())
          for q in Analyst.SHAPE_GRID["2comp"]]
    print(f"    terminal half-life spanned by the box: {min(hl):.3f} d to "
          f"{max(hl):.0f} d")

    print("\n  CHECK 5. This script's fit against i1_estimator's own, on the "
          "well-specified")
    print("  case, so that 'the repository's estimator' is a claim and not a "
          "label.")
    ref = committed_ref()
    A = make_analysts((tn, wn), ref, "peak", which=("fixed 3.5 d",))["fixed 3.5 d"]
    phi = A.phi(np.zeros(0), MONTHLY)
    worst_t, worst_k = 0.0, 0.0
    for sd in range(6):
        m, s, _ = make_record(phi, sd)
        f = A.fit(m, s)
        r = I1.varpro_exact(m, phi, s, pen=True)
        worst_t = max(worst_t, abs(f["theta"] / r["theta"] - 1.0))
        worst_k = max(worst_k, abs(f["Km"] / r["Km"] - 1.0))
    print(f"    over 6 records, worst relative difference in theta "
          f"{worst_t:.3e}, in K_m {worst_k:.3e}")

    print("\n  CHECK 6. The Footing object reproduces the committed "
          "constants at theta = 3000,")
    print("  which is what lets Part 4 move the capacity without any constant "
          "being edited.")
    gaps = {"sigma": abs(FOOT.sigma - I1.SIG), "prior on b": abs(FOOT.s_b - S_B),
            "L": abs(FOOT.L - LB), "b_loc": abs(FOOT.bloc - BLOC),
            "m_bar": abs(FOOT.m_bar - M_BAR)}
    print("    " + ",  ".join(f"{k} off by {v:.3e}" for k, v in gaps.items()))
    print(f"    at 5% of nominal it becomes: {Footing(TH * 0.05)}")

    print("\n  CHECK 7. The dose train the mg/kg scale is anchored on.")
    mean_mg, n_doses = published_mean_mgkg()
    print(f"    published per-regimen counts total {n_doses} doses, "
          f"dose-weighted mean {mean_mg:.5f} mg/kg.")
    print(f"    so a train carrying mg/kg weights is put on the committed "
          f"scale by")
    print(f"    dividing by {mean_mg:.5f}, and D_AMP = {D_AMP:.0f} then "
          f"stands for {mean_mg:.3f} mg/kg.")
    for kind in ("C5 de-escalated", "C1 Q3W then Q2W", "composite"):
        tt, aa = irregular_train(kind)
        gaps = np.diff(tt[tt <= 200.0])
        print(f"    {kind:18s} first 200 d: {len(tt[tt <= 200])} doses, "
              f"gaps {gaps.min():.0f} to {gaps.max():.0f} d, "
              f"mg/kg {aa[tt <= 200].min():.2f} to {aa[tt <= 200].max():.2f}")
    return dict(committed_rel=rel, twocomp_rel=rel2, fit_rel=worst_t,
                mean_mgkg=mean_mg)


# ==========================================================================
#  Part 1. The grid, on the nominal uniform train, at normal capacity.
# ==========================================================================
def part_1(convention="peak", seeds=SEEDS, coverage=True):
    banner(f"(1) TRUTH AGAINST ANALYST, NOMINAL 14 d TRAIN, theta = {TH:.0f} "
           f"({convention}-matched)")
    print("  The record is the paper's own monthly design, n = 25, load CV "
          "10%, and the")
    print("  repository's noise model. Only the SHAPE of the increment "
          "differs between")
    print("  the generator and the fit. `theta_hat/theta` above 1 is capacity "
          "read too")
    print("  high. `cover` is the fraction of seeds whose 95% profile "
          "interval contains")
    print(f"  the truth, against a nominal 95%; with {len(seeds)} seeds its "
          f"own standard error")
    print(f"  at 95% is {100 * np.sqrt(0.95 * 0.05 / len(seeds)):.1f} points, "
          f"so this separates 95% from 50% and not")
    print("  95% from 88%.")
    ref = committed_ref()
    tn, wn = nominal_train()
    analysts = make_analysts((tn, wn), ref, convention)
    store = {}
    for tname in KERNELS:
        inc = truth_increment(tname, (tn, wn), convention, ref, MONTHLY)
        print(f"\n  TRUTH: {tname}   (increment over the 25 visits "
              f"{inc.min():.0f} to {inc.max():.0f})")
        first = True
        for aname, A in analysts.items():
            c = cell(A, inc, seeds=seeds, coverage=coverage)
            store[(tname, aname)] = c
            if first:
                print(cell_header())
                first = False
            print(cell_row(aname, c)
                  + f"   p={A.n_par}, expected resid {c['expected_resid']:.3f}")
            if A.shape_kind == "exp":
                print(f"      half-life recovered: "
                      f"{spread(c['hl'], '{:.2f}')} d")
    return store


# ==========================================================================
#  Part 2. The trial's real, irregular dose train.
# ==========================================================================
def part_2(convention="peak", seeds=SEEDS, delay=7.0):
    banner("(2) THE TRIAL'S IRREGULAR DOSE TRAIN")
    print("  Two questions, and they are different. (a) With the analyst "
          "HOLDING the real")
    print("  train, does an irregular train change what a wrong SHAPE costs? "
          "An uneven,")
    print("  amplitude-varying train is better excitation, so it might help. "
          "(b) With the")
    print("  analyst assuming the NOMINAL Q2W train while the truth ran the "
          "real one, what")
    print("  does getting the train itself wrong cost? That is a record "
          "error, not a")
    print("  model error, and it is the one a clinic can actually avoid.")
    ref = committed_ref()
    tn, wn = nominal_train()
    mean_mg, _ = published_mean_mgkg()
    truths = ("exp 3.5 d", "exp 7.0 d", "Attarwala 2-comp")
    keep = ("fixed 3.5 d", "half-life free")
    out = {}
    for kind in ("nominal Q2W", "C5 de-escalated", "C1 Q3W then Q2W",
                 "composite"):
        if kind == "nominal Q2W":
            tt, aa = tn, wn
        else:
            tt, aa = irregular_train(kind, delay=delay)
            aa = aa / mean_mg          # mg/kg onto the committed D_AMP scale
        print(f"\n  TRAIN: {kind}")
        anal_true = make_analysts((tt, aa), ref, convention, which=keep)
        anal_nom = make_analysts((tn, wn), ref, convention, which=keep)
        for tname in truths:
            inc = truth_increment(tname, (tt, aa), convention, ref, MONTHLY)
            print(f"    truth {tname:18s} (increment {inc.min():.0f} to "
                  f"{inc.max():.0f})")
            for lab, anal in (("analyst holds the real train", anal_true),
                              ("analyst assumes nominal Q2W", anal_nom)):
                if kind == "nominal Q2W" and lab.endswith("Q2W"):
                    continue
                for aname, A in anal.items():
                    c = cell(A, inc, seeds=seeds)
                    out[(kind, tname, lab, aname)] = c
                    print("    " + cell_row(f"{aname}, {lab}", c, width=48))
    return out


# ==========================================================================
#  Part 3. What the amplitude convention costs, which is a footing and not a
#  result.
# ==========================================================================
def part_3(seeds=SEEDS):
    banner("(3) THE AMPLITUDE CONVENTION: THE SAME GRID UNDER ALL THREE")
    print("  Comparing two shapes needs a statement of what is held equal, "
          "and there is no")
    print("  right answer. `peak` holds the steady peak increment, so the "
          "capacity ceiling")
    print("  of Section 7 does not move. `dose` holds the single-dose peak, "
          "so each")
    print("  administration deposits the same and a slower kernel "
          "accumulates more.")
    print("  `mean` holds the cycle-mean exposure. Coverage is not run here; "
          "this part is")
    print("  the bias column only, and it is the sensitivity of Part 1 to a "
          "choice.")
    tn, wn = nominal_train()
    ref = committed_ref()
    out = {}
    print(f"\n    {'truth':18s} {'analyst':18s}" +
          "".join(f"{c:>26s}" for c in REF_KEYS))
    for tname in KERNELS:
        for aname in ("fixed 3.5 d", "half-life free", "two-comp free"):
            cells = []
            for conv in REF_KEYS:
                A = make_analysts((tn, wn), ref, conv, which=(aname,))[aname]
                inc = truth_increment(tname, (tn, wn), conv, ref, MONTHLY)
                c = cell(A, inc, seeds=seeds, coverage=False)
                out[(tname, aname, conv)] = c
                r = c["ratio"]
                cells.append(f"{np.median(r):8.3f} [{r.min():6.3f},"
                             f"{r.max():7.3f}]")
            print(f"    {tname:18s} {aname:18s}" + "".join(cells))
    print("\n    The peak-matched column is Part 1. A convention that lets a "
          "slower kernel")
    print("    deliver more enzyme moves the bias with it, so the SIZE of the "
          "bias is a")
    print("    statement about the convention as much as about the shape; the "
          "DIRECTION is")
    print("    the part that should be read as a finding.")
    return out


# ==========================================================================
#  Part 4. The same question at a disease-range capacity.
# ==========================================================================
def part_4(convention="peak", seeds=SEEDS):
    th = TH * 0.05
    banner(f"(4) THE SAME GRID AT theta = {th:.0f}, 5% OF NOMINAL")
    print("  Section 7 already establishes that at this capacity the record "
          "does not bound")
    print("  theta from below at any schedule, so a bias figure here is not "
          "the operative")
    print("  statement and is reported only so the two ranges are not "
          "confused. What is")
    print("  worth reading is the residual column: if a wrong shape is "
          "inaudible at normal")
    print("  capacity, is it more so here, where the enzymatic flux is a "
          "twentieth of the")
    print("  load? The footing is re-solved at the new capacity in "
          "n4_theta_scope's")
    print("  convention A, the physical one: the recorded load stays put and "
          "the marker")
    print("  rises.")
    F = Footing(th, convention="A")
    print(f"    {F}")
    ref = committed_ref()
    tn, wn = nominal_train()
    analysts = make_analysts((tn, wn), ref, convention,
                             which=("fixed 3.5 d", "half-life free"))
    out = {}
    for tname in KERNELS:
        inc = truth_increment(tname, (tn, wn), convention, ref, MONTHLY)
        print(f"\n  TRUTH: {tname}")
        print(cell_header())
        for aname, A in analysts.items():
            c = cell(A, inc, seeds=seeds, F=F)
            out[(tname, aname)] = c
            print(cell_row(aname, c))
    return out


# ==========================================================================
#  Part 5. What measurement would pin the shape.
# ==========================================================================
#: relative (log-scale) measurement error on one increment reading. NOTHING in
#: the record supplies this: brief R2 finds no human hepatic PCC protein or
#: activity time course at all, and REVIEW_LOG.md R4 already retracts the one
#: in-vivo repeatability figure the reviews offered. So it is swept and named.
INC_CV = (0.05, 0.10, 0.20, 0.30, 0.50)


def steady_profile(kernel, ref, convention, tau, period=T_DOSE):
    """The steady increment at cycle phase tau, on the committed footing."""
    A = amplitude(kernel, convention, ref)
    return A * kernel.steady(np.asarray(tau, float), period)


def shape_columns(tau, ref, convention, t_half=TH_HALF0, period=T_DOSE,
                  readout="delta", theta=TH):
    """
    The two sensitivity columns of a set of increment measurements.

    The parameter vector is (ln A, ln t_half): the amplitude has to be carried
    and profiled out, because nothing independently fixes how much increment
    one dose delivers, which is why a single phase can never pin a half-life.

    `readout` is what the assay returns.
      delta      the increment itself, e.g. a hepatic protein or activity, so
                 the amplitude column is flat at 1.
      capacity   theta + Delta, which is what any FUNCTIONAL assay returns, so
                 both columns are damped by Delta/(theta + Delta).
    """
    tau = np.atleast_1d(np.asarray(tau, float))
    h = 1e-5

    def d(lh, la):
        return np.exp(la) * steady_profile(kernel_exp(np.exp(lh)), ref,
                                           convention, tau, period)

    def lny(lh, la):
        v = d(lh, la)
        return np.log(v if readout == "delta" else theta + v)

    lh0 = np.log(t_half)
    c_half = (lny(lh0 + h, 0.0) - lny(lh0 - h, 0.0)) / (2 * h)
    c_amp = (lny(lh0, h) - lny(lh0, -h)) / (2 * h)
    return c_amp, c_half


def separation(tau, ref, convention, period=T_DOSE, readout="delta",
               theta=TH):
    """
    The per-phase log departure between the committed 3.5 d shape and the
    published two-compartment one, AFTER each has been given its best free
    amplitude.

    The free amplitude matters: two shapes differing only by a scale are not
    distinguishable at all, so what discriminates is the log difference
    measured against its own mean, and a design is a PAIR of phases whose
    departures point opposite ways.
    """
    tau = np.atleast_1d(np.asarray(tau, float))

    def ly(k):
        d = steady_profile(k, ref, convention, tau, period)
        return np.log(d if readout == "delta" else theta + d)

    d = ly(kernel_attarwala()) - ly(kernel_exp(TH_HALF0))
    return d - d.mean()


def best_two_phase(tau, ref, convention, readout="delta", **kw):
    """
    The best pair of phases for sd(ln t_half), in closed form.

    For two measurements with rows (a_i, a_i s_i) the 2x2 information matrix
    has determinant a_1^2 a_2^2 (s_1 - s_2)^2, so the variance of ln t_half is
    (a_1^2 + a_2^2) / (a_1 a_2 (s_1 - s_2))^2. That is a scalar in every pair
    at once, so the search is a broadcast rather than a loop, and for the
    `delta` readout, where a is flat, it reduces to the familiar answer: the
    two phases whose sensitivity to the half-life is furthest apart.
    """
    tau = np.atleast_1d(np.asarray(tau, float))
    a, c = shape_columns(tau, ref, convention, readout=readout, **kw)
    s = np.divide(c, a, out=np.zeros_like(c), where=np.abs(a) > 0)
    A2 = a ** 2
    num = A2[:, None] + A2[None, :]
    den = (A2[:, None] * A2[None, :]) * (s[:, None] - s[None, :]) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        var = np.where(den > 0, num / den, np.inf)
    k = int(np.nanargmin(var))
    i, j = k // tau.size, k % tau.size
    return float(np.sqrt(var[i, j])), float(tau[i]), float(tau[j])


def part_5(convention="peak"):
    banner("(5) WHAT MEASUREMENT WOULD PIN THE SHAPE")
    print("  Part 1 says the shape decides the answer and the fit does not "
          "reliably")
    print("  announce a wrong one. So the shape has to come from outside the "
          "marker")
    print("  record. This part prices that measurement: how many direct "
          "readings of the")
    print("  increment, at which phases of the dose cycle, and at what "
          "measurement error.")
    print("\n  THREE THINGS FIXED FIRST, because each changes the answer.")
    print("    (i) THE READOUT.")
    print("        delta     the increment itself: a hepatic PCC protein or "
          "activity. Brief")
    print("                  R2 reports this has NEVER been measured in a "
          "human at any")
    print("                  time point, so this row prices an assay that "
          "does not exist.")
    print("        capacity  theta + Delta, which is what a FUNCTIONAL assay "
          "returns, so the")
    print("                  sensitivity is damped by Delta/(theta + Delta).")
    print("    (ii) THE ACCESSIBLE WINDOW. Phase 0 is the instant of the "
          "infusion and is")
    print("        where the two families differ most, because one jumps and "
          "the other")
    print("        rises. It is also the phase a clinic cannot book. Every "
          "design below is")
    print("        therefore reported against the earliest phase it is "
          "allowed to use.")
    print("    (iii) THE MEASUREMENT CV, swept because nothing supplies it: "
          "brief R2 finds")
    print("        no human hepatic time course, and REVIEW_LOG.md R4 "
          "retracts the one")
    print("        in-vivo repeatability figure the September reviews "
          "offered.")
    ref = committed_ref()
    #  the gate: a truth of 5.0 d read as 3.5 d is already a bias Part 1
    #  measures, so a useful measurement separates them at two sd.
    gate = float(np.log(5.0 / TH_HALF0) / 2.0)

    print("\n  (a) WHERE THE INFORMATION IS, by phase.")
    print(f"      {'phase (d)':>10s} {'Delta(3.5 d)':>13s} "
          f"{'d lnDelta/d ln t_half':>22s} {'ln(2comp/exp)':>15s}")
    taus = np.array([0.0, 0.04, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0,
                     13.0, 13.75, 13.96])
    d35 = steady_profile(kernel_exp(TH_HALF0), ref, convention, taus)
    datt = steady_profile(kernel_attarwala(), ref, convention, taus)
    h = 1e-5
    sens = (np.log(steady_profile(kernel_exp(np.exp(np.log(TH_HALF0) + h)),
                                  ref, convention, taus))
            - np.log(steady_profile(kernel_exp(np.exp(np.log(TH_HALF0) - h)),
                                    ref, convention, taus))) / (2 * h)
    for a, b, c, e in zip(taus, d35, sens, np.log(datt / d35)):
        print(f"      {a:10.2f} {b:13.0f} {c:22.3f} {e:15.3f}")
    print("      The sensitivity to the half-life is ZERO at the phase where "
          "the amplitude")
    print("      alone sets the level and grows with phase, so the two "
          "questions want")
    print("      different draws: a half-life wants the ENDS of the cycle, a "
          "family wants")
    print("      the phase where one shape is still rising and the other "
          "already falling.")

    def need(sd_per_cv, cv, target):
        """Even N >= 2 reaching `target`, from the two-point sd per unit CV."""
        return int(2 * max(1, np.ceil((sd_per_cv * cv / target) ** 2)))

    def best_pair(lo, hi, readout, step=0.02):
        return best_two_phase(np.round(np.arange(lo, hi + 1e-9, step), 3),
                              ref, convention, readout=readout)

    print("\n  (b) ESTIMATING THE HALF-LIFE inside the single-exponential "
          "family. The")
    print("      amplitude is profiled out, so ONE phase identifies nothing "
          "at any N: two")
    print("      is the minimum, and the model being linear in the two log "
          "parameters it")
    print("      is also the optimum, so N is split evenly between two "
          "phases.")
    print(f"      The gate is separating 3.5 d from 5.0 d at two sd, "
          f"sd(ln t_half) <= {gate:.3f}.")
    for readout in ("delta", "capacity"):
        print(f"\n      readout = {readout}")
        print(f"        {'earliest phase':>15s} {'best pair (d)':>20s} "
              f"{'sd per unit CV':>15s}" +
              "".join(f"{f'N at {int(100*c)}%':>11s}" for c in INC_CV))
        for lo in (0.04, 0.25, 1.0, 2.0):
            v, a, b = best_pair(lo, T_DOSE - lo, readout)
            row = "".join(f"{need(v, c, gate):11d}" for c in INC_CV)
            print(f"        {lo:15.2f} {f'{a:.2f} and {b:.2f}':>20s} "
                  f"{v:15.4f}" + row)
    print("\n      N is the total number of increment measurements, half at "
          "each phase.")

    print("\n  (c) TELLING THE TWO FAMILIES APART, which is the harder "
          "question and the one")
    print("      Part 1 says actually decides the answer. Both shapes are "
          "given their best")
    print("      free amplitude first, so what is left is shape alone and the "
          "statistic is")
    print("      the DIFFERENCE of two log readings, whose sd is the "
          "measurement CV times")
    print("      sqrt(2). A 3-sigma call is the target.")
    for readout in ("delta", "capacity"):
        print(f"\n      readout = {readout}")
        print(f"        {'earliest phase':>15s} {'best pair (d)':>20s} "
              f"{'separation/CV':>14s}" +
              "".join(f"{f'N at {int(100*c)}%':>11s}" for c in INC_CV))
        for lo in (0.04, 0.25, 1.0, 2.0):
            grid = np.round(np.arange(lo, T_DOSE - lo + 1e-9, 0.02), 3)
            d = separation(grid, ref, convention, readout=readout)
            j, k = int(np.argmax(d)), int(np.argmin(d))
            sep = float(abs(d[j] - d[k]) / np.sqrt(2.0))
            row = "".join(f"{int(2 * max(1, np.ceil((3.0 * c / sep) ** 2))):11d}"
                          for c in INC_CV)
            print(f"        {lo:15.2f} "
                  f"{f'{grid[j]:.2f} and {grid[k]:.2f}':>20s} {sep:14.4f}"
                  + row)
    print("\n      The phases that separate the FAMILIES are not the phases "
          "that estimate a")
    print("      half-life WITHIN one, and how much the window costs is the "
          "practical")
    print("      content of this part: the information both questions want "
          "sits in the")
    print("      hours after an infusion, which is the one draw the trial's "
          "own schedule")
    print("      already brings a patient in for.")
    return dict(gate=gate)


# ==========================================================================
def main(argv=()):
    quick = "--quick" in argv
    seeds = tuple(range(4)) if quick else SEEDS
    print("=" * 96)
    print("  S1. THE INCREMENT TIME COURSE: ANALYST AGAINST TRUTH")
    print("=" * 96)
    if quick:
        print("  --quick: 4 seeds, for development. Not a reportable run.")
    part_0()
    s1 = part_1(seeds=seeds)
    null = s1[("exp 3.5 d", "fixed 3.5 d")]
    print("\n  THE RESIDUAL REFERENCE, stated once. sqrt((n-p)/n) is the "
          "UNPENALISED value")
    print("  and the paper's fit is penalised, so the well-specified fit sits "
          "above it. The")
    print("  reference to read the residual column against is the measured "
          "null cell, truth")
    print(f"  exp 3.5 d against analyst fixed 3.5 d, which is "
          f"{spread(null['resid'], '{:.3f}')}")
    print(f"  against the sqrt((n-p)/n) = {null['expected_resid']:.3f} an "
          f"unpenalised fit would give.")
    part_2(seeds=seeds)
    part_3(seeds=seeds)
    part_4(seeds=seeds)
    part_5()
    print("\n" + "=" * 96)
    print("  Independent re-check: audit/q9_recheck.py")
    print("=" * 96)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(tuple(sys.argv[1:])))
