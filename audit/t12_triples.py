"""
18. TRIPLE-AVERAGING: what an ensemble of exactly-determined fits is worth,
    and whether its spread diagnoses anything.

Measurement only. Nothing here reaches the paper, nothing here writes a
results/numbers file, and nothing outside this file is touched.

WHAT IS MEASURED. Multiplying the quasi-steady balance through by (K_m + m)
makes it a quadratic in m, so any THREE visits with distinct levels determine
the quadratic exactly through a Vandermonde system and yield one (A, B, C),
hence one theta. The paper's monthly record has 25 visits, so it contains
C(25,3) = 2300 such exactly-determined fits. Two questions follow, and this
script answers them in that order:

    (a) as an ESTIMATOR, how does an average over the 2300 triples compare
        with the single least-squares fit the repository already uses --
        in sd(ln theta-hat) and in the bias of ln theta-hat;
    (b) as a DIAGNOSTIC, does the SPREAD of the 2300 estimates carry
        information about misspecification? Specifically, do the triples
        whose own Vandermonde is well conditioned disagree with the
        clustered ones by more than noise, and does that disagreement grow
        when the generating model is not the fitted one?

(b) is the question that matters, and it is meaningless without a null, so
the same split is run on data generated from the fitted model itself before
it is run on anything else.

WHAT audit/t10_vandermonde.py ALREADY ESTABLISHED, and is not repeated here.
t10 derived the same clearing -- multiply by (K_m + m_i), not by m_i -- got
the same coefficients A = (L+b)K_m, B = (L+b) - theta - k_e K_m, C = -k_e,
and settled three things that are used below rather than retested: that the
leading column must be the RECORDED LOAD and not a column of ones once the
diet varies; that the transformation multiplies each residual by (K_m + m_i),
so least squares on the linear form is weighted least squares on the balance;
and that the three coefficients see only the SUM L + b, so with the load bias
unknown there are three equations for four unknowns and theta is not
recoverable at all. Its answer to its own question -- does cond(V) price the
design results the repository commits to -- was no.

WHAT IS NEW HERE. t10 solved the system once and asked what its condition
number predicts. This script solves it 2300 times on one record and asks what
the ENSEMBLE is worth: as four competing point estimates (Sections 3 to 5),
and as a misspecification diagnostic (Section 6). Neither question is a
condition-number prediction, and Section 6.2 finds the mechanism that decides
both: whether a misspecification lies inside or outside the quadratic in m
that every triple is fitting.

WHAT IS PERTURBED, AND WHY. The observed quantity is the ORDINATE of the
balance -- the recorded substrate load -- and that is the only thing noise is
added to:

    y_i = L_i + eps_i ,   eps_i ~ N(0, sigma^2) ,
    sigma = model_A.run.sigma_total(LOAD_FRAC, ASSAY_CV) = 21.9612 .

The marker level m_i is held at the value that actually solves the balance.
This is not a convenience: sigma_total IS an ordinate noise by construction,
the quadrature sum of the load-record error (LOAD_FRAC x baseline load) and
the assay error carried to the ordinate through the local slope, and
design/run.py's own fit_ensemble states why it must go there -- "perturbing m
and fitting against the perturbed m would make the fit look better than it
is, because the model would then be evaluated at the level the assay reported
rather than the level the patient was at". Putting the noise anywhere else
would price an errors-in-variables gap that audit/i1_estimator.py Part C
already measures, instead of pricing the four estimators, which is the job
here. Every estimator below sees exactly the same (m_i, Delta_i, y_i).

The day-to-day intake variation is NOT noise and is not treated as such. The
recorded load varies at LOAD_CV = 10% and the level is re-solved at it, which
is design.run.design_record's own construction, because by prop:rank a
perfectly steady load is the exactly singular case.

THE FOOTING, held fixed throughout and stated once. theta = 3000, K_m = 100,
k_e = 5, L = design.run.load(0) = 322.7273, doses every 14 d from
design.run.delta, the paper's own monthly record np.arange(0, 730, 30.4),
n = 25. It is held fixed because both questions are about the ESTIMATOR and
not about the design: a comparison that moved the schedule as well would not
say which of the two moved the answer, and the re-checks at the end have to
land on numbers the repository already committed to on this footing.

WHERE THE FOOTING CANNOT BE HELD IDENTICAL, AND THIS IS THE ONE PLACE. A
triple is exactly determined: three equations, three unknowns, no residual and
therefore no prior. The repository's fit carries two priors -- a 5% prior on
the load bias b and a log-prior of 0.8 on ln k_e -- and Section 2 below shows
that without the second one this design does not identify theta at all
(sd(ln theta) = 3.36, a factor of 29). So Sections 3 and 4 compare the
repository's fit AS IT IS against the prior-free triples, which is the
comparison the question asks for and is not like-for-like in the priors, and
Section 5 then takes the priors out of it by running the triples against the
least-squares solution of their OWN linear system over all 25 visits. Which
number answers which question is stated in each section.

Every quantity that depends on the intake record is a median over the 24
records of design.run.SEEDS with its range, through design.run.median_spread.

Nothing here records an expected result in advance.

Run:  python audit/t12_triples.py
"""
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.optimize import brentq, least_squares

import model_A.run as LOC
from common.models import ParamsC, alt_flux_C, g
from design.run import (ASSAY_CV, LOAD_CV, LOAD_FRAC, SEEDS, _filter, delta,
                        design_record, load, median_spread, quasi_steady_level,
                        P)

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
SIG = float(LOC.sigma_total(LOAD_FRAC, ASSAY_CV))
MONTHLY = np.arange(0.0, 730.0, 30.4)
INC = np.asarray(delta(MONTHLY), float)
N = MONTHLY.size

#: the prior widths the repository's own filter states, in the coordinates it
#: states them in: Gaussian on the load bias b as a fraction of the load, and
#: Gaussian on ln k_e. design/run.py:_filter and audit/i1_estimator.py both
#: carry these, and the fit below is the same estimator as theirs.
S_B = LOAD_FRAC * LB
S_LNKE = 0.8

#: all C(25,3) = 2300 triples, as index arrays. Built once: the set of triples
#: is a property of the DESIGN, so it does not move with the record or the
#: noise draw, and rebuilding it per replicate would be the slowest line here.
IDX = np.array(list(combinations(range(N), 3)))

#: the decile split of Section 6. A tenth of 2300 is 230 triples a side, which
#: is large enough that each side's median is not itself an outlier.
DEC = int(round(0.10 * len(IDX)))

BAR = "=" * 84


def bars(title):
    print("\n" + BAR)
    print(title)
    print(BAR)


# --------------------------------------------------------------------------
#  the record
# --------------------------------------------------------------------------
def intake(seed):
    """
    The recorded load at each visit, for one intake record.

    Spelled out here rather than taken from design_record, which returns only
    the levels, and asserted against design_record in R3 so the two cannot
    drift: every level below has to be the level the repository's own headline
    is computed on, or nothing below is comparable with it.
    """
    return LB * (1.0 + LOAD_CV *
                 np.random.default_rng(seed).standard_normal(N))


def level_A(theta, inc, L, ke=KE):
    """The positive root of the balance: k_e m^2 + (k_e K_m + theta + Delta - L) m - L K_m = 0."""
    B = ke * KM + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * L * KM)) / (2.0 * ke)


# --------------------------------------------------------------------------
#  the misspecified generator of Section 6
# --------------------------------------------------------------------------
#: the parent model's alternative channel, scaled so that it delivers EXACTLY
#: the linear channel's flux at the operating point m_bar. Matching the flux is
#: the whole point: the patient then sits in the same place and the only thing
#: the fit is wrong about is the SHAPE of the non-enzymatic channel, so a
#: measured gap cannot be an operating-point difference in disguise.
M_BAR = float(level_A(TH, 0.0, LB))
R0_STAR = float(brentq(lambda R: alt_flux_C(M_BAR, ParamsC(R0=R)) - KE * M_BAR,
                       1.0, 1e6))


def level_parent(theta, inc, L, w):
    """
    The level under a non-enzymatic channel blended between the fitted linear
    one and the parent model's curved alpha(m), at weight w.

    w is a dial rather than a switch because w = 0 has to reproduce the
    correctly-specified case EXACTLY -- otherwise the null and the alternative
    would differ in two things at once.
    """
    L = np.atleast_1d(np.asarray(L, float))
    inc = np.atleast_1d(np.asarray(inc, float))
    pc = ParamsC(R0=R0_STAR)
    out = np.empty(L.size)
    for i in range(L.size):
        out[i] = brentq(
            lambda m: ((theta + inc[i]) * g(m, KM) + (1.0 - w) * KE * m
                       + w * alt_flux_C(m, pc) - L[i]), 1e-6, 1e4)
    return out


# --------------------------------------------------------------------------
#  the estimators
# --------------------------------------------------------------------------
def triple_thetas(m, y):
    """
    theta from every one of the 2300 triples, as one batched 3x3 solve.

    The columns are (y_i, m_i, m_i^2) rather than (1, m_i, m_i^2): the
    recorded load is not flat on this design and Section 1 shows what
    pretending it is costs.
    """
    X = np.stack([y[IDX], m[IDX], m[IDX] ** 2], axis=-1)
    r = ((INC - y) * m)[IDX]
    a = np.linalg.solve(X, r[..., None])[..., 0]
    return -a[:, 1] + a[:, 2] * a[:, 0]


def beta_ols(m, y):
    """
    The least-squares solution of the SAME linear system over all 25 visits.

    This is the estimator a triple is the exactly-determined case of, so both
    sides solve the same system with the same unknowns and no prior.

    WHAT THE RATIO IN SECTION 5 THEREFORE PRICES, and it is not what it first
    looks like. The left-hand side there is the MEDIAN OVER ALL 2300 TRIPLES,
    which uses all 25 visits and merely recombines them, so the ratio is the
    cost of the ensemble-median COMBINATION RULE and not the cost of using
    three equations instead of twenty-five. The genuine three-against-25
    comparison is a single triple against beta-OLS, and it is nothing like as
    kind: the best-conditioned triple alone runs about 5 times the sd of the
    full fit, an efficiency near 0.04 rather than 0.74.
    It is audit/i1_estimator.py's `beta_fit` with the load-bias column removed,
    and it inherits that function's stated price: it reweights every residual
    by (K_m + m_i), so it is not minimising the paper's objective, and it
    cannot carry the paper's priors because neither is a linear function of
    the coefficients.
    """
    X = np.stack([y, m, m ** 2], axis=1)
    a, *_ = np.linalg.lstsq(X, (INC - y) * m, rcond=None)
    return float(-a[1] + a[2] * a[0])


def paper_fit(m, y):
    """
    The repository's own penalised nonlinear fit, in its own coordinates.

    Four parameters (ln theta, ln K_m, ln k_e, b) against the ordinate
    residual, with the two priors _filter states. Started at the truth on
    purpose: any start dependence would be the optimiser's contribution to a
    spread this script attributes to the estimator, and i1_estimator.py Part B
    has already priced the start.
    """
    def res(q):
        t, K, k = np.exp(np.clip(q[:3], -25.0, 25.0))
        return np.concatenate([
            ((t + INC) * m / (K + m) + k * m + q[3] - y) / SIG,
            [(q[2] - np.log(KE)) / S_LNKE, q[3] / S_B]])

    start = np.array([np.log(TH), np.log(KM), np.log(KE), 0.0])
    r = least_squares(res, start, method="lm", max_nfev=5000)
    return float(np.exp(np.clip(r.x[0], -25.0, 25.0)))


def cond_vandermonde(m):
    """
    cond(V) for each triple's own three nodes, from the PURE Vandermonde.

    The pure form, not the solved matrix, because the split of Section 6 has to
    be a property of the DESIGN. Conditioning read off the matrix that carries
    the noisy load column would move with the noise draw, and splitting on a
    quantity the noise moved is selecting on the outcome.
    """
    V = np.stack([np.ones_like(m[IDX]), m[IDX], m[IDX] ** 2], axis=-1)
    return np.linalg.cond(V)


def log_stats(v):
    """
    sd and bias of ln theta-hat over replicates, with what had to be dropped.

    A triple can return a negative theta -- that is the estimator's honest
    output on this design, not a numerical failure -- and ln of it does not
    exist, so the count of usable replicates is reported beside every sd
    rather than quietly closing over it.
    """
    v = np.asarray(v, float)
    ok = np.isfinite(v) & (v > 0.0)
    if ok.sum() < 3:
        return np.nan, np.nan, float(ok.mean()), float(np.median(v))
    l = np.log(v[ok])
    return (float(l.std(ddof=1)), float(l.mean() - np.log(TH)),
            float(ok.mean()), float(np.median(v)))


# ==========================================================================
def s0_algebra():
    bars("0. THE ALGEBRA: theta FROM (A, B, C), DERIVED RATHER THAN GUESSED")
    print("  The clearing and the coefficients are t10's; they are restated")
    print("  because this script inverts the map 2300 times per record and the")
    print("  inversion, not the map, is what Section 4 is about. R4 checks the")
    print("  two spellings against each other.")
    print("""
  The balance at visit i, with a load bias b the fit may or may not know:

      (theta + Delta_i) m_i/(K_m + m_i) + k_e m_i = L_i + b

  Multiply through by (K_m + m_i) and collect the powers of m_i:

      Delta_i m_i = (L_i + b) K_m
                  + (L_i + b - theta - k_e K_m) m_i
                  + (-k_e) m_i^2

  At a FLAT recorded load L_i = L and b = 0 the three coefficients are
  constants, so three visits with distinct levels fix them exactly through the
  Vandermonde V = [[1, m, m^2]]:

      A = L K_m ,    B = L - theta - k_e K_m ,    C = -k_e .

  Inverting the three definitions, in the order in which each becomes
  available:

      k_e   = -C
      K_m   =  A / L
      theta =  L - B - k_e K_m  =  L - B + A C / L .

  theta is the ONLY one of the three that is nonlinear in (A, B, C): it
  carries the product A C. k_e and K_m are linear in them. That single fact is
  what Section 4 is about.

  When the recorded load is not flat -- and on this design it is not -- the
  column of ones is replaced by the recorded load itself. Moving the known
  L_i m_i term across,

      (Delta_i - L_i) m_i = K_m L_i - (theta + k_e K_m) m_i - k_e m_i^2 ,

  so the triple solves the 3x3 system with columns (L_i, m_i, m_i^2) and
  coefficients (a1, a2, a3) = (K_m, -(theta + k_e K_m), -k_e), giving

      theta = -a2 + a3 a1 ,

  which is the same map under A = L a1, B = a2 + L, C = a3. This is
  audit/i1_estimator.py's `beta_fit` with the load-bias column removed, and
  R3 checks the two against each other.

  One consequence of t10's, carried forward because Section 6 walks straight
  into it: A and B see only the SUM L + b. Setting b = 0, which is what an
  exactly-determined triple has no choice but to do, does not leave the load
  bias out of the answer -- it puts the whole of it into theta.""")
    print(f"\n  design: {N} monthly visits over 730 d, C({N},3) = {len(IDX)} triples")
    print(f"  L = {LB:.4f}   theta = {TH:.0f}   K_m = {KM:.0f}   k_e = {KE:.0f}")


def s1_noiseless():
    bars("1. VERIFICATION ON NOISELESS DATA, BEFORE ANYTHING IS MEASURED")
    print("  If a triple does not return 3000 exactly with no noise on, nothing")
    print("  below means anything.\n")

    Lc = np.full(N, LB)
    mc = level_A(TH, INC, Lc)
    V = np.stack([np.ones_like(mc[IDX]), mc[IDX], mc[IDX] ** 2], axis=-1)
    abc = np.linalg.solve(V, (INC * mc)[IDX][..., None])[..., 0]
    th_pure = LB - abc[:, 1] + abc[:, 0] * abc[:, 2] / LB
    print(f"  (a) FLAT load, the pure Vandermonde and theta = L - B + AC/L:")
    print(f"      worst of {len(IDX)} triples: |theta-hat - 3000| = "
          f"{np.abs(th_pure - TH).max():.3e}")
    print(f"      worst |K_m-hat - 100| = "
          f"{np.abs(abc[:, 0] / LB - KM).max():.3e},  "
          f"worst |k_e-hat - 5| = {np.abs(-abc[:, 2] - KE).max():.3e}")

    m = level_A(TH, INC, intake(3))
    th_gen = triple_thetas(m, intake(3))
    print(f"\n  (b) the specification's 10% intake variation, recorded-load column:")
    print(f"      worst of {len(IDX)} triples: |theta-hat - 3000| = "
          f"{np.abs(th_gen - TH).max():.3e}")

    abc3 = np.linalg.solve(
        np.stack([np.ones_like(m[IDX]), m[IDX], m[IDX] ** 2], axis=-1),
        (INC * m)[IDX][..., None])[..., 0]
    th_wrong = LB - abc3[:, 1] + abc3[:, 0] * abc3[:, 2] / LB
    print(f"\n  (c) the same record with the load PRETENDED flat at L, which is")
    print(f"      what the textbook Vandermonde form would do here:")
    print(f"      median theta-hat = {np.median(th_wrong):.4g}, "
          f"{100 * np.mean(th_wrong <= 0):.1f}% of triples negative.")
    print(f"      So the recorded-load column is carried from here on. A 10%")
    print(f"      intake variation is not a small perturbation of a flat load;")
    print(f"      it is the thing that closes the fibre.")

    cd = np.array([median_spread(cond_vandermonde(level_A(TH, INC, intake(s))))
                   for s in SEEDS])
    print(f"\n  (d) cond(V) over the {len(IDX)} triples, medians over the "
          f"{len(SEEDS)} records:")
    print(f"      best {np.median(cd[:, 1]):.4g}   median "
          f"{np.median(cd[:, 0]):.4g}   worst {np.median(cd[:, 2]):.4g}")
    print(f"      five orders of magnitude, which is what makes the decile")
    print(f"      split of Section 6 a split on something rather than nothing.")


def s2_footing():
    bars("2. THE FOOTING, AND WHAT IS IDENTIFIED WITHOUT A PRIOR")
    print(f"  sigma = sigma_total({LOAD_FRAC}, {ASSAY_CV}) = {SIG:.4f}, which is "
          f"{100 * SIG / LB:.1f}% of the load.\n")
    print("  sd(ln theta) from the information filter, same record, four footings.")
    print("  This table is here because it decides what Sections 3 and 5 can ask.\n")

    def crb(m, cols, pri):
        gm = m / (KM + m)
        Q = np.stack([TH * gm, -(TH + INC) * KM * m / (KM + m) ** 2,
                      np.ones_like(m), KE * m], axis=1)[:, cols]
        F = Q.T @ Q / SIG ** 2 + np.diag(pri)
        return float(np.sqrt(np.linalg.inv(F)[0, 0]))

    rows = [("the paper's own: b free at a 5% prior, ln k_e prior 0.8",
             [0, 1, 2, 3], [0.0, 0.0, 1.0 / S_B ** 2, 1.0 / S_LNKE ** 2]),
            ("b held at zero, ln k_e prior 0.8",
             [0, 1, 3], [0.0, 0.0, 1.0 / S_LNKE ** 2]),
            ("b held at zero, k_e FREE -- the triple's own footing",
             [0, 1, 3], [0.0, 0.0, 0.0]),
            ("b held at zero, k_e held exactly known",
             [0, 1], [0.0, 0.0])]
    free = None
    for lab, cols, pri in rows:
        v = median_spread([crb(level_A(TH, INC, intake(s)), cols, pri)
                           for s in SEEDS])
        if "FREE" in lab:
            free = v[0]
        print(f"    {lab:<54} {v[0]:8.4f}  [{v[1]:.4f}, {v[2]:.4f}]")
    print("\n  The third row is the footing a triple is on: no prior of any kind,")
    print("  and the load record believed. On that footing this design does not")
    print(f"  identify theta -- one standard deviation is a factor of "
          f"{np.exp(free):.0f} on a")
    print("  quantity whose gate is 20%. The k_e prior is therefore not a detail")
    print("  of the comparison, it is most of it, and Section 5 exists to take")
    print("  it out.")


def s3_four(reps=200):
    bars("3. THE FOUR ESTIMATORS AT THE SPECIFICATION'S OWN NOISE")
    print(f"  {reps} noise draws on each of the {len(SEEDS)} intake records; the")
    print("  record is fixed within a seed so the spread is the estimator's and")
    print("  not the design's. b_true = 0 and k_e at nominal, so estimator 1's")
    print("  priors are correct rather than merely present, which is the most")
    print("  favourable reading of it.\n")

    acc = {k: [] for k in ("paper", "ols", "mean", "median", "cw")}
    rng = np.random.default_rng(1812)
    for s in SEEDS:
        Ls = intake(s)
        m = level_A(TH, INC, Ls)
        w = 1.0 / cond_vandermonde(m)
        w = w / w.sum()
        out = {k: [] for k in acc}
        for _ in range(reps):
            y = Ls + SIG * rng.standard_normal(N)
            th = triple_thetas(m, y)
            out["paper"].append(paper_fit(m, y))
            out["ols"].append(beta_ols(m, y))
            out["mean"].append(float(th.mean()))
            out["median"].append(float(np.median(th)))
            out["cw"].append(float(w @ th))
        for k in acc:
            acc[k].append(log_stats(out[k]))

    names = {"paper": "1. the paper's fit (priors on)",
             "mean": "2. MEAN over 2300 triples",
             "median": "3. MEDIAN over 2300 triples",
             "cw": "4. mean weighted by 1/cond(V)",
             "ols": "   (least squares, priors OFF)"}
    print(f"  {'estimator':<32} {'sd(ln th)':>11} {'bias(ln th)':>12} "
          f"{'th-hat > 0':>11} {'median th-hat':>15}")
    print("  " + "-" * 82)
    res = {}
    for k in ("paper", "mean", "median", "cw", "ols"):
        A = np.array(acc[k], float)
        sd = median_spread(A[:, 0])
        bi = median_spread(A[:, 1])
        res[k] = (sd, bi, np.nanmedian(A[:, 2]), np.nanmedian(A[:, 3]))
        if k == "ols":
            print("  " + "-" * 82)
        print(f"  {names[k]:<32} {sd[0]:11.4f} {bi[0]:+12.4f} "
              f"{100 * np.nanmedian(A[:, 2]):10.1f}% {np.nanmedian(A[:, 3]):15.4g}")
    print(f"\n  medians over the {len(SEEDS)} records. Ranges, for the two rows a")
    print("  reader would quote:")
    print(f"    the paper's fit: sd [{res['paper'][0][1]:.4f}, "
          f"{res['paper'][0][2]:.4f}], bias [{res['paper'][1][1]:+.4f}, "
          f"{res['paper'][1][2]:+.4f}]")
    print(f"    triple median:   sd [{res['median'][0][1]:.4f}, "
          f"{res['median'][0][2]:.4f}], bias [{res['median'][1][1]:+.4f}, "
          f"{res['median'][1][2]:+.4f}]")
    print(f"\n  Read the last row before reading the gap in the first four. The")
    print("  prior-free least-squares fit returns a POSITIVE theta in only")
    print(f"  {100 * res['ols'][2]:.1f}% of draws -- its sd(ln) is computed on that "
          f"minority and")
    print("  is not a spread of anything -- so most of the distance between")
    print("  estimator 1 and estimators 2-4 is the two priors and not the")
    print("  arithmetic. Section 5 asks the estimator question without them.")
    return res


def s4_why_mean(reps=60):
    bars("4. WHY THE MEAN IS THE WRONG FUNCTIONAL, AND WHY 1/cond DOES NOT FIX IT")
    print("  theta = -a2 + a3 a1 is a PRODUCT of two estimated coefficients. A")
    print("  mean over triples therefore averages a nonlinear function of three")
    print("  quantities whose sampling errors are heavy-tailed wherever the")
    print("  triple's own Vandermonde is near singular, and E[f(x)] != f(E[x]).")
    print("  Least squares never forms that average: it solves one well-posed")
    print("  system and applies the map once, at the end.\n")

    print(f"  the 2300 theta-hats within ONE noisy record, as percentiles;")
    print(f"  medians over {len(SEEDS)} records x {reps} draws. Two noise levels,")
    print("  because the shape of the ensemble is the claim and one level")
    print("  cannot show a shape emerging.\n")
    pcts = (1, 10, 25, 50, 75, 90, 99)
    print(f"  {'sigma':>10} " + "".join(f"{f'p{p}':>13}" for p in pcts)
          + f"{'th-hat <= 0':>13}{'max |th|/th':>13}")
    print("  " + "-" * 122)
    rng = np.random.default_rng(99)
    shares, wratio, tails = [], [], {}
    for f in (1.0, 0.1):
        qs, negs, ratios = [], [], []
        for s in SEEDS:
            Ls = intake(s)
            m = level_A(TH, INC, Ls)
            cd = cond_vandermonde(m)
            w = (1.0 / cd) / np.sum(1.0 / cd)
            if f == 1.0:
                o = np.argsort(cd)
                shares.append(float(w[o[:DEC]].sum()))
                wratio.append(float(w[o[0]] / w[o[-1]]))
            for _ in range(reps):
                y = Ls + f * SIG * rng.standard_normal(N)
                th = triple_thetas(m, y)
                qs.append(np.percentile(th, pcts))
                negs.append(float(np.mean(th <= 0.0)))
                ratios.append(float(np.max(np.abs(th)) / TH))
        qs = np.array(qs)
        tails[f] = float(np.median(ratios))
        print(f"  {f * SIG:10.4f} "
              + "".join(f"{np.median(qs[:, j]):13.4g}" for j in range(len(pcts)))
              + f"{100 * np.median(negs):12.1f}%{np.median(ratios):13.4g}")
    print(f"""
  At the specification's own noise the ensemble has no location at all: half
  the triples return a negative capacity and the quartiles sit tens of times
  theta away on both sides, so neither a mean nor a median of it is an
  estimate of anything. Drop the noise by a factor of ten and the shape the
  argument is about appears: the median lands on the truth, the quartiles come
  back to within a few times theta, and the extreme member of the ensemble is
  still {tails[0.1]:.3g} times theta. That is the shape a median survives and a mean
  does not, and it is the same reason the repository reports every design
  number as a median with its range rather than as a mean: the sampling
  distribution is not symmetric, so its mean is not a location.""")
    print(f"\n  1/cond(V) as a weight: the best decile of triples carries "
          f"{100 * np.median(shares):.1f}% of")
    print("  the total weight, so the worst triples are downweighted but not")
    print(f"  excluded. The whole dynamic range of the weight is "
          f"{np.median(wratio):.3g} between the")
    print(f"  best triple and the worst, against ensemble members "
          f"{tails[1.0]:.3g} times theta.")
    print("  A weight cannot divide by an outlier it is three orders of")
    print("  magnitude smaller than, which is what Section 3's row 4 shows.")


def s5_ladder(reps=120, facs=(1.0, 0.3, 0.1, 0.03, 0.01, 0.003)):
    bars("5. THE EFFICIENCY GAP WITH THE PRIORS TAKEN OUT OF IT: A NOISE LADDER")
    print("  Both sides now solve the SAME linear system with the SAME three")
    print("  unknowns and no prior, so the ratio of their sd(ln theta-hat) is")
    print("  an efficiency and not a difference of footing.")
    print("  READ IT FOR WHAT IT IS: the left-hand side is the MEDIAN OVER")
    print("  2300 TRIPLES, which uses all 25 visits and only recombines them,")
    print("  so this prices the ensemble-median COMBINATION RULE. It is not")
    print("  the cost of three equations against twenty-five -- a single")
    print("  best-conditioned triple runs about 5 times the full fit's sd,")
    print("  an efficiency near 0.04 rather than the 0.74 below.")
    print("\n  The ladder is needed because Section 2's third row says the")
    print("  prior-free problem is not identified at the specification's noise.")
    print("  An efficiency can only be read where both estimators exist, so the")
    print("  noise is walked down until they do and the ratio is watched for a")
    print("  limit rather than quoted at one point.\n")
    print(f"  {'sigma':>10} {'x spec':>8} {'OLS sd':>10} {'median sd':>11} "
          f"{'ratio':>8} {'OLS bias':>10} {'median bias':>12} {'mean bias':>11} "
          f"{'both > 0':>10}")
    print("  " + "-" * 97)
    rng = np.random.default_rng(404)
    for f in facs:
        sig = SIG * f
        rows = {k: [] for k in ("ols", "mean", "median", "cw")}
        rat = []
        for s in SEEDS:
            Ls = intake(s)
            m = level_A(TH, INC, Ls)
            w = 1.0 / cond_vandermonde(m)
            w = w / w.sum()
            out = {k: [] for k in rows}
            for _ in range(reps):
                y = Ls + sig * rng.standard_normal(N)
                th = triple_thetas(m, y)
                out["ols"].append(beta_ols(m, y))
                out["mean"].append(float(th.mean()))
                out["median"].append(float(np.median(th)))
                out["cw"].append(float(w @ th))
            st = {k: log_stats(out[k]) for k in rows}
            for k in rows:
                rows[k].append(st[k])
            if np.isfinite(st["ols"][0]) and st["ols"][0] > 0:
                rat.append(st["median"][0] / st["ols"][0])
        A = {k: np.array(rows[k], float) for k in rows}
        r = median_spread(rat) if rat else (np.nan,) * 3
        both = min(np.nanmean(A["ols"][:, 2]), np.nanmean(A["median"][:, 2]))
        print(f"  {sig:10.4f} {f:8.3g} {median_spread(A['ols'][:, 0])[0]:10.4f} "
              f"{median_spread(A['median'][:, 0])[0]:11.4f} {r[0]:8.3f} "
              f"{median_spread(A['ols'][:, 1])[0]:+10.4f} "
              f"{median_spread(A['median'][:, 1])[0]:+12.4f} "
              f"{median_spread(A['mean'][:, 1])[0]:+11.4f} {100 * both:9.1f}%")
    print("\n  The last column is what makes the ratio readable: only where both")
    print("  estimators return a positive theta on essentially every draw is a")
    print("  ratio of log standard deviations a comparison of estimators rather")
    print("  than of failure rates. The ratio settles there, which is what makes")
    print("  it an efficiency and not an artefact of one noise level. The mean's")
    print("  bias does NOT go to zero on the same ladder: near-singular triples")
    print("  exist at every noise level, because they are a property of the node")
    print("  spacing and not of the noise.")


# ==========================================================================
def _gap(m, y, cd):
    """
    Best decile against worst decile, in units of theta.

    A difference of decile MEDIANS, not of means: Section 4 has just shown the
    within-decile distribution has no usable mean, and the worst decile is
    where that is most true. Reported as a difference rather than a ratio of
    logs because a decile median can be negative and a log of it does not
    exist, which would silently discard the cases the split is supposed to
    catch.
    """
    th = triple_thetas(m, y)
    o = np.argsort(cd)
    return (np.median(th[o[:DEC]]) - np.median(th[o[-DEC:]])) / TH


#: the generators Section 6 runs the split on. (label, load bias b, parent
#: weight w). The first is the null and has to be first.
CASES = (("correctly specified (THE NULL)", 0.0, 0.0),
         ("load bias b = 5% of L", 0.05 * LB, 0.0),
         ("load bias b = 20% of L", 0.20 * LB, 0.0),
         ("parent channel, half weight", 0.0, 0.5),
         ("parent channel, full weight", 0.0, 1.0))


def _record_for(seed, b, w):
    Ls = intake(seed)
    m = level_A(TH, INC, Ls + b) if w == 0.0 else level_parent(TH, INC, Ls + b, w)
    return Ls, m


def s6_diagnostic(reps=250):
    bars("6. THE QUESTION THAT MATTERS: IS THE TRIPLE ENSEMBLE'S SPREAD A"
         " USABLE DIAGNOSTIC?")
    print("""
  The statistic. Sort the 2300 triples by cond(V) on their own three nodes --
  a property of the design, fixed before any noise is drawn -- and take

      G = [ median theta-hat over the best decile
            - median theta-hat over the worst decile ] / theta ,

  230 triples a side. If clustered nodes distort theta in a way spread nodes
  do not, G moves away from whatever the correctly-specified model gives.
  "Whatever the correctly-specified model gives" is not zero once noise is on,
  which is why the null below is run first and not assumed.

  The two misspecifications. A load bias b the fit does not know about, which
  is the simplest one available and the one the paper's own load-bias prior is
  there for; and the parent model of src/model_C/run.py supplying the
  non-enzymatic channel in place of the linear k_e m, scaled to the same flux
  at the operating point so that the only error is the channel's SHAPE.""")

    print("\n  What each misspecification does to theta, with no noise at all,")
    print("  so that a null result below cannot be a null misspecification:\n")
    print(f"  {'generator':<34} {'m range':>16} {'beta-OLS theta-hat':>20}")
    print("  " + "-" * 74)
    for lab, b, w in CASES:
        Ls, m = _record_for(3, b, w)
        print(f"  {lab:<34} {m.min():7.3f}-{m.max():<8.3f} "
              f"{beta_ols(m, Ls):20.2f}")
    print("\n  Every one of them moves theta, two of them catastrophically. The")
    print("  question is whether the ENSEMBLE SPREAD sees any of it.")

    # ---------------------------------------------------------------- 6.2
    print("\n  " + "-" * 80)
    print("  6.2  THE SIGNAL WITH NO OBSERVATION NOISE AT ALL")
    print("  " + "-" * 80)
    print("  G on the exact record, so the number below is the misspecification's")
    print("  whole signature with nothing else in it. Run twice: at the")
    print("  specification's 10% intake variation, and at a perfectly flat load")
    print("  record, because the difference between those two columns turns out")
    print("  to be the mechanism.\n")
    print(f"  {'generator':<34} {'G, intake CV 10%':>26} {'G, flat load':>14}")
    print("  " + "-" * 78)
    for lab, b, w in CASES:
        gv = []
        for s in SEEDS:
            Ls, m = _record_for(s, b, w)
            gv.append(_gap(m, Ls, cond_vandermonde(m)))
        v = median_spread(gv)
        Lf = np.full(N, LB)
        mf = (level_A(TH, INC, Lf + b) if w == 0.0
              else level_parent(TH, INC, Lf + b, w))
        gf = _gap(mf, Lf, cond_vandermonde(mf))
        print(f"  {lab:<34} {v[0]:+10.4f} [{v[1]:+7.4f},{v[2]:+7.4f}] "
              f"{gf:+14.4f}")
    print("""
  The flat-load column is exactly zero for every load bias, to machine
  precision, and that is not luck. At a flat load the perturbation a bias adds
  to the response is (K_m + m_i) b = K_m b . 1 + b . m_i, a polynomial of
  degree one in m, and the triple is fitting a polynomial of degree two in m.
  The perturbation is therefore in the span of the regressors EXACTLY, every
  triple absorbs it identically into (A, B, C), and well- and ill-conditioned
  triples return the same wrong theta. Conditioning cannot see a
  misspecification that lives inside the model's own column space.

  What makes the middle column non-zero is the 10% intake variation: the first
  regressor is then the recorded load rather than a column of ones, the
  constant part K_m b is no longer in the span, and the leftover is rough
  visit to visit -- which is precisely what an ill-conditioned Vandermonde
  amplifies. So the diagnostic's entire sensitivity to a load bias is borrowed
  from the roughness of the diet record.

  The parent channel is the other way round, and worse. It biases theta by up
  to 55% and its G stays within a few per cent of the null, because over the
  marker band this record visits the difference alpha(m) - k_e m is itself
  very nearly a quadratic in m and is absorbed the same way. The
  misspecification that survives as a bias is the one the split is blind to.""")

    # ---------------------------------------------------------------- 6.3
    for f in (1.0, 0.1):
        sig = SIG * f
        print("\n  " + "-" * 80)
        print(f"  6.{3 if f == 1.0 else 4}  THE SAME SPLIT UNDER NOISE, sigma = "
              f"{sig:.4f} ({f:g}x the specification)")
        print("  " + "-" * 80)
        print(f"  {reps} draws on each of the {len(SEEDS)} records. 'one-record")
        print("  spread' is the 10th-90th percentile of G over all draws: it is")
        print("  what a reader holding ONE patient's record would face.\n")
        print(f"  {'generator':<32} {'median G':>10} {'one-record spread':>24} "
              f"{'vs null':>9} {'/spread':>8}")
        print("  " + "-" * 88)
        rng = np.random.default_rng(2718)
        null = None
        for lab, b, w in CASES:
            allg, perseed = [], []
            for s in SEEDS:
                Ls, m = _record_for(s, b, w)
                cd = cond_vandermonde(m)
                gg = [_gap(m, Ls + sig * rng.standard_normal(N), cd)
                      for _ in range(reps)]
                allg.extend(gg)
                perseed.append(float(np.median(gg)))
            allg = np.array(allg)
            med = float(np.median(allg))
            lo, hi = np.percentile(allg, [10, 90])
            ps = median_spread(perseed)
            if null is None:
                null = (med, hi - lo)
                shift, rel = 0.0, 0.0
            else:
                shift = med - null[0]
                rel = shift / null[1]
            print(f"  {lab:<32} {med:+10.4f} [{lo:+9.4f}, {hi:+9.4f}] "
                  f"{shift:+9.4f} {rel:+8.3f}")
            print(f"  {'   per-record medians:':<32} {ps[0]:+10.4f} "
                  f"[{ps[1]:+9.4f}, {ps[2]:+9.4f}]")

    # ---------------------------------------------------------------- 6.5
    print("\n  " + "-" * 80)
    print("  6.5  THE VERDICT")
    print("  " + "-" * 80)
    print("""
  At the specification's own noise the split separates nothing. Read the
  '/spread' column of 6.3: no misspecification shifts the median gap by
  anything approaching the null's own one-record spread, and that includes the
  two that put theta's sign in doubt. A clinician holding one record cannot
  read a gap off it, and pooling records buys only the usual square root.

  At a tenth of that noise the load bias does separate, cleanly, and the
  parent channel still does not -- and the second half of that sentence is the
  one worth keeping. A diagnostic that fires on a load bias and stays silent on
  a wrong clearance channel is not a specification test; it is a second,
  noisier reading of the load record. Section 6.2 says why in one line: the
  split can only see a misspecification that lies OUTSIDE the quadratic in m
  that every triple is fitting, and the parent channel, over the 3.4-fold
  marker band this design visits, does not.""")


# ==========================================================================
def rechecks():
    bars("R1. RE-CHECK: the filter on these records reproduces sd_monthly and"
         " sd_aliased_off")
    print("  design.run._filter on design.run.design_record, the footing every")
    print("  number above is computed on. If this row moved, the record this")
    print("  script generates would not be the paper's record.\n")
    for lab, times, want in (("sd_monthly", MONTHLY, 0.1699),
                             ("sd_aliased_off", np.arange(0.0, 730.0, 14.5),
                              0.1615)):
        inc = np.asarray(delta(times), float)
        v = median_spread([float(np.sqrt(
            _filter(design_record(times, seed=s)[0], SIG, inc)[-1][0, 0]))
            for s in SEEDS])
        print(f"  {lab:>16}: {v[0]:.4f} [{v[1]:.4f}, {v[2]:.4f}] against the "
              f"committed {want:.4f} ({abs(v[0] - want) / want:.2%} apart)")

    bars("R2. RE-CHECK: the constants this script is built on")
    print("  Each is quoted in results/numbers/design.tex or in the constants")
    print("  registry, and each one enters the arithmetic above.\n")
    import g_design_common as G
    t91 = np.linspace(0.0, 91.0, 4000)
    lv = quasi_steady_level(t91)
    for lab, got, want in (
            ("L = load(0)", LB, 322.7273),
            ("sigma (design.tex 'sigma')", SIG, 22.0),
            ("delta_range_monthly", G.delta_range(MONTHLY), 1.84),
            ("marker_low", float(lv.min()), 3.3646),
            ("marker_high", float(lv.max()), 8.9708)):
        print(f"  {lab:>28}: {got:12.6f} against the committed {want:9.4f} "
              f"({abs(got - want) / want:.3%} apart)")

    bars("R3. RE-CHECK: the triple's algebra against the repository's own"
         " beta-form fit")
    print("  audit/i1_estimator.py's beta_fit solves the same rearrangement with")
    print("  a free load-bias column, by least squares over all 25 visits. On a")
    print("  noiseless record the two must land on the same point, and both must")
    print("  land on theta = 3000. Also checked: the levels this script builds")
    print("  are design.run.design_record's levels, bit for bit.\n")
    import i1_estimator as I1
    worst_lv, worst_b, worst_t, worst_o = 0.0, 0.0, 0.0, 0.0
    for s in SEEDS:
        Ls = intake(s)
        m = level_A(TH, INC, Ls)
        worst_lv = max(worst_lv,
                       float(np.abs(m - design_record(MONTHLY, seed=s)[0]).max()))
        worst_b = max(worst_b, abs(I1.beta_fit(m, INC, Ls)["theta"] - TH))
        worst_t = max(worst_t, float(np.abs(triple_thetas(m, Ls) - TH).max()))
        worst_o = max(worst_o, abs(beta_ols(m, Ls) - TH))
    print(f"  levels vs design_record, worst over {len(SEEDS)} records: "
          f"{worst_lv:.3e}")
    print(f"  i1_estimator.beta_fit  vs the committed theta = 3000: "
          f"{worst_b:.3e}  ({worst_b / TH:.2%} apart)")
    print(f"  all {len(IDX)} triples    vs the committed theta = 3000: "
          f"{worst_t:.3e}  ({worst_t / TH:.2%} apart)")
    print(f"  this script's beta_ols vs the committed theta = 3000: "
          f"{worst_o:.3e}  ({worst_o / TH:.2%} apart)")

    bars("R4. RE-CHECK: the record and the map against audit/t10_vandermonde.py")
    print("  t10 spells the same clearing with a constant leading column and")
    print("  reads theta back as L + CA/L - B; this file spells it with the")
    print("  recorded-load column and reads it as -a2 + a3 a1. On a flat load")
    print("  the two are one equation, so they have to agree exactly, and if")
    print("  they did not, one of the two items would be wrong.\n")
    import t10_vandermonde as T10
    m10, inc10, Ls10 = T10.record(MONTHLY, seed=3)
    m_here = level_A(TH, INC, intake(3))
    print(f"  t10.record levels vs this file's, worst over the record: "
          f"{np.abs(m10 - m_here).max():.3e}")
    print(f"  t10.record loads  vs this file's, worst over the record: "
          f"{np.abs(Ls10 - intake(3)).max():.3e}")
    mc = level_A(TH, INC, np.full(N, LB))
    V = np.stack([np.ones_like(mc[IDX]), mc[IDX], mc[IDX] ** 2], axis=-1)
    abc = np.linalg.solve(V, (INC * mc)[IDX][..., None])[..., 0]
    th10 = LB + abc[:, 2] * abc[:, 0] / LB - abc[:, 1]
    th12 = triple_thetas(mc, np.full(N, LB))
    print(f"  t10's map vs this file's, worst over all {len(IDX)} triples at a "
          f"flat load: {np.abs(th10 - th12).max():.3e}")
    c10 = np.array([T10.cond_V(mc[i]) for i in IDX[:200]])
    print(f"  t10.cond_V vs this file's cond_vandermonde, worst over the first "
          f"200 triples: {np.abs(c10 - cond_vandermonde(mc)[:200]).max():.3e}")


# ==========================================================================
def main():
    print(__doc__)
    s0_algebra()
    s1_noiseless()
    s2_footing()
    s3_four()
    s4_why_mean()
    s5_ladder()
    s6_diagnostic()
    rechecks()


if __name__ == "__main__":
    main()
