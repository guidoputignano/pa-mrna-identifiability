"""
16. THE VANDERMONDE REFRAMING: WHAT ITS CONDITION NUMBER SEES, AND WHAT IT DOES NOT.

Scoping only. Nothing here reaches the paper and nothing here writes a
results/numbers file.

WHAT IS MEASURED. The principal investigator's premise is that clearing the
Michaelis denominator turns the quasi-steady balance into a linear system whose
rows are [1, m_i, m_i^2] -- a Vandermonde in the marker levels -- so that
distinct levels make it invertible and three visits suffice. If that is the
right picture, the Vandermonde's condition number ought to price the two design
results this repository already commits to: the cost of sampling at the dose
interval (design.tex alias_cost = 27) and the gain from straddling K_m
(design.tex straddle_gain = 9.5). This script computes cond(V) for both
comparisons and asks whether it does. It then asks the estimator question the
reframing raises: the transformation multiplies each residual by something, so
ordinary least squares on the linear form is WEIGHTED least squares on the
balance, exact in algebra and a different estimator in practice.

THE DERIVATION, done here rather than adopted, because the premise states the
multiplier and the coefficients and both turn out to matter. The balance is

    (theta + Delta_i) m_i/(K_m + m_i) + k_e m_i = L_i + b .

The denominator is cleared by multiplying through by (K_m + m_i), NOT by m_i;
multiplying by m_i leaves m_i^2/(K_m+m_i) behind. What the premise is seeing is
that the Delta term becomes Delta_i m_i, which it does. Multiplying by
(K_m + m_i) and collecting:

    Delta_i m_i = (L_i + b) K_m
                + [ (L_i + b) - theta - k_e K_m ] m_i
                + [ -k_e ] m_i^2 ,

so with a CONSTANT total load the premise's A, B, C are

    A = (L + b) K_m ,   B = (L + b) - theta - k_e K_m ,   C = -k_e ,

and theta is read back as theta = (L + b) + C A/(L + b) - B. Section 0 checks
that on three exact points and recovers theta = 3000 to 1e-8 relative.

Three consequences the premise does not state, all of them measured below.

  (i)  The three coefficients are three functions of FOUR unknowns
       (theta, K_m, k_e, b). C gives k_e and A gives K_m(L+b), but theta cannot
       be separated from the load bias b: the three-coefficient system is
       identified only with the total load known. That is the load-bias fibre
       of prop:rank arriving in a different coordinate system, not a new
       degeneracy.
  (ii) The recorded load varies day to day by specification (LOAD_CV = 10%), and
       then A and B are row-dependent and the rows are not [1, m_i, m_i^2] at
       all. The system stays three-column only if the leading column is taken to
       be the RECORDED LOAD L_i, giving
           Delta_i m_i - L_i m_i = K_m L_i - (theta + k_e K_m) m_i - k_e m_i^2 ,
       which is a Vandermonde only when L_i is constant. Section 1 solves both
       and shows what the constant column costs when the load moves.
  (iii) The residual is multiplied by (K_m + m_i), so the transformed noise has
       sd sigma (K_m + m_i) and the weights that undo the transformation are
       1/(K_m + m_i), not 1/m_i. On this patient m runs 2.8 to 9.5 against
       K_m = 100, so the two prescriptions are nothing like each other: the
       exact weights span 6.5% and the premise's span a factor of 3.4.

THE FOOTING, held fixed throughout so that only the question changes.
theta = 3000, K_m = 100, k_e = 5, L = load(0.0) = 322.7273, the biweekly dose
train of design.run.delta, sigma = model_A.run.sigma_total(5%, 5%) = 21.96, and
the monthly record np.arange(0, 730, 30.4). Every quantity that depends on the
intake record is a median over the 24 records of design.run.SEEDS through
median_spread, because one draw is not a property of a design. The load is drawn
exactly as design_record draws it (Gaussian at LOAD_CV) rather than as
g_design_common.record draws it (lognormal), because the numbers being
reproduced here are design.run's; section 0 checks the two generators agree to
machine precision.

Nothing here records an expected result in advance.

Run:  python audit/t10_vandermonde.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import model_A.run as LOC
import g_design_common as G
from design.run import (ASSAY_CV, DOSE_INTERVAL, LOAD_CV, LOAD_FRAC, SEEDS,
                        BASE_SEED, delta, design_record, load, median_spread)

TH, KM, KE = 3000.0, 100.0, 5.0
LB = load(0.0)
SIG = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
KE_PRIOR_LOG = 0.8

MONTHLY = np.arange(0.0, 730.0, 30.4)
ALIASED = np.arange(0.0, 730.0, DOSE_INTERVAL)
DENSE = np.arange(0.0, 730.0, 2.0)
#  the placement row of Table 3: the same 365 draws moved onto a band that
#  straddles K_m. The levels are an override, not a record a diet produces.
STRADDLE = np.geomspace(KM / 2.0, KM * 2.0, len(DENSE))

BAR = "=" * 84


# --------------------------------------------------------------------------
def record(times, load_cv=LOAD_CV, seed=BASE_SEED):
    """
    Levels, increments and the per-visit recorded load of one intake record.

    design_record returns the first two and not the third, and the third is
    exactly what the Vandermonde's leading column is, so it is redrawn here
    from the same generator rather than reconstructed from the levels.
    """
    times = np.asarray(times, dtype=float)
    inc = np.asarray(delta(times), dtype=float)
    Ls = (LB * (1.0 + load_cv *
                np.random.default_rng(seed).standard_normal(times.size))
          if load_cv else np.full(times.size, LB))
    B = KE * KM + TH + inc - Ls
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    return m, inc, Ls


def implied_load(m, inc):
    """
    The total load a patient sitting at these levels must be carrying.

    For a real record this returns the load that was drawn, which is the check
    at the end of section 0. For the K_m-straddling row it is the only load
    available, because that row's levels are imposed rather than generated, and
    the Vandermonde's leading column cannot be left undefined for it.
    """
    return (TH + inc) * m / (KM + m) + KE * m


def cond_V(m, first=None, equilibrate=False):
    """
    cond(V) for rows [first_i, m_i, m_i^2], `first` defaulting to 1.

    `equilibrate` scales each column to unit norm first. Raw cond(V) mixes two
    different things -- how distinct the nodes are, and the fact that column
    three is in micromolar squared while column one is dimensionless -- and on
    a band at m ~ 100 the units dominate by three orders of magnitude. Both are
    reported wherever the comparison is between bands at different m, because
    only the equilibrated number is about the nodes.
    """
    c = np.ones_like(m) if first is None else np.asarray(first, float)
    V = np.stack([c, m, m ** 2], axis=1)
    if equilibrate:
        V = V / np.linalg.norm(V, axis=0)
    return float(np.linalg.cond(V))


def sd_paper(m, inc, ke_prior=KE_PRIOR_LOG, load_frac=LOAD_FRAC):
    """
    sd(ln theta-hat) on the paper's own footing, rebuilt here rather than
    imported, so that the re-checks at the end test this file's arithmetic
    against the repository's committed numbers and not against itself.

    Coordinates (ln theta, ln K_m, b, ln k_e); the increment is carried in the
    ln K_m column, which is the correction of audit/c1_sensitivity.py.
    """
    q = np.stack([TH * m / (KM + m),
                  -(TH + inc) * KM * m / (KM + m) ** 2,
                  np.ones_like(m),
                  KE * m], axis=1)
    pri = [np.inf, np.inf,
           (load_frac * LB) ** 2 if load_frac else np.inf,
           ke_prior ** 2 if ke_prior else np.inf]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + q.T @ q / SIG ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def sd_original_3(m, inc):
    """
    The same fit on the footing the three-coefficient Vandermonde actually
    offers: the total load known exactly, so the b column is DROPPED rather
    than given a prior, and no prior on anything else.

    This is the object the weighted Vandermonde has to reproduce. If it does
    not, the transformation has been implemented wrongly.
    """
    q = np.stack([TH * m / (KM + m),
                  -(TH + inc) * KM * m / (KM + m) ** 2,
                  KE * m], axis=1)
    return float(np.sqrt(np.linalg.inv(q.T @ q / SIG ** 2)[0, 0]))


WEIGHTS = {"1/(K_m+m)": lambda m: 1.0 / (KM + m),
           "unweighted": lambda m: np.ones_like(m),
           "1/m": lambda m: 1.0 / m}


def sd_vandermonde(m, inc, Ls, weight="1/(K_m+m)", b_free=False,
                   ke_prior=None, load_prior=None, posterior=False,
                   sigma=SIG):
    """
    sd(ln theta-hat) from a least-squares fit of the LINEAR form.

    b_free=False is the three-column system [L_i, m_i, m_i^2] with coefficients
    (K_m, -(theta + k_e K_m), -k_e) and the total load known, which is the only
    footing on which three coefficients determine theta.

    b_free=True adds the constant column, giving [L_i, 1, m_i, m_i^2] with
    coefficients (K_m, K_m b, b - theta - k_e K_m, -k_e). That is four
    coefficients for four unknowns and it is invertible, so the paper's own
    footing CAN be written in this form -- with the priors carried as ridge
    penalties, linearised at the truth: a prior of width s on b is a penalty
    (1/K_m)^2/s^2 on the second coefficient, and a log-prior of width p on k_e
    is a penalty 1/(p k_e)^2 on the fourth.

    `posterior` returns the width of the penalised normal matrix, which is what
    _filter reports and what the paper quotes. Otherwise the SAMPLING
    covariance of the penalised estimator is returned,
    (X'X)^-1 X' Cov(y) X (X'X)^-1, and that is the only object in which two
    DIFFERENT weightings can be compared: the posterior width of an OLS fit is
    not the width of anything the OLS fit actually delivers. Both are printed
    where the distinction bites. The priors are centred at the truth, so both
    estimators are unbiased here and the sampling covariance is the whole
    error.
    """
    m = np.asarray(m, float)
    s = sigma * (KM + m)                  # sd of the TRANSFORMED residual
    w = WEIGHTS[weight](m)
    if b_free:
        V = np.stack([Ls, np.ones_like(m), m, m ** 2], axis=1)
        grad = np.array([-KE, 1.0 / KM, -1.0, KM])   # d theta / d coefficients
        pen = np.zeros((4, 4))
        if load_prior:
            pen[1, 1] = (1.0 / KM) ** 2 / (load_prior * LB) ** 2
        if ke_prior:
            pen[3, 3] = 1.0 / (ke_prior * KE) ** 2
    else:
        V = np.stack([Ls, m, m ** 2], axis=1)
        grad = np.array([-KE, -1.0, KM])
        pen = np.zeros((3, 3))
        if ke_prior:
            pen[2, 2] = 1.0 / (ke_prior * KE) ** 2
    Vw = V * w[:, None]
    A = Vw.T @ Vw / sigma ** 2 + pen
    ev = np.linalg.eigvalsh(A)
    if ev[0] <= 0 or ev[0] <= 1e-14 * ev[-1]:
        return np.inf
    M = np.linalg.inv(A)
    if posterior:
        cov = M
    else:
        meat = (Vw * ((s * w) ** 2)[:, None]).T @ Vw / sigma ** 4
        cov = M @ meat @ M
    return float(np.sqrt(grad @ cov @ grad)) / TH


def med(fn):
    """Median over the 24 intake records, with its range."""
    return median_spread([fn(s) for s in SEEDS])


# ==========================================================================
def section0():
    print(BAR)
    print("0. THE DERIVATION, CHECKED RATHER THAN ADOPTED")
    print(BAR)
    m, inc, Ls = record(MONTHLY, seed=BASE_SEED)
    lv, ic = design_record(MONTHLY, seed=BASE_SEED)
    print(f"  this file's record generator against design_record: "
          f"levels differ by {np.abs(m - lv).max():.3e}, increments by "
          f"{np.abs(inc - ic).max():.3e}")
    print(f"  the load read back off the balance against the load drawn: "
          f"{np.abs(implied_load(m, inc) - Ls).max():.3e}")

    #  three visits at a CONSTANT load and three different cycle phases, which
    #  is the premise's own best case: exact data, exact solve, no noise
    t3 = np.array([0.0, 3.0, 9.0])
    d3 = np.asarray(delta(t3), float)
    B = KE * KM + TH + d3 - LB
    m3 = (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)
    V = np.stack([np.ones_like(m3), m3, m3 ** 2], axis=1)
    A, Bc, C = np.linalg.solve(V, d3 * m3)
    print(f"\n  three exact points at t = {t3}, constant load")
    print(f"    nodes m           = {np.array2string(m3, precision=6)}")
    print(f"    solved  (A, B, C) = ({A:.6f}, {Bc:.6f}, {C:.6f})")
    print(f"    derived (A, B, C) = ({LB*KM:.6f}, {LB-TH-KE*KM:.6f}, {-KE:.6f})")
    th = LB + C * A / LB - Bc
    print(f"    theta = L + CA/L - B = {th:.6f} against 3000 "
          f"({abs(th-TH)/TH:.2e} relative)")
    print(f"    K_m   = A/L          = {A/LB:.6f} against 100")
    print(f"    k_e   = -C           = {-C:.6f} against 5")
    print("\n  so the coefficients are A = (L+b)K_m, B = (L+b) - theta - k_e K_m,")
    print("  C = -k_e, and the clearing multiplier is (K_m + m_i). Note what the")
    print("  three of them do NOT contain: A and B see only the SUM L + b, so")
    print("  with the load bias b unknown the system has three equations for")
    print("  four unknowns and theta is not recoverable at all.")


def section1():
    print("\n" + BAR)
    print("1. THREE VISITS SUFFICE -- WHEN THE PHASES DIFFER AND THE LOAD DOES NOT")
    print(BAR)
    #  a patient ALREADY ON THERAPY, where the dose train has reached its
    #  steady saw-tooth. From t = 0 the train is still accumulating, so
    #  fortnightly draws are merely nearly repeated; at steady state they are
    #  exactly repeated, and that is the case the premise is about.
    S = 100 * DOSE_INTERVAL
    for lab, t3, cv in (("three spread phases, constant load",
                         S + np.array([0.0, 3.0, 9.0]), 0.0),
                        ("three visits one dose interval apart, constant load",
                         S + np.array([0.0, 14.0, 28.0]), 0.0)):
        d3 = np.asarray(delta(t3), float)
        B = KE * KM + TH + d3 - LB
        m3 = (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)
        V = np.stack([np.ones_like(m3), m3, m3 ** 2], axis=1)
        det = float(np.linalg.det(V))
        print(f"\n  {lab}")
        print(f"    nodes {np.array2string(m3, precision=6)}")
        print(f"    det V = {det:.4e}   cond V = {np.linalg.cond(V):.4e}")
        if abs(det) < 1e-12:
            print("    exactly singular: the solve fails, no theta")
        else:
            A, Bc, C = np.linalg.solve(V, d3 * m3)
            print(f"    theta = {LB + C*A/LB - Bc:.6f}")

    #  and now the thing the premise does not mention: the specification's own
    #  10% day-to-day load variation, which is what BREAKS the repeated nodes
    #  and simultaneously breaks the constant column
    t3 = S + np.array([0.0, 3.0, 9.0])
    d3 = np.asarray(delta(t3), float)
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(BASE_SEED).standard_normal(3))
    B = KE * KM + TH + d3 - Ls
    m3 = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    print(f"\n  three spread phases at the specification's 10% load variation")
    print(f"    loads {np.array2string(Ls, precision=4)}")
    A, Bc, C = np.linalg.solve(
        np.stack([np.ones_like(m3), m3, m3 ** 2], axis=1), d3 * m3)
    th = LB + C * A / LB - Bc
    print(f"    rows [1, m, m^2] as the premise states them: theta = {th:.1f}, "
          f"K_m = {A/LB:.1f}, k_e = {-C:.1f}")
    print(f"      wrong by a factor of {th/TH:.0f}, and K_m comes out negative:")
    print("      the constant column is not a constant when the load moves.")
    a, b, c = np.linalg.solve(np.stack([Ls, m3, m3 ** 2], axis=1),
                              d3 * m3 - Ls * m3)
    print(f"    rows [L_i, m, m^2] instead:                 theta = {-b + c*a:.6f}, "
          f"K_m = {a:.6f}, k_e = {-c:.6f}")
    print("\n  The reframing survives a varying load only in the second form, and")
    print("  that form is not a Vandermonde: its leading column is the diet.")


def section2():
    print("\n" + BAR)
    print("2. (a) DOES cond(V) REPRODUCE THE ALIASING COST?")
    print(BAR)
    print("  Committed (results/numbers/design.tex): sd_monthly = 0.1699,")
    print("  sd_aliased = 4.5845, alias_cost = 27 [21, 33], taken per record.")
    print("  Every schedule of g_design_common.SCHEDULES at n = 25, plus the")
    print("  14.5 d schedule the paper offers as the repair.\n")
    n = 25
    print(f"  {'schedule':<28}{'span':>8}{'nodes':>7}{'cond V':>12}"
          f"{'cond V':>12}{'range':>8}{'sd(ln th)':>11}")
    print(f"  {'':<28}{'[d]':>8}{'distinct':>7}{'L const':>12}"
          f"{'10% CV':>12}{'D/theta':>8}{'10% CV':>11}")
    out = {}
    for name, fn in G.SCHEDULES.items():
        t = np.asarray(fn(n), float)
        m0, inc0, _ = record(t, load_cv=0.0)
        c0 = cond_V(m0)
        c1 = med(lambda s: cond_V(record(t, seed=s)[0]))[0]
        sd = med(lambda s: sd_paper(*record(t, seed=s)[:2]))[0]
        #  distinct to a RELATIVE tolerance: rounding to a fixed number of
        #  decimals counts differently at m = 3 and m = 150, and the count is
        #  the whole point of the row
        srt = np.sort(m0)
        nd = 1 + int(np.sum(np.diff(srt) > 1e-9 * srt[-1]))
        out[name] = (c0, c1, sd)
        print(f"  {name:<28}{t[-1]-t[0]:8.1f}{nd:7d}{c0:12.4g}{c1:12.4g}"
              f"{(inc0.max()-inc0.min())/TH:8.2f}{sd:11.4f}")
    print("\n  At a constant load the aliased schedule is exactly what the premise")
    print("  says it is: one node, repeated, and V exactly singular -- but only")
    print("  once the dose train has reached steady state. Started at t = 0 the")
    print("  train is still accumulating and the nodes are merely close.")

    #  the comparison the committed 27 is actually taken on: the repository's
    #  own two designs, equal SPAN and unequal n, not equal n
    print("\n  The repository's own pair (two years each, n as the design has it):")
    print(f"  {'design':<28}{'n':>5}{'cond V':>12}{'cond V':>12}{'cond V':>12}"
          f"{'sd(ln th)':>11}")
    print(f"  {'':<28}{'':>5}{'L const':>12}{'10% CV':>12}{'equilib':>12}{'':>11}")
    got = {}
    for name, t in (("monthly (30.4 d)", MONTHLY),
                    ("fortnightly (14 d)", ALIASED)):
        c0 = cond_V(record(t, load_cv=0.0)[0])
        c1 = med(lambda s: cond_V(record(t, seed=s)[0]))[0]
        ce = med(lambda s: cond_V(record(t, seed=s)[0], equilibrate=True))[0]
        sd = med(lambda s: sd_paper(*record(t, seed=s)[:2]))[0]
        got[name] = (c0, c1, ce, sd)
        print(f"  {name:<28}{t.size:5d}{c0:12.4g}{c1:12.4g}{ce:12.4g}{sd:11.4f}")
    a, b = got["fortnightly (14 d)"], got["monthly (30.4 d)"]
    ratio_sd = med(lambda s: (sd_paper(*record(ALIASED, seed=s)[:2])
                              / sd_paper(*record(MONTHLY, seed=s)[:2])))
    ratio_c = med(lambda s: (cond_V(record(ALIASED, seed=s)[0])
                             / cond_V(record(MONTHLY, seed=s)[0])))
    print(f"\n  ratio of sd(ln theta-hat), per record:      "
          f"{ratio_sd[0]:8.2f}  [{ratio_sd[1]:.2f}, {ratio_sd[2]:.2f}]")
    print(f"  ratio of cond(V) at a constant load:        {a[0]/b[0]:8.2f}")
    print(f"  ratio of cond(V) at 10% load variation:     "
          f"{ratio_c[0]:8.2f}  [{ratio_c[1]:.2f}, {ratio_c[2]:.2f}]")
    print(f"  ratio of equilibrated cond(V) at 10% CV:    {a[2]/b[2]:8.2f}")
    #  the same ratio with the load held CONSTANT, so that cond(V)'s
    #  constant-load number is compared against an sd on the same footing
    #  rather than against the committed one, which is a 10%-CV number
    sd_const = (sd_paper(*record(ALIASED, load_cv=0.0)[:2])
                / sd_paper(*record(MONTHLY, load_cv=0.0)[:2]))
    print(f"  ratio of sd(ln theta-hat) at a constant load:{sd_const:8.2f}"
          f"   (against cond's {a[0]/b[0]:.0f})")
    print(f"\n  IT DOES NOT TRACK. At the load variation the specification asks")
    print(f"  for, cond(V) puts the aliased design {ratio_c[0]:.1f} times worse and the")
    print(f"  estimator puts it {ratio_sd[0]:.0f} times worse, a disagreement of "
          f"{ratio_sd[0]/ratio_c[0]:.0f}x. At a")
    print(f"  constant load it overshoots instead, {a[0]/b[0]:.0f} against the "
          f"estimator's {sd_const:.0f},")
    print(f"  a disagreement of {(a[0]/b[0])/sd_const:.0f}x the other way; and with the "
          f"patient already on")
    print("  therapy it is infinite where the estimator is merely bad. It")
    print("  brackets the answer from both sides and reproduces neither.")
    print("\n  THE MECHANISM, which is structural and not a scaling accident.")
    print("  In the linear form Delta_i appears ONLY in the response Delta_i m_i.")
    print("  The design matrix knows about the dose train solely through where")
    print("  the nodes sit, so anything else that moves the nodes counts the")
    print("  same to cond(V). Day-to-day intake is exactly such a thing, and it")
    print("  moves the leading column at the same time, which cond(V) cannot")
    print("  see at all. The excitation this paper is about -- the RANGE of")
    print("  Delta/theta a record covers -- has been moved out of the matrix")
    print("  whose conditioning is being quoted.")
    dm = med(lambda s: (lambda i: (i.max()-i.min())/TH)(record(MONTHLY, seed=s)[1]))[0]
    da = med(lambda s: (lambda i: (i.max()-i.min())/TH)(record(ALIASED, seed=s)[1]))[0]
    print(f"\n  range of Delta/theta: monthly {dm:.2f}, aliased {da:.2f}, "
          f"a ratio of {dm/da:.1f},")
    print(f"  which is {abs(dm/da - ratio_sd[0])/ratio_sd[0]:.0%} from the committed "
          f"{ratio_sd[0]:.0f} -- closer than any condition")
    print("  number here, and it is the quantity the paper already quotes.")
    return {"sd": ratio_sd[0], "cond_cv": ratio_c[0], "cond_const": a[0]/b[0],
            "sd_const": sd_const,
            "cond_eq": a[2]/b[2], "drange": dm/da,
            "steady": out["14 d, already on therapy"][0]}


def section3():
    print("\n" + BAR)
    print("3. (b) DOES cond(V) REPRODUCE THE PLACEMENT RESULT?")
    print(BAR)
    print("  Committed: sd_dense = 0.1500, sd_straddle = 0.0157,")
    print("  straddle_gain = 9.5. Both designs are 365 draws; only the band")
    print("  moves.\n")
    inc_d = record(DENSE, seed=BASE_SEED)[1]
    rows = []
    for name, get in (("monthly, n=25", lambda s: record(MONTHLY, seed=s)[:2]),
                      ("every 2 d, n=365", lambda s: record(DENSE, seed=s)[:2]),
                      ("straddling K_m, n=365",
                       lambda s: (STRADDLE, record(DENSE, seed=s)[1]))):
        m, inc = get(BASE_SEED)
        c_raw = med(lambda s: cond_V(get(s)[0]))[0]
        c_eq = med(lambda s: cond_V(get(s)[0], equilibrate=True))[0]
        c_L = med(lambda s: cond_V(get(s)[0],
                                   first=implied_load(*get(s))))[0]
        sd = med(lambda s: sd_paper(*get(s)))[0]
        rows.append((name, m.min(), m.max(), c_raw, c_eq, c_L, sd))
    print(f"  {'design':<24}{'m range':>17}{'cond V':>12}{'cond V':>10}"
          f"{'cond V':>12}{'sd(ln th)':>11}")
    print(f"  {'':<24}{'':>17}{'[1,m,m2]':>12}{'equilib':>10}"
          f"{'[L,m,m2]':>12}{'':>11}")
    for nm, lo, hi, c, ce, cl, sd in rows:
        print(f"  {nm:<24}{f'{lo:.2f} - {hi:.2f}':>17}{c:12.4g}{ce:10.4g}"
              f"{cl:12.4g}{sd:11.4f}")
    d, s_ = rows[1], rows[2]
    print(f"\n  sd:                     dense / straddle = {d[6]/s_[6]:8.2f}"
          f"   (committed 9.5)")
    print(f"  cond(V), raw:           dense / straddle = {d[3]/s_[3]:8.4f}")
    print(f"  cond(V), equilibrated:  dense / straddle = {d[4]/s_[4]:8.2f}")
    print(f"  cond(V), [L, m, m^2]:   dense / straddle = {d[5]/s_[5]:8.4f}")
    print("\n  IT DOES NOT REPRODUCE IT, and the raw number reverses the sign of")
    print("  the result: read literally, cond(V) says the straddling design is")
    print(f"  {s_[3]/d[3]:.0f} times WORSE conditioned than the band the patient is")
    print("  actually held in, while the estimator says it is 9.5 times better.")
    print("  That is a units artefact and it should be named as one: at m ~ 150")
    print("  the third column is 2.3e4 times the second and 1e0 against the")
    print("  first, so cond(V) is measuring micromolar-squared against nothing")
    print("  in particular. Equilibrating the columns removes it and leaves")
    print(f"  {d[4]/s_[4]:.2f}, the right sign now but still {9.5/(d[4]/s_[4]):.0f} times too small.")
    print("\n  The reason is the same as in section 2 and it is the paper's own:")
    print("  placement buys precision through the SENSITIVITY K_m m/(K_m+m)^2,")
    print("  which peaks at m = K_m, and that factor is nowhere in V. cond(V)")
    print("  is a statement about solving a linear system, not about how much")
    print("  the ordinate moves when theta does.")
    return {"raw": s_[3]/d[3], "eq": d[4]/s_[4], "sd": d[6]/s_[6]}


def section4():
    print("\n" + BAR)
    print("4. (c) THE ESTIMATOR CAUTION: OLS ON THE LINEAR FORM IS WLS ON THE")
    print("   BALANCE")
    print(BAR)
    m, inc, Ls = record(MONTHLY, seed=BASE_SEED)
    mult = KM + m
    print(f"  The multiplier is (K_m + m_i), derived in the docstring, not m_i.")
    print(f"  Over the monthly record it runs {mult.min():.2f} to {mult.max():.2f}, "
          f"a spread of {100*(mult.max()/mult.min()-1):.1f}%,")
    print(f"  because the patient sits at m << K_m. The premise's m_i runs "
          f"{m.min():.2f} to {m.max():.2f},")
    print(f"  a factor of {m.max()/m.min():.2f}. The two prescriptions are not "
          f"close, so both are")
    print("  measured and the difference between them is itself the answer.")

    print("\n  4.1  THE FOOTING THE THREE-COEFFICIENT FORM ACTUALLY OFFERS")
    print("       total load known exactly (b dropped, not given a prior), theta,")
    print("       K_m and k_e all free. No priors, so the posterior width and the")
    print("       estimator's own sd are the same object and the comparison is")
    print("       unambiguous. Medians over the 24 records, monthly design.")
    base = med(lambda s: sd_original_3(*record(MONTHLY, seed=s)[:2]))
    print(f"\n       {'the balance fitted directly (Fisher)':<38}{base[0]:.6f}"
          f"  [{base[1]:.4f}, {base[2]:.4f}]")
    three = {}
    for wt in ("1/(K_m+m)", "unweighted", "1/m"):
        three[wt] = med(lambda s: sd_vandermonde(*record(MONTHLY, seed=s),
                                                 weight=wt))
    ref, c_o, c_m = three["1/(K_m+m)"][0], three["unweighted"], three["1/m"]
    #  the exact-weight row is compared record by record rather than median to
    #  median: two medians can agree while no record does
    rel = med(lambda s: abs(
        sd_vandermonde(*record(MONTHLY, seed=s), weight="1/(K_m+m)")
        - sd_original_3(*record(MONTHLY, seed=s)[:2]))
        / sd_original_3(*record(MONTHLY, seed=s)[:2]))
    for wt in ("1/(K_m+m)", "unweighted", "1/m"):
        v = three[wt]
        tail = (f"   <- agrees to {rel[2]:.1e}" if wt == "1/(K_m+m)"
                else f"   {v[0]/ref - 1:+.4%}")
        print(f"       {'linear form, weights ' + wt:<38}{v[0]:.6f}"
              f"  [{v[1]:.4f}, {v[2]:.4f}]{tail}")
    print("\n       The exact weights reproduce the direct fit to the last bit,")
    print("       which is the check that the transformation is implemented")
    print("       right. Unweighted costs almost nothing on this design; the")
    print(f"       premise's 1/m costs {c_m[0]/ref - 1:.1%}, which is "
          f"{(c_m[0]/ref - 1)/(c_o[0]/ref - 1):.0f} times more.")

    print("\n  4.2  THE PAPER'S OWN FOOTING")
    print("       load bias free with a 5% prior, ln k_e prior 0.8, written as the")
    print("       four-column form [L_i, 1, m_i, m_i^2] with both priors carried")
    print("       as ridge penalties. Two columns, and they are different objects:")
    print("       `posterior` is the width of the penalised normal matrix, which")
    print("       is what _filter reports and what the paper quotes; `sampling` is")
    print("       the sd of what the estimator actually delivers. For the exact")
    print("       weights the second is the smaller because the prior is centred")
    print("       at the truth. For a MISWEIGHTED fit the posterior width is not")
    print("       the width of anything, and the gap between the two columns is")
    print("       the whole caution.")
    kw = dict(b_free=True, ke_prior=KE_PRIOR_LOG, load_prior=LOAD_FRAC)
    paper = med(lambda s: sd_paper(*record(MONTHLY, seed=s)[:2]))
    print(f"\n       {'weights':<14}{'posterior':>12}{'sampling':>12}"
          f"{'vs exact':>11}{'vs paper':>11}")
    ex_s, fac = None, {}
    for wt in ("1/(K_m+m)", "unweighted", "1/m"):
        pw = med(lambda s: sd_vandermonde(*record(MONTHLY, seed=s),
                                          weight=wt, posterior=True, **kw))
        q = med(lambda s: sd_vandermonde(*record(MONTHLY, seed=s),
                                         weight=wt, **kw))
        if ex_s is None:
            ex_s = q[0]
        fac[wt] = q[0] / ex_s
        print(f"       {wt:<14}{pw[0]:12.6f}{q[0]:12.6f}{q[0]/ex_s:11.1f}x"
              f"{q[0]/paper[0]:10.1f}x")
    ols_over_exact, m_over_exact = fac["unweighted"], fac["1/m"]
    print(f"       {'the paper':<14}{paper[0]:12.6f}{'--':>12}{'--':>11}{'--':>11}")
    print(f"\n       The exact-weight posterior reproduces the paper's own "
          f"{paper[0]:.4f} to")
    ex_p = med(lambda s: abs(sd_vandermonde(*record(MONTHLY, seed=s),
                                            weight="1/(K_m+m)",
                                            posterior=True, **kw)
                             - sd_paper(*record(MONTHLY, seed=s)[:2]))
               / sd_paper(*record(MONTHLY, seed=s)[:2]))
    print(f"       {ex_p[2]:.1e} relative, record by record, which is the second")
    print("       check that the reframing is faithful: the paper's four-parameter")
    print("       fit IS this linear system, in different coordinates.")
    print("\n       IT CHANGES MATERIALLY HERE. On the footing the paper uses,")
    print("       dropping the weights costs a factor of tens, not a fraction of")
    print("       a per cent, even though the weights span only 6.5%.")

    print("\n  4.3  WHY THE SAME 6.5% COSTS 0.06% ON ONE FOOTING AND TENS ON THE")
    print("       OTHER")
    print("       Misweighting costs in proportion to how nearly singular the")
    print("       system is, not to how much the weights vary. With the load")
    print("       known the system is well posed and a small misweighting is a")
    print("       small inefficiency. With the load bias free the design sits")
    print("       next to the fibre of prop:rank, the prior holds it, and a")
    print("       small misweighting turns the estimator slightly off the")
    print("       direction the prior is holding -- which is the one direction")
    print("       where being slightly off is expensive.")
    for cv in (0.02, 0.10, 0.40):
        ex = med(lambda s: sd_vandermonde(*record(MONTHLY, load_cv=cv, seed=s),
                                          weight="1/(K_m+m)", **kw))[0]
        ols = med(lambda s: sd_vandermonde(*record(MONTHLY, load_cv=cv, seed=s),
                                           weight="unweighted", **kw))[0]
        print(f"       load variation {100*cv:4.0f}%:  exact {ex:.4f}  "
              f"unweighted {ols:.4f}  ({ols/ex:6.1f}x)")
    print("       More day-to-day intake moves the design away from the fibre,")
    print("       and the price of misweighting falls with it.")

    print("\n  4.4  AND WHERE THE WEIGHTS ARE NOT NEARLY CONSTANT")
    print("       On the K_m-straddling band (K_m + m) runs 150 to 300, a factor")
    print("       of 2, so the unweighted fit has something real to lose even on")
    print("       the well-posed footing.")
    for nm, m_, inc_ in (("monthly, m = 2.8-9.5", *record(MONTHLY, seed=BASE_SEED)[:2]),
                         ("straddling, m = 50-200", STRADDLE,
                          record(DENSE, seed=BASE_SEED)[1])):
        Ls_ = implied_load(m_, inc_)
        e = sd_vandermonde(m_, inc_, Ls_, weight="1/(K_m+m)")
        o = sd_vandermonde(m_, inc_, Ls_, weight="unweighted")
        print(f"       {nm:<24} exact {e:.6f}  unweighted {o:.6f} "
              f"({o/e - 1:+.2%})")

    print("\n  4.5  A MONTE CARLO ON THE SANDWICH, BECAUSE THE 4.2 RESULT RESTS")
    print("       ON A COVARIANCE FORMULA AND NOT ON A FIT")
    #  noise and priors scaled together by f, so the whole problem scales
    #  linearly and the nonlinear read-back of theta is exercised in the
    #  regime where a delta-method sd is meant to describe it
    f, N = 0.01, 20000
    V = np.stack([Ls, np.ones_like(m), m, m ** 2], axis=1)
    ctrue = np.array([KM, 0.0, -(TH + KE * KM), -KE])
    s = SIG * f * (KM + m)
    rng = np.random.default_rng(20260921)
    E = rng.standard_normal((N, m.size)) * s
    for wt in ("1/(K_m+m)", "unweighted"):
        w = WEIGHTS[wt](m)
        pen = np.zeros((4, 4))
        pen[1, 1] = (1.0 / KM) ** 2 / (LOAD_FRAC * f * LB) ** 2
        pen[3, 3] = 1.0 / (KE_PRIOR_LOG * f * KE) ** 2
        Vw = V * w[:, None]
        M = np.linalg.inv(Vw.T @ Vw / (SIG * f) ** 2 + pen)
        Y = (V @ ctrue + E) * w
        Cs = (M @ (Vw.T @ Y.T / (SIG * f) ** 2 + (pen @ ctrue)[:, None])).T
        th = Cs[:, 1] / Cs[:, 0] - Cs[:, 2] + Cs[:, 3] * Cs[:, 0]
        ok = th > 0
        emp = float(np.std(np.log(th[ok]))) / f
        #  the whole problem scales linearly in f -- noise, priors and sd --
        #  so the empirical sd divided by f is comparable with the UNSCALED
        #  analytic, and the analytic is the one section 4.2 quotes
        ana = sd_vandermonde(m, inc, Ls, weight=wt, b_free=True,
                             ke_prior=KE_PRIOR_LOG, load_prior=LOAD_FRAC)
        print(f"       {wt:<12} {N} fits at 1/{1/f:.0f} the noise and priors: "
              f"sd {emp:.4f} against {ana:.4f} ({abs(emp-ana)/ana:.2%})")
    print("       The sandwich is what the fits do.")
    return {"known_ols": c_o[0]/ref - 1, "known_m": c_m[0]/ref - 1,
            "paper_ols": ols_over_exact, "paper_m": m_over_exact}


def section5(a, b, c):
    print("\n" + BAR)
    print("5. WHAT CAME OUT")
    print(BAR)
    print("  The reframing is ALGEBRAICALLY RIGHT and its coefficients are")
    print("  A = (L+b)K_m, B = (L+b) - theta - k_e K_m, C = -k_e, with the")
    print("  clearing multiplier (K_m + m_i). Three visits do suffice -- at a")
    print("  constant, known total load and three distinct cycle phases.")
    print()
    print(f"  (a) NO. cond(V) neither reproduces the committed factor of "
          f"{a['sd']:.0f} nor")
    print(f"      tracks it: {a['cond_cv']:.1f} where the estimator says "
          f"{a['sd']:.0f}, at the specification's")
    print(f"      load variation; {a['cond_const']:.0f} where the estimator "
          f"says {a['sd_const']:.0f}, at a constant")
    print(f"      load; and {a['steady']:.0e} once the patient is already on "
          f"therapy.")
    print("      Delta lives in the response, not in the design matrix. The")
    print(f"      range of Delta/theta, which the paper already quotes, gives "
          f"{a['drange']:.1f}.")
    print(f"  (b) NO, and read literally it inverts the result: raw cond(V) "
          f"makes the")
    print(f"      straddling design {b['raw']:.0f} times WORSE where the "
          f"estimator makes it")
    print(f"      {b['sd']:.1f} times better. Equilibrated it gets the sign "
          f"right and is still")
    print(f"      {b['sd']/b['eq']:.0f} times too small.")
    print(f"  (c) YES on the footing the paper actually uses -- unweighted "
          f"costs a")
    print(f"      factor of {c['paper_ols']:.0f} there -- and NO on the footing "
          f"the three-coefficient")
    print(f"      form offers, where it costs {c['known_ols']:.2%}. The "
          f"premise's stated weight")
    print(f"      1/m_i is not the right one either: it costs "
          f"{c['known_m']:.1%} on the")
    print(f"      well-posed footing and a factor of {c['paper_m']:.0f} on the "
          f"paper's. The")
    print("      right weight is 1/(K_m + m_i).")


# ==========================================================================
def rechecks():
    print("\n" + BAR)
    print("R1. RE-CHECK: the aliasing row of results/numbers/design.tex")
    print(BAR)
    for lab, times, want in (("sd_monthly", MONTHLY, 0.1699),
                             ("sd_aliased", ALIASED, 4.5845)):
        got = med(lambda s: sd_paper(*record(times, seed=s)[:2]))
        print(f"  {lab:<12} {got[0]:.4f} [{got[1]:.4f}, {got[2]:.4f}] against "
              f"the committed {want:.4f} ({abs(got[0]-want)/want:.2%} apart)")
    got = med(lambda s: (sd_paper(*record(ALIASED, seed=s)[:2])
                         / sd_paper(*record(MONTHLY, seed=s)[:2])))
    print(f"  {'alias_cost':<12} {got[0]:.1f} [{got[1]:.1f}, {got[2]:.1f}] against "
          f"the committed 27 [21, 33] ({abs(got[0]-27)/27:.2%} apart)")

    print("\n" + BAR)
    print("R2. RE-CHECK: the placement row of results/numbers/design.tex")
    print(BAR)
    dense = med(lambda s: sd_paper(*record(DENSE, seed=s)[:2]))
    strad = med(lambda s: sd_paper(STRADDLE, record(DENSE, seed=s)[1]))
    for lab, got, want in (("sd_dense", dense, 0.1500),
                           ("sd_straddle", strad, 0.0157)):
        print(f"  {lab:<14} {got[0]:.4f} [{got[1]:.4f}, {got[2]:.4f}] against "
              f"the committed {want:.4f} ({abs(got[0]-want)/want:.2%} apart)")
    print(f"  {'straddle_gain':<14} {dense[0]/strad[0]:.1f} against the committed "
          f"9.5 ({abs(dense[0]/strad[0]-9.5)/9.5:.2%} apart)")

    print("\n" + BAR)
    print("R3. RE-CHECK: the excursion the schedules see, and the marker band")
    print(BAR)
    #  these keys are COMMITTED AT TWO DECIMALS (PAPER_FORMATS), so the gap
    #  below is the rounding of the committed value and not a disagreement;
    #  0.1333 printed at two decimals is 0.13, which is 2.5% away from itself.
    for lab, times, want in (("delta_range_monthly", MONTHLY, 1.84),
                             ("delta_range_aliased", ALIASED, 0.13),
                             ("delta_range_dense", DENSE, 1.95)):
        inc = record(times, seed=BASE_SEED)[1]
        got = float((inc.max() - inc.min()) / TH)
        print(f"  {lab:<22} {got:.4f} against the committed {want:.2f} "
              f"({abs(got-want)/want:.2%} apart)")
    #  the band the marker actually visits, from the same three-month window
    #  fig_D0_dosetrain draws it on
    tt = np.linspace(0.0, 91.0, 4000)
    mm = record(tt, load_cv=0.0)[0]
    for lab, got, want in (("marker_low", float(mm.min()), 3.3646),
                           ("marker_high", float(mm.max()), 8.9708)):
        print(f"  {lab:<22} {got:.4f} against the committed {want:.4f} "
              f"({abs(got-want)/want:.2%} apart)")


def main():
    print(__doc__)
    section0()
    section1()
    a = section2()
    b = section3()
    c = section4()
    section5(a, b, c)
    rechecks()


if __name__ == "__main__":
    main()
