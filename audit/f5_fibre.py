"""
P5 -- is the rank-one degeneracy "recovered through curvature", or is it an
exact fibre?

The paper (Section 6, and rem:second_order) says three things about the null
direction of prop:rank:

    (i)   the design is flat along it only to FIRST order;
    (ii)  the noiseless misfit therefore grows as eps^2, chi^2 as eps^4, so the
          combination is "recovered through curvature rather than slope";
    (iii) the recoverable displacement therefore scales as sqrt(noise).

A second reader proposes a different object. Multiplying the quasi-steady
balance through by (K_m + m) makes it linear in

    beta = (K_m, K_m b, b - theta - k_e K_m, -k_e)

against the regressors (L_k, 1, m_k, m_k^2).  At a CONSTANT recorded load the
first two columns collapse onto each other, so only the combination K_m (L + b)
survives, and the claim is that the resulting one-parameter family

    K_m(b) = K_0 L/(L + b),   theta(b) = theta_0 + b + k_e (K_0 - K_m(b)),  k_e fixed

fits the noiseless record EXACTLY at every b -- not to second order, exactly.
If that holds, (i)-(iii) are all statements about a straight line drawn through
a curved exact fibre, and the eps^4 is the fibre's curvature rather than any
information the record contains.

This script decides between the two readings. Nothing below records an expected
result in advance; the only constants quoted from outside are the four numbers
already printed in the paper, and they are labelled where they are used.

Run:  python audit/f5_fibre.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

from common.models import ParamsA, g, clearance_A
from design.run import P, delta, load, LOAD_FRAC, ASSAY_CV, DOSE_INTERVAL
import model_A.run as LOC

TH0, KM0, KE0 = P.theta, P.Km, P.ke
LB = load(0.0)                                   # the recorded load, constant
SIG = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)       # ordinate sd, load units

# numbers the paper prints, quoted here only so the script can say whether it
# reproduces them.  paper/sections/06_excitation.tex:502-513
PAPER_PRIOR_HALFWIDTH = 31.6
PAPER_PROFILE_HALFWIDTH = 36.9
PAPER_CHI2_OVER_EPS4 = 2.08e-6

MONTHLY = np.arange(0.0, 730.0, 30.4)
TWO_DAILY = np.arange(0.0, 730.0, 2.0)


def banner(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


# --------------------------------------------------------------------------
def levels(times, th=TH0, km=KM0, ke=KE0, b=0.0, Ls=None):
    """Quasi-steady marker at each visit: the positive root of the balance."""
    inc = np.asarray(delta(times), dtype=float)
    tot = (np.full(np.shape(times), LB) if Ls is None else np.asarray(Ls, float)) + b
    B = ke * km + th + inc - tot
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * tot * km)) / (2.0 * ke), inc


def load_residual(m, inc, th, km, ke, b, Ls=None):
    """Balance residual on the load scale, which is the fit's own ordinate."""
    tot = (np.full(np.shape(m), LB) if Ls is None else np.asarray(Ls, float)) + b
    return (th + inc) * g(m, km) + ke * m - tot


def chi2(m, inc, th, km, ke, b, Ls=None, sigma=SIG):
    r = load_residual(m, inc, th, km, ke, b, Ls)
    return float(np.sum(r ** 2) / sigma ** 2)


def fibre(b, th0=TH0, km0=KM0, ke=KE0, L=LB):
    """The candidate exact fibre, parametrised by the load bias b."""
    km = km0 * L / (L + b)
    th = th0 + b + ke * (km0 - km)
    return th, km, ke, b


# --------------------------------------------------------------------------
def s1_symbolic():
    banner("S1  symbolic: what a constant-load record identifies")

    th, km, ke, b, L, m, dl = sp.symbols(
        "theta K_m k_e b L m Delta", positive=True)

    bal = sp.Eq((th + dl) * m / (km + m) + ke * m, L + b)
    dl_of_m = sp.simplify(sp.solve(bal, dl)[0])
    print("  Solving the balance for the increment that produces level m:")
    print(f"    Delta(m) = {sp.simplify(sp.expand(dl_of_m))}")

    # collect it against the basis 1/m, 1, m
    expr = sp.expand(sp.simplify(dl_of_m))
    A = sp.simplify(sp.limit(expr * m, m, 0))
    C = sp.simplify(sp.limit(expr / m, m, sp.oo))
    B = sp.simplify(sp.expand(expr - A / m - C * m))
    print()
    print("  In the basis {1/m, 1, m} the increment-to-level map has coefficients")
    print(f"    A (on 1/m) = {sp.factor(A)}")
    print(f"    B (on 1)   = {sp.factor(B)}")
    print(f"    C (on m)   = {sp.factor(C)}")
    resid = sp.simplify(expr - (A / m + B + C * m))
    print(f"    remainder  = {resid}")
    assert resid == 0, "the three-term form is not exact"
    print()
    print("  Those three numbers are EVERYTHING a record of (m_k, Delta_k) pairs")
    print("  at a constant load can contain: {1/m, 1, m} are independent as")
    print("  functions of m, so two parameter vectors are observationally")
    print("  identical iff they agree on (A, B, C).")

    # the fibre: fix (A, B, C) and solve for (theta, K_m) as functions of b
    th0, km0 = sp.symbols("theta_0 K_0", positive=True)
    A0 = A.subs({th: th0, km: km0, b: 0})
    B0 = B.subs({th: th0, km: km0, b: 0})
    sol = sp.solve([sp.Eq(A, A0), sp.Eq(B, B0)], [th, km], dict=True)
    print()
    print("  Fixing (A, B, C) at the truth (theta_0, K_0, k_e, b = 0) and solving")
    print("  for (theta, K_m) as functions of b:")
    for s in sol:
        km_b = sp.simplify(s[km])
        th_b = sp.simplify(sp.expand(s[th]))
        print(f"    K_m(b)   = {km_b}")
        print(f"    theta(b) = {sp.collect(th_b, ke)}")
        print(f"    k_e      = unconstrained by A and B, so k_e(b) = k_e  (C fixes it)")
    print()
    print("  That is a ONE-PARAMETER family of exact observational twins, which")
    print("  is the rank-one deficiency prop:rank reports -- but as an exact")
    print("  curve, not as a first-order flat direction.")

    # the tangent at b = 0, against prop:rank's null vector
    print()
    print("  Tangent of the fibre at b = 0, against prop:rank's n = (L + k_e K_m, -K_m, 0, L):")
    km_b, th_b = sp.simplify(sol[0][km]), sp.simplify(sol[0][th])
    tang = [sp.simplify(sp.diff(th_b, b).subs(b, 0)),
            sp.simplify(sp.diff(km_b, b).subs(b, 0)), sp.Integer(0), sp.Integer(1)]
    print(f"    d/db (theta, K_m, k_e, b) at b=0 = "
          f"({tang[0]}, {tang[1]}, {tang[2]}, {tang[3]})")
    n = [(L + ke * km0) / L, -km0 / L, sp.Integer(0), sp.Integer(1)]
    diff = [sp.simplify(a - c) for a, c in zip(tang, n)]
    print(f"    n normalised to unit b-component     = "
          f"({sp.expand(n[0])}, {n[1]}, {n[2]}, {n[3]})")
    print(f"    difference = {diff}")
    assert all(d == 0 for d in diff), "tangent does not match prop:rank"
    print("  => prop:rank's null direction IS the tangent to this fibre at b = 0.")
    return sol[0]


# --------------------------------------------------------------------------
def s2_exact_or_second_order(times=MONTHLY, tag="monthly x 2 y"):
    banner(f"S2  misfit along the straight line vs along the curve ({tag})")

    m, inc = levels(times)
    print(f"  n = {len(times)} visits, marker span "
          f"[{m.min():.4f}, {m.max():.4f}], Delta/theta span "
          f"[{inc.min()/TH0:.4f}, {inc.max()/TH0:.4f}]")
    print(f"  ordinate sd sigma = {SIG:.4f} (load units), recorded load L = {LB:.4f}")

    # the straight null direction, normalised to unit b-component
    nvec = np.array([(LB + KE0 * KM0) / LB, -KM0 / LB, 0.0, 1.0])
    print(f"\n  straight direction n (unit b-component) = "
          f"{np.array2string(nvec, precision=6)}")

    print(f"\n  {'b displacement':>15} | {'rms resid, LINE':>17} {'chi2, LINE':>13}"
          f" | {'rms resid, CURVE':>18} {'chi2, CURVE':>13}")
    print("  " + "-" * 88)
    for eps in (1e-2, 1e-1, 1.0, 3.0, 10.0, 30.0, 36.9, 100.0, 200.0):
        pl = np.array([TH0, KM0, KE0, 0.0]) + eps * nvec
        rl = load_residual(m, inc, *pl)
        pc = fibre(eps)
        rc = load_residual(m, inc, *pc)
        print(f"  {eps:15.4f} | {np.sqrt(np.mean(rl**2)):17.6e} "
              f"{np.sum(rl**2)/SIG**2:13.6e} | {np.sqrt(np.mean(rc**2)):18.6e} "
              f"{np.sum(rc**2)/SIG**2:13.6e}")

    # the power law on the line
    print("\n  power of the LINE's chi^2 in the displacement, and chi^2/eps^4:")
    print(f"    {'eps':>10} {'chi2':>16} {'local slope d ln chi2/d ln eps':>32}"
          f" {'chi2/eps^4':>14}")
    grid = np.geomspace(1e-3, 1e2, 21)
    c2 = np.array([chi2(m, inc, *(np.array([TH0, KM0, KE0, 0.0]) + e * nvec))
                   for e in grid])
    for i, (e, c) in enumerate(zip(grid, c2)):
        sl = (np.nan if i == 0 else
              np.log(c2[i] / c2[i - 1]) / np.log(grid[i] / grid[i - 1]))
        print(f"    {e:10.4g} {c:16.6e} {sl:32.4f} {c/e**4:14.6e}")

    print(f"\n  the paper prints chi^2/eps^4 = {PAPER_CHI2_OVER_EPS4:.3g} "
          f"(06_excitation.tex:513)")
    print(f"  Delta chi^2 = 3.84 on the LINE at eps = "
          f"{(3.84 / np.median(c2 / grid**4))**0.25:.4f}")
    print(f"  the paper prints +/- {PAPER_PROFILE_HALFWIDTH} for that quantity "
          f"(06_excitation.tex:503)")
    print(f"  the load prior's own 95% half-width, 1.96 * {LOAD_FRAC} * L = "
          f"{1.96 * LOAD_FRAC * LB:.4f}; the paper prints "
          f"+/- {PAPER_PRIOR_HALFWIDTH} (06_excitation.tex:502)")
    return m, inc, nvec


# --------------------------------------------------------------------------
def s3_true_profile(times=MONTHLY):
    banner("S3  a TRUE profile in b: re-optimise the other three at each b")

    m, inc = levels(times)
    nvec = np.array([(LB + KE0 * KM0) / LB, -KM0 / LB, 0.0, 1.0])

    def profile(b, x0):
        """min over (theta, K_m, k_e) of chi^2 at this fixed b."""
        def res(x):
            th, lk, ke = x[0], x[1], x[2]
            return load_residual(m, inc, th, np.exp(lk), ke, b) / SIG
        s = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        return float(np.sum(s.fun ** 2)), s.x

    print(f"  {'b':>10} | {'chi2 along LINE':>17} {'chi2 TRUE profile':>19}"
          f" | {'theta_hat':>12} {'K_m_hat':>10} {'k_e_hat':>9}")
    print("  " + "-" * 88)
    x0 = np.array([TH0, np.log(KM0), KE0])
    for b in (0.0, 1.0, 10.0, 30.0, 36.9, 100.0, 200.0, 300.0, -100.0, -200.0):
        pl = np.array([TH0, KM0, KE0, 0.0]) + b * nvec
        cl = chi2(m, inc, *pl)
        cp, xh = profile(b, x0)
        print(f"  {b:10.2f} | {cl:17.6e} {cp:19.6e} | {xh[0]:12.4f} "
              f"{np.exp(xh[1]):10.4f} {xh[2]:9.5f}")

    print("\n  and against the closed-form fibre point at the same b:")
    print(f"  {'b':>10} | {'theta fibre':>13} {'K_m fibre':>11} {'k_e fibre':>10}")
    print("  " + "-" * 54)
    for b in (30.0, 100.0, 200.0, -100.0):
        th, km, ke, _ = fibre(b)
        print(f"  {b:10.2f} | {th:13.4f} {km:11.4f} {ke:10.5f}")

    # the fibre is unbounded above; what bounds it is admissibility in b
    print("\n  The fibre runs to b -> +infinity with theta, K_m and k_e all")
    print("  positive, so nothing in the model closes it. What closes it is a")
    print("  bound on the load bias. theta implied by each plausible bound:")
    # theta(b) is increasing in b and runs to -infinity as b -> -L, so the
    # admissible end is its root. brentq rather than a hand-rolled bisection:
    # an earlier version of this line inverted the bracket update and printed
    # -215.1515, which is where K_m reaches 300 and where theta is still 1785.
    lo = float(brentq(lambda x: fibre(x)[0], -LB + 1e-9, 0.0, xtol=1e-12))
    print(f"    theta(b) reaches zero at b = {lo:.4f} "
          f"({100*lo/LB:+.2f}% of the recorded load)")
    print(f"\n  {'bound on |b|':>28} {'theta interval':>34} {'fold':>9}")
    print("  " + "-" * 74)
    for frac, lab in ((0.05, "the paper's load prior, 5%"),
                      (0.10, "10% of the recorded load"),
                      (0.20, "20% of the recorded load"),
                      (0.50, "50% of the recorded load")):
        blo, bhi = max(-frac * LB, lo), frac * LB
        t1, t2 = fibre(blo)[0], fibre(bhi)[0]
        print(f"  {lab:>28} [{min(t1,t2):14.4f}, {max(t1,t2):14.4f}] "
              f"{max(t1,t2)/max(min(t1,t2), 1e-12):9.3f}")
    print(f"\n  and with b >= 0 imposed (unrecorded load cannot be negative):")
    for frac in (0.05, 0.20, 0.50):
        t1, t2 = fibre(0.0)[0], fibre(frac * LB)[0]
        print(f"    0 <= b <= {frac:.0%} of L:  theta in "
              f"[{t1:.4f}, {t2:.4f}], fold {t2/t1:.4f}")
    return m, inc


# --------------------------------------------------------------------------
def s3b_which_design_gave_36_9():
    banner("S3b  which construction produced the paper's +/-36.9 and 2.08e-6?")

    print("  chi^2 along a straight line is a SECOND-ORDER quantity, so it is")
    print("  parametrisation dependent: the same null direction walked in")
    print("  natural (theta, K_m) and in log (ln theta, ln K_m) coordinates")
    print("  gives different answers. prop:rank is stated in log parameters")
    print("  (06_excitation.tex:332, 'in (ln theta, ln K_m, L_0)') and the")
    print("  information filter's columns are derivatives in those same")
    print("  coordinates (src/design/run.py:324-329), so the LOG column is the")
    print("  paper's own construction and the natural one is not.\n")
    c_th = (1.0 + KM0 * KE0 / LB) / TH0        # eq:rank
    c_km = -1.0 / LB
    print(f"  c_theta = {c_th:.6e}, c_Km = {c_km:.6e} (eq:rank, log parameters)")
    print(f"\n  {'design':>22} {'n':>5} | {'NATURAL chi2/eps^4':>19}"
          f" {'half-width':>11} | {'LOG chi2/eps^4':>15} {'half-width':>11}")
    print("  " + "-" * 94)
    designs = [("two-daily x 2 y", np.arange(0.0, 730.0, 2.0)),
               ("weekly x 2 y", np.arange(0.0, 730.0, 7.0)),
               ("monthly x 2 y", np.arange(0.0, 730.0, 30.4)),
               ("two-monthly x 2 y", np.arange(0.0, 730.0, 60.0)),
               ("daily x 2 y", np.arange(0.0, 730.0, 1.0)),
               ("two-daily x 1 y", np.arange(0.0, 365.0, 2.0)),
               ("weekly x 1 y", np.arange(0.0, 365.0, 7.0)),
               ("14.5 d x 2 y", np.arange(0.0, 730.0, 14.5))]
    for lab, t in designs:
        m, inc = levels(t)
        row = []
        for mode in ("natural", "log"):
            e = 1.0
            th = TH0 * (1 + e * c_th) if mode == "natural" else TH0 * np.exp(e * c_th)
            km = KM0 * (1 + e * c_km) if mode == "natural" else KM0 * np.exp(e * c_km)
            k = chi2(m, inc, th, km, KE0, e)
            # the quartic extrapolation, and the exact crossing, which differ
            # because chi^2 drifts off the quartic by the time it reaches 3.84
            f = lambda x, mo=mode: chi2(
                m, inc,
                TH0 * (1 + x * c_th) if mo == "natural" else TH0 * np.exp(x * c_th),
                KM0 * (1 + x * c_km) if mo == "natural" else KM0 * np.exp(x * c_km),
                KE0, x) - 3.84
            row += [k, float(brentq(f, 1e-3, 1e4))]
        print(f"  {lab:>22} {len(t):5d} | {row[0]:19.6e} {row[1]:11.4f}"
              f" | {row[2]:15.6e} {row[3]:11.4f}")

    print(f"\n  the paper prints chi^2/eps^4 = {PAPER_CHI2_OVER_EPS4:.3g} and "
          f"+/-{PAPER_PROFILE_HALFWIDTH} (06_excitation.tex:502-513).")
    print(f"  Note the printed pair is internally a QUARTIC EXTRAPOLATION:")
    print(f"    (3.84/{PAPER_CHI2_OVER_EPS4:.3g})^(1/4) = "
          f"{(3.84/PAPER_CHI2_OVER_EPS4)**0.25:.3f} against the printed "
          f"{PAPER_PROFILE_HALFWIDTH}.")
    print("  A true profile does not produce a constant chi^2/eps^4 to read a")
    print("  half-width off, so that consistency is itself evidence the number")
    print("  was taken along a ray.")
    print(f"\n  the load prior's own 95% half-width, 1.96 * {LOAD_FRAC} * L = "
          f"{1.96 * LOAD_FRAC * LB:.4f}; the paper prints "
          f"+/-{PAPER_PRIOR_HALFWIDTH} (06_excitation.tex:502)")
    print("\n  The half-width is strongly design dependent in either")
    print("  parametrisation, so no single number defends 'every covariance")
    print("  this note reports', which is what rem:second_order claims.")


# --------------------------------------------------------------------------
def s3c_how_well_is_k_e_pinned():
    banner("S3c  k_e is EXACTLY identified and still barely estimable")

    print("  C = -k_e is an exact invariant of the record, so k_e is identified")
    print("  with no fibre at all. That is a statement about rank, not about")
    print("  precision. The three regressors {1/m, 1, m} are evaluated over the")
    print("  marker band the design actually visits:\n")
    for lab, t in (("monthly x 2 y", MONTHLY), ("two-daily x 2 y", TWO_DAILY)):
        m, inc = levels(t)
        X = np.column_stack([1.0 / m, np.ones_like(m), m])
        Xn = X / np.linalg.norm(X, axis=0)
        print(f"    {lab:>16}: marker band [{m.min():.4f}, {m.max():.4f}], "
              f"fold {m.max()/m.min():.3f}, cond(X) = {np.linalg.cond(Xn):.4e}")
    print("\n  So 'exactly identified' and 'usefully estimated' are as far apart")
    print("  here as prop:rank's own rem:rankvsprec says they are, and this")
    print("  time the gap runs the other way: an exact invariant that the band")
    print("  is too narrow to resolve.")


# --------------------------------------------------------------------------
def s4_with_noise(times=MONTHLY, reps=200, seeds=range(8)):
    banner("S4  with noise: does the likelihood pick a point on the fibre?")

    m, inc = levels(times)
    print(f"  {reps} noisy replicates at each of {len(list(seeds))} seeds, ordinate")
    print(f"  noise sd {SIG:.4f}, monthly design, no priors, constant load.")
    print("  Every summary below is quoted with its spread ACROSS SEEDS, because")
    print("  an earlier version of this section quoted one seed and two of its")
    print("  numbers (a 23% negative-theta rate and a k_e median of -0.879) are")
    print("  seed artefacts: the sign of the k_e median is close to a coin toss.\n")

    def fit(y, b_free=True):
        def res(x):
            th, lk, ke = x[0], np.exp(x[1]), x[2]
            b = x[3] if b_free else 0.0
            return ((th + inc) * g(m, lk) + ke * m - (LB + b) - y) / SIG
        s = least_squares(res, np.array([TH0, np.log(KM0), KE0, 0.0]),
                          xtol=1e-14, ftol=1e-14)
        return s.x, float(np.sum(s.fun ** 2))

    per_seed = []
    for sd in seeds:
        rng = np.random.default_rng(sd)
        bs, ths, kms, kes = [], [], [], []
        for _ in range(reps):
            x, _ = fit(rng.normal(0.0, SIG, m.size))
            bs.append(x[3]); ths.append(x[0])
            kms.append(np.exp(x[1])); kes.append(x[2])
        per_seed.append(tuple(map(np.asarray, (bs, ths, kms, kes))))

    def spread(fn):
        v = np.array([fn(*p) for p in per_seed], float)
        return v.min(), np.median(v), v.max()

    print(f"    {'quantity':>28} {'truth':>10} | {'min':>12} {'median':>12} "
          f"{'max':>12}   (over seeds)")
    print("  " + "-" * 92)
    for nm, i, tr in (("median b_hat", 0, 0.0), ("median theta_hat", 1, TH0),
                      ("median K_m_hat", 2, KM0), ("median k_e_hat", 3, KE0)):
        lo, md, hi = spread(lambda *p, i=i: np.median(p[i]))
        print(f"    {nm:>28} {tr:10.4f} | {lo:12.4f} {md:12.4f} {hi:12.4f}")
    lo, md, hi = spread(lambda *p: float(np.mean(p[3] < 0)))
    print(f"    {'fraction with k_e_hat < 0':>28} {'':>10} | "
          f"{lo:12.4f} {md:12.4f} {hi:12.4f}")
    lo, md, hi = spread(lambda *p: float(np.mean(p[1] <= 0)))
    print(f"    {'fraction with theta_hat <= 0':>28} {'':>10} | "
          f"{lo:12.4f} {md:12.4f} {hi:12.4f}")
    print("\n    The sign of the fitted k_e is close to a coin toss, so no single")
    print("    seed's median may be quoted. What is stable across every seed is")
    print("    the collapse itself, and the invariants below.")
    bs, ths, kms, kes = per_seed[0]

    # the fibre invariants are what the data actually determine
    print("\n  The invariants (A, B, C) = (K_m(L+b), L+b-theta-k_e K_m, -k_e)")
    print("  are what a constant-load record determines. How well are those")
    print("  estimated, against how badly b, theta and K_m are?")
    A_t, B_t, C_t = KM0 * LB, LB - TH0 - KE0 * KM0, -KE0
    A_h = kms * (LB + bs)
    B_h = (LB + bs) - ths - kes * kms
    C_h = -kes
    for nm, v, tr in (("A", A_h, A_t), ("B", B_h, B_t), ("C", C_h, C_t)):
        print(f"    {nm}: truth {tr:14.4f}   median {np.median(v):14.4f}   "
              f"relative IQR {np.subtract(*np.percentile(v,[75,25]))/abs(tr):.4f}")
    for nm, v, tr in (("b", bs, LB), ("theta", ths, TH0), ("K_m", kms, KM0)):
        print(f"    {nm:>5}: relative IQR "
              f"{np.subtract(*np.percentile(v,[75,25]))/abs(tr):.4f}")

    # b fixed at 0 costs how much?
    y = rng.normal(0.0, SIG, m.size)
    _, c_free = fit(y, True)
    _, c_fix = fit(y, False)
    print(f"\n    one replicate: chi2 with b free {c_free:.4f}, b fixed at 0 "
          f"{c_fix:.4f}, difference {c_fix - c_free:.4e}")
    return bs


# --------------------------------------------------------------------------
def s5_varying_load(times=MONTHLY, cvs=(0.0, 0.005, 0.02, 0.10), seed=3):
    banner("S5  the fibre under a VARYING load: does it close, and how fast?")

    print("  The fibre exists because L_k and 1 are the same column when L is")
    print("  constant. With L varying they separate. Below: the exact-fibre")
    print("  residual at b = 30 (which is zero at CV 0) as the CV grows, and a")
    print("  true profile half-width in b at each CV.\n")
    rng = np.random.default_rng(seed)
    print(f"  {'load CV':>9} {'rms resid on fibre(b=30)':>26} {'chi2':>13} "
          f"{'95% profile half-width in b':>29}")
    print("  " + "-" * 82)
    for cv in cvs:
        Ls = (LB * np.exp(rng.normal(-0.5 * np.log1p(cv ** 2),
                                     np.sqrt(np.log1p(cv ** 2)), times.size))
              if cv else np.full(times.size, LB))
        m, inc = levels(times, Ls=Ls)
        th, km, ke, b = fibre(30.0)
        r = load_residual(m, inc, th, km, ke, b, Ls=Ls)
        c = float(np.sum(r ** 2) / SIG ** 2)

        def prof(bv):
            """Warm-started at the closed-form fibre point, and bounded, so a
            flat direction returns a flat chi^2 rather than an optimiser
            failure."""
            th_s, km_s, _, _ = fibre(bv)
            th_s = max(th_s, 1e-3)

            def res(x):
                return load_residual(m, inc, np.exp(x[0]), np.exp(x[1]),
                                     np.exp(x[2]), bv, Ls=Ls) / SIG
            best = np.inf
            for start in (np.array([np.log(th_s), np.log(km_s), np.log(KE0)]),
                          np.array([np.log(TH0), np.log(KM0), np.log(KE0)])):
                s = least_squares(res, start, xtol=1e-15, ftol=1e-15,
                                  gtol=1e-15,
                                  bounds=(np.array([-20., -20., -20.]),
                                          np.array([30., 30., 30.])))
                best = min(best, float(np.sum(s.fun ** 2)))
            return best

        CAP = 1e5
        lo, hi = 0.0, CAP
        if prof(CAP) < 3.84:
            hw = f"> {CAP:.0e}"
        else:
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                (lo, hi) = (mid, hi) if prof(mid) < 3.84 else (lo, mid)
            hw = f"{lo:.4f}"
        print(f"  {cv*100:8.2f}% {np.sqrt(np.mean(r**2)):26.6e} {c:13.6e} {hw:>29}")


# --------------------------------------------------------------------------
def s6_full_ode(times=MONTHLY, bvals=(30.0, 100.0, 200.0)):
    banner("S6  is the fibre exact for the ODE, or only for the algebraic locus?")

    print("  The fibre is derived on the quasi-steady balance. The record is")
    print("  generated by the ODE. Below: integrate dm/dt = L + b - (theta +")
    print("  Delta(t)) g(m) - k_e m at the truth and at fibre points, sample at")
    print("  the visits, and compare against the assay sd on the marker.\n")

    def run(th, km, ke, b, T=760.0):
        f = lambda t, y: [LB + b - (th + float(delta(t))) * g(y[0], km) - ke * y[0]]
        m0 = levels(np.array([0.0]), th, km, ke, b)[0][0]
        s = solve_ivp(f, (0.0, T), [m0], dense_output=True, rtol=1e-11,
                      atol=1e-13, method="LSODA")
        return s.sol(times)[0]

    base = run(TH0, KM0, KE0, 0.0)
    assay_sd = (ASSAY_CV) * float(np.mean(base))
    print(f"    marker at the visits: mean {base.mean():.6f}, "
          f"assay sd at CV {ASSAY_CV:.0%} = {assay_sd:.6f}")
    print(f"\n  {'b':>8} {'max |m_fibre - m_true|':>24} {'as a fraction of assay sd':>27}")
    print("  " + "-" * 62)
    for b in bvals:
        mm = run(*fibre(b))
        d = np.max(np.abs(mm - base))
        print(f"  {b:8.2f} {d:24.6e} {d/assay_sd:27.3e}")


# --------------------------------------------------------------------------
def s7_reconcile():
    banner("S7  reconciling with StructuralIdentifiability.jl")

    print("  The Julia verdict is about the TRANSIENT from a generic initial")
    print("  condition, where every time derivative of the output is available.")
    print("  Below: how many independent equations a transient supplies against")
    print("  what the quasi-steady locus supplies, at one visit.")

    th, km, ke, b, L, m = sp.symbols("theta K_m k_e b L m", positive=True)
    dl = sp.Symbol("Delta", nonnegative=True)
    rhs = L + b - (th + dl) * m / (km + m) - ke * m
    print(f"\n    dm/dt        = {rhs}")
    d2 = sp.simplify(sp.diff(rhs, m) * rhs)
    print(f"    d2m/dt2      = (dm/dt) * d/dm (dm/dt), which brings in")
    print(f"                   d/dm = {sp.simplify(sp.diff(rhs, m))}")
    print()
    print("  On the locus dm/dt = 0, so every higher derivative vanishes too and")
    print("  the second equation is 0 = 0. The transient's extra equations are")
    print("  exactly the ones the clinic does not sample.")
    print()
    print("  Time constant of the marker at the operating point:")
    tau = 1.0 / (TH0 * KM0 / (KM0 + 10.0) ** 2 + KE0)
    print(f"    tau = 1/(theta K_m/(K_m+m)^2 + k_e) = {tau:.6f} d = {tau*24:.4f} h")
    print(f"    against a therapeutic half-life of 3.5 d: separation "
          f"{3.5/tau:.1f}x")
    print()
    print("  So there is no contradiction: globally identifiable FROM THE")
    print("  TRANSIENT, and a one-parameter exact fibre ON THE LOCUS at constant")
    print("  load. The two statements are about different records.")


# --------------------------------------------------------------------------
def main():
    print(__doc__)
    s1_symbolic()
    s2_exact_or_second_order(MONTHLY, "monthly x 2 y")
    s2_exact_or_second_order(TWO_DAILY, "two-daily x 2 y")
    s3_true_profile()
    s3b_which_design_gave_36_9()
    s3c_how_well_is_k_e_pinned()
    s4_with_noise()
    s5_varying_load()
    s6_full_ode()
    s7_reconcile()


if __name__ == "__main__":
    main()
