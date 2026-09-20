"""
What the design costs at a disease-range capacity.

Phase 7 priority item 1, route (a). The author's decision of September 2026 is
to relabel the local block's `theta = 3000` as a NORMAL-CAPACITY illustration
rather than re-anchor it, and to measure what re-anchoring would do before
submission rather than leave a reviewer to find it.

THE CONTRADICTION BEING MEASURED. `src/common/models.py` annotates the local
`theta = 3000` as `V_max`, the maximal enzymatic velocity, and the parent
block's `theta = 300` as "~10% of normal". Those two annotations put normal at
3000 and the local patient at 100% of it. Residual PCC activity in propionic
acidemia is under 5% of control (`pcc_residual_activity`), which on the same
scale is `theta <= 150`. So the illustrative point is not a patient.

TWO CONVENTIONS, and they answer different questions.

    A  "same diet". The recorded load L stays at its theta = 3000 value and the
       marker rises. This is the physical reading: losing capacity does not
       change what a patient eats.
    B  "managed patient". The marker is held at m_bar = 10 and L is re-solved,
       which is the convention `n1_km_cost.py` uses for K_m because there the
       question was about curvature at a fixed operating point.

Convention A is the paper's own, and this script verifies rather than assumes
that: Section 11 reports the elasticity |d ln m / d ln theta| falling from
0.915 at nominal to 0.173 at 5%, and only convention A reproduces 0.173.
Convention B gives 0.219 at the same theta. Both are reported below.

WHAT MOVES WITH THETA, and all of it is recomputed rather than held:
the operating level, the elasticity, sigma (through `fit_ab`'s slope at the new
level), the load-bias prior `S_B = LOAD_FRAC * L`, the excitation ratio
Delta/theta, and the capacity ceiling.

THE CEILING IS THE SURPRISE AND IT IS REPORTED FIRST. `audit/h1_admissibility.py`
check 5 reads `normal = theta / frac`, so it takes theta to be the RESIDUAL. At
theta = 3000 and a 5% reference that puts normal at 60000 and the booked
D = 6000 well inside it. At a disease-range theta the same reading puts normal
at 3000 and the booked increment above it, so D = 6000 is not admissible for
the patient the disease describes. Every precision below is therefore reported
twice: at the booked D, which is the comparable number, and at the largest D
the ceiling allows at that theta, which is the physical one.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

from common.models import ParamsA
from design.run import delta, D_AMP, LAMBDA, DOSE_INTERVAL, LOAD_FRAC, ASSAY_CV
import model_A.run as LOC

KM0, KE, M_BAR0, TH0 = 100.0, 5.0, 10.0, 3000.0
D0 = D_AMP
GATE = 0.20
SEEDS = range(24)

#: `delta` sums 200 decaying doses per call, so it is evaluated ONCE on a long
#: monthly grid and sliced thereafter. Every schedule below is a prefix of this
#: grid, and `delta_at` only rescales, so nothing is approximated by doing so.
GRID = np.arange(0.0, 600 * 30.4, 30.4)
GRID_INC = np.asarray(delta(GRID), float)
N_MONTHLY = int(np.sum(GRID < 730.0))
MONTHLY = GRID[:N_MONTHLY]

#: the densest schedule Section 7 lists, kept so the disease-range question is
#: not answered only on the schedule a clinic would actually run.
DENSE = np.arange(0.0, 730.0, 2.0)
DENSE_INC = np.asarray(delta(DENSE), float)

#: the residual-activity reference settled by the author at 5%. See
#: audit/h1_admissibility.py, which adopts the same figure.
REFERENCE_FRAC = 0.05
#: theta at the reference, on the scale models.py's own two annotations imply
TH_DISEASE = TH0 * REFERENCE_FRAC


def delta_at(n_or_t, D):
    """
    The same dose train, scaled. `delta` is linear in its amplitude.

    Takes a visit COUNT and slices the precomputed grid, or an explicit time
    array, in which case it pays for a fresh evaluation.
    """
    if np.isscalar(n_or_t):
        inc = GRID_INC[:int(n_or_t)]
    else:
        tt = np.asarray(n_or_t, float)
        if tt.size == DENSE.size and tt[1] == DENSE[1]:
            inc = DENSE_INC
        elif tt.size <= GRID.size and np.allclose(tt, GRID[:tt.size]):
            inc = GRID_INC[:tt.size]
        else:
            inc = np.asarray(delta(tt), float)
    return inc * (D / D0)


def peak_increment(D):
    return D / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))


def level(theta, Km, ke, L, inc=0.0):
    """Positive root of (theta + inc) m/(Km+m) + ke m = L."""
    B = ke * Km + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * L * Km)) / (2.0 * ke)


def elasticity(theta, Km, ke, m):
    """|d ln m / d ln theta| at the untreated balance, in closed form."""
    g = m / (Km + m)
    gp = Km / (Km + m) ** 2
    return float(theta * g / (m * (theta * gp + ke)))


def anchor(theta, Km=KM0, ke=KE, convention="A"):
    """(L, m_bar) at this theta under the named convention."""
    if convention == "A":
        L = TH0 * M_BAR0 / (Km + M_BAR0) + ke * M_BAR0
        return L, level(theta, Km, ke, L)
    L = theta * M_BAR0 / (Km + M_BAR0) + ke * M_BAR0
    return L, M_BAR0


def footing(theta, Km=KM0, ke=KE, convention="A"):
    L, m_bar = anchor(theta, Km, ke, convention)
    p = ParamsA(theta=theta, Km=Km, ke=ke)
    sigma = float(LOC.sigma_total(LOAD_FRAC, ASSAY_CV, p=p, m_bar=m_bar))
    return L, m_bar, sigma, LOAD_FRAC * L


def sd_ln_theta(theta, D=D0, Km=KM0, ke=KE, convention="A", load_cv=0.10,
                seed=0, ke_prior=0.8, load_frac=None, times=MONTHLY):
    """sd(ln theta-hat) from the information matrix in (ln th, ln Km, b, ln ke)."""
    lf = LOAD_FRAC if load_frac is None else load_frac
    L, m_bar, _, _ = footing(theta, Km, ke, convention)
    p = ParamsA(theta=theta, Km=Km, ke=ke)
    sigma = float(LOC.sigma_total(lf, ASSAY_CV, p=p, m_bar=m_bar))
    S_B = lf * L
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        loads = np.full(times.size, L)
    inc = delta_at(times, D)
    B = ke * Km + theta + inc - loads
    m = (-B + np.sqrt(B ** 2 + 4.0 * ke * loads * Km)) / (2.0 * ke)
    g = m / (Km + m)
    Q = np.stack([theta * g,
                  -(theta + inc) * Km * m / (Km + m) ** 2,
                  np.ones_like(m),
                  ke * m], axis=1)
    F = (np.diag([0.0, 0.0, 1.0 / S_B ** 2, 1.0 / ke_prior ** 2])
         + Q.T @ Q / sigma ** 2)
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over_seeds(**kw):
    v = np.array([sd_ln_theta(seed=s, **kw) for s in SEEDS], float)
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


def visits_to_target(theta, D=D0, cap=None, **kw):
    """
    Monthly visits needed to reach the gate, median over seeds.

    Doubling then bisection, not a linear scan: precision improves with the
    record, so the predicate is monotone in n up to seed noise, and the
    bisection is confirmed by evaluating both sides of the answer.
    """
    cap = (GRID.size if cap is None else cap)

    def met(n):
        return over_seeds(theta=theta, D=D, times=GRID[:n], **kw)[0] <= GATE

    n = 4
    if met(n):
        return n
    while n < cap:
        n2 = min(2 * n, cap)
        if met(n2):
            lo, hi = n, n2
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if met(mid):
                    hi = mid
                else:
                    lo = mid
            return hi
        n = n2
    return None


def crossing(f, lo, hi, iters=60):
    """Bisect on a monotone-in-practice quantity. Returns None if no crossing."""
    flo, fhi = f(lo), f(hi)
    if (flo - GATE) * (fhi - GATE) > 0:
        return None
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        fmid = f(mid)
        if (flo - GATE) * (fmid - GATE) <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return float(np.sqrt(lo * hi))


def d_ceiling(theta, frac=REFERENCE_FRAC):
    """Largest D with theta + Delta_peak <= normal, normal = theta/frac."""
    return float(theta * (1.0 / frac - 1.0) * (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL)))


# ------------------------------------------------------------------ re-checks
def recheck(rows):
    """
    Four independent routes, each sharing no algebra with what it checks.

    R1  the closed-form operating level against a bracketed root solve of the
        balance itself.
    R2  the closed-form elasticity against a central difference in ln theta,
        Richardson-extrapolated.
    R3  the ordinate is L identically along the balance, so its TOTAL
        derivative in ln theta must vanish. A first version of this check
        compared the measured value against eps|f|/h and failed, and the
        threshold was the thing that was wrong: the closed-form root subtracts
        two numbers of order 3200 to get one of order 100, so the ordinate
        carries about thirty times the round-off that |f| alone suggests, and
        no hand-derived constant for that is trustworthy. The check is
        therefore made scale-free instead of loosened: a round-off artefact
        grows as 1/h when h is halved, and a genuine nonzero derivative does
        not move at all. Measuring the ratio settles which it is with no
        threshold to choose.
    R4  Section 11 of the paper prints the marker and the elasticity at seven
        residual fractions. Those were computed elsewhere, before this script
        existed, so reproducing them is an outside check on the whole footing.
    """
    print("\n  RE-CHECKS")
    worst_lvl = worst_el = 0.0
    ratios, mags, ord_scale = [], [], 1.0
    for theta, conv in rows:
        L, m_bar, _, _ = footing(theta, convention=conv)
        root = brentq(lambda m: theta * m / (KM0 + m) + KE * m - L,
                      1e-9, 1e9, xtol=1e-14, rtol=8.9e-16)
        worst_lvl = max(worst_lvl, abs(root - m_bar) / m_bar)

        def lnm(lt):
            return np.log(level(np.exp(lt), KM0, KE, L))
        lt = np.log(theta)
        d1 = [(lnm(lt + h) - lnm(lt - h)) / (2 * h) for h in (1e-4, 5e-5)]
        rich = (4 * d1[1] - d1[0]) / 3.0
        worst_el = max(worst_el, abs(abs(rich) - elasticity(theta, KM0, KE, m_bar)))

        inc = delta_at(MONTHLY, D0)
        loads = np.full(MONTHLY.size, L)

        def ordinate(lt):
            th = np.exp(lt)
            mm = level(th, KM0, KE, loads, inc)
            return (th + inc) * mm / (KM0 + mm) + KE * mm
        tot = []
        for h in (2e-5, 1e-5, 5e-6):
            tot.append(float(np.abs(
                (ordinate(lt + h) - ordinate(lt - h)) / (2 * h)).max()))
        ratios += [tot[1] / max(tot[0], 1e-300), tot[2] / max(tot[1], 1e-300)]
        mags.append(tot[0])
        ord_scale = float(np.abs(ordinate(lt)).max())

    print(f"    R1 closed-form level against brentq: worst relative "
          f"{worst_lvl:.3e} against 1e-12")
    print(f"    R2 closed-form elasticity against Richardson central "
          f"difference: worst {worst_el:.3e} against 1e-8")
    print(f"    R3 the ordinate is L identically along the balance. Halving h "
          f"scales its apparent")
    print(f"       total derivative by {np.median(ratios):.2f} (median over "
          f"{len(ratios)} halvings, range "
          f"{min(ratios):.2f} to {max(ratios):.2f}; the individual ratios of a "
          f"round-off")
    print(f"       quantity are themselves noise, which is why the median is "
          f"the statistic).")
    print(f"       Round-off predicts 2, a genuine derivative predicts 1. "
          f"It is round-off, and")
    print(f"       the magnitude agrees: the largest apparent derivative is "
          f"{max(mags):.2e} on an")
    print(f"       ordinate of {ord_scale:.1f}, a relative "
          f"{max(mags)/ord_scale:.1e}.")

    # R4: Section 11's own table, computed before this script existed
    s11 = {1.00: (1.00, 0.915), 0.50: (1.82, 0.807), 0.20: (3.41, 0.537),
           0.10: (4.57, 0.321), 0.05: (5.40, 0.173), 0.02: (6.00, 0.072),
           0.01: (6.22, 0.036)}
    worst_m = worst_e = 0.0
    print(f"    R4 Section 11's table, recomputed:")
    print(f"       {'% nominal':>10} {'m/managed':>20} {'elasticity':>20}")
    for f, (mm, el) in sorted(s11.items(), reverse=True):
        _, m_bar, _, _ = footing(TH0 * f, convention="A")
        e = elasticity(TH0 * f, KM0, KE, m_bar)
        worst_m = max(worst_m, abs(m_bar / M_BAR0 - mm))
        worst_e = max(worst_e, abs(e - el))
        print(f"       {100*f:9.0f}% {mm:9.2f} vs {m_bar/M_BAR0:8.2f} "
              f"{el:10.3f} vs {e:9.3f}")
    print(f"       worst disagreement {worst_m:.4f} on the marker and "
          f"{worst_e:.4f} on the elasticity,")
    print(f"       against the half-unit each is printed to, 0.005 and 0.0005.")
    ok = (worst_lvl < 1e-12 and worst_el < 1e-8
          and 1.7 < np.median(ratios) < 2.3
          and worst_m < 0.005 and worst_e < 0.0005)
    print(f"    RE-CHECKS {'PASS' if ok else 'FAIL'}")
    return ok


# --------------------------------------------------------------- the tables
#  Section 6 typed these two tables into the .tex and Section 5 typed the
#  third. They are generated here because this is where the convention lives:
#  the recorded load is held at the nominal patient's and the operating level,
#  sigma, the load-bias prior and the excitation are re-solved at each
#  capacity. Reproducing that convention anywhere else would duplicate it.
def scope_numbers():
    """
    The scope result as generated numbers, so the three sections that quote it
    cannot drift apart.

    Section 2.5, Section 6.2 and the abstract all carry the band and the
    enzymatic share, and all three carried them as hand-typed literals. The
    sweep is information-matrix work rather than simulation, so it is cheap
    enough to run wherever the tables are written.
    """
    from common.numbers import write_numbers
    f_m = lambda th: over_seeds(theta=th, D=D0, convention="A",
                                times=MONTHLY)[0]
    lo_edge = crossing(f_m, 3000.0, 300.0)
    hi_edge = crossing(f_m, 6000.0, 30000.0)
    lo_dense = crossing(lambda th: over_seeds(theta=th, D=D0, convention="A",
                                              times=DENSE)[0], 3000.0, 100.0)
    scan = [(th, f_m(th)) for th in np.geomspace(1000.0, 1.0e5, 22)]
    best = min(scan, key=lambda r: r[1])

    def share(th):
        """theta's own share of the load: the flux identity of (16)."""
        L, m, _, _ = footing(th, convention="A")
        enz = th * m / (KM0 + m)
        return enz / (enz + KE * m)

    N = {
        "band_lo": lo_edge, "band_hi": hi_edge,
        "band_lo_pct": 100.0 * lo_edge / TH0,
        "band_hi_pct": 100.0 * hi_edge / TH0,
        "band_fold": hi_edge / lo_edge,
        "band_best": best[1], "band_best_pct": 100.0 * best[0] / TH0,
        "band_lo_dense_pct": 100.0 * lo_dense / TH0,
        "band_dense_gain_pct": 100.0 * (1.0 - lo_dense / lo_edge),
        "disease_factor": lo_edge / TH_DISEASE,
        "share_nominal_pct": 100.0 * share(TH0),
        "share_10pct": 100.0 * share(0.10 * TH0),
        "share_5pct": 100.0 * share(REFERENCE_FRAC * TH0),
        "sd_10pct": over_seeds(theta=0.10 * TH0, D=D0, convention="A",
                               times=MONTHLY)[0],
        "sd_5pct": over_seeds(theta=REFERENCE_FRAC * TH0, D=D0,
                              convention="A", times=MONTHLY)[0],
        #  the residual activity below which the capacity-ceiling check
        #  cannot report a breach, because theta cancels from its verdict
        "ceiling_blind_pct": 100.0 / (1.0 + peak_increment(D0) / TH0),
    }
    fmt = {"band_lo": "{:.0f}", "band_hi": "{:.0f}", "band_lo_pct": "{:.0f}",
           "band_hi_pct": "{:.0f}", "band_fold": "{:.1f}",
           "band_best": "{:.3f}", "band_best_pct": "{:.0f}",
           "band_lo_dense_pct": "{:.0f}", "band_dense_gain_pct": "{:.0f}",
           "disease_factor": "{:.0f}", "share_nominal_pct": "{:.1f}",
           "share_10pct": "{:.1f}", "share_5pct": "{:.1f}",
           "sd_10pct": "{:.4f}", "sd_5pct": "{:.4f}",
           "ceiling_blind_pct": "{:.0f}"}
    write_numbers("scope", N, fmt)
    return N


def table_elasticity():
    """How far the marker moves per unit of capacity, capacity by capacity."""
    from common.style import write_table
    rows = []
    for pct in (100, 50, 20, 10, 5, 2, 1):
        th = TH0 * pct / 100.0
        _, m = anchor(th, convention="A")
        rows.append([f"${pct}\\%$", f"${m / M_BAR0:.2f}$",
                     f"${elasticity(th, KM0, KE, m):.3f}$"])
    write_table("tab_S_elasticity",
                [r"$\theta$, as \% of nominal", r"marker, $\times$ managed",
                 r"$\lvert\partial\ln m/\partial\ln\theta\rvert$"],
                rows, colspec="rrr",
                note="closed form at the untreated balance, R2 of this script")
    print("  wrote tab_S_elasticity.tex")


def table_scope():
    """What the monthly design delivers at each capacity, and the visits to the gate."""
    from common.style import write_table
    rows = []
    for pct in (100, 20, 10, 5, 2):
        th = TH0 * pct / 100.0
        _, m = anchor(th, convention="A")
        med, lo, hi = over_seeds(theta=th, D=D0, convention="A")
        n = visits_to_target(th, D=D0, convention="A")
        sd = (r"$\mathbf{%.4f}$" % med) if med < GATE else f"${med:.4f}$"
        rows.append([f"${pct}\\%$", f"${m / M_BAR0:.2f}$",
                     f"${elasticity(th, KM0, KE, m):.3f}$", sd,
                     ("---" if n is None or not np.isfinite(n) else f"${int(n)}$")])
    write_table("tab_S_scope",
                [r"$\theta$, as \% of nominal", r"marker, $\times$ managed",
                 r"$\lvert\partial\ln m/\partial\ln\theta\rvert$",
                 r"$\operatorname{sd}(\ln\hat\theta)$",
                 r"monthly visits to $%.2f$" % GATE],
                rows, colspec="rrrrr",
                note=("medians over %d intake records at a 10%% day-to-day "
                      "variability; a dash means the gate is not reached at "
                      "any number of monthly visits up to 600" % len(list(SEEDS))))
    print("  wrote tab_S_scope.tex")


def table_ke_prior():
    """What a tighter prior on the alternative clearance is worth, at two capacities."""
    from common.style import write_table
    rows = []
    for s, label in ((0.8, r"$0.8$ (as published here)"), (0.2, r"$0.2$"),
                     (0.05, r"$0.05$"), (0.001, r"exact")):
        nom = over_seeds(theta=TH0, D=D0, convention="A", ke_prior=s)[0]
        dis = over_seeds(theta=0.05 * TH0, D=D0, convention="A", ke_prior=s)[0]
        rows.append([label, f"${np.exp(s):.2f}$", f"${nom:.4f}$", f"${dis:.4f}$"])
    write_table("tab_S_keprior",
                [r"prior on $\ln k_e$", "factor per sd",
                 r"$\operatorname{sd}(\ln\hat\theta)$ at nominal",
                 r"at $5\%$ of nominal"],
                rows, colspec="lrrr",
                note=("medians over %d intake records; the last row is a prior "
                      "of 0.001, which is k_e known to a tenth of a per cent"
                      % len(list(SEEDS))))
    print("  wrote tab_S_keprior.tex")


def main():
    scope_numbers()
    table_elasticity()
    table_scope()
    table_ke_prior()
    print("=" * 96)
    print("  WHAT THE DESIGN COSTS AT A DISEASE-RANGE CAPACITY")
    print("=" * 96)
    print(f"\n  normal on the scale models.py's annotations imply: {TH0:.0f}")
    print(f"  residual reference {REFERENCE_FRAC:.0%}  ->  disease-range "
          f"theta <= {TH_DISEASE:.0f}")
    print(f"  booked D = {D0:.0f}, peak increment {peak_increment(D0):.0f}, "
          f"monthly record of {MONTHLY.size} visits\n")

    print("  1. THE OPERATING POINT MOVES, AND THE ELASTICITY WITH IT")
    print(f"  {'theta':>8} {'conv':>5} {'L':>9} {'m_bar':>9} {'|dlnm/dlnth|':>13}"
          f" {'Delta_pk/theta':>15} {'sigma':>9}")
    thetas = [TH0, 600.0, 300.0, 150.0, 100.0, 60.0, 30.0]
    for th in thetas:
        for conv in ("A", "B"):
            L, m_bar, sigma, _ = footing(th, convention=conv)
            print(f"  {th:8.0f} {conv:>5} {L:9.2f} {m_bar:9.3f} "
                  f"{elasticity(th, KM0, KE, m_bar):13.4f} "
                  f"{peak_increment(D0)/th:15.1f} {sigma:9.3f}")
    print("\n     Section 11 reports 0.915 at nominal and 0.173 at 5%. "
          "Convention A reproduces")
    a150 = elasticity(150.0, KM0, KE, footing(150.0, convention='A')[1])
    b150 = elasticity(150.0, KM0, KE, footing(150.0, convention='B')[1])
    print(f"     both ({elasticity(TH0, KM0, KE, M_BAR0):.3f} and {a150:.3f}); "
          f"convention B gives {b150:.3f} at the same theta.")
    print("     So the paper's own scope section already uses convention A, "
          "and the operating")
    print(f"     level it implies at theta = 150 is m = "
          f"{footing(150.0, convention='A')[1]:.1f} uM, which sits in the "
          "untreated disease")
    print("     range the m_bar row quotes rather than at the managed 10.")

    print("\n  2. THE CAPACITY CEILING, AND WHY THE CHECK CANNOT FAIL AS CODED")
    print("     h1 check 5 evaluates used = (theta + Delta_peak)/normal with")
    print("     normal = theta/frac, and scales Delta with theta across blocks so")
    print("     that Delta_peak/theta is 2.1333 in both. Those two together give")
    print(f"     used = frac x (1 + {peak_increment(D0)/TH0:.4f}) = "
          f"{1+peak_increment(D0)/TH0:.4f} x frac, with theta cancelled out, so the")
    print(f"     check can only fail above frac = "
          f"{1/(1+peak_increment(D0)/TH0):.3f}, which is a residual activity")
    print("     no one claims for this disease. As coded it cannot report a breach.")
    print("     What makes it bite is holding D at its booked value instead of")
    print("     scaling it, which is what a mg/kg dose is:\n")
    print(f"  {'theta':>8} {'% nominal':>10} {'normal=th/frac':>15} "
          f"{'th+Delta_pk':>13} {'used':>9} {'admissible?':>12} {'largest D':>11}")
    for th in thetas:
        nrm = th / REFERENCE_FRAC
        tot = th + peak_increment(D0)
        print(f"  {th:8.0f} {100*th/TH0:9.0f}% {nrm:15.0f} {tot:13.0f} "
              f"{tot/nrm:9.4f} {'yes' if tot <= nrm else 'NO':>12} "
              f"{d_ceiling(th):11.0f}")
    print(f"\n     At theta = {TH_DISEASE:.0f} the booked D = {D0:.0f} drives "
          f"the enzyme to "
          f"{(TH_DISEASE+peak_increment(D0))/(TH_DISEASE/REFERENCE_FRAC):.2f}"
          f" times normal.")
    print(f"     The largest admissible increment there is "
          f"D = {d_ceiling(TH_DISEASE):.0f}, a factor of "
          f"{D0/d_ceiling(TH_DISEASE):.2f} below the booked value.")
    print("     Under the relabel the local block is a NORMAL-capacity")
    print("     illustration, so on it the ceiling reads theta + Delta_peak <= theta,")
    print("     which no positive dose satisfies. The constraint belongs to the")
    print("     parent block, where it is a statement about a patient.")

    print("\n  3. PRECISION AND HEADROOM, ON THE PAPER'S OWN CONVENTION A")
    print(f"  {'theta':>7} {'%nom':>5} {'conv':>5} {'D':>6} {'schedule':>9} "
          f"{'n':>4} {'sd(ln th)':>10} {'spread':>19} {'visits':>7} "
          f"{'K_m gate':>9} {'head':>6} {'load gate':>10} {'head':>6}")
    cases = []
    for th in (TH0, 600.0, 300.0, TH_DISEASE, 60.0):
        cases.append((th, "A", D0, "monthly", MONTHLY))
    for th in (300.0, TH_DISEASE):
        cases.append((th, "A", d_ceiling(th), "monthly", MONTHLY))
    cases += [(TH0, "A", D0, "2-daily", DENSE),
              (TH_DISEASE, "A", D0, "2-daily", DENSE),
              (TH_DISEASE, "A", d_ceiling(TH_DISEASE), "2-daily", DENSE),
              (TH_DISEASE, "B", D0, "monthly", MONTHLY)]
    for th, conv, D, nm, tt in cases:
        med, lo, hi = over_seeds(theta=th, D=D, convention=conv, times=tt)
        n = (visits_to_target(th, D=D, convention=conv) if nm == "monthly"
             else None)
        kmg = crossing(lambda K: over_seeds(theta=th, D=D, Km=K,
                                            convention=conv, times=tt)[0],
                       5.0, 3000.0)
        lfg = crossing(lambda f: over_seeds(theta=th, D=D, convention=conv,
                                            times=tt, load_frac=f)[0],
                       0.002, 0.80)
        print(f"  {th:7.0f} {100*th/TH0:4.0f}% {conv:>5} {D:6.0f} {nm:>9} "
              f"{tt.size:4d} {med:10.4f} [{lo:8.4f},{hi:8.4f}] "
              f"{(str(n) if nm == 'monthly' else '  --'):>7} "
              f"{('%.0f' % kmg) if kmg else '  none':>9} "
              f"{('%.2f' % (kmg/KM0)) if kmg else '   --':>6} "
              f"{('%.1f%%' % (100*lfg)) if lfg else ' none':>10} "
              f"{('%.2f' % (lfg/LOAD_FRAC)) if lfg else '   --':>6}")
    print(f"\n     `visits` is monthly visits to reach sd = {GATE}, median over "
          f"{len(list(SEEDS))} records, searched to {GRID.size}.")
    print("     `none` in the gate columns means the design misses the target at "
          "every value swept,")
    print("     so there is no headroom to report rather than a large one.")

    print("\n  3b. THE BAND OF CAPACITIES THE DESIGN ACTUALLY WORKS IN")
    print("      sd(ln theta) is NOT monotone in theta, so the gate is crossed")
    print("      twice and each crossing is bracketed separately. Below the band")
    print("      the marker stops responding to capacity (Section 11). Above it")
    print("      the theta column saturates at L while the assay term keeps")
    print("      growing, so the information per draw stops improving.")
    scan = [(th, over_seeds(theta=th, D=D0, convention="A", times=MONTHLY)[0])
            for th in np.geomspace(1000.0, 1.0e5, 22)]
    print(f"  {'theta':>9} {'% nominal':>10} {'sd(ln theta)':>13}")
    for th, v in scan:
        mark = "  <= gate" if v <= GATE else ""
        print(f"  {th:9.0f} {100*th/TH0:9.0f}% {v:13.4f}{mark}")
    f_m = lambda th: over_seeds(theta=th, D=D0, convention="A",
                                times=MONTHLY)[0]
    lo_edge = crossing(f_m, 3000.0, 300.0)
    hi_edge = crossing(f_m, 6000.0, 30000.0)
    best_th = min(scan, key=lambda r: r[1])
    lo_dense = crossing(lambda th: over_seeds(theta=th, D=D0, convention="A",
                                              times=DENSE)[0], 3000.0, 100.0)
    print(f"\n     monthly, {MONTHLY.size} draws: the gate holds for theta in "
          f"[{lo_edge:.0f}, {hi_edge:.0f}],")
    print(f"       that is [{100*lo_edge/TH0:.0f}%, {100*hi_edge/TH0:.0f}%] of "
          f"nominal, a band of {hi_edge/lo_edge:.1f} fold.")
    print(f"       Best precision {best_th[1]:.4f} at theta = {best_th[0]:.0f}, "
          f"{100*best_th[0]/TH0:.0f}% of nominal.")
    print(f"     every two days for two years, {DENSE.size} draws: the lower "
          f"edge moves only to")
    print(f"       theta = {lo_dense:.0f}, {100*lo_dense/TH0:.0f}% of nominal. "
          f"Fifteen times the sampling buys")
    print(f"       {100*(1-lo_dense/lo_edge):.0f}% on the edge, which is the "
          f"quantitative form of Section 11's")
    print(f"       `this is not a statement about sampling effort`.")
    print(f"\n     THE FINDING. The disease is defined at under "
          f"{100*REFERENCE_FRAC:.0f}% of nominal. The lower edge of")
    print(f"     the working band is {100*lo_edge/TH0:.0f}%. The design "
          f"therefore works on the illustrative")
    print(f"     point and nowhere near the patient, by a factor of "
          f"{lo_edge/TH_DISEASE:.0f} in capacity.")

    print("\n  4. WHAT IT WOULD TAKE AT A DISEASE-RANGE CAPACITY")
    base = over_seeds(theta=TH_DISEASE, D=D0, convention="A", times=MONTHLY)[0]
    print(f"     monthly, theta = {TH_DISEASE:.0f}: sd = {base:.4f} against a "
          f"gate of {GATE}, a factor of {base/GATE:.1f}.")
    print(f"     Precision scales as the square root of the visit count, so "
          f"closing that factor")
    print(f"     by sampling alone needs about "
          f"{MONTHLY.size*(base/GATE)**2:,.0f} visits, which at a monthly "
          f"cadence is")
    print(f"     {MONTHLY.size*(base/GATE)**2/12:,.0f} years.")
    dense_sd = over_seeds(theta=TH_DISEASE, D=D0, convention="A",
                          times=DENSE)[0]
    print(f"     Every two days for two years ({DENSE.size} draws) gives "
          f"{dense_sd:.4f}, still {dense_sd/GATE:.0f} times the gate.")
    sig_needed = crossing(
        lambda f: over_seeds(theta=TH_DISEASE, D=D0, convention="A",
                             load_frac=f)[0], 1e-5, 0.80)
    print(f"     A perfect load record does not close it either: the sweep "
          f"finds no recording")
    print(f"     error in [0.001%, 80%] that reaches the gate "
          f"({'none' if sig_needed is None else '%.4f' % sig_needed}).")
    print(f"     This is Section 11's own conclusion, which says no schedule, "
          f"no excursion and no")
    print(f"     precision recovers a quantity the observable has stopped "
          f"responding to. What was")
    print(f"     missing there was the number, and the number is "
          f"{base:.2f} against {GATE}.")

    ok = recheck([(TH0, "A"), (300.0, "A"), (TH_DISEASE, "A"),
                  (TH_DISEASE, "B"), (60.0, "A")])
    print("\n" + "=" * 96)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
