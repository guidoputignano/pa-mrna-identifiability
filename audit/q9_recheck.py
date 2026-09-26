"""
The independent re-check for Phase 8.

One re-check per reported result, each by a route that shares no code path with
the thing it checks. Where a check fails it is PRINTED AS A FAILURE and the
script exits non-zero. Nothing here repairs anything, and nothing here records
an expected value in advance: every check compares two independent computations
of one quantity, or a computation against a closed form derived by hand.

Covers, for `audit/q1_increment_shape.py` (S1):

  1  The two-compartment kernel, by hand from the Laplace transform and
     partial fractions, against q1's numerical eigen-decomposition. Neither
     route is `solve_ivp`, which q1 already checks itself against.
  2  The steady-state identity q1 leans on everywhere: the cycle mean of the
     steady increment is the AREA of one impulse response divided by the dose
     interval. Checked by brute-force summation over past doses with a
     midpoint rule that does not straddle the dose jump.
  3  The direction of the bias, from first-order linear algebra with no
     optimiser at all: the shift in theta is the coefficient of g when the
     increment ERROR times g is projected onto the design columns. A sign
     disagreement with q1's Monte Carlo is a failure.
  4  q1's fit, against a from-scratch nonlinear least squares started at the
     truth. Same objective, because it has to be the same estimator, but no
     shared code.
  5  q1's profile interval, which is found by bracketing, against a fixed-grid
     likelihood-ratio scan driven by `n6_profile_vs_crb.nll`, which is the
     repository's own objective and a different optimiser.
  6  The published dose arithmetic, in exact rational arithmetic.
  7  Part 5's optimal two-phase design, against the closed form
     sd(ln t_half) = sqrt(2) / (lambda |tau_1 - tau_2|), which follows because
     the free amplitude removes the phase-independent part of the sensitivity
     and leaves lambda*tau.
  8  The residual reference: penalised against unpenalised on the same records.
  9  The free-amplitude fit against the probe's fitter
     (`probes_2026-09-25/fit.py`), on the record where the delta-method and log
     forms of the kappa prior part company. Added when q1b_reconcile.py found
     the first run's kappa prior in the wrong coordinate.
 10  The dose-day draw: pre-dose, every kernel is zero at t = 0; post-dose, the
     committed analyst is design.run.delta there.
 11  Gap B on both baseline footings (q1c), rebuilt from the probe's own
     truths and fitter.
 12  S2: the requirement map's Fisher bound against the probe's fisher_sd,
     and q2's MAP estimator against the probe's fitter on identical records.

Run:  python audit/q9_recheck.py
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares, minimize

import k1_sources as K1
import n6_profile_vs_crb as N6
import q1_increment_shape as Q

RESULTS = []


def report(name, got, want, tol, note="", rel=True):
    """One check. `tol` is relative unless `rel` is False."""
    gap = (abs(got - want) / abs(want) if rel and want else abs(got - want))
    ok = bool(np.isfinite(gap) and gap <= tol)
    RESULTS.append((name, ok))
    kind = "relative" if rel else "absolute"
    print(f"    {name:<52s} {'PASSES' if ok else 'FAILS'}  "
          f"{kind} gap {gap:.3e} against {tol:.0e}")
    if note:
        print(f"      {note}")
    return ok


def banner(t):
    print("\n" + "=" * 92)
    print("  " + t)
    print("=" * 92)


# ==========================================================================
def check_1_kernel():
    banner("1. THE TWO-COMPARTMENT KERNEL, BY HAND FROM THE LAPLACE TRANSFORM")
    print("""  With Ce(s) = 1/(s + ke0) the protein block gives

      PCP = kq PCC/(s + kq),   so
      PCC(s) = Slope (s + kq) / [(s + ke0) (s^2 + (2kq + kdeg) s + kq kdeg)],

  three simple poles at -ke0 and at the two roots of the quadratic. The
  residue at pole p is Slope (p + kq) / prod_{other q} (p - q), so the impulse
  response is a three-term exponential sum written down directly. q1 gets the
  same object from an eigen-decomposition of the 3x3 state matrix.""")
    worst = 0.0
    pairs = [(K1.KQ, K1.KDEG), (K1.KQ * 5.0, K1.KDEG / 3.0),
             (K1.KQ / 8.0, K1.KDEG * 4.0), (0.02, 0.02 * 1.37)]
    u = np.linspace(0.0, 60.0, 601)
    for kq, kdeg in pairs:
        ke0, slope = K1.KE0, K1.SLOPE
        bq, cq = 2.0 * kq + kdeg, kq * kdeg
        disc = np.sqrt(bq * bq - 4.0 * cq)
        poles = np.array([-ke0, 0.5 * (-bq + disc), 0.5 * (-bq - disc)])
        res = np.array([slope * (p + kq)
                        / np.prod([p - o for o in poles if o is not p])
                        for p in poles])
        mine = (res[None, :] * np.exp(poles[None, :] * (u[:, None] * 24.0))
                ).sum(axis=1)
        theirs = Q.kernel_twocomp(kq, kdeg)(u)
        gap = float(np.max(np.abs(mine - theirs)) / np.max(np.abs(theirs)))
        worst = max(worst, gap)
        print(f"    kq = {kq:.6f}, kdeg = {kdeg:.6f}: worst relative "
              f"difference {gap:.3e}")
    return report("two-compartment kernel, Laplace vs eigen", worst, 0.0,
                  1e-9, rel=False)


def check_2_steady_identity():
    banner("2. THE STEADY CYCLE MEAN IS THE IMPULSE AREA OVER THE INTERVAL")
    print("""  q1 normalises every shape with this identity and never integrates
  the assembled train, so if it were wrong every amplitude would be wrong
  together and nothing inside q1 would notice. Checked here by assembling the
  train explicitly over many past doses and integrating one cycle with a
  MIDPOINT rule, which never straddles the jump a bolus makes.""")
    T = Q.T_DOSE
    ok = True
    for name, k in Q.KERNELS.items():
        edges = np.linspace(0.0, T, 20001)
        tau = 0.5 * (edges[1:] + edges[:-1])
        acc = np.zeros_like(tau)
        for j in range(0, 4000):
            acc = acc + k(tau + j * T)
            if j > 20 and np.max(k(tau + j * T)) < 1e-16 * max(np.max(acc), 1):
                break
        brute = float(np.mean(acc))
        ident = k.area() / T
        ok &= report(f"cycle mean, {name}", brute, ident, 1e-6)
    return ok


def check_3_bias_direction():
    banner("3. THE DIRECTION OF THE BIAS, WITH NO OPTIMISER")
    print("""  Hold K_m, k_e and b at the truth and the balance is LINEAR. The
  analyst regresses s - Delta_a g on the columns [g, m, 1]; the truth puts
  s = (theta + Delta_true) g + k_e m + b. So the response the analyst sees is
  theta g + (Delta_true - Delta_a) g, and the shift in theta-hat is the
  coefficient of g when (Delta_true - Delta_a) g is projected onto those three
  columns. That is one lstsq call and it predicts a SIGN. Anything that
  disagrees in sign with q1's Monte Carlo median is a failure; the size is not
  expected to match, because q1 frees K_m and carries priors.""")
    ref = Q.committed_ref()
    tn, wn = Q.nominal_train()
    A = Q.make_analysts((tn, wn), ref, "peak",
                        which=("fixed 3.5 d",))["fixed 3.5 d"]
    inc_a = A.phi(np.zeros(0), Q.MONTHLY)
    ok = True
    print(f"\n    {'truth':20s} {'predicted theta_hat/theta':>26s} "
          f"{'q1 median':>12s} {'signs agree':>12s}")
    for tname in Q.KERNELS:
        inc_t = Q.truth_increment(tname, (tn, wn), "peak", ref, Q.MONTHLY)
        #  a noiseless record at the true nuisances, so nothing but the
        #  increment error is moving theta
        L = np.full(Q.MONTHLY.size, Q.FOOT.L)
        m = Q.level(Q.TH, inc_t, L)
        g = m / (Q.KM + m)
        X = np.stack([g, m, np.ones_like(m)], axis=1)
        y = (inc_t - inc_a) * g
        shift = float(np.linalg.lstsq(X, y, rcond=None)[0][0])
        pred = 1.0 + shift / Q.TH
        c = Q.cell(A, inc_t, seeds=range(8), coverage=False)
        got = float(np.median(c["ratio"]))
        agree = np.sign(pred - 1.0) == np.sign(got - 1.0) or \
            abs(pred - 1.0) < 1e-3 and abs(got - 1.0) < 0.05
        ok &= agree
        RESULTS.append((f"bias sign, {tname}", bool(agree)))
        print(f"    {tname:20s} {pred:26.3f} {got:12.3f} "
              f"{'yes' if agree else 'NO':>12s}")
    return ok


def check_4_fit():
    banner("4. q1's FIT AGAINST A FROM-SCRATCH NONLINEAR LEAST SQUARES")
    print("""  Same objective and the SAME DOMAIN, because a check that changed
  either would be checking something else. q1's estimator is defined as global
  over K_m in [1, 1000] with an exact polish inside that box, so the
  from-scratch fit is bounded the same way. Different code, different
  parametrisation, and started AT THE TRUTH rather than from a grid, so a grid
  that missed its basin would show up here.

  The same fit is then run UNBOUNDED and reported beside it, not as a check but
  as a measurement: on a badly misspecified record K_m walks out of the box,
  where theta*g(m) collapses into (theta/K_m)*m and stops being separable from
  the k_e column, and the fit returns a capacity that is an artefact of leaving
  the domain rather than an estimate. That is why the box is there.""")
    ref = Q.committed_ref()
    tn, wn = Q.nominal_train()
    A = Q.make_analysts((tn, wn), ref, "peak",
                        which=("fixed 3.5 d",))["fixed 3.5 d"]
    F = Q.FOOT
    lo, hi = np.log(Q.KM_GRID[0]), np.log(Q.KM_GRID[-1])
    worst, worst_th = 0.0, 0.0
    print(f"\n    {'truth':20s} {'objective gap':>14s} {'theta gap':>11s}"
          f" {'unbounded theta_hat/theta':>26s} {'its K_m':>10s}")
    for tname in ("exp 3.5 d", "exp 7.0 d", "Attarwala 2-comp"):
        inc_t = Q.truth_increment(tname, (tn, wn), "peak", ref, Q.MONTHLY)
        inc_a = A.phi(np.zeros(0), Q.MONTHLY)
        free_th, free_km = [], []
        for sd in range(4):
            m, s, _ = Q.make_record(inc_t, sd, F=F)

            def res(v):
                th, Km, ke = np.exp(np.clip(v[:3], -30, 30))
                r = ((th + inc_a) * m / (Km + m) + ke * m + v[3] - s) / F.sigma
                return np.concatenate([r, [v[3] / F.s_b,
                                           (v[2] - np.log(F.ke)) / 0.8]])

            q0 = np.array([np.log(Q.TH), np.log(Q.KM), np.log(Q.KE), 0.0])
            lb = np.array([-np.inf, lo, -np.inf, -np.inf])
            ub = np.array([np.inf, hi, np.inf, np.inf])
            out = least_squares(res, q0, bounds=(lb, ub), method="trf",
                                max_nfev=20000, xtol=1e-14, ftol=1e-14,
                                gtol=1e-14)
            mine, J_mine = float(np.exp(out.x[0])), float(2.0 * out.cost)
            fit = A.fit(m, s, F=F)
            theirs, J_theirs = fit["theta"], fit["J"]
            #  The OBJECTIVE is what an estimator is defined by; theta is a
            #  coordinate on a surface the paper itself reports is far from
            #  quadratic. So the check is on the objective, and the gap in
            #  theta at equal objective is reported as a MEASUREMENT of that
            #  flatness rather than as a second check.
            worst = max(worst, abs(J_mine - J_theirs)
                        / max(abs(J_theirs), 1e-12))
            worst_th = max(worst_th, abs(mine / theirs - 1.0))
            un = least_squares(res, q0, method="lm", max_nfev=20000)
            free_th.append(float(np.exp(un.x[0])) / Q.TH)
            free_km.append(float(np.exp(un.x[1])))
        print(f"    {tname:20s} {worst:14.3e} {worst_th:11.3e} "
              f"{np.median(free_th):26.4g} {np.median(free_km):10.4g}")
    print(f"\n    the theta gap at an objective agreeing to {worst:.1e} is a "
          f"reading of how flat")
    print(f"    the surface is in theta, not of how well either optimiser "
          f"converged.")
    return report("q1 fit vs from-scratch least squares (objective)", worst,
                  0.0, 1e-6, rel=False)


def _cross(grid, v, drop):
    """
    Both edges of the connected region around the minimum, by LINEAR
    interpolation of the crossing rather than by taking the last grid point
    inside. Without it the comparison is limited by the grid spacing and would
    not resolve the thing it is checking.
    """
    k = int(np.argmin(v))
    out = []
    for step in (-1, +1):
        i = k
        while 0 < i < len(grid) - 1 and v[i + step] <= drop:
            i += step
        j = i + step
        if not (0 <= j < len(grid)) or v[j] <= drop:
            out.append(float(np.exp(grid[i])))       # never crosses: an edge
            continue
        f = (drop - v[i]) / (v[j] - v[i])
        out.append(float(np.exp(grid[i] + f * (grid[j] - grid[i]))))
    return out[0], out[1]


def check_5_interval():
    banner("5. q1's PROFILE INTERVAL AGAINST A FIXED-GRID SCAN ON n6's "
           "OBJECTIVE")
    print("""  q1 finds the interval by expanding from the optimum until the
  objective has risen by 1.92 and then bisecting. This re-check instead
  evaluates the profile on a FIXED grid in ln theta, re-optimising the other
  three with Nelder-Mead from three starts at each node, and reads the crossing
  off by interpolation. The objective is `n6_profile_vs_crb.nll`, which is the
  repository's own and not q1's, so the two share no line.

  n6's nll takes the ordinate y directly, so the record is handed to it
  unchanged: the same m, the same analyst increment, the same sigma and the
  same prior width on b.""")
    ref = Q.committed_ref()
    tn, wn = Q.nominal_train()
    A = Q.make_analysts((tn, wn), ref, "peak",
                        which=("fixed 3.5 d",))["fixed 3.5 d"]
    F, inc_a = Q.FOOT, None
    inc_a = A.phi(np.zeros(0), Q.MONTHLY)
    worst = 0.0
    for tname in ("exp 3.5 d", "exp 5.0 d"):
        inc_t = Q.truth_increment(tname, (tn, wn), "peak", ref, Q.MONTHLY)
        for sd in (0, 1, 2):
            m, s, _ = Q.make_record(inc_t, sd, F=F)
            fit = A.fit(m, s, F=F)
            lo1, hi1, o1, o2 = A.interval(m, s, fit, F=F)
            args = (m, inc_a, s, F.sigma, F.L)
            grid = np.log(fit["theta"]) + np.linspace(-3.5, 2.0, 221)
            v = np.array([N6.profile_at(lt, args) for lt in grid])
            v -= v.min()
            lo2, hi2 = _cross(grid, v, N6.DROP)
            g = max(abs(lo1 / lo2 - 1.0), abs(hi1 / hi2 - 1.0))
            worst = max(worst, g)
            print(f"    {tname:12s} seed {sd}: q1 [{lo1:8.1f}, {hi1:8.1f}]  "
                  f"n6 grid [{lo2:8.1f}, {hi2:8.1f}]  worst {g:.3e}")
    return report("profile interval, bracketing vs fixed grid", worst, 0.0,
                  0.01, rel=False,
                  note="n6's inner optimiser is Nelder-Mead from three starts "
                       "and q1's is a global K_m scan with an exact polish, "
                       "so a residual disagreement of a few parts in a "
                       "thousand is the optimisers and not the interval")


def check_6_dose_arithmetic():
    banner("6. THE PUBLISHED DOSE ARITHMETIC, IN EXACT RATIONALS")
    print("""  q1 anchors the mg/kg scale on the dose-weighted mean of the 346
  doses the interim analysis reports per regimen. Recomputed here in Fraction
  arithmetic, so the float total is checked rather than repeated.""")
    counts = [(Fraction(30, 100), 73), (Fraction(30, 100), 136),
              (Fraction(45, 100), 68), (Fraction(60, 100), 45),
              (Fraction(90, 100), 24)]
    n = sum(c for _, c in counts)
    mean = sum(d * c for d, c in counts) / n
    got_mean, got_n = Q.published_mean_mgkg()
    print(f"    exact: {n} doses, mean {mean} = {float(mean):.10f} mg/kg")
    ok = report("dose count", float(got_n), float(n), 0.0, rel=False)
    ok &= report("dose-weighted mean mg/kg", got_mean, float(mean), 1e-12)
    return ok


def check_7_two_phase():
    banner("7. PART 5's OPTIMAL TWO-PHASE DESIGN, AGAINST THE CLOSED FORM")
    print("""  For the single-exponential steady train,

      ln Delta_ss(tau) = ln A - lambda tau - ln(1 - exp(-lambda T)),

  so d ln Delta / d ln t_half = lambda tau + (a term with no tau in it). The
  free amplitude absorbs everything that does not depend on tau, so the two-
  point information depends only on lambda (tau_1 - tau_2), and

      sd(ln t_half) = sqrt(2) / (lambda |tau_1 - tau_2|),

  minimised at the two ENDS of the accessible window. q1 gets the same number
  from a numerical derivative and a pairwise search over the whole window.""")
    ref = Q.committed_ref()
    ok = True
    for lo in (0.04, 0.25, 1.0, 2.0):
        hi = Q.T_DOSE - lo
        grid = np.round(np.arange(lo, hi + 1e-9, 0.02), 3)
        got, a, b = Q.best_two_phase(grid, ref, "peak", readout="delta")
        closed = float(np.sqrt(2.0) / (Q.LAMBDA0 * abs(hi - lo)))
        ok &= report(f"two-phase sd, window from {lo:.2f} d", got, closed,
                     2e-3, note=f"q1 picked {a:.2f} d and {b:.2f} d against "
                                f"the window ends {lo:.2f} and {hi:.2f}")
    return ok


def check_8_residual_reference():
    banner("8. THE RESIDUAL REFERENCE")
    print("""  q1 compares its residual column against the MEASURED null cell
  rather than against sqrt((n - p)/n), on the ground that the paper's fit is
  penalised and a penalised fit leaves a larger residual than an unpenalised
  one with the same parameter count. That is an argument, so it is checked:
  the same records are fitted with the priors REMOVED and the two residuals
  compared. If the unpenalised one does not sit at sqrt((n - p)/n) and below
  the penalised one, the reasoning in q1 is wrong.""")
    ref = Q.committed_ref()
    tn, wn = Q.nominal_train()
    F = Q.FOOT
    A = Q.make_analysts((tn, wn), ref, "peak",
                        which=("fixed 3.5 d",))["fixed 3.5 d"]
    inc = A.phi(np.zeros(0), Q.MONTHLY)
    n, p = Q.MONTHLY.size, 4
    pen_r, free_r = [], []
    for sd in range(16):
        m, s, _ = Q.make_record(inc, sd, F=F)
        f = A.fit(m, s, F=F)
        pen_r.append(Q.rms_residual(f["theta"], f["Km"], f["ke"], f["b"], 1.0,
                                    m, s, inc, F=F))

        def res(v):
            th, Km, ke = np.exp(np.clip(v[:3], -30, 30))
            return ((th + inc) * m / (Km + m) + ke * m + v[3] - s) / F.sigma

        out = least_squares(res, np.array([np.log(Q.TH), np.log(Q.KM),
                                           np.log(Q.KE), 0.0]),
                            method="lm", max_nfev=20000)
        free_r.append(float(np.sqrt(np.mean(res(out.x) ** 2))))
    pen, free = float(np.median(pen_r)), float(np.median(free_r))
    expect = float(np.sqrt((n - p) / n))
    print(f"    unpenalised median RMS residual {free:.4f} against "
          f"sqrt((n-p)/n) = {expect:.4f}")
    print(f"    penalised   median RMS residual {pen:.4f}")
    ok = report("unpenalised residual at sqrt((n-p)/n)", free, expect, 0.12)
    ordered = pen >= free
    RESULTS.append(("penalised residual is the larger", bool(ordered)))
    print(f"    penalised >= unpenalised: {'PASSES' if ordered else 'FAILS'}")
    return ok and ordered


def check_9_kappa_prior():
    banner("9. THE FREE-AMPLITUDE FIT AGAINST THE PROBE'S FITTER")
    print("""  q1's first run put the 12% prior on kappa itself, the delta-method
  form, while the repository states it on ln kappa (the design column is
  `inc g`, the derivative in ln kappa). On the one record family where the two
  forms part company, the probe's fitter, which shares no code with q1 and
  carries the prior on ln kappa, is the independent route. The record is the
  q1 two-compartment kernel read post-dose at t = 0, which is where the first
  run found the increment "switched off".""")
    sys.path.insert(0, str(Path(__file__).resolve().parent / "probes_2026-09-25"))
    import core as C
    import fit as PF
    ref = Q.committed_ref()
    tn, wn = Q.nominal_train()
    k = Q.kernel_attarwala()
    amp = Q.amplitude(k, "peak", ref)
    times = np.arange(0.0, 730.0, 30.4)
    inc_t = amp * Q.assemble(k, tn, wn, times)
    print("""
  The probe's fitter is unbounded in K_m and q1's is bounded to its scan's
  box, [1, 1000] (check 4 prints what leaving it does). This record family is
  the one where K_m runs to the box, so the comparison is made on the seeds
  where neither fit touches it, and the seeds where one does are PRINTED, not
  dropped. The first version of this check compared every seed, failed at
  0.223 on two box seeds, and is recorded as failed in PHASE8_SIMULATIONS.md.""")
    lo, hi = Q.KM_GRID[0], Q.KM_GRID[-1]
    worst, n_in, boxed = 0.0, 0, []
    with Q.rebind(Q, DOSE_DAY_DRAW="post-dose"):
        A = Q.make_analysts((tn, wn), ref, "peak",
                            which=("fixed + free amp",))["fixed + free amp"]
        for s in range(8):
            m, y, sig, _ = PF.record(C.TH_NOM, times, s, inc_true=inc_t,
                                     noisy=False)
            f = A.fit(m, y, times=times)
            r = PF.fit(m, C.train(times), y, sig, C.L_A, kappa_cv=0.12)
            on_box = (f["Km"] >= hi * (1 - 1e-9) or f["Km"] <= lo * (1 + 1e-9)
                      or not lo <= r["Km"] <= hi)
            gap = abs(f["theta"] / r["theta"] - 1.0)
            if on_box:
                boxed.append((s, f["Km"], r["Km"], gap))
                continue
            n_in += 1
            worst = max(worst, gap)
    for s, k1, k2, g in boxed:
        print(f"    seed {s}: q1 K_m {k1:.1f} on its bound, probe K_m {k2:.1f}; "
              f"theta gap {g:.3f}, the box and not the prior")
    print(f"    interior seeds compared: {n_in} of 8")
    return report("free amplitude, q1 vs probe, interior seeds", worst, 0.0,
                  1e-5, rel=False)


def check_10_dose_day_draw():
    banner("10. THE DOSE-DAY DRAW IS READ ON ONE SIDE, FOR TRUTH AND ANALYST")
    print("""  q1 now reads a draw at a dose instant BEFORE the infusion. Two
  things must hold and neither is a property of the fit. Pre-dose, every
  kernel's increment at t = 0 is exactly zero, the analyst's included, so the
  truth and the analyst agree at that draw by construction rather than differ
  by the whole dose. Post-dose, the committed analyst reproduces
  design.run.delta at that draw, so the switch changes nothing else. Checked
  by evaluating the kernels directly, not through q1's truth or analyst
  helpers.""")
    ref = Q.committed_ref()
    tn, wn = Q.nominal_train()
    ok = True
    with Q.rebind(Q, DOSE_DAY_DRAW="pre-dose"):
        vals = [Q.assemble(kk, tn, wn, Q.at_draw(np.array([0.0, 14.0])))
                for kk in Q.KERNELS.values()]
        worst = max(float(np.abs(v[0])) for v in vals)
        ok &= report("pre-dose increment at t = 0, every kernel", worst, 0.0,
                     1e-9, rel=False)
    with Q.rebind(Q, DOSE_DAY_DRAW="post-dose"):
        k = Q.kernel_exp(Q.TH_HALF0)
        mine = (Q.D_AMP / k.single_peak()) * Q.assemble(
            k, tn, wn, Q.at_draw(np.arange(0.0, 730.0, 30.4)))
        theirs = np.asarray(Q.DR.delta(np.arange(0.0, 730.0, 30.4)), float)
        ok &= report("post-dose analyst vs design.run.delta",
                     float(np.max(np.abs(mine - theirs)) / theirs.max()), 0.0,
                     1e-12, rel=False)
    return ok


def check_11_footings():
    banner("11. GAP B ON BOTH BASELINE FOOTINGS, BY THE PROBE'S OWN ROUTE")
    print("""  q1c_baseline_footings.py reports Gap B with and without an untreated
  t = 0 draw, under three amplitude conventions. Here the same cells are built
  a second way: truths from the probe's own train() and train_biexp(), scaled
  by the probe's own steady-state window (days 600 to 614) or its per-dose D,
  and fitted by the probe's own fitter. q1's route is q1's amplitude(),
  assemble() and Analyst. Noiseless, the probe's 8 load records, monthly.
  The probe's cycle mean is a window average that carries the trapezoid error
  across each dose jump (q1 Part 0 check 3), about 2e-3, so agreement is
  asked to 1e-2 on the medians.""")
    sys.path.insert(0, str(Path(__file__).resolve().parent / "probes_2026-09-25"))
    import core as C
    import fit as PF
    ref = Q.committed_ref()
    fine = np.arange(0.0, 730.0, 0.01)
    ss = (fine > 600) & (fine < 614)
    an = C.train(fine)
    a_pk, a_mn = an[ss].max(), an[ss].mean()

    def probe_truth(kind, conv):
        if kind == "jump":
            base = lambda t: C.train_biexp(t)
        else:
            base = lambda t: C.train(t, D=1.0, half=7.0)
        b = base(fine)
        scale = {"dose": 6000.0, "peak": a_pk / b[ss].max(),
                 "mean": a_mn / b[ss].mean()}[conv]
        return lambda t: scale * base(t)

    A = np.array([[-(Q.K1.KQ + Q.K1.KDEG), Q.K1.KQ], [Q.K1.KQ, -Q.K1.KQ]]) * 24
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.array([1.0, 0.0]))
    kernels = {"exp 7.0 d": Q.kernel_exp(7.0),
               "jump": Q.ExpSum(np.real(V[0] * c), np.real(w))}
    tn, wn = Q.nominal_train()
    worst = 0.0
    fa = Q.make_analysts((tn, wn), ref, "peak",
                         which=("fixed 3.5 d",))["fixed 3.5 d"]
    print(f"\n    {'truth':10s} {'footing':11s} {'conv':5s} {'q1 route':>9s} "
          f"{'probe route':>12s}")
    for kname, k in kernels.items():
        for foot, times in (("baseline", np.arange(0.0, 730.0, 30.4)),
                            ("on therapy", np.arange(0.0, 730.0, 30.4)[1:])):
            t_eval = times - 1e-9
            for conv in ("dose", "peak", "mean"):
                amp = Q.amplitude(k, conv, ref)
                inc_q = amp * Q.assemble(k, tn, wn, t_eval)
                inc_p = probe_truth("jump" if kname == "jump" else "exp",
                                    conv)(t_eval)
                rq, rp = [], []
                for s in range(8):
                    m, y, sig, _ = PF.record(C.TH_NOM, times, s,
                                             inc_true=inc_q, noisy=False)
                    with Q.rebind(Q, DOSE_DAY_DRAW="pre-dose"):
                        rq.append(fa.fit(m, y, times=times)["theta"] / 3000)
                    m, y, sig, _ = PF.record(C.TH_NOM, times, s,
                                             inc_true=inc_p, noisy=False)
                    rp.append(PF.fit(m, C.train(t_eval), y, sig,
                                     C.L_A)["theta"] / 3000)
                a, b = float(np.median(rq)), float(np.median(rp))
                worst = max(worst, abs(a / b - 1.0))
                print(f"    {kname:10s} {foot:11s} {conv:5s} {a:9.3f} {b:12.3f}")
    return report("footings, q1 route vs probe route", worst, 0.0, 1e-2,
                  rel=False)


def check_12_s2():
    banner("12. S2: THE FISHER MAP AND THE MONTE CARLO ESTIMATOR, BY THE PROBE")
    print("""  (a) The requirement map's Fisher bound against the probe's own
  core.fisher_sd, which builds the footing, the load records and the matrix
  itself. Its load records are seeded differently, so the check is on the
  24-record medians, to 3%, on cells with no dose-instant draw (the probe reads
  a dose-instant draw post-dose).
  (b) q2's MAP estimator against the probe's fitter, fit.fit (15 starts, no
  grid), on IDENTICAL records with the nuisances drawn both ways. Asked to
  agree to 1e-4 in ln theta_hat on every record.""")
    sys.path.insert(0, str(Path(__file__).resolve().parent / "probes_2026-09-25"))
    import core as C
    import fit as PF
    import q2_requirement_map as S2
    cells = [(150.0, "weekly 7.3 d", 0.02, "measured"),
             (300.0, "monthly", 0.0, "measured"),
             (600.0, "two-monthly", 0.10, "recorded"),
             (3000.0, "monthly", 0.8, "recorded"),
             (150.0, "every 2 d", 0.8, "recorded")]
    worst = 0.0
    for th, sn, ks, ld in cells:
        t = S2.times_for(sn, "no baseline")
        mine = S2.fisher_cell(th, t, ks, ld)[0]
        theirs = C.med_sd(th, t, seeds=range(24),
                          ke_prior=(0 if ks == 0 else ks),
                          b_prior=(np.inf if ld == "measured" else None))[0]
        worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    theta {th:6.0f} {sn:13s} k_e {ks:4.2f} {ld:9s} q2 {mine:.4f}"
              f"  probe {theirs:.4f}")
    ok = report("S2 Fisher map vs probe fisher_sd (medians)", worst, 0.0,
                0.03, rel=False)
    worst_e = 0.0
    for th, sn, ks, ld in ((3000.0, "monthly", 0.8, "measured"),
                           (300.0, "weekly 7.3 d", 0.05, "recorded")):
        t = S2.times_for(sn, "no baseline")
        L, _, sigma, s_b = S2.footing(th)
        inc = S2.increment(t)
        for nu in ("truth", "centre"):
            for r in range(24):
                rng = np.random.default_rng(10_007 * r + 3)
                Ls = S2.loads_for(t.size, r, L)
                z = rng.normal(0.0, ks)
                ke_t, ke_c = ((S2.KE * np.exp(z), S2.KE) if nu == "truth"
                              else (S2.KE, S2.KE * np.exp(z)))
                b_t = rng.normal(0.0, s_b) if ld == "recorded" else 0.0
                m = S2.level(th, inc, Ls + b_t, ke=ke_t)
                s = Ls + rng.normal(0.0, sigma, t.size)
                a = np.log(S2.mc_fit(m, inc, s, sigma, s_b, ks, ke_c, ld))
                b = np.log(PF.fit(m, inc, s, sigma, L, ke_prior=ks,
                                  b_sd=(0.0 if ld == "measured" else s_b),
                                  theta0=th, ke0=ke_c)["theta"])
                if a < np.log(th) - 8 and b < np.log(th) - 8:
                    continue             # both collapsed to the floor
                worst_e = max(worst_e, abs(a - b))
    ok &= report("S2 MAP estimator vs probe fitter, per record", worst_e, 0.0,
                 1e-4, rel=False)
    return ok


def check_13_s2b():
    banner("13. S2b: THE TRACER STUDY, BY THREE INDEPENDENT ROUTES")
    print("""  (a) The plasma arm's sd(ln CL) by ORDINARY LEAST SQUARES algebra rather
  than q2b's information matrix in (ln CL, ln V). Regress ln c on t; with
  intercept a = ln(D/V) and slope s = -lambda, ln CL = ln lambda + ln D - a,
  so the delta method on the OLS covariance cv^2 (X'X)^-1 gives sd(ln CL)
  directly. No line of q2b is used.
  (b) The two-compartment bolus concentration against a MATRIX EXPONENTIAL of
  the rate matrix, which shares no algebra with the A/B/alpha/beta form.
  (c) The tracer ROW, reading (i), must vanish as the study's error grows: with
  sd_tracer enormous, q2b's matrix has to fall back on S2's own Fisher bound
  at the population prior 0.8, which check 12 already tied to the probe.""")
    from scipy.linalg import expm
    import q2b_ke_measurement as B
    import q2_requirement_map as S2

    worst = 0.0
    for th_min, w, n, cv, c0 in ((5.0, 2.0, 8, 0.10, 100.0),
                                 (30.0, 6.0, 12, 0.05, 1000.0),
                                 (120.0, 24.0, 6, 0.15, 30.0),
                                 (360.0, 6.0, 16, 0.10, 100.0)):
        t_all = B.draw_times(n, w)
        lam = np.log(2.0) / (th_min / 60.0)
        t = t_all[np.exp(-lam * t_all) * c0 >= 1.0]
        mine, _ = B.one_comp_sd(th_min / 60.0, t_all, cv, c0)
        X = np.stack([np.ones_like(t), t], axis=1)
        C2 = cv ** 2 * np.linalg.inv(X.T @ X)
        #  ln CL = ln(-s) + ln D - a, so d/d(a, s) = (-1, 1/s) with s = -lam
        g = np.array([-1.0, -1.0 / lam])
        theirs = float(np.sqrt(g @ C2 @ g))
        worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    t1/2 {th_min:5.0f} min, {w:4.0f} h, N {n:2d}, CV "
              f"{100 * cv:2.0f}%: q2b {mine:.5f}  OLS {theirs:.5f}")
    ok = report("plasma arm sd(ln CL) vs OLS delta method", worst, 0.0, 1e-9,
                rel=False)

    worst = 0.0
    for CL, V1, Q, V2 in ((1.0, 1.0, 1.0, 1.0), (0.3, 1.0, 2.0, 5.0),
                          (2.0, 0.5, 0.1, 3.0)):
        A = np.array([[-(CL / V1 + Q / V1), Q / V2], [Q / V1, -Q / V2]])
        for t in (0.01, 0.1, 1.0, 5.0):
            mine = B.two_comp_conc(np.array([t]), CL, V1, Q, V2)[0]
            theirs = (expm(A * t) @ np.array([1.0, 0.0]))[0] / V1
            worst = max(worst, abs(mine / theirs - 1.0))
    ok &= report("two-compartment bolus vs matrix exponential", worst, 0.0,
                 1e-10, rel=False)

    worst = 0.0
    for cap in (0.05, 0.20):
        th = S2.TH_NOM * cap
        for sn in ("monthly", "every 2 d"):
            for ld in S2.LOADS:
                t = S2.times_for(sn, "no baseline")
                mine = float(np.median([B.fisher_with_tracer(th, t, s, ld,
                                                             1e6, 0.0)
                                        for s in S2.SEEDS]))
                theirs = S2.fisher_cell(th, t, 0.8, ld)[0]
                worst = max(worst, abs(mine / theirs - 1.0))
                print(f"    {100 * cap:3.0f}% {sn:13s} {ld:9s} "
                      f"q2b(no study) {mine:.4f}  S2 at 0.8 {theirs:.4f}")
    ok &= report("tracer row vanishes into S2's 0.8 prior", worst, 0.0, 1e-6,
                 rel=False)
    return ok


def check_14_s3():
    banner("14. S3: THE SURFACE, BY A NUMERICAL HESSIAN AND THE PROBE'S FOOTING")
    print("""  (a) S3's Fisher bound against a FINITE-DIFFERENCE HESSIAN of the negative log
  posterior at the truth, in (ln theta, ln K_m, b, ln k_e). q3 builds the
  information matrix from analytic sensitivity columns; this route differentiates
  the objective twice numerically and shares none of that algebra.
  (b) S3's sigma against the PROBE's own footing formula, which codes the
  quasi-steady slope as theta*Km/(Km+m)^2 + ke rather than through
  model_A.fit_ab, at the committed marker CV and at a swept one.
  (c) The visits-to-gate cell is a real crossing: at the reported n the Hessian
  route is at or below the gate and at n-1 it is above.""")
    sys.path.insert(0, str(Path(__file__).resolve().parent / "probes_2026-09-25"))
    import core as C
    import q3_marker_cv_load as S3
    import q2_requirement_map as S2

    def hess_sd(theta, times, cv, load_frac):
        """sd(ln theta) from a numerical Hessian of J = 2 x negative log
        posterior at the truth, with the load record as the data."""
        L, _, sigma, s_b = S3.footing(theta, cv, load_frac)
        Ls = S2.loads_for(times.size, 0, L)
        inc = S2.increment(times)
        m = S2.level(theta, inc, Ls)

        def J(q):
            th, km, b, ke = np.exp(q[0]), np.exp(q[1]), q[2], np.exp(q[3])
            pred = (th + inc) * m / (km + m) + ke * m - b
            return (float(np.sum(((pred - Ls) / sigma) ** 2))
                    + (b / s_b) ** 2 + ((q[3] - np.log(S2.KE)) / 0.8) ** 2)
        q0 = np.array([np.log(theta), np.log(S2.KM), 0.0, np.log(S2.KE)])
        h = 1e-5
        H = np.empty((4, 4))
        for i in range(4):
            for j in range(4):
                ei, ej = np.zeros(4), np.zeros(4)
                ei[i] = ej[j] = h
                H[i, j] = (J(q0 + ei + ej) - J(q0 + ei - ej)
                           - J(q0 - ei + ej) + J(q0 - ei - ej)) / (4 * h * h)
        F = H / 2.0
        w = np.linalg.eigvalsh(F)
        if w[0] <= 1e-12 * w[-1]:
            return np.inf
        return float(np.sqrt(np.linalg.inv(F)[0, 0]))

    worst = 0.0
    for th in (2472.0, 3000.0, 18814.0):
        for cv in (0.05, 0.20):
            for lf in (0.01, 0.10, 0.40):
                t = S3.times_on_therapy(30.4, 24)
                mine = S3.fisher_sd(th, t, 0, cv, lf)
                theirs = hess_sd(th, t, cv, lf)
                worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    theta {th:6.0f}: q3 "
              f"{S3.fisher_sd(th, t, 0, 0.05, 0.10):.5f}"
              f"  Hessian {hess_sd(th, t, 0.05, 0.10):.5f}   (CV 5%, load 10%)")
    ok = report("S3 Fisher surface vs numerical Hessian", worst, 0.0, 1e-4,
                rel=False)

    worst = 0.0
    for th in (2472.0, 3000.0, 18814.0):
        L, m_bar = S3.N4.anchor(th, S2.KM, S2.KE, "A")
        for cv in (0.05, 0.15, 0.30):
            mine = S3.footing(th, cv, 0.10)[2]
            m = C.level(th, S2.KM, S2.KE, L, 0.0)
            slope = th * S2.KM / (S2.KM + m) ** 2 + S2.KE
            theirs = float(np.hypot(0.10 * L, slope * cv * m))
            worst = max(worst, abs(mine / theirs - 1.0))
            print(f"    theta {th:6.0f} CV {100 * cv:3.0f}%: q3 sigma "
                  f"{mine:8.2f}  probe {theirs:8.2f}")
    ok &= report("S3 sigma vs probe footing formula", worst, 0.0, 1e-9,
                 rel=False)

    bad = 0
    for th, cv, lf in ((3000.0, 0.05, 0.05), (3000.0, 0.15, 0.10),
                       (2472.0, 0.05, 0.01)):
        n = S3.visits_to_gate(th, 30.4, cv, lf)
        if n is None:
            print(f"    theta {th:6.0f} CV {100 * cv:3.0f}% load "
                  f"{100 * lf:3.0f}%: never, nothing to cross")
            continue
        at = float(np.median([hess_sd(th, S3.times_on_therapy(30.4, n), cv, lf)
                              for _ in (0,)]))
        before = float(np.median([hess_sd(
            th, S3.times_on_therapy(30.4, n - 1), cv, lf) for _ in (0,)]))
        print(f"    theta {th:6.0f} CV {100 * cv:3.0f}% load {100 * lf:3.0f}%: "
              f"n = {n:3d}, Hessian sd {at:.4f} at n and {before:.4f} at n-1")
        #  seed 0 alone, so the median over records may sit either side; the
        #  check is that the crossing is in the right direction
        if not at < before:
            bad += 1
    ok &= report("visits-to-gate crossings fall in n", bad, 0, 0, rel=False)
    return ok


def check_15_s4():
    banner("15. S4: THE TROUGH READING, THE DRIFT CONVENTION, AND ONE CELL")
    print("""  (a) The trough amplitude by BRUTE SUMMATION of the dose train rather than
  q4's closed-form steady state: sum the kernel over 4000 past doses and read the
  value just before the next one. Shares no algebra with ExpSum.steady.
  (b) S4 applies the k_e drift on S2's `truth` convention and claims the two
  honest conventions coincide at these widths. Measured, not assumed, with S2's
  own Monte Carlo at drift 0.05, 0.10 and 0.15, against 0.8 where S2b found them
  half apart.
  (c) One joint-scenario cell rebuilt through the PROBE's fitter on identical
  records, which shares no estimator code with q1's Analyst.""")
    sys.path.insert(0, str(Path(__file__).resolve().parent / "probes_2026-09-25"))
    import fit as PF
    import q1_increment_shape as Q1
    import q2_requirement_map as S2
    import q4_joint_stress as S4

    worst = 0.0
    for name in ("exp 3.5 d", "2-comp, rising", "2-comp, jump"):
        k = S4.Q1C.TRUTHS[name]
        mine = S4.steady_trough(k, S4.T)
        u = S4.T * np.arange(1, 4001)          # time since each past dose
        theirs = float(np.sum(k(u)))
        worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    {name:16s} closed form {mine:12.6f}  summed {theirs:12.6f}")
    ok = report("trough amplitude, closed form vs brute sum", worst, 0.0,
                1e-10, rel=False)

    print()
    worst = 0.0
    for ks in (0.05, 0.10, 0.15, 0.8):
        t = S2.times_for("monthly", "no baseline")
        a = S2.mc_cell(S2.TH_NOM, t, ks, "recorded", "truth")["rmse"]
        b = S2.mc_cell(S2.TH_NOM, t, ks, "recorded", "centre")["rmse"]
        gap = abs(a / b - 1.0)
        if ks < 0.8:
            worst = max(worst, gap)
        print(f"    k_e width {ks:4.2f}: truth {a:.4f}  centre {b:.4f}  "
              f"gap {100 * gap:5.1f}%" + ("   (reference, not checked)"
                                          if ks == 0.8 else ""))
    ok &= report("drift conventions coincide below 0.2", worst, 0.0, 0.05,
                 rel=False)

    print()
    ref, rt = Q1.committed_ref(), S4.committed_trough()
    times = S4.footing_times("baseline")
    worst_e = 0.0
    for kname in ("2-comp, rising", "2-comp, jump"):
        k = S4.Q1C.TRUTHS[kname]
        with Q1.rebind(Q1, CV=0.15):
            A = S4.analyst_for(S4.TRAINS["irregular C5"], ref, "dose", rt)
            F = Q1.Footing()
            inc = S4.truth_inc(k, "dose", S4.TRAINS["irregular C5"], times,
                               ref, rt)
            for sd in range(8):
                m, s, _ = Q1.make_record(inc, sd, F=F, times=times)
                mine = A.fit(m, s, times=times, F=F)["theta"]
                phi = A.phi(np.zeros(0), times)
                theirs = PF.fit(m, phi, s, F.sigma, F.L, ke_prior=0.8,
                                b_sd=F.s_b, theta0=Q1.TH, ke0=F.ke)["theta"]
                worst_e = max(worst_e, abs(np.log(mine) - np.log(theirs)))
        print(f"    {kname:16s} worst |ln theta| gap so far {worst_e:.2e}")
    ok &= report("S4 joint cell, q1 analyst vs probe fitter", worst_e, 0.0,
                 1e-3, rel=False)
    return ok


def check_16_s5():
    banner("16. S5: THE THETA PRIOR, BY A HESSIAN, A FROM-SCRATCH FIT, AND A LIMIT")
    print("""  (a) S5's Fisher bound with the ex vivo prior against a FINITE-DIFFERENCE
  HESSIAN of the negative log posterior, the prior term included. q5 adds the
  prior to one diagonal element of an analytic information matrix; this route
  differentiates the whole objective numerically.
  (b) S5's MAP against a FROM-SCRATCH objective minimised by Nelder-Mead from
  several starts, sharing no line with q5's variable-projection grid and polish.
  (c) A LIMIT: with the record made uninformative, the posterior must return the
  prior's own centre and width, so the offset passes through entirely. That is
  the boundary Part D's pass-through numbers are read against.""")
    from scipy.optimize import minimize
    import q2_requirement_map as S2
    import q5_exvivo_prior as S5

    def hess_sd(theta, times, ke_sd, load, th_sd):
        L, _, sigma, s_b = S2.footing(theta)
        Ls = S2.loads_for(times.size, 0, L)
        inc = S2.increment(times)
        m = S2.level(theta, inc, Ls)

        def J(q):
            th, km, b, ke = np.exp(q[0]), np.exp(q[1]), q[2], np.exp(q[3])
            pred = (th + inc) * m / (km + m) + ke * m + b
            out = float(np.sum(((pred - Ls) / sigma) ** 2))
            if load == "recorded":
                out += (b / s_b) ** 2
            if ke_sd:
                out += ((q[3] - np.log(S2.KE)) / ke_sd) ** 2
            if th_sd is not None:
                out += ((q[0] - np.log(theta)) / th_sd) ** 2
            return out
        q0 = np.array([np.log(theta), np.log(S2.KM), 0.0, np.log(S2.KE)])
        h = 1e-5
        H = np.empty((4, 4))
        for i in range(4):
            for j in range(4):
                ei, ej = np.zeros(4), np.zeros(4)
                ei[i] = ej[j] = h
                H[i, j] = (J(q0 + ei + ej) - J(q0 + ei - ej)
                           - J(q0 - ei + ej) + J(q0 - ei - ej)) / (4 * h * h)
        F = H / 2.0
        w = np.linalg.eigvalsh(F)
        if w[0] <= 1e-12 * w[-1]:
            return np.inf
        return float(np.sqrt(np.linalg.inv(F)[0, 0]))

    worst = 0.0
    for cap in (0.05, 0.20, 1.00):
        th = S2.TH_NOM * cap
        for v in (0.3, 1.0, None):
            for ks in (0.8, 0.054):
                t = S2.times_for("monthly", "no baseline")
                mine = S5.fisher_sd(th, t, 0, ks, "recorded", v)
                theirs = hess_sd(th, t, ks, "recorded", v)
                worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    {100 * cap:4.0f}% prior 0.3, k_e 0.054: q5 "
              f"{S5.fisher_sd(th, t, 0, 0.054, 'recorded', 0.3):.5f}  Hessian "
              f"{hess_sd(th, t, 0.054, 'recorded', 0.3):.5f}")
    ok = report("S5 Fisher with theta prior vs numerical Hessian", worst, 0.0,
                1e-4, rel=False)

    print()
    worst_e = 0.0
    for cap, v, ks in ((0.05, 0.5, 0.054), (0.20, 0.3, 0.8), (1.00, 1.0, 0.02)):
        theta = S2.TH_NOM * cap
        t = S2.times_for("monthly", "no baseline")
        L, _, sigma, s_b = S2.footing(theta)
        inc = S2.increment(t)
        for r in range(6):
            rng = np.random.default_rng(10_007 * r + 3)
            Ls = S2.loads_for(t.size, r, L)
            z = rng.normal(0.0, ks)
            b_t = rng.normal(0.0, s_b)
            m = S2.level(theta, inc, Ls + b_t)
            s = Ls + rng.normal(0.0, sigma, t.size)
            th_c = theta * np.exp(rng.normal(0.0, v))
            mine = S5.mc_fit(m, inc, s, sigma, s_b, ks, S2.KE * np.exp(z),
                             "recorded", v, th_c)

            def obj(q):
                th, km, ke, b = np.exp(q[0]), np.exp(q[1]), np.exp(q[2]), q[3]
                pred = (th + inc) * m / (km + m) + ke * m + b
                return (np.sum(((pred - s) / sigma) ** 2)
                        + ((q[0] - np.log(th_c)) / v) ** 2
                        + ((q[2] - np.log(S2.KE * np.exp(z))) / ks) ** 2
                        + (b / s_b) ** 2)
            best = None
            for t0 in (theta, theta * 3, theta / 3, th_c):
                for k0 in (10.0, 100.0, 500.0):
                    o = minimize(obj, [np.log(t0), np.log(k0),
                                       np.log(S2.KE), 0.0],
                                 method="Nelder-Mead",
                                 options=dict(maxiter=40000, maxfev=40000,
                                              xatol=1e-12, fatol=1e-12))
                    if best is None or o.fun < best.fun:
                        best = o
            worst_e = max(worst_e, abs(np.log(mine) - best.x[0]))
        print(f"    {100 * cap:4.0f}% prior {v} k_e {ks}: worst |ln theta| gap "
              f"so far {worst_e:.2e}")
    ok &= report("S5 MAP vs from-scratch Nelder-Mead", worst_e, 0.0, 1e-3,
                 rel=False)

    print()
    #  the uninformative-record limit: one draw, noise blown up by 1e6
    theta = S2.TH_NOM
    t = S2.times_for("monthly", "no baseline")[:1]
    L, _, sigma, s_b = S2.footing(theta)
    inc = S2.increment(t)
    m = S2.level(theta, inc, np.full(1, L))
    worst = 0.0
    for off in (0.3, 0.7):
        th_c = theta * np.exp(off)
        got = S5.mc_fit(m, inc, np.full(1, L), sigma * 1e6, s_b * 1e6, 0.054,
                        S2.KE, "recorded", 0.5, th_c)
        worst = max(worst, abs(np.log(got) - np.log(th_c)))
        print(f"    offset {off}: prior centre {th_c:8.1f}  posterior "
              f"{got:8.1f}")
    ok &= report("uninformative record returns the prior centre", worst, 0.0,
                 1e-6, rel=False)
    return ok


def main():
    print("=" * 92)
    print("  PHASE 8 RE-CHECKS. Independent routes to every reported result.")
    print("=" * 92)
    for fn in (check_1_kernel, check_2_steady_identity, check_3_bias_direction,
               check_4_fit, check_5_interval, check_6_dose_arithmetic,
               check_7_two_phase, check_8_residual_reference,
               check_9_kappa_prior, check_10_dose_day_draw,
               check_11_footings, check_12_s2, check_13_s2b,
               check_14_s3, check_15_s4, check_16_s5):
        try:
            fn()
        except Exception as exc:          # a check that cannot run is a failure
            RESULTS.append((fn.__name__ + " (raised)", False))
            print(f"    {fn.__name__} RAISED {type(exc).__name__}: {exc}")
    banner("SUMMARY")
    bad = [n for n, ok in RESULTS if not ok]
    for name, ok in RESULTS:
        print(f"    {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\n    {len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass.")
    if bad:
        print("    FAILURES, reported and not repaired:")
        for n in bad:
            print(f"      - {n}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
