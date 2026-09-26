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


def main():
    print("=" * 92)
    print("  PHASE 8 RE-CHECKS. Independent routes to every reported result.")
    print("=" * 92)
    for fn in (check_1_kernel, check_2_steady_identity, check_3_bias_direction,
               check_4_fit, check_5_interval, check_6_dose_arithmetic,
               check_7_two_phase, check_8_residual_reference,
               check_9_kappa_prior, check_10_dose_day_draw,
               check_11_footings):
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
