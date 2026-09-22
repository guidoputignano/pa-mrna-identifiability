"""
What the paper can still claim at a disease-range capacity, and why.

`audit/n4_theta_scope.py` established that the monthly design meets its
sd(ln theta) = 0.20 gate only for theta between 2472 and 18814, that is 82% to
627% of nominal, and that propionic acidemia sits a factor of 16 below the
lower edge. That is a scope limit on ONE estimand. This script asks what is
left, because a paper that says what can be recovered in the disease is
stronger than one that only says theta cannot.

Six questions, in the order they have to be answered.

  S1  WHY, structurally. Is the failure a one-line consequence of the balance
      rather than a fact discovered by sweeping? If it is, the reason belongs
      in the main text and the sweep is only its confirmation.
  S2  WHAT IS LEFT. The full covariance at disease-range capacity, not just
      its theta corner: K_m, k_e, the load bias, and the best-determined
      direction in the three log coordinates.
  S3  THE DESIGN RULES. Aliasing, phase coverage and recording the load are
      statements about the observable, not about theta. Do they still hold
      where theta does not?
  S4  THE DOSE. The design fails partly because Delta swamps theta. Does a
      SMALLER dose open a window, and where is the trade-off between too
      little excitation and too much swamping?
  S5  THE ALTERNATIVE CHANNEL. If S1 is right, theta is recovered as the
      residual of L minus k_e m, so a tighter k_e should matter far more in
      the disease than at nominal. That is a testable prediction and it says
      what an external measurement would have to measure.
  S6  THE COHORT STUDY. The abstract reports a median of 0.19 with just over
      half meeting the target. Where does that cohort sit relative to the band?

Footing throughout is n4's, imported rather than restated so the two cannot
drift: convention A (the recorded load held, the marker allowed to rise), which
n4 verifies is the convention Section 11's own table uses.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import n4_theta_scope as N4
from n4_theta_scope import (KM0, KE, TH0, D0, GATE, SEEDS, MONTHLY, DENSE,
                            TH_DISEASE, footing, level, elasticity, delta_at)
from design.run import LOAD_FRAC, ASSAY_CV, DOSE_INTERVAL
from common.models import ParamsA
import model_A.run as LOC

COORDS = ("ln theta", "ln K_m", "b", "ln k_e")


def covariance(theta, D=D0, Km=KM0, ke=KE, times=MONTHLY, seed=0,
               load_cv=0.10, ke_prior=0.8, load_frac=None):
    """The full 4x4 posterior covariance, not only its theta corner."""
    lf = LOAD_FRAC if load_frac is None else load_frac
    L, m_bar, _, _ = footing(theta, Km, ke, "A")
    p = ParamsA(theta=theta, Km=Km, ke=ke)
    sigma = float(LOC.sigma_total(lf, ASSAY_CV, p=p, m_bar=m_bar))
    S_B = lf * L
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    loads = (L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
             if load_cv else np.full(times.size, L))
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
        return None, S_B, sigma
    return np.linalg.inv(F), S_B, sigma


def marginals(theta, **kw):
    """Median marginal sd of each coordinate over the seed sweep."""
    out = []
    for s in SEEDS:
        P, S_B, sigma = covariance(theta, seed=s, **kw)
        if P is None:
            continue
        out.append(np.sqrt(np.diag(P)))
    if not out:
        return None, None, None
    A = np.array(out)
    P, S_B, sigma = covariance(theta, seed=SEEDS[0], **kw)
    return np.median(A, axis=0), S_B, sigma


def tight_direction(theta, **kw):
    """
    The best-determined direction among the three LOG coordinates.

    b is in flux units, so it is left out: an eigenvector mixing a log
    coordinate with a dimensional one has no units and cannot be read.
    """
    Ps = []
    for s in SEEDS:
        P, _, _ = covariance(theta, seed=s, **kw)
        if P is not None:
            Ps.append(P[np.ix_([0, 1, 3], [0, 1, 3])])
    if not Ps:
        return None, None
    P = np.median(np.array(Ps), axis=0)
    w, V = np.linalg.eigh(P)
    return float(np.sqrt(w[0])), V[:, 0]


def s1_structure():
    print("\n  S1. THE REASON, IN ONE LINE FROM THE BALANCE")
    print("      At the untreated balance theta g(m) + k_e m = L, so the column")
    print("      the information matrix carries for ln theta,")
    print("")
    print("          d(ordinate)/d(ln theta) = theta g(m) = L - k_e m ,")
    print("")
    print("      IS the flux through the enzyme. The sensitivity of the record to")
    print("      capacity is the capacity's own share of the load. In the disease")
    print("      the alternative channel carries almost all of it, so there is")
    print("      almost nothing left for theta to be read from.\n")
    print(f"  {'theta':>8} {'% nominal':>10} {'m_bar':>8} {'theta g = L - k_e m':>20} "
          f"{'enzyme share':>13} {'vs k_e column':>14}")
    worst = 0.0
    for th in (TH0, 2472.0, 1500.0, 600.0, 300.0, TH_DISEASE, 60.0):
        L, m_bar, _, _ = footing(th, convention="A")
        col = th * m_bar / (KM0 + m_bar)
        identity = L - KE * m_bar
        worst = max(worst, abs(col - identity) / L)
        print(f"  {th:8.0f} {100*th/TH0:9.0f}% {m_bar:8.2f} "
              f"{col:12.2f} = {identity:5.2f} {col/L:12.1%} {col/(KE*m_bar):14.3f}")
    print(f"\n      The identity theta g = L - k_e m holds to {worst:.2e} "
          f"relative, which is the")
    print("      point: it is an identity, not a fit. The last column is the "
          "enzymatic flux")
    print("      over the alternative flux, and it is what the estimator has to "
          "separate.")
    print("      At nominal the enzyme carries five times the alternative "
          "channel; at 5% of")
    print("      nominal it carries a fifth of it, so theta is a small "
          "difference of two")
    print("      large numbers, and k_e is itself known only through a prior.")
    return worst


def s2_whats_left():
    print("\n  S2. WHAT IS STILL IDENTIFIED AT A DISEASE-RANGE CAPACITY")
    print(f"  {'theta':>8} {'% nom':>6} " + " ".join(f"{c:>12}" for c in COORDS)
          + f" {'b prior':>9} {'b learned':>10}")
    for th in (TH0, 600.0, 300.0, TH_DISEASE, 60.0):
        md, S_B, sigma = marginals(th)
        if md is None:
            print(f"  {th:8.0f} {100*th/TH0:5.0f}%  singular")
            continue
        print(f"  {th:8.0f} {100*th/TH0:5.0f}% "
              + " ".join(f"{x:12.4f}" for x in md)
              + f" {S_B:9.2f} {S_B/md[2]:10.2f}x")
    print("\n      `b prior` is the load-record prior's own sd and `b learned` is")
    print("      how much tighter the posterior on b is than that prior, so a")
    print("      number above 1 means the data said something about the load bias.")
    print(f"\n  {'theta':>8} {'% nom':>6} {'tight sd':>9}   best-determined direction "
          f"in (ln theta, ln K_m, ln k_e)")
    for th in (TH0, 600.0, 300.0, TH_DISEASE, 60.0):
        sd, v = tight_direction(th)
        if sd is None:
            continue
        v = v / np.abs(v).max()
        print(f"  {th:8.0f} {100*th/TH0:5.0f}% {sd:9.4f}   "
              f"({v[0]:+.3f}, {v[1]:+.3f}, {v[2]:+.3f})")


def s3_design_rules():
    print("\n  S3. DO THE DESIGN RULES SURVIVE WHERE THETA DOES NOT")
    T = 730.0
    scheds = [("monthly", np.arange(0.0, T, 30.4)),
              ("two-monthly", np.arange(0.0, T, 60.0)),
              ("fortnightly = the dose interval", np.arange(0.0, T, DOSE_INTERVAL)),
              ("14.5-daily, offset half a day", np.arange(0.0, T, DOSE_INTERVAL + 0.5)),
              ("six-monthly", np.arange(180.0, T + 1.0, 180.0)),
              ("every two days", DENSE)]
    for th in (TH0, TH_DISEASE):
        print(f"\n    theta = {th:.0f} ({100*th/TH0:.0f}% of nominal)")
        print(f"    {'schedule':>32} {'n':>4} " +
              " ".join(f"{c:>11}" for c in COORDS) + f" {'tight':>8}")
        base = {}
        for nm, tt in scheds:
            md, _, _ = marginals(th, times=tt)
            sd, _ = tight_direction(th, times=tt)
            if md is None:
                print(f"    {nm:>32} {tt.size:4d}   singular")
                continue
            base[nm] = md
            print(f"    {nm:>32} {tt.size:4d} "
                  + " ".join(f"{x:11.4f}" for x in md)
                  + f" {sd:8.4f}")
        if "fortnightly = the dose interval" in base:
            a = base["fortnightly = the dose interval"]
            o = base["14.5-daily, offset half a day"]
            print(f"      aliasing cost, coordinate by coordinate: "
                  + ", ".join(f"{c} {a[i]/o[i]:.2f}x"
                              for i, c in enumerate(COORDS)))


def s4_dose_sweep():
    print("\n  S4. DOES A SMALLER DOSE OPEN A WINDOW AT DISEASE-RANGE CAPACITY")
    print("      Too small and the therapy channel supplies no excitation; too")
    print("      large and Delta swamps theta in the (theta + Delta) g term.")
    for th in (TH_DISEASE, 300.0, TH0):
        print(f"\n    theta = {th:.0f} ({100*th/TH0:.0f}% of nominal), "
              f"ceiling-admissible D <= {N4.d_ceiling(th):.0f}")
        print(f"    {'D':>9} {'D/theta':>9} {'Delta_pk/theta':>15} "
              f"{'sd(ln theta) med':>17} {'spread':>21} {'admissible':>11}")
        best = (np.inf, None)
        for D in (10.0, 30.0, 100.0, 300.0, 600.0, 1000.0, 2000.0, 2672.0,
                  4000.0, 6000.0, 12000.0, 30000.0):
            med, lo, hi = N4.over_seeds(theta=th, D=D, convention="A")
            pk = N4.peak_increment(D)
            adm = "yes" if th + pk <= th / N4.REFERENCE_FRAC else "no"
            if med < best[0]:
                best = (med, D)
            print(f"    {D:9.0f} {D/th:9.2f} {pk/th:15.2f} {med:17.4f} "
                  f"[{lo:9.4f},{hi:9.4f}] {adm:>11}")
        print(f"      best over this sweep: {best[0]:.4f} at D = {best[1]:.0f}"
              f"{'  <= REACHES THE GATE' if best[0] <= GATE else ''}")


def s5_ke_prior():
    print("\n  S5. WOULD KNOWING THE ALTERNATIVE CLEARANCE RESCUE CAPACITY")
    print("      S1 predicts this matters far more in the disease than at nominal,")
    print("      because there theta is the small residual of L minus k_e m.")
    print(f"    {'s(ln k_e)':>10} {'factor/sd':>10} " +
          "  ".join(f"theta={th:.0f}".rjust(16) for th in (TH0, 300.0, TH_DISEASE)))
    for s in (0.8, 0.4, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001):
        cells = []
        for th in (TH0, 300.0, TH_DISEASE):
            med, _, _ = N4.over_seeds(theta=th, D=D0, convention="A", ke_prior=s)
            cells.append(f"{med:16.4f}")
        print(f"    {s:10.3f} {np.exp(s):10.2f} " + "  ".join(cells))
    print("\n      The last row is k_e known to a tenth of a per cent, which no")
    print("      measurement supplies; it is there to show the limit rather than")
    print("      to propose anything.")


def s6_cohort():
    print("\n  S6. WHERE THE COHORT STUDY SITS RELATIVE TO THE BAND")
    lo_mult, hi_mult = 1.0 / 3.0, 3.0
    lo, hi = TH0 * lo_mult, TH0 * hi_mult
    band_lo, band_hi = 2472.0, 18814.0
    print(f"      audit/e8_cohort_unblocked.py draws theta log-uniform over a "
          f"9-fold span,")
    print(f"      theta in [{lo:.0f}, {hi:.0f}], that is "
          f"[{100*lo/TH0:.0f}%, {100*hi/TH0:.0f}%] of nominal.")
    print(f"      The band where the monthly design meets its gate is "
          f"[{band_lo:.0f}, {band_hi:.0f}].")
    frac_below = (np.log(band_lo) - np.log(lo)) / (np.log(hi) - np.log(lo))
    print(f"      Log-uniform, the fraction of that cohort BELOW the band's "
          f"lower edge is {frac_below:.1%}.")
    print(f"      The abstract reports just over half the cohort meeting the "
          f"target, and")
    print(f"      attributes it to the dose being fixed while capacity is not. "
          f"On these")
    print(f"      numbers the band edge accounts for {frac_below:.0%} of the "
          f"misses on its own.")
    print(f"      The disease range, theta <= {TH_DISEASE:.0f}, is outside the "
          f"cohort's span entirely.")
    return frac_below


def recheck():
    """
    Two routes that share no code with the above.

    R1  the covariance assembled here against n4's scalar sd(ln theta), which
        builds its own information matrix. They agree or one of them is wrong.
    R2  the S1 identity theta g(m) = L - k_e m re-derived from a bracketed root
        solve of the balance rather than from the closed form.
    """
    from scipy.optimize import brentq
    print("\n  RE-CHECKS")
    worst_cov = 0.0
    for th in (TH0, 600.0, TH_DISEASE):
        for s in (0, 7, 19):
            P, _, _ = covariance(th, seed=s)
            mine = np.sqrt(P[0, 0])
            theirs = N4.sd_ln_theta(theta=th, D=D0, convention="A", seed=s)
            worst_cov = max(worst_cov, abs(mine - theirs) / theirs)
    worst_id = 0.0
    for th in (TH0, 600.0, 300.0, TH_DISEASE, 60.0):
        L, _, _, _ = footing(th, convention="A")
        m = brentq(lambda x: th * x / (KM0 + x) + KE * x - L, 1e-9, 1e9,
                   xtol=1e-14, rtol=8.9e-16)
        worst_id = max(worst_id, abs(th * m / (KM0 + m) - (L - KE * m)) / L)
    print(f"    R1 this file's covariance against n4's independent build: "
          f"worst relative {worst_cov:.2e} against 1e-12")
    print(f"    R2 the S1 identity from a bracketed root rather than the closed "
          f"form: worst {worst_id:.2e} against 1e-14")
    ok = worst_cov < 1e-12 and worst_id < 1e-14
    print(f"    RE-CHECKS {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print("=" * 100)
    print("  WHAT THE PAPER CAN STILL CLAIM AT A DISEASE-RANGE CAPACITY")
    print("=" * 100)
    s1_structure()
    s2_whats_left()
    s3_design_rules()
    s4_dose_sweep()
    s5_ke_prior()
    s6_cohort()
    ok = recheck()
    print("\n" + "=" * 100)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
