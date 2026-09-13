"""
2D.6 / 2D.8 -- the decision characteristic and the cohort number, on the
therapy-channel design.

The design is the one the specification now asks for: routine marker series
sampled across dose cycles at outpatient density, with the dosing record, and
the delivered increment characterised to `kappa_cv`. Everything is simulated
through fit_ensemble, so the estimator here is the same one 1.9 always used --
what has changed is that it now carries the increment.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D

P = D.P
base = D.clearance_A(D.M_BAR, P)


def design(kind):
    if kind == "fixed point (the old specification)":
        return np.geomspace(D.M_BAR, 5 * D.M_BAR, 9), None
    t = {"across dose cycles, monthly (n=24)": np.arange(0.0, 730.0, 30.4),
         "across dose cycles, two-daily (n=365)": np.arange(0.0, 730.0, 2.0),
         "one dose cycle, every two days (n=7)": np.arange(0.0, D.DOSE_INTERVAL, 2.0)}[kind]
    return D.quasi_steady_level(t), np.asarray(D.delta(t), dtype=float)


def main():
    print(__doc__)
    kinds = ["fixed point (the old specification)",
             "one dose cycle, every two days (n=7)",
             "across dose cycles, monthly (n=24)",
             "across dose cycles, two-daily (n=365)"]

    print("=" * 96)
    print("SIGN ERROR at the decision boundary, and the cohort it implies")
    print("=" * 96)
    print(f"{'design':>40} {'kappa':>7} | {'p(err) near':>12} {'treat':>8} {'fail to':>8} | "
          f"{'cohort':>9}")
    print(f"{'':>40} {'':>7} | {'boundary':>12} {'unneeded':>8} {'treat':>8} | {'N':>9}")
    for kind in kinds:
        lv, inc = design(kind)
        for kcv in (0.0, 0.159):
            if inc is None and kcv:
                continue
            p_bad = D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=21,
                                 s_peak=3.6 * base, increments=inc,
                                 kappa_cv=kcv)["p"]
            e_safe = D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=21,
                                  s_peak=2.98 * base, increments=inc,
                                  kappa_cv=kcv)["p"]
            e_uns = D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=21,
                                 s_peak=3.58 * base, increments=inc,
                                 kappa_cv=kcv)["p"]
            N = D.cohort_size(p_bad)
            lab = kind if kcv == 0.0 else ""
            print(f"{lab:>40} {kcv:7.3f} | {p_bad:11.1%} {e_safe:8.1%} {e_uns:8.1%} | "
                  f"{('%d' % N) if np.isfinite(N) else 'infinite':>9}")

    print()
    print("=" * 96)
    print("THE COHORT AGAINST PA PREVALENCE")
    print("=" * 96)
    lv, inc = design("across dose cycles, two-daily (n=365)")
    for kcv, lab in ((0.0, "Delta known exactly"), (0.159, "Delta to 15.9%, the gate"),
                     (0.30, "Delta to 30%")):
        p = D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=21,
                         s_peak=3.6 * base, increments=inc, kappa_cv=kcv)["p"]
        N = D.cohort_size(p)
        print(f"  {lab:>26}: p(err) = {p:5.1%}   N = "
              f"{('%d' % N) if np.isfinite(N) else 'infinite'}")
    print()
    print("  Reference denominators:")
    print("    phase 1/2 mRNA-3927 trial enrolled 16 (trial_phase12)")
    print("    PA newborn incidence ~1:105,000 in the US / Europe; the highest reported")
    print("    is Saudi Arabia at roughly 1:14,000 (nbs_* rows in refs.bib)")
    for N in (16, 50, 100, 250):
        pass
    print()
    print("  What error rate a feasible cohort can tolerate, at +/-10 points on the")
    print("  fraction with mu < 0:")
    from scipy.optimize import brentq
    for N in (16, 30, 60, 120, 250, 500):
        try:
            pm = brentq(lambda p: D.cohort_size(p) - N, 1e-6, 0.499, xtol=1e-6)
            print(f"    N = {N:4d} patients  ->  per-patient sign error must be below {pm:5.1%}")
        except ValueError:
            print(f"    N = {N:4d} patients  ->  unattainable at any error rate")


if __name__ == "__main__":
    main()
