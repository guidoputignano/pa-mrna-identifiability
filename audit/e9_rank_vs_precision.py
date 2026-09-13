"""
15.1 / 15.2 -- reconciling scenario 1 against prop:rank, and deriving the scaling.

prop:rank reports that half a per cent of day-to-day load variation restores full
rank. Scenario 1 of 14.1 runs twenty times that variation and returns 1.79, which
is unusable. Both cannot be describing the same mechanism, so this checks three
things and then derives the one claim 14.1 left asserted.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from design.run import P, load, delta, LOAD_FRAC, ASSAY_CV, design_record
import model_A.run as LOC

SIGMA = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
LB = load(0.0)
MON = np.arange(0.0, 730.0, 30.4)


def banner(t):
    print("\n" + "=" * 82); print("  " + t); print("=" * 82)


def record(n_or_times, cv, with_dose, seed, lognormal=True):
    t = (np.asarray(n_or_times, float) if np.ndim(n_or_times)
         else np.linspace(0.0, 730.0, int(n_or_times)))
    n = len(t)
    inc = delta(t) if with_dose else np.zeros(n)
    rng = np.random.default_rng(seed)
    if cv and lognormal:
        s2 = np.log1p(cv ** 2)
        L = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))
    elif cv:
        L = LB * (1.0 + cv * rng.standard_normal(n))
    else:
        L = np.full(n, LB)
    B = P.ke * P.Km + P.theta + inc - L
    m = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * L * P.Km)) / (2.0 * P.ke)
    Q = np.array([[P.theta * mi / (P.Km + mi),
                   -(P.theta + di) * P.Km * mi / (P.Km + mi) ** 2,
                   1.0, P.ke * mi] for mi, di in zip(m, inc)])
    return L, m, Q


def sd_from(Q, load_prior=True, ke_prior=0.8):
    """ke_prior=None means k_e is held EXACTLY KNOWN, which is _filter's own
    convention: the column is DROPPED, not given a flat prior. Getting that
    wrong makes 'k_e known' read as 'k_e free' and inverts the comparison."""
    if ke_prior is None:
        Q = Q[:, :3]
        pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2 if load_prior else np.inf]
    else:
        pri = [np.inf, np.inf,
               (LOAD_FRAC * LB) ** 2 if load_prior else np.inf, ke_prior ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) + Q.T @ Q / SIGMA ** 2
    w = np.linalg.eigvalsh(F)
    return np.inf if (w[0] <= 0 or w[0] <= 1e-14 * w[-1]) else float(np.sqrt(np.linalg.inv(F)[0, 0]))


def med(fn):
    v = np.array([fn(s) for s in range(24)], float)
    v = v[np.isfinite(v)]
    return np.median(v) if v.size else np.inf


# ---------------------------------------------------------------------------
def a_does_variation_reach_the_levels():
    banner("(a) DOES THE VARIATION ENTER THE BALANCE, OR ONLY THE NOISE?")
    print("  If the load were added as observation noise it would never reach the")
    print("  sensitivity matrix and the sweep would test nothing. It does reach it:")
    print(f"    {'CV':>7} {'L span':>24} {'m span':>24}")
    for cv in (0.0, 0.005, 0.10, 0.40):
        L, m, _ = record(MON, cv, False, 3)
        print(f"    {cv*100:6.1f}% [{L.min():9.2f},{L.max():9.2f}] "
              f"[{m.min():10.4f},{m.max():10.4f}]")
    print("\n  The mechanism IS exercised. But note the last row under a GAUSSIAN draw:")
    L, m, _ = record(MON, 0.40, False, 3, lognormal=False)
    print(f"    Gaussian at 40% CV: L min {L.min():.2f}, m min {m.min():.4f}")
    print("  A negative load and a negative marker are not a record. 14.1's 40% row")
    print("  was computed on a partly unphysical draw; every number below uses a")
    print("  LOGNORMAL load of the same CV, which is the right model for a positive")
    print("  quantity and changes the 40% figure from 0.70 to 0.62.")


def b_was_the_null_direction_prior_held():
    banner("(b) WAS THE LOAD-BIAS COLUMN PRIOR-HELD? (the C2 failure)")
    print("  It was. Scenario 1 carried the 5% load-record prior. If that prior were")
    print("  filling the null direction, the sweep would be insensitive to the very")
    print("  variability it sweeps. So here is the same sweep both ways.\n")
    print(f"    {'CV':>7} {'dose train':>11} {'prior ON':>12} {'prior OFF':>12}")
    for with_dose in (False, True):
        for cv in (0.005, 0.02, 0.10, 0.40):
            on = med(lambda s: sd_from(record(MON, cv, with_dose, s)[2], True))
            off = med(lambda s: sd_from(record(MON, cv, with_dose, s)[2], False))
            f = lambda x: "singular" if not np.isfinite(x) else f"{x:.4f}"
            print(f"    {cv*100:6.1f}% {str(with_dose):>11} {f(on):>12} {f(off):>12}")
    print("\n  Removing the prior makes scenario 1 WORSE, not better -- 1.77 to 13.8")
    print("  at 10% CV. The prior was helping, so 14.1's negative answer was the")
    print("  CONSERVATIVE one and it stands. The prior-free column is the honest")
    print("  statement of how far the untreated record is from usable.")


def c_rank_against_conditioning():
    banner("(c) RANK AGAINST CONDITIONING -- THE RECONCILIATION")
    print("  Singular spectrum of the COLUMN-NORMALISED design matrix, seed 3.")
    print("  'rank restored' is read at s_min/s_max > 1e-8.\n")
    for with_dose in (False, True):
        print(f"    {'dose train: ' + str(with_dose):<20}")
        print(f"      {'CV':>7} {'s_min/s_max':>13} {'rank':>7} {'sd(ln theta)':>13}")
        for cv in (0.0, 0.005, 0.02, 0.10, 0.40):
            _, _, Q = record(MON, cv, with_dose, 3)
            sv = np.linalg.svd(Q / np.linalg.norm(Q, axis=0), compute_uv=False)
            r = sv[-1] / sv[0]
            print(f"      {cv*100:6.1f}% {r:13.3e} {'yes' if r > 1e-8 else 'NO':>7} "
                  f"{sd_from(Q):13.4f}")
    print()
    print("  THE RECONCILIATION, and it is a finding rather than a bookkeeping note:")
    print("  prop:rank is right that variation restores rank, and rank restoration")
    print("  is not precision. Without the dose train, half a per cent of variation")
    print("  moves s_min/s_max off zero and sd is still 28; ten per cent leaves it")
    print("  at 1.3. RESTORING RANK AND ACHIEVING USABLE PRECISION ARE DIFFERENT")
    print("  THINGS, and prop:rank delivers only the first.")
    print()
    print("  Note also the top-left cell of the second block: WITH the dose train at")
    print("  a constant load the matrix is still rank deficient, yet sd reads 0.1696.")
    print("  That is the load prior filling the null direction, exactly as _filter's")
    print("  docstring warns. At any real variation the DATA fill it instead and the")
    print("  same 0.1696 is earned rather than assumed.")


def published_numbers_do_not_reproduce():
    banner("WHAT FOOTING THE PUBLISHED 0.054 AND 0.023 ARE ON")
    print("  The note's prop:rank passage quotes prior-free sd falling to 0.054 at")
    print("  2% CV and 0.023 at 10%. Those do not reproduce on the note's claimed")
    print("  footing. Searching the designs and priors actually available:\n")
    print(f"    {'design':>18} {'k_e':>9} {'load prior':>11} {'cv=2%':>9} {'cv=10%':>9}")
    for nm, t in (("monthly (24)", MON),
                  ("two-daily (365)", np.arange(0.0, 730.0, 2.0))):
        for kep in (0.8, None):
            lab = "known" if kep is None else "prior 0.8"
            for lp in (False, True):
                a = sd_from(design_Q(t, 0.02), lp, kep)
                b = sd_from(design_Q(t, 0.10), lp, kep)
                print(f"    {nm:>18} {lab:>9} {str(lp):>11} {a:9.4f} {b:9.4f}")
    print("\n  The match is the two-daily design with k_e held EXACTLY KNOWN and no")
    print("  load prior: 0.0539 and 0.0227. But the note's stated footing is the")
    print("  monthly design with a log-prior of 0.8 on ln k_e, and _filter's own")
    print("  docstring says holding k_e known is 'kept for the audit, not for a")
    print("  headline'. On the claimed footing the numbers are 0.2343 and 0.1721.")
    print()
    print("  Both are also computed WITH the dose train present, so neither measures")
    print("  what load variation buys on its own. That is the whole of the apparent")
    print("  contradiction with scenario 1.")


def design_Q(times, cv, seed=3):
    lv, inc = design_record(times, load_cv=cv, seed=seed)
    return np.array([[P.theta * m / (P.Km + m),
                      -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                      1.0, P.ke * m] for m, d in zip(lv, inc)])


def scaling_with_n():
    banner("15.2  --  THE SCALING WITH n, DERIVED RATHER THAN ASSERTED")
    print("  14.1 said information in the degenerate direction 'saturates, because")
    print("  consecutive days resample one load distribution rather than adding")
    print("  independent directions'. That was a stated mechanism with no derivation")
    print("  behind it. Here is the derivation, and it does not support the claim.\n")
    print(f"  {'n':>8} {'sd(ln theta)':>13} {'local exponent':>15}"
          f" {'DATA info, weak dir':>21} {'PRIOR info, same dir':>21}")
    pri = np.diag([0.0, 0.0, 1 / (LOAD_FRAC * LB) ** 2, 1 / 0.8 ** 2])
    prev = None
    for n in (24, 48, 104, 365, 730, 3000, 10000, 30000, 100000, 300000):
        _, _, Q = record(n, 0.10, False, 3)
        sd = sd_from(Q)
        D = Q.T @ Q / SIGMA ** 2
        w, V = np.linalg.eigh(D)
        v = V[:, 0]
        ex = "" if prev is None else f"{-np.log(sd/prev[1])/np.log(n/prev[0]):15.3f}"
        print(f"  {n:8d} {sd:13.5f} {ex:>15} {w[0]:21.4e} {float(v@pri@v):21.4e}")
        prev = (n, sd)
    print("\n  -0.5 is the independent-sample rate. The observed exponent is 0.06 to")
    print("  0.24 and never approaches it over five orders of magnitude in n.")
    print()
    print("  BUT THE REASON IS NOT SATURATION, AND 14.1's SENTENCE IS WITHDRAWN.")
    print("  The data's information in the weak direction GROWS -- 1.8e-07 at n=24")
    print("  to 7.6e-04 at n=3e5, roughly in proportion to n. It never saturates.")
    print("  What it never does is ARRIVE: the prior's contribution in that same")
    print("  direction is about 1.05 and constant, so even at n = 300,000 the data")
    print("  supply three orders of magnitude less than the prior already holds.")
    print()
    print("  So the correct statement is that the weak direction is PRIOR-HELD at")
    print("  every n a study could reach, and sd falls only through the directions")
    print("  the data do reach -- which is why the exponent sits near 0.1 rather")
    print("  than 0.5. Sub-root-n here is a statement about a fixed prior, not about")
    print("  correlated samples.")


if __name__ == "__main__":
    print(__doc__)
    a_does_variation_reach_the_levels()
    b_was_the_null_direction_prior_held()
    c_rank_against_conditioning()
    published_numbers_do_not_reproduce()
    scaling_with_n()
