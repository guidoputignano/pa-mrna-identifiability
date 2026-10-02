"""
Item 1 (the ratio route): does kappa need estimating at all? If the estimand is
capacity as a FRACTION of normal (which is what the figure's x-axis and the gate
both are), then theta_frac = A_patient / mean(A_healthy) cancels the population
mapping kappa by construction -- no known-theta anchor anywhere, only healthy
controls, who are obtainable.

Item 3 (the assumption): the ratio rests on kappa being the same in patients and
controls. How much would kappa have to differ to break 0.20?

Footing: the repo's theta = 3000 is V_max (absolute, n4_theta_scope.py:10), but
every result is reported as '% of normal' and sd(ln theta) is scale-invariant, so
all prior numbers already apply to the fraction. The ratio route supplies the
direct fraction reading WITHOUT the kappa calibration the absolute route needed.
"""
import numpy as np

import core as C


def sd_ratio(Y, sig_healthy, N, nread=1, seeds=C.SEEDS):
    """sd(ln theta_frac) at 5%: record + an ex vivo RATIO reading of theta/theta_normal.
    The reading's precision is hypot(patient repeatability Y/sqrt(nread), healthy-mean
    SE sig_healthy/sqrt(N)); kappa cancels in the ratio, so there is no kappa term."""
    P = np.hypot(Y / np.sqrt(nread), sig_healthy / np.sqrt(N))

    def one(s):
        r = C.record(C.TH_5, seed=s)
        row = np.array([1.0, 0, 0, 0])          # direct reading of ln theta (fraction)
        return C.sd_th(C.fim(r, noise="perdraw") + np.outer(row, row) / P ** 2)
    return C.over_seeds(one, seeds)[0]


def main():
    print("=" * 92)
    print("  ITEM 1. THE RATIO ROUTE: theta_frac = A_patient / mean(A_healthy), kappa CANCELS")
    print("=" * 92)
    print("""
  The estimand is capacity as a FRACTION of normal (the figure's x-axis; the gate
  sd(ln theta)=0.20 is scale-free). A_patient = kappa*theta_patient,
  mean(A_healthy) = kappa*theta_normal: the ratio is theta_patient/theta_normal and
  kappa drops out. No known-theta anchor is needed -- only N healthy controls,
  whose mean carries an error sig_healthy/sqrt(N) that falls with N. sig_healthy is
  the between-healthy-subject spread of ln A (biological PCC-activity spread + any
  per-subject kappa + repeatability); the repository has NO figure for it (the 51%
  it carries is the 3-HP marker baseline between mice, not PCC activity), so it is
  swept below and flagged.
""")
    print("  sd(ln theta_hat) at 5% of normal, ratio route (no kappa anchor):")
    for Y, yl in ((0.20, "patient reading 20% (one read)"),
                  (0.10, "patient reading 10% (or 20% x4 reads)")):
        print(f"\n  {yl}:")
        hdr = "sig_healthy \\ N controls"
        print(f"    {hdr:>26}  {'N=10':>7}  {'N=30':>7}  {'N=100':>7}")
        for sig in (0.30, 0.50, 0.70):
            cells = "".join(f"{sd_ratio(Y, sig, N):8.3f}" + ("*" if sd_ratio(Y, sig, N) <= 0.20 else " ")
                            for N in (10, 30, 100))
            print(f"    healthy spread {int(sig*100):>3}%        {cells}")
    print("""
  READING: with the patient read to ~10% (achievable by averaging a few assays) and
  ~30-100 obtainable healthy controls, the ratio route reaches ~0.11-0.19 at 5%
  capacity WITH NO known-theta anchor. kappa never has to be estimated: it cancels.
  The 'calibrate kappa from known-theta patients' framing (t27) answered a question
  that need not be asked. The one real cost is a healthy control panel, which exists.
  (At 20% patient repeatability it is borderline, ~0.21-0.26; the lever is averaging
  the patient's own reads, not more controls, once N >~ 30.)
""")

    print("=" * 92)
    print("  ITEM 3. THE ASSUMPTION: kappa IDENTICAL IN PATIENTS AND CONTROLS")
    print("=" * 92)
    print("""
  The ratio is unbiased only if the cell-to-whole-body mapping kappa is the same in
  patients and controls. Under chronic overload (different expression, different
  tissue contribution) it may not be. A systematic difference
  b = ln(kappa_patient / kappa_control) passes straight into ln theta_hat as a BIAS,
  so the total error is rmse = hypot(sd, b).
""")
    # at an achievable sd (Y=10%, sig=0.5, N=30) what b breaks 0.20?
    for Y, sig, N in ((0.10, 0.50, 30), (0.10, 0.50, 100), (0.20, 0.50, 30)):
        sd = sd_ratio(Y, sig, N)
        b_break = np.sqrt(max(0.20 ** 2 - sd ** 2, 0.0))
        print(f"    at sd = {sd:.3f} (read {int(Y*100)}%, {N} controls, spread {int(sig*100)}%): "
              f"a kappa difference of ln-ratio {b_break:.3f} "
              f"(kappa_patient {100*(np.exp(b_break)-1):+.0f}%) makes rmse = 0.20")
    print("""
  So the ratio tolerates a systematic kappa difference of only ~15-20% between
  chronically-overloaded patients and healthy controls before the bias alone reaches
  the target. Whether chronic propionate overload changes the PCC cell-to-whole-body
  mapping by less than ~15% is a BIOLOGICAL question the statistics cannot settle,
  and it is the assumption on which the ratio route stands or falls. The repository
  says nothing about it.
""")


if __name__ == "__main__":
    main()
