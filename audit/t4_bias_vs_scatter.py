"""
C2. THE TWO THINGS CALLED "THE LOAD RECORD": BIAS AGAINST PER-VISIT SCATTER.

Scoping only. Nothing here reaches the paper.

The paper asks for intake to be recorded "to better than about 10%", and that
single figure is doing two jobs. Proposition 3 says the fibre runs along a
SYSTEMATIC offset b: a record that is wrong by a constant fraction moves the
estimate along the degeneracy, and a prior on b is what closes it. Zero-mean
PER-VISIT error is a different thing entirely and enters only through sigma,
the ordinate noise of the fit.

Section 5.3 collapses the two. Its sweep moves a single "recording error"
through both channels at once, which is why the accuracy matters so much more
than the variability, and why the requirement is as tight as it is. This
separates them.

    bias       the prior width on the constant offset, s_b = bias * L
    scatter    the per-visit relative error, entering sigma as scatter * L

The paper's own configuration is bias = scatter = 5%, and its 20% row is
bias = scatter = 20%.

WHY IT MATTERS OPERATIONALLY. These two are bought by different studies. The
per-visit scatter needs every patient to keep a better diary, for ever. The
bias needs one weighed-intake calibration substudy, in a subset, once: the
PKU figure the paper cites, actual natural protein at 167% of prescribed with
a median 40% above it, is exactly a bias of that kind and exactly what such a
substudy would measure. If the bias term dominates, the cheap study buys the
gate and the expensive one is not needed.

Run:  python audit/t4_bias_vs_scatter.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

import model_A.run as LOC
from common.models import ParamsA
from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, SEEDS, load, median_spread)

KM, KE = 100.0, 5.0
LB = load(0.0)
KE_PRIOR_LOG = 0.8
GATE = 0.20
MONTH = 30.4
N_VISITS = 25
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))


def op_level(theta):
    B = KE * KM + theta - LB
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)


def sd_of(theta, seed, bias, scatter, n_visits=N_VISITS):
    """
    sd(ln theta) with the two halves of the load record set independently.

    `bias` is the prior width on the constant offset; `scatter` is the
    per-visit relative error, which reaches the fit only through sigma.
    bias = scatter reproduces the paper's own sweep.
    """
    t = np.arange(0.0, n_visits * MONTH, MONTH)[:n_visits]
    inc = PEAK * np.exp(-LAMBDA * np.mod(t - T_DOSE, DOSE_INTERVAL))
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    sig = LOC.sigma_total(scatter, ASSAY_CV,
                          p=ParamsA(theta=theta, Km=KM, ke=KE),
                          m_bar=op_level(theta))
    pri_b = np.inf if bias is None else (bias * LB) ** 2
    F = np.diag([0.0, 0.0, 0.0 if np.isinf(pri_b) else 1.0 / pri_b,
                 1.0 / KE_PRIOR_LOG ** 2])
    for m, d in zip(lv, inc):
        gm = m / (KM + m)
        q = np.array([theta * gm, -(theta + d) * KM * m / (KM + m) ** 2,
                      1.0, KE * m])
        F = F + np.outer(q, q) / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over(theta, bias, scatter):
    return median_spread([sd_of(theta, s, bias, scatter) for s in SEEDS])


def bias_numbers():
    """
    The load decomposition Sections 6 and 8 quote, generated rather than typed.

    WHY THIS FUNCTION NOW EXISTS. This script's docstring said "Scoping only.
    Nothing here reaches the paper." It owns the paper's third-tightest
    requirement, stated as prose literals in TWO sections and compared with
    nothing (REVIEW_LOG.md R81).
    """
    from common.numbers import write_numbers
    TH = 3000.0
    base = over(TH, LOAD_FRAC, LOAD_FRAC)[0]
    both = over(TH, 0.20, 0.20)[0]
    bias = over(TH, 0.20, LOAD_FRAC)[0]
    scat = over(TH, LOAD_FRAC, 0.20)[0]
    tot = (bias - base) + (scat - base)
    N = {"bias_base": base, "bias_both20": both,
         "bias_alone20": bias, "bias_scatter20": scat,
         "bias_scatter_share": 100.0 * (scat - base) / tot,
         "bias_offset_share": 100.0 * (bias - base) / tot,
         "bias_tight": over(TH, 0.01, 0.01)[0],
         "bias_loose": both,
         "bias_sweep_lo": over(TH, 0.01, LOAD_FRAC)[0],
         "bias_sweep_hi": over(TH, 5.00, LOAD_FRAC)[0]}
    #  the format is the precision the PAPER prints, so a macro renders exactly
    #  the literal it replaces and no number changes in the conversion
    fmt = {k: "{:.3f}" for k in N}
    fmt.update({"bias_scatter_share": "{:.0f}", "bias_offset_share": "{:.0f}",
                "bias_sweep_lo": "{:.4f}", "bias_sweep_hi": "{:.4f}"})
    write_numbers("bias", N, fmt)
    return N


def main():
    print(__doc__)
    bar = "=" * 78
    TH = 3000.0

    print(bar)
    print("1. THE PAPER'S SWEEP, AND WHAT IT IS ACTUALLY MOVING")
    print(bar)
    print(f"  {N_VISITS} monthly draws at nominal capacity, medians over "
          f"{len(SEEDS)} records\n")
    print(f"  {'recording error':>16} {'both (the paper)':>24} "
          f"{'bias only':>24} {'scatter only':>24}")
    for e in (0.01, 0.05, 0.10, 0.20, 0.40):
        both = over(TH, e, e)[0]
        b_only = over(TH, e, LOAD_FRAC)[0]      # scatter held at the paper's 5%
        s_only = over(TH, LOAD_FRAC, e)[0]      # bias held at the paper's 5%
        print(f"  {100*e:15.0f}% {both:11.4f}{'  *' if both < GATE else '   '}"
              f"            {b_only:11.4f}             {s_only:11.4f}")

    print("\n" + bar)
    print("2. THE DECOMPOSITION AT THE 20% ROW")
    print(bar)
    base = over(TH, LOAD_FRAC, LOAD_FRAC)[0]
    both = over(TH, 0.20, 0.20)[0]
    b_only = over(TH, 0.20, LOAD_FRAC)[0]
    s_only = over(TH, LOAD_FRAC, 0.20)[0]
    print(f"  the paper's 5%/5% baseline          {base:.4f}")
    print(f"  at 20% in both channels             {both:.4f}"
          f"   ({both/base:.2f}x)")
    print(f"  at 20% in the BIAS alone            {b_only:.4f}"
          f"   ({b_only/base:.2f}x)")
    print(f"  at 20% in the SCATTER alone         {s_only:.4f}"
          f"   ({s_only/base:.2f}x)")
    lb, ls = b_only - base, s_only - base
    print(f"\n  of the {both-base:.4f} the paper's 20% row adds to the "
          f"baseline,")
    print(f"    the bias channel accounts for    {lb:.4f}  "
          f"({100*lb/(lb+ls):.0f}% of the two taken separately)")
    print(f"    the scatter channel accounts for {ls:.4f}  "
          f"({100*ls/(lb+ls):.0f}%)")

    print("\n" + bar)
    print("3. WHAT A BIAS-CALIBRATION SUBSTUDY WOULD BUY")
    print(bar)
    print("  The per-visit scatter is left where it is and only the mean")
    print("  offset is tightened, which is what one weighed-intake substudy")
    print("  in a subset measures.\n")
    for scat, lab in ((LOAD_FRAC, "scatter at the paper's 5%"),
                      (0.20, "scatter at 20%"),
                      (0.40, "scatter at 40%, the PKU reading")):
        print(f"  {lab}")
        for bias in (None, 0.40, 0.20, 0.10, 0.05, 0.02, 0.01):
            m, lo, hi = over(TH, bias, scat)
            nm = "no prior at all" if bias is None else f"bias to {100*bias:.0f}%"
            mark = " *" if m < GATE else "  "
            print(f"     {nm:>18}: {m:9.4f} [{lo:.4f}, {hi:.4f}]{mark}")
        print()

    print(bar)
    print("4. THE TRADE, STATED AS A REQUIREMENT ON EACH CHANNEL ALONE")
    print(bar)
    lo_b, hi_b = over(TH, 0.01, LOAD_FRAC)[0], over(TH, 5.0, LOAD_FRAC)[0]
    if max(lo_b, hi_b) < GATE:
        print(f"  with the scatter at the paper's 5%, the gate holds at EVERY "
              f"bias from 1% to 500%")
        print(f"    ({lo_b:.4f} at 1%, {hi_b:.4f} at 500%): the bias channel "
              f"never binds at all")
    else:
        b_star = brentq(lambda b: over(TH, b, LOAD_FRAC)[0] - GATE, 0.01, 5.0,
                        xtol=1e-4)
        print(f"  with the scatter at the paper's 5%, the bias may be as loose "
              f"as {100*b_star:.0f}%")
    try:
        s_star = brentq(lambda s: over(TH, LOAD_FRAC, s)[0] - GATE, 0.01, 2.0,
                        xtol=1e-4)
        print(f"  with the bias at the paper's 5%, the scatter must be tighter "
              f"than {100*s_star:.1f}%")
    except ValueError:
        print("  with the bias at 5%, no scatter in [1%, 200%] reaches the gate")
    try:
        both_star = brentq(lambda e: over(TH, e, e)[0] - GATE, 0.01, 2.0,
                           xtol=1e-4)
        print(f"  moved together, as the paper moves them, the requirement is "
              f"{100*both_star:.1f}%")
    except ValueError:
        print("  moved together, no value in [1%, 200%] reaches the gate")

    print("\n" + bar)
    print("5. DOES THE SPLIT SURVIVE INTO THE DISEASE RANGE?")
    print(bar)
    for th, lab in ((3000.0, "100%"), (300.0, "10%"), (150.0, "5%")):
        row = [f"  theta = {th:6.0f} ({lab:>4}):"]
        for bias, scat, nm in ((LOAD_FRAC, LOAD_FRAC, "5/5"),
                               (0.01, LOAD_FRAC, "bias 1%"),
                               (LOAD_FRAC, 0.01, "scatter 1%"),
                               (0.01, 0.01, "both 1%")):
            row.append(f"{nm} {over(th, bias, scat)[0]:8.4f}")
        print("   ".join(row))

    print("\n" + bar)
    print("R1. RE-CHECK: bias = scatter reproduces the paper's own sweep")
    print(bar)
    for e, want in ((0.01, 0.1590), (0.05, 0.1699), (0.10, 0.1991),
                    (0.20, 0.2876), (0.40, 0.5049)):
        got = over(TH, e, e)[0]
        print(f"  {100*e:4.0f}%: {got:.4f} against tab_D_loadprec's {want:.4f}"
              f"   ({abs(got-want)/want:.2%} apart)")

    print("\nR2. RE-CHECK: dropping the bias prior is not the same as "
          "loosening it")
    loose = over(TH, 5.0, LOAD_FRAC)[0]
    none = over(TH, None, LOAD_FRAC)[0]
    print(f"  a 500% bias prior gives {loose:.4f}; none at all gives "
          f"{none:.4f}")
    print(f"  they agree to {abs(loose-none)/none:.2%}, so the prior is inert "
          f"once it is wide")


    #  the numbers the paper cites, written where the build can see them
    bias_numbers()

if __name__ == "__main__":
    main()
