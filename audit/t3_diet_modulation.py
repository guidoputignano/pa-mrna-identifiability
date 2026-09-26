"""
C1. A DESIGNED, RECORDED MODULATION OF DIETARY PROTEIN.

Scoping only. Nothing here reaches the paper.

Section 5.4 prices spontaneous load variability and finds it insufficient:
1.79 at a 10% day-to-day coefficient of variation, flooring near 0.376 however
variable the intake is made. But that is variation nobody planned, nobody
phased and nobody recorded exactly. The question here is different: what is a
DELIBERATE modulation worth, prescribed on a known cycle, recorded as
prescribed, and carried as a known regressor alongside the dose train.

WHAT CHANGES AND WHAT DOES NOT. A recorded modulation adds no parameter: the
estimator knows L_i, so the load-bias column stays a single unknown constant.
What it changes is which levels the record visits. A bigger, phase-controlled
swing in L moves m further across phi(m) = K_m/(K_m+m), which is the direction
that separates theta from K_m. That is the same mechanism the dose train uses,
applied to the other side of the balance.

TWO READINGS OF THE RECORDING ERROR, and the honest one is reported first.

  conservative   the recording error applies to total intake however it
                 arose, so sigma is unchanged and the only gain is the
                 larger excursion. This is what is reported unless stated.
  optimistic     the prescribed component is known exactly and only the
                 spontaneous part carries recording error, so sigma falls
                 as well. Reported alongside, because the truth is between.

THE ALIASING QUESTION IS ASKED TOO. Section 5.3 found that an unrecorded
antibiotic cycle becomes harmless when its period equals the sampling
interval, because every draw lands at the same phase. A designed modulation
has the same geometry with the sign reversed: at the wrong period it delivers
nothing, and that is a protocol trap worth naming.

CLINICAL SAFETY IS ASSUMED, NOT ESTABLISHED. Every amplitude below is a
fraction of prescribed natural protein. Whether a swing of that size is safe
in a patient whose decompensations are driven by protein load is a question
for a metabolic physician and nothing here bears on it.

Run:  python audit/t3_diet_modulation.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

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
#: the modulation period Part 1 finds best, hoisted to module scope so the
#: generated-numbers function uses the same value main() prints.
P_GOOD = 21.0
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))


def op_level(theta, Lb=LB):
    B = KE * KM + theta - Lb
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * Lb * KM)) / (2.0 * KE)


def sigma_of(theta, load_frac=LOAD_FRAC):
    return LOC.sigma_total(load_frac, ASSAY_CV,
                           p=ParamsA(theta=theta, Km=KM, ke=KE),
                           m_bar=op_level(theta))


def sd_of(theta, seed, amp=0.0, period_d=0.0, dose_train=True,
          n_visits=N_VISITS, optimistic=False):
    """
    sd(ln theta) with a prescribed square-wave modulation of the load.

    `amp` is the fractional swing about the prescribed mean, so amp = 0.15 is
    a diet alternating between 115% and 85% of the usual natural protein.
    `period_d` is the full cycle in days. amp = 0 is the paper's own design.
    """
    t = np.arange(0.0, n_visits * MONTH, MONTH)[:n_visits]
    inc = (PEAK * np.exp(-LAMBDA * np.mod(t - T_DOSE, DOSE_INTERVAL))
           if dose_train else np.zeros_like(t))
    mod = (1.0 + amp * np.sign(np.sin(2.0 * np.pi * t / period_d))
           if amp and period_d else np.ones_like(t))
    Ls = LB * mod * (1.0 + LOAD_CV *
                     np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)

    #  the recording error: over total intake, or over the spontaneous part
    #  alone when the prescribed swing is taken to be known exactly
    frac = LOAD_FRAC * (LOAD_CV / np.hypot(LOAD_CV, amp)) if (optimistic and amp) \
        else LOAD_FRAC
    sig = sigma_of(theta, frac)

    F = np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * LB) ** 2, 1.0 / KE_PRIOR_LOG ** 2])
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


def over(theta, **kw):
    return median_spread([sd_of(theta, s, **kw) for s in SEEDS])


def diet_numbers():
    """
    The prescribed-modulation figures Sections 8 and 11 quote.

    WHY THIS FUNCTION NOW EXISTS. See REVIEW_LOG.md R80 and R81: this script
    is absent from `make figures` and its results reach the paper only as
    prose literals.
    """
    from common.numbers import write_numbers
    TH_DIS = 150.0
    base = over(TH_DIS)[0]
    #  ONLY the value that reproduces is generated. This script's own output
    #  gives a 6.2% gain at a 15% swing and 24.1% at 50%, against the paper's
    #  "7%" and "32%" at 08_specification.tex:372 and 11_scope.tex:204. The
    #  baseline is not in doubt: 2.3198 at a 50% swing rounds to the paper's
    #  own 2.32, so the capacity and period are the ones the paper used. The
    #  two percentages therefore DISAGREE with their generator, and neither is
    #  adopted here: wiring the computed value would silently change a
    #  published number and wiring the published one would falsify the
    #  generator. Reported instead (REVIEW_LOG.md R88).
    N = {"diet_sd50_disease": over(TH_DIS, amp=0.50, period_d=P_GOOD)[0],
         "diet_gain15_pct_UNRECONCILED": 100.0 * (
             1.0 - over(TH_DIS, amp=0.15, period_d=P_GOOD)[0] / base),
         "diet_gain50_pct_UNRECONCILED": 100.0 * (
             1.0 - over(TH_DIS, amp=0.50, period_d=P_GOOD)[0] / base)}
    fmt = {"diet_sd50_disease": "{:.2f}",
           "diet_gain15_pct_UNRECONCILED": "{:.1f}",
           "diet_gain50_pct_UNRECONCILED": "{:.1f}"}
    write_numbers("diet", N, fmt)
    return N


def main():
    print(__doc__)
    bar = "=" * 78
    caps = [(3000.0, "100%"), (300.0, "10%"), (150.0, "5%")]
    #  a period that is not a multiple of the 30.4 d sampling interval and not
    #  a multiple of the 14 d dose interval either

    print(bar)
    print("1. WITH AND WITHOUT THE DOSE TRAIN, AT A 15% SWING ON 21 DAYS")
    print(bar)
    print(f"  {N_VISITS} monthly draws, medians over {len(SEEDS)} intake "
          f"records; the modulation is prescribed and recorded\n")
    print(f"  {'theta':>8} {'dose train':>11} {'no modulation':>24} "
          f"{'+/-15% on 21 d':>24} {'gain':>7}")
    for th, lab in caps:
        for train in (True, False):
            a = over(th, dose_train=train)
            b = over(th, amp=0.15, period_d=P_GOOD, dose_train=train)
            gain = a[0] / b[0] if np.isfinite(b[0]) and b[0] > 0 else np.nan
            print(f"  {th:8.0f} {str(train):>11} "
                  f"{a[0]:9.4f} [{a[1]:7.4f},{a[2]:7.4f}] "
                  f"{b[0]:9.4f} [{b[1]:7.4f},{b[2]:7.4f}] {gain:6.2f}x")

    print("\n" + bar)
    print("2. THE AMPLITUDE SWEEP, WITH THE DOSE TRAIN, PERIOD 21 DAYS")
    print(bar)
    print(f"  {'amplitude':>10}" + "".join(f"{l:>26}" for _, l in caps))
    for amp in (0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50):
        row = [f"  {100*amp:9.0f}%"]
        for th, _ in caps:
            m, lo, hi = over(th, amp=amp, period_d=P_GOOD)
            mark = "*" if m < GATE else " "
            row.append(f"{m:12.4f} [{lo:6.4f},{hi:6.4f}]{mark}")
        print("".join(row))

    print("\n" + bar)
    print("3. THE PERIOD SWEEP AT A 15% SWING: THE ALIASING TRAP")
    print(bar)
    print("  The sampling interval is 30.4 d and the dose interval is 14 d.")
    print("  A modulation at either, or at a multiple, lands every draw at the")
    print("  same phase of it and delivers nothing.\n")
    print(f"  {'period':>10}" + "".join(f"{l:>26}" for _, l in caps))
    for per in (7.0, 14.0, 15.2, 21.0, 28.0, 30.4, 45.0, 60.8, 91.0):
        row = [f"  {per:8.1f} d"]
        for th, _ in caps:
            m, lo, hi = over(th, amp=0.15, period_d=per)
            row.append(f"{m:12.4f} [{lo:6.4f},{hi:6.4f}] ")
        tag = ""
        if abs(per - 30.4) < 1e-9 or abs(per - 60.8) < 1e-9:
            tag = "  <- a multiple of the sampling interval"
        if abs(per - 14.0) < 1e-9 or abs(per - 28.0) < 1e-9:
            tag = "  <- a multiple of the dose interval"
        print("".join(row) + tag)

    print("\n" + bar)
    print("4. A RECORDED MODULATION AGAINST THE SAME SWING ARRIVING BY ITSELF")
    print(bar)
    print("  Spontaneous variation of the same size is already in the record")
    print("  as a day-to-day CV. The comparison holds the total swing fixed")
    print("  and moves it between the two channels.\n")
    for th, lab in caps:
        designed = over(th, amp=0.15, period_d=P_GOOD)[0]
        #  the same 15% arriving spontaneously: raise the day-to-day CV so the
        #  record carries the same total variation, with nothing prescribed
        spont = median_spread([sd_of_cv(th, s, cv=np.hypot(LOAD_CV, 0.15))
                               for s in SEEDS])[0]
        print(f"  theta = {th:6.0f} ({lab:>4}): designed {designed:8.4f}   "
              f"spontaneous {spont:8.4f}   "
              f"ratio {spont/designed:5.3f}")

    print("\n" + bar)
    print("5. THE OPTIMISTIC READING: THE PRESCRIBED SWING KNOWN EXACTLY")
    print(bar)
    for th, lab in caps:
        a = over(th, amp=0.15, period_d=P_GOOD)
        b = over(th, amp=0.15, period_d=P_GOOD, optimistic=True)
        print(f"  theta = {th:6.0f} ({lab:>4}): conservative {a[0]:8.4f}   "
              f"optimistic {b[0]:8.4f}")

    print("\n" + bar)
    print("R1. RE-CHECK: zero amplitude reproduces the paper")
    print(bar)
    m, lo, hi = over(3000.0)
    print(f"  {m:.4f} [{lo:.4f}, {hi:.4f}] against Table 3's "
          f"0.1699 [0.1684, 0.1709]")
    print("\nR2. WHICH WAY THE AMPLITUDE RUNS AT NOMINAL CAPACITY")
    seq = [over(3000.0, amp=a, period_d=P_GOOD)[0]
           for a in (0.0, 0.05, 0.10, 0.15, 0.20, 0.30)]
    print(f"  {[round(x, 4) for x in seq]}")
    print(f"  monotonically better: "
          f"{all(b <= a for a, b in zip(seq, seq[1:]))}")
    print("  It is not. At nominal capacity a larger prescribed swing is")
    print("  slightly WORSE, because the dose train already supplies the")
    print("  excitation and the extra load movement buys nothing against the")
    print("  load-bias column it also feeds. The sign is the result here.")
    print("\nR3. RE-CHECK: at the sampling period the modulation is inert")
    at = over(3000.0, amp=0.30, period_d=30.4)[0]
    none = over(3000.0)[0]
    print(f"  a 30% swing at 30.4 d gives {at:.4f} against {none:.4f} with no "
          f"modulation at all: {abs(at-none)/none:.2%} apart")


    #  the numbers the paper cites, written where the build can see them
    diet_numbers()

def sd_of_cv(theta, seed, cv, n_visits=N_VISITS):
    """The paper's design with the day-to-day CV raised, nothing prescribed."""
    t = np.arange(0.0, n_visits * MONTH, MONTH)[:n_visits]
    inc = PEAK * np.exp(-LAMBDA * np.mod(t - T_DOSE, DOSE_INTERVAL))
    Ls = LB * (1.0 + cv * np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    sig = sigma_of(theta)
    F = np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * LB) ** 2, 1.0 / KE_PRIOR_LOG ** 2])
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


if __name__ == "__main__":
    main()
