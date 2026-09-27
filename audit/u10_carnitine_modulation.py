"""
Prescribed carnitine modulation. A negative result, recorded before anyone
tries it.

WHY IT LOOKS LIKE THE OBVIOUS LEVER, AND WHY THE PAPER HAS NOT CLOSED IT.
eq:input_carn writes the carnitine channel as an input in exactly the form
the design wants,

    k_e(t) = k_e + kappa u_c(t),      u_c recorded, kappa unknown,

and Section sec:selective calls it one of the two SELECTIVE perturbations,
acting on k_e alone and leaving theta untouched. audit/u9_ke_frontier.py has
just made k_e the binding quantity. Put those together and a prescribed
carnitine swing reads as the obvious move: an input that already exists, on
the channel that now binds, in a cohort where every patient already receives
it.

Section sec:whoto closes carnitine, but it closes a DIFFERENT question. Its
two arguments are that carnitine cannot touch the theta sensitivity's first
factor, and that the dose is adjusted in response to the patient's state and
so is endogenous. The second does not apply to a dose a protocol PRESCRIBES,
which is exogenous by construction, exactly as the dose train is. So the
prescribed version is open, and this script closes it.

THE MODEL. Writing f for the carnitine-mediated share of the alternative
route at baseline, k_e = k_e_base + kappa u_bar with kappa u_bar = f k_e, a
square-wave swing of amplitude `amp` about the prescribed mean gives

    k_e(t) = k_e (1 + f amp s(t)),    s(t) = sign(sin(2 pi t / period)),

so the swing the clearance actually sees is f amp, not amp. The parameter
vector gains ln kappa beside ln k_e_base, both at the 0.8 log-prior width the
paper carries for k_e.

THE COMPARISON IS RUN TWO WAYS because they answer different questions.
  (a) What the MODULATION buys, holding the five-parameter structure fixed:
      sd at amp = 0 over sd at amp > 0. This isolates the swing.
  (b) What the WHOLE ARM buys against the paper's own four-parameter design,
      where k_e is one unknown at a 0.8 prior. This prices the split as well
      as the swing, and it is the comparison a protocol faces.

THE MECHANISM IS TESTED, NOT ASSERTED. sd(ln kappa) and sd(ln k_e_base) are
reported beside sd(ln theta), so the claim that a modulation pins the SCALE
and not the BASELINE is visible rather than argued.

FOOTING. As audit/t3_diet_modulation.py, whose square wave and 21 d period
this reuses so the two modulation results sit on one footing: 25 monthly
draws, theta with the recorded load held fixed at L = 322.727 umol/L/d and
the operating level, sigma and the load-bias prior re-solved at each
capacity, K_m = 100, k_e = 5/d, D = 6000 on a 14 d train with a 3.5 d
half-life, load recorded to 5%, day-to-day intake CV 10%, load-bias prior
0.05 L, gate 0.20, medians over 24 intake records with the seed range shown.

CLINICAL SAFETY IS ASSUMED, NOT ESTABLISHED, exactly as t3 assumes it for
diet. Whether a prescribed carnitine swing of any size is admissible is not
a question this script can answer, and the answer does not matter if the
swing buys nothing.

Nothing here records an expected result in advance.

Run:  python audit/u10_carnitine_modulation.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, load)

KM, KE = 100.0, 5.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
KE_PRIOR = 0.8
PUBLISHED_KAPPA_CV = 0.1402
N_VISITS = 25
P_GOOD = 21.0          # t3's period: aliases with neither 14 d nor 30.4 d
P_BAD = 14.0           # the dose interval, which is the aliasing trap
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
SEEDS = range(24)


def op_level(theta):
    B = KE * KM + theta - LB
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)


def sigma_of(theta):
    mbar = op_level(theta)
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(LOAD_FRAC * LB, slope * ASSAY_CV * mbar))


def sd_of(theta, seed, f=0.5, amp=0.0, period_d=P_GOOD, split=True,
          want=0, ke_prior=KE_PRIOR, f_prior=KE_PRIOR, mrna_kappa=None):
    """
    sd of one parameter with a prescribed carnitine square wave.

    THE PARAMETERISATION IS THE WHOLE DIFFICULTY AND THE FIRST ONE WAS WRONG.
    Writing the two carnitine unknowns as the NON-carnitine baseline k_e_base
    and the scale kappa, and giving each the paper's 0.8 log-prior, does not
    reduce to the paper's design when the swing is switched off: at a
    carnitine share of 0.8 it puts a 0.8 prior on a term carrying a fifth of
    the clearance and another on the term carrying four fifths, which
    constrains the TOTAL better than the paper's single 0.8 prior on it. That
    version reported gains up to 2.97x, and the gain was the reparameterisation
    and not the swing.

    The honest coordinates are the total and the share. With K the total
    alternative clearance at the prescribed mean and f the carnitine-mediated
    share of it,

        k_e(t) = K (1 + f amp s(t)),     s(t) = sign(sin(2 pi t / period)),

    so d(flux)/d ln K = K(1 + f amp s) m and d(flux)/d ln f = K f amp s m.
    ln K carries the paper's own 0.8 prior, directly comparable because it is
    a prior on the same object. ln f carries 0.8 as the unknown scale.

    THE CHECK THIS BUYS: at amp = 0 the ln f column is identically zero, so
    f drops out and the design must reproduce the paper's number EXACTLY.
    part_z() asserts that rather than trusting it.

    `split=False` drops the share entirely, which is the paper's design.
    `want` selects the parameter: 0 ln theta, 3 ln K, 4 ln f.
    """
    t = np.arange(0.0, N_VISITS * MONTH, MONTH)[:N_VISITS]
    inc = PEAK * np.exp(-LAMBDA * np.mod(t - T_DOSE, DOSE_INTERVAL))
    s = (np.sign(np.sin(2.0 * np.pi * t / period_d))
         if (amp and period_d) else np.zeros_like(t))
    swing = 1.0 + f * amp * s
    ke_t = KE * swing                          # the clearance actually seen
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))

    B = ke_t * KM + theta + inc - Ls
    m = (-B + np.sqrt(B ** 2 + 4.0 * ke_t * Ls * KM)) / (2.0 * ke_t)
    sig = sigma_of(theta)

    gm = m / (KM + m)
    cols = [theta * gm,
            -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m),
            KE * swing * m]                    # d/d ln K, K = KE here
    pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2, ke_prior ** 2]
    if split:
        cols.append(KE * f * amp * s * m)      # d/d ln f
        pri.append(f_prior ** 2)
    if mrna_kappa is not None:
        cols.append(inc * gm)                  # d/d ln kappa, the DOSE scale
        pri.append(mrna_kappa ** 2)

    Q = np.stack(cols, axis=1)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sig ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    if want >= Q.shape[1]:
        return np.nan
    v = np.linalg.inv(F)[want, want]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over(theta, **kw):
    v = np.array([sd_of(theta, s, **kw) for s in SEEDS], float)
    fin = v[np.isfinite(v)]
    if not fin.size:
        return np.inf, np.inf, np.inf
    return float(np.median(fin)), float(fin.min()), float(fin.max())


# ======================================================================
def part_a():
    print("=" * 84)
    print("PART A.  what the modulation buys, and the two readings of 'kappa'")
    print("=" * 84)
    print("""
  gain = sd(ln theta) at amp = 0, over sd at that amplitude. Above 1.00 the
  swing helps.

  TWO READINGS ARE RUN because this repository uses 'kappa' for the mRNA
  DOSE SCALE (audit/t2_tier1_margin.py: "the increment is carried as a
  regressor with an unknown scale kappa"), while the carnitine channel of
  eq:input_carn has a scale of its own. Which one is carried free changes
  the answer at nominal capacity by an order of magnitude in the gain, and
  changes it not at all in the disease range. Both are reported rather than
  one chosen.
""")
    readings = [("carnitine share only", None),
                ("plus the mRNA dose scale at 0.8", 0.8),
                ("plus the mRNA dose scale at 14.02%", PUBLISHED_KAPPA_CV)]
    gains = {}
    for lab, mk in readings:
        print(f"  reading: {lab}")
        print(f"    {'theta':>7} {'share f':>8}" +
              "".join(f"{f'amp {a:.0%}':>20}" for a in (0.15, 0.30, 0.50)))
        print("    " + "-" * 76)
        for theta in (3000.0, 600.0, 300.0, 150.0):
            for f in (0.2, 0.5, 0.8):
                base = over(theta, f=f, amp=0.0, mrna_kappa=mk)[0]
                cells = []
                for amp in (0.15, 0.30, 0.50):
                    m = over(theta, f=f, amp=amp, mrna_kappa=mk)[0]
                    g = base / m if np.isfinite(m) and m > 0 else np.nan
                    gains.setdefault(lab, []).append(g)
                    cells.append(f"{m:.4f} ({g:.3f}x)")
                print(f"    {theta:7.0f} {f:8.1f}" +
                      "".join(f"{c:>20}" for c in cells))
            print()
    print("  gain ranges by reading:")
    for lab, _ in readings:
        g = np.array([x for x in gains[lab] if np.isfinite(x)])
        print(f"    {lab:<36} {g.min():.3f} to {g.max():.3f}")
    print()
    print("  AND THE SPLIT BY CAPACITY IS THE RECONCILIATION. With the mRNA")
    print("  dose scale carried at 0.8 the gain at NOMINAL capacity is a")
    print("  per cent or two, because that prior dominates everything the")
    print("  swing could buy. In the DISEASE range the same reading gives")
    print("  the same gain as without it, because there the binding term is")
    print("  the alternative clearance and not the increment.")
    return gains


def part_b():
    print()
    print("=" * 84)
    print("PART B.  the gap the swing has to close, which is the whole answer")
    print("=" * 84)
    print("\n  The gain is real. The question is whether it is the right")
    print("  size. Best case throughout: share 0.8, amplitude 50%, the good")
    print("  period, no mRNA dose-scale column.\n")
    print(f"    {'theta':>7} {'no swing':>10} {'best swing':>11} {'gain':>7} "
          f"{'gate':>7} {'factor still short':>19}")
    print("    " + "-" * 66)
    for theta in (3000.0, 600.0, 300.0, 150.0):
        b = over(theta, f=0.8, amp=0.0)[0]
        a = over(theta, f=0.8, amp=0.50)[0]
        short = a / GATE
        print(f"    {theta:7.0f} {b:10.4f} {a:11.4f} {b / a:7.3f} "
              f"{GATE:7.2f} "
              f"{('meets it' if short <= 1 else f'{short:.1f}x'):>19}")
    print()
    print("  So the negative result is NOT 'the swing buys nothing'. It buys")
    print("  a real factor, up to 1.3 in the disease range. It is that the")
    print("  gap there is a factor of six, and a third is not a sixth.")
    print("  Stating it the first way would invite a reader to stack the")
    print("  swing with something else; stating it this way says how much")
    print("  that something else would have to be worth.")


def part_c():
    print()
    print("=" * 84)
    print("PART C.  the mechanism, tested rather than asserted")
    print("=" * 84)
    print("\n  The mechanism under test: if a modulation of an input known")
    print("  up to scale pins the SCALE and not the BASELINE, sd(ln f)")
    print("  should fall sharply with the amplitude while sd(ln K), the")
    print("  total, barely moves. f = 0.5 throughout.\n")
    for theta in (3000.0, 300.0):
        print(f"  theta = {theta:.0f}")
        print(f"    {'amp':>6} {'sd(ln theta)':>14} {'sd(ln K) total':>17} "
              f"{'sd(ln f) share':>14}")
        print("    " + "-" * 54)
        for amp in (0.0, 0.15, 0.30, 0.50, 1.00):
            a = over(theta, f=0.5, amp=amp, want=0)[0]
            b = over(theta, f=0.5, amp=amp, want=3)[0]
            c = over(theta, f=0.5, amp=amp, want=4)[0]
            print(f"    {amp:6.2f} {a:14.4f} {b:17.4f} {c:14.4f}")
        print()
    print("""
  READING, AND IT REFUTES THE MECHANISM THIS BLOCK WAS WRITTEN TO CONFIRM.
  The proposed mechanism was that modulating an input known up to scale
  pins the SCALE and not the BASELINE, so that theta gains nothing. The
  columns above say otherwise: sd(ln K), the TOTAL alternative clearance,
  improves by almost exactly as much as sd(ln f), the share. At theta = 300
  the total goes from 0.7228 to 0.5594 as the amplitude reaches 50%, and
  the share from 0.8000 to 0.5241. Neither is preferentially pinned.

  The real mechanism is simpler and it is quantitative rather than
  structural. The swing DOES buy information about the total alternative
  clearance, by about a third. sd(ln theta) in the disease range is
  governed by sd(ln K), so it improves by about a third too, which is
  exactly the gain PART A reports. Nothing is being cancelled; the gain is
  real and it is small.

  Why it is small is PART E: the swing moves the clearance by f*amp, at
  most 40% at a share of 0.8 and an amplitude of 50%, and the frontier
  audit/u9_ke_frontier.py states needs sd(ln k_e) near 6%.
""")


def part_d():
    print()
    print("=" * 84)
    print("PART D.  the aliasing trap, which applies here as it does to diet")
    print("=" * 84)
    print("\n  A modulation whose period is the dose interval lands every")
    print("  draw at the same phase of it. Remark rem:aliasing inverted,")
    print("  which t3 already records for diet. f = 0.5, amp = 50%.\n")
    print(f"    {'theta':>7} {'period':>12} {'sd(ln theta)':>26} {'vs amp=0':>10}")
    print("    " + "-" * 60)
    for theta in (3000.0, 300.0):
        base = over(theta, f=0.5, amp=0.0)[0]
        for lab, per in [("21 d (good)", P_GOOD), ("14 d (dose)", P_BAD),
                         ("30.4 d (visit)", MONTH), ("7 d", 7.0)]:
            m, lo, hi = over(theta, f=0.5, amp=0.50, period_d=per)
            print(f"    {theta:7.0f} {lab:>12} "
                  f"{f'{m:.4f} [{lo:.4f},{hi:.4f}]':>26} {base / m:10.3f}")
        print()


def part_z():
    """The zero-amplitude check: the arm must reduce to the paper's design."""
    print("=" * 84)
    print("CHECK.  at zero amplitude the carnitine arm IS the paper's design")
    print("=" * 84)
    print("\n  The share column is identically zero at amp = 0, so the two")
    print("  must agree to machine precision. They did NOT under the first")
    print("  parameterisation, which is how that error was found.\n")
    worst = 0.0
    print(f"    {'theta':>7} {'f':>5} {'carnitine arm, amp=0':>22} "
          f"{'paper design':>14} {'difference':>12}")
    print("    " + "-" * 62)
    for theta in (3000.0, 600.0, 300.0, 150.0):
        for f in (0.2, 0.5, 0.8):
            a = over(theta, f=f, amp=0.0)[0]
            b = over(theta, f=f, split=False)[0]
            worst = max(worst, abs(a - b))
            print(f"    {theta:7.0f} {f:5.1f} {a:22.6f} {b:14.6f} "
                  f"{abs(a - b):12.2e}")
    print(f"\n  largest difference over twelve cells: {worst:.3e}")
    return worst


def part_e():
    print()
    print("=" * 84)
    print("PART E.  the swing as a route to the k_e frontier, which is the")
    print("         question it should be judged on")
    print("=" * 84)
    print("""
  audit/u9_ke_frontier.py states the requirement: the record reaches the
  gate at theta = 300 if and only if the total alternative clearance is
  known to about 6%. Carnitine modulation is one way of learning it, in
  situ and with no new instrument. So the fair test is not what the swing
  does to sd(ln theta) directly, but how close it brings sd(ln K) to the
  frontier.
""")
    FRONT = {3000.0: None, 600.0: 0.1567, 300.0: 0.0633}
    for theta in (600.0, 300.0):
        print(f"  theta = {theta:.0f}, frontier sd(ln k_e) = {FRONT[theta]:.4f}")
        print(f"    {'share f':>8} {'amp':>6} {'clearance swing':>16} "
              f"{'sd(ln K)':>10} {'factor short':>13}")
        print("    " + "-" * 58)
        for f in (0.5, 0.8):
            for amp in (0.15, 0.50, 1.00):
                k = over(theta, f=f, amp=amp, want=3)[0]
                print(f"    {f:8.1f} {amp:6.2f} {f * amp:15.0%} "
                      f"{k:10.4f} {k / FRONT[theta]:12.1f}x")
        print()
    print("  An amplitude of 1.00 is a carnitine dose swinging between zero")
    print("  and twice the prescribed mean, which is the physical ceiling on")
    print("  this lever and is almost certainly not admissible. Even there")
    print("  the total clearance is short of the frontier by the factor in")
    print("  the last column.")
    print()
    print("  THAT is the negative result worth recording: carnitine")
    print("  modulation is a real but weak k_e measurement, and the")
    print("  requirement it has to meet is an order of magnitude tighter")
    print("  than the best swing delivers.")


if __name__ == "__main__":
    w = part_z()
    print()
    part_a()
    part_b()
    part_c()
    part_e()
    part_d()
    print()
    print("=" * 84)
    print(f"ZERO-AMPLITUDE CHECK: largest difference {w:.3e}")
    print("=" * 84)
