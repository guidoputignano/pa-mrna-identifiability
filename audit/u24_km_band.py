"""
Derive or withdraw `06_excitation.tex:103`: does sd(ln K_m) fall from 0.587 to
0.209 at m/K_m = 0.01 when the increment is carried, with sd(ln theta) not
improving at all?

WHY THIS SCRIPT EXISTS. R141 established that NEITHER literal has a generator.
Their only trace in the repository is a sentence in `d0c_mechanism.py`'s
docstring, and reading that script confirms it computes three constructions,
none of which is the one the paper's sentence describes: the real dose train,
the same levels with the increment reassigned, and a single marker level at the
MEDIAN of the record. **There is no idealised band at m/K_m = 0.01 anywhere in
the repository.** A number whose only source is a comment cannot be wired,
because the guard would then vouch for a value someone typed by hand.

THE HONEST DIFFICULTY, stated before any number. "An idealised band with Delta
sweeping monotonically with m" does not name a construction. It leaves the
band's width, its number of points and the increment's range free, and
different choices give different answers. So a script that searched those
choices until 0.587 and 0.209 appeared would not be deriving the claim, it
would be reverse-engineering it, and it would certify whatever it found.

WHAT THIS DOES INSTEAD. It fixes the construction by the most natural reading
available, states every free choice in the open, computes what that gives, and
then SWEEPS each free choice to show whether the claim's CONTENT survives them.
The content is two comparative statements, and they are what a referee would
check:
    (1) carrying the increment improves sd(ln K_m) SUBSTANTIALLY;
    (2) carrying the increment does NOT improve sd(ln theta).
If those hold across the sweep, the sentence is sound and its two literals are
restated at whatever the stated construction gives. If they do not, the
sentence is withdrawn. Either way the literals stop being orphans.

The Fisher construction, priors and noise are taken from `d0c_mechanism.py`
UNCHANGED, so this is comparable with the script the sentence was written
beside rather than a second convention.

Deterministic: this is a Fisher computation with no sampling, so there is no
seed spread to report. That is stated rather than omitted, because a missing
range is otherwise indistinguishable from a suppressed one.

Run:  python audit/u24_km_band.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import design.run as D
import model_A.run as LOC

P = D.P
Lb = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def sd(levels, incs, carry=True):
    """
    sd(ln theta) and sd(ln K_m), IDENTICAL to d0c_mechanism.py's `sd`.

    Copied rather than imported so that a later edit to d0c cannot silently
    move this result, and so the two can be compared as written.
    """
    pri = [0.8 ** 2, 0.8 ** 2, (D.LOAD_FRAC * Lb) ** 2, 0.8 ** 2]
    F = np.linalg.inv(np.diag(pri))
    for m, d in zip(levels, incs):
        eff = P.theta + (d if carry else 0.0)
        q = np.array([P.theta * m / (P.Km + m),
                      -eff * P.Km * m / (P.Km + m) ** 2, 1.0, P.ke * m])
        F = F + np.outer(q, q) / SIG ** 2
    C = np.linalg.inv(F)
    return float(np.sqrt(C[0, 0])), float(np.sqrt(C[1, 1]))


def band(centre_ratio=0.01, half_decade=0.5, n=365, dmax_mult=2.0):
    """
    An idealised band centred at m/K_m = `centre_ratio`, with the increment
    sweeping monotonically with m.

    EVERY FREE CHOICE IS AN ARGUMENT, and the sweep below moves each one.
      centre_ratio  the operating point the sentence names, m/K_m
      half_decade   the band's half-width in decades of m
      n             how many points, i.e. the record's length
      dmax_mult     the increment sweeps 0 to dmax_mult * theta
    Monotone means the smallest increment sits at the smallest level, which is
    the reading "sweeping monotonically with m" admits.
    """
    c = centre_ratio * P.Km
    levels = np.logspace(np.log10(c) - half_decade,
                         np.log10(c) + half_decade, n)
    incs = np.linspace(0.0, dmax_mult * P.theta, n)
    return levels, incs


def row(**kw):
    lv, dd = band(**kw)
    t_off, k_off = sd(lv, dd, carry=False)
    t_on, k_on = sd(lv, dd, carry=True)
    return t_off, t_on, k_off, k_on


PAPER_KM = (0.587, 0.209)

if __name__ == "__main__":
    print("=" * 94)
    print("06_excitation.tex:103, DERIVE OR WITHDRAW: sd(ln K_m) 0.587 -> 0.209 "
          "at m/K_m = 0.01")
    print("=" * 94)
    print("""
  The claim has TWO parts and they are checked separately:
    (1) carrying the increment improves sd(ln K_m) substantially
    (2) carrying the increment does NOT improve sd(ln theta)
  The paper states part (1) with the literals 0.587 and 0.209.
""")

    print("-" * 94)
    print("PART A. The stated construction, every free choice named")
    print("-" * 94)
    base = dict(centre_ratio=0.01, half_decade=0.5, n=365, dmax_mult=2.0)
    for k, v in base.items():
        print(f"    {k:<14} {v}")
    t_off, t_on, k_off, k_on = row(**base)
    print()
    print(f"    sd(ln K_m)    dropped {k_off:8.4f}   carried {k_on:8.4f}   "
          f"ratio {k_on / k_off:6.3f}")
    print(f"    sd(ln theta)  dropped {t_off:8.4f}   carried {t_on:8.4f}   "
          f"ratio {t_on / t_off:6.3f}")
    print()
    print(f"    the paper's pair for sd(ln K_m): {PAPER_KM[0]} -> {PAPER_KM[1]}")
    d0 = abs(k_off - PAPER_KM[0])
    d1 = abs(k_on - PAPER_KM[1])
    print(f"    departure from it: {d0:.4f} and {d1:.4f}")
    print(f"    reproduces the paper's literals to 3 decimals: "
          f"{'YES' if max(d0, d1) < 5e-4 else 'NO'}")

    print()
    print("-" * 94)
    print("PART B. Does the CONTENT survive the free choices? Each row moves one.")
    print("-" * 94)
    print(f"  {'variation':>34} {'sdKm off':>10} {'sdKm on':>9} {'ratio':>7}"
          f" {'sdTh off':>10} {'sdTh on':>9} {'ratio':>7}")
    sweeps = [("half_decade", (0.25, 0.5, 1.0, 1.5)),
              ("n", (50, 120, 365, 1000)),
              ("dmax_mult", (0.5, 1.0, 2.0, 4.0)),
              ("centre_ratio", (0.003, 0.01, 0.03, 0.1))]
    verdicts = []
    for name, vals in sweeps:
        for v in vals:
            kw = dict(base)
            kw[name] = v
            t_off, t_on, k_off, k_on = row(**kw)
            rk, rt = k_on / k_off, t_on / t_off
            verdicts.append((rk, rt))
            print(f"  {name + ' = ' + str(v):>34} {k_off:10.4f} {k_on:9.4f} "
                  f"{rk:7.3f} {t_off:10.4f} {t_on:9.4f} {rt:7.3f}")

    print()
    print("-" * 94)
    print("PART C. The verdict on each half of the claim")
    print("-" * 94)
    worst_km = max(r for r, _ in verdicts)
    best_th = min(r for _, r in verdicts)
    print(f"    (1) sd(ln K_m) improves: the WORST ratio over the whole sweep is "
          f"{worst_km:.3f}")
    print(f"        {'HOLDS' if worst_km < 0.9 else 'FAILS'} "
          f"(an improvement at every setting means every ratio < 1)")
    print(f"    (2) sd(ln theta) does not improve: the BEST ratio is {best_th:.3f}")
    print(f"        {'HOLDS' if best_th > 0.95 else 'FAILS'} "
          f"(no improvement means no ratio meaningfully < 1)")
    print()
    print("    Deterministic Fisher computation. No sampling, so no seed spread "
          "exists to report.")
