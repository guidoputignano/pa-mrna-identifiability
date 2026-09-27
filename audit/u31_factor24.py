"""
WHETHER THE FACTOR OF 24 IS DERIVABLE.

THE CLAIM. `01_question.tex:69` and `A1_open.tex:126` both say that the
critical capacity is "an operational definition built on two chosen
constants, and across their defensible range it moves by a factor of 24".
`S1.20` calls the claim it supports "the second of the paper's four claims
and the reason it ends in a data specification rather than a dosing rule".

WHY THIS FILE EXISTS. R178 found that the factor appears exactly twice, that
`01_question` sources it to a 49-page section, and that the section derives
nothing of the kind. The derived threshold sensitivities in the registry are
1.56 over a 32-fold range of `K_c`, 1.37 over an eightfold range of `M_i`
and 1.6 per cent for the calibration, so nothing in print is within an order
of magnitude of 24.

The author set three outcomes. Derivable and never written down, so derive
and wire it. Derivable and far smaller, so restate the claim on the real
number. Or no derivation, so withdraw the claim and put `S1.16.6`'s window
in its place.

THE FOOTING. `theta_crit` is called through the functions that own it,
`common.models.theta_crit` for the Model A threshold of equation S34 and
`model_C.run.theta_crit_intrinsic` for Proposition 2's, never reimplemented.
`s_peak` is built as `model_A.run.s_peak_default` builds it, as a fold of
`baseline_load()`, because the registry states the illness amplitude as a
fold and not absolutely.

WHAT "DEFENSIBLE RANGE" MEANS HERE, and it is the whole question. Neither
occurrence of the 24 says what range it used, so this sweeps every range the
repository itself states and reports the factor each one gives, rather than
choosing one and calling the answer derived.

Run:  python audit/u31_factor24.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import model_A.run as A
import model_C.run as C
from common.models import ParamsA, theta_crit
from common.params import _A as P, M_TOL, M_OP, ILLNESS_FOLD, S_BASE_C

#: The ranges the repository states, each with where it states it. A range
#: invented here would make the answer a choice rather than a derivation, so
#: every one of these is sourced.
MTOL_RANGES = {
    "15 to 160, the two ends the m_tol registry row quotes":
        (15.0, 160.0),
    "12 to 320, the sweep fig_A7_margin panel (c) draws":
        (12.0, 320.0),
    "10 to 40, m_bar to the paper's convention":
        (M_OP, M_TOL),
}
FOLD_RANGES = {
    "1.18 to 1.50, Model C's own fold to the paper's convention":
        (1.18, ILLNESS_FOLD),
    "1.18 to 2.27, Model C's fold to a 227% rise read as 2.27x":
        (1.18, 2.27),
    "1.18 to 3.25, Model C's fold to a 325% rise read as 3.25x":
        (1.18, 3.25),
    "1.18 to 4.25, Model C's fold to a rise BY 325%":
        (1.18, 4.25),
    "1.0 to 6.0, the sweep fig_A7_margin panel (b) draws":
        (1.0, 6.0),
}
TARGET = 24.0


def thc_A(m_tol, fold):
    """Equation S34's threshold, through the function that owns it."""
    return float(theta_crit(fold * A.baseline_load(), m_tol, P))


def factor(vals):
    v = np.array([x for x in vals if np.isfinite(x) and x > 0])
    if v.size < 2:
        return float("nan"), float("nan"), float("nan")
    return float(v.max() / v.min()), float(v.min()), float(v.max())


def part_a():
    print("PART A. theta_crit^A over m_tol alone, at the paper's fold of "
          f"{ILLNESS_FOLD:.2f}")
    print(f"  the anchor: theta_crit at m_tol = {M_TOL:.0f} is "
          f"{thc_A(M_TOL, ILLNESS_FOLD):.1f}")
    for name, (lo, hi) in MTOL_RANGES.items():
        grid = np.linspace(lo, hi, 2001)
        f, mn, mx = factor([thc_A(x, ILLNESS_FOLD) for x in grid])
        print(f"    {name}")
        print(f"      theta_crit in [{mn:8.1f}, {mx:8.1f}]   factor {f:6.2f}")
    print()


def part_b():
    print("PART B. theta_crit^A over the illness fold alone, at "
          f"m_tol = {M_TOL:.0f}")
    print("  NOTE, and it is the reason the fold cannot supply a large factor:")
    print("  s_peak enters S34 affinely with a positive offset, so the fold")
    print("  moves theta_crit LESS than proportionally.")
    for name, (lo, hi) in FOLD_RANGES.items():
        grid = np.linspace(lo, hi, 2001)
        f, mn, mx = factor([thc_A(M_TOL, x) for x in grid])
        print(f"    {name}")
        print(f"      theta_crit in [{mn:8.1f}, {mx:8.1f}]   factor {f:6.2f}")
    print()


def part_c():
    print("PART C. the JOINT range, which is what 'two chosen constants' means")
    best = None
    for mn_name, (ml, mh) in MTOL_RANGES.items():
        for fn_name, (fl, fh) in FOLD_RANGES.items():
            mg = np.linspace(ml, mh, 301)
            fg = np.linspace(fl, fh, 301)
            vals = [thc_A(x, y) for x in mg for y in fg]
            f, mn, mx = factor(vals)
            flag = "  <-- nearest 24" if best is None or abs(f - TARGET) < abs(best[0] - TARGET) else ""
            if best is None or abs(f - TARGET) < abs(best[0] - TARGET):
                best = (f, mn_name, fn_name, mn, mx)
            print(f"    m_tol {ml:6.1f}-{mh:6.1f}  fold {fl:4.2f}-{fh:4.2f}"
                  f"   theta_crit [{mn:8.1f}, {mx:8.1f}]   factor {f:7.2f}{flag}")
    print()
    print(f"  the widest defensible box gives a factor of {best[0]:.2f}")
    print(f"    m_tol: {best[1]}")
    print(f"    fold : {best[2]}")
    print(f"  against the {TARGET:.0f} in print: "
          f"{'REPRODUCED' if abs(best[0] - TARGET) < 0.5 else 'NOT REPRODUCED'}")
    print()
    return best


def part_d():
    print("PART D. Proposition 2's INTRINSIC threshold, which drops m_tol")
    print("  the claim is about a definition resting on TWO conventions, so")
    print("  the intrinsic one is the control: it rests on s_peak alone.")
    base = S_BASE_C
    vals = {}
    for name, (lo, hi) in FOLD_RANGES.items():
        out = []
        for fold in np.linspace(lo, hi, 60):
            try:
                out.append(float(C.theta_crit_intrinsic(C.PC, base * fold / ILLNESS_FOLD)))
            except Exception:
                pass
        f, mn, mx = factor(out)
        vals[name] = f
        print(f"    {name}")
        print(f"      theta_crit in [{mn:8.1f}, {mx:8.1f}]   factor {f:6.2f}")
    print()
    return vals


def part_e(best):
    print("PART E. what range WOULD give 24, stated so the claim can be")
    print("  checked against the literature rather than against my choices.")
    #  solve for the m_tol upper end that gives exactly 24 at each fold range
    for fn_name, (fl, fh) in FOLD_RANGES.items():
        lo_m = 15.0
        f_at = lambda mh: factor([thc_A(x, y)
                                  for x in np.linspace(lo_m, mh, 200)
                                  for y in np.linspace(fl, fh, 60)])[0]
        hi = None
        for cand in np.linspace(20.0, 5000.0, 400):
            if f_at(cand) >= TARGET:
                hi = cand
                break
        if hi is None:
            print(f"    fold {fl:4.2f}-{fh:4.2f}: m_tol would have to exceed "
                  f"5000 uM to reach {TARGET:.0f}")
        else:
            print(f"    fold {fl:4.2f}-{fh:4.2f}: reaches {TARGET:.0f} only at "
                  f"m_tol up to {hi:.0f} uM, against a convention of "
                  f"{M_TOL:.0f} and m_bar of {M_OP:.0f}")
    print()


def part_f():
    print("PART F. the WINDOW, which is what S1.16.6 actually supports")
    a0, b0 = A.fit_ab(P)
    base = A.baseline_load()
    for x in (15.0, M_TOL, 160.0):
        lo, hi = A.margin_window(a0, b0, x)
        w = 100.0 * (hi - lo) / (0.5 * (lo + hi))
        print(f"    m_tol = {x:6.1f}:  s_peak undetermined on "
              f"({lo:7.1f}, {hi:7.1f}) = ({lo/base:4.2f}x, {hi/base:4.2f}x) "
              f"baseline,  width {w:5.1f}% of centre")
    print()


if __name__ == "__main__":
    print(__doc__.split("Run:")[0].strip()[:0] or "", end="")
    print("=" * 72)
    part_a()
    part_b()
    best = part_c()
    part_d()
    part_e(best)
    part_f()
    print("=" * 72)
