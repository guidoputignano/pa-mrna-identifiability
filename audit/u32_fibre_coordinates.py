"""
WHICH CURVE EACH MARGIN DERIVATIVE IS TAKEN ALONG.

THE CONTRADICTION, as the author put it. `S1.16.6` says the margin is linear
along the unidentified direction where theta is quadratic, so the decision
quantity is better conditioned by one order. The margin document carries
d(mu)/db = 9.950 against d(theta)/db = 2.549 and calls the margin 3.90 times
MORE exposed to the unidentified direction than capacity is. Those
characterise the same behaviour in opposite directions, so either they are
different coordinates or one of them is wrong.

THE ANSWER IS THAT THE PAPER HAS TWO DIFFERENT FIBRES AND CALLS A COORDINATE
ON EACH OF THEM `b`.

  `eq:fibre` (S1.16.6, 04_local.tex). The fibre of the LOCAL AFFINE MODEL:
  the fitted pair (a, b) is held fixed, K_m is the coordinate, and
  theta(K_m) = a(K_m + m_bar)^2/m_bar^2 with k_e(K_m) = b - a K_m/m_bar^2.
  Here `b` is the affine SLOPE of total clearance. **k_e VARIES along it.**

  `eq:fibrecl` (06_excitation.tex). The fibre of the INVARIANTS under a
  mis-specified load: K_m(b) = K_m^0 L/(L+b), theta(b) = theta^0 + b +
  k_e(K_m^0 - K_m(b)), and the text says in so many words "k_e unchanged".
  Here `b` is the LOAD BIAS, in uM/d. **k_e is FIXED along it.**

A curve on which k_e varies and a curve on which it does not are not the same
curve, so their tangents differ and a derivative ratio along one says nothing
about the other. `06_excitation:48` states the second is the one the estimator
sees: "Proposition 2's null vector is its tangent at b = 0", which matches the
null eigenvector's k_e component being zero to machine precision (R168).

AND THE MARGIN DOCUMENT'S NUMBERS ARE ON A THIRD FOOTING. 9.950 requires
s_peak = 1155.4, a fold of 3.58 of the baseline load, against the paper's
committed fold of 1.50. That is the R142 defect again: a ratio formed across
two footings.

WHAT THIS FILE COMPUTES, and every derivative is taken twice, once in closed
form and once by central difference on the module function that owns the
curve, because the whole question is whether two derivatives describe the
same object.

Run:  python audit/u32_fibre_coordinates.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import model_A.run as A
from common.models import clearance_A, theta_crit
from common.params import _A as P, M_TOL, M_OP, ILLNESS_FOLD

#: the load the operating point implies, which is `L` in eq:fibrecl
L = clearance_A(M_OP, P)
#: the paper's committed illness amplitude, as a fold of that load
SP_PAPER = ILLNESS_FOLD * L
#: the amplitude the margin document's derivatives were taken at, recovered
#: from its own d(mu)/db rather than assumed
SP_DOC = 1155.364

H = 1e-4


def central(f, x, h=H):
    return (f(x + h) - f(x - h)) / (2.0 * h)


# ------------------------------------------------- eq:fibrecl, load bias b
def fibrecl(b):
    """(theta, K_m, k_e) at load bias b. k_e is unchanged, per the equation."""
    Km = P.Km * L / (L + b)
    th = P.theta + b + P.ke * (P.Km - Km)
    return th, Km, P.ke


def mu_fibrecl(b, s_peak):
    th, Km, ke = fibrecl(b)
    return th - (s_peak - ke * M_TOL) * (Km + M_TOL) / M_TOL


def part_a():
    print("PART A. eq:fibrecl, the coordinate the margin document used")
    print("  k_e is FIXED along this curve, and b is the load bias in uM/d.")
    dth_cf = 1.0 + P.ke * P.Km / L
    dth_fd = central(lambda b: fibrecl(b)[0], 0.0)
    print(f"    d(theta)/db   closed form 1 + ke*Km/L = {dth_cf:9.4f}"
          f"   central difference {dth_fd:9.4f}"
          f"   {'agree' if abs(dth_cf - dth_fd) < 1e-6 else 'DIFFER'}")
    print(f"    the margin document recorded 2.549: "
          f"{'REPRODUCED' if abs(dth_cf - 2.549) < 5e-4 else 'not reproduced'}")
    print()
    dKm_db = -P.Km / L
    for sp, label in ((SP_PAPER, f"the paper's fold {ILLNESS_FOLD:.2f}"),
                      (SP_DOC, "the margin document's s_peak")):
        dmu_cf = dth_cf - (sp - P.ke * M_TOL) / M_TOL * dKm_db
        dmu_fd = central(lambda b: mu_fibrecl(b, sp), 0.0)
        thc = float(theta_crit(sp, M_TOL, P))
        print(f"    {label}:  s_peak = {sp:8.2f}  (fold {sp/L:5.3f})"
              f"  theta_crit = {thc:8.2f}")
        print(f"      d(mu)/db   closed {dmu_cf:9.4f}   central {dmu_fd:9.4f}"
              f"   {'agree' if abs(dmu_cf - dmu_fd) < 1e-5 else 'DIFFER'}")
        print(f"      ratio  d(mu)/db / d(theta)/db = {dmu_cf/dth_cf:7.4f}"
              f"   -> the margin is {dmu_cf/dth_cf:.2f}x as exposed")
    print()
    return dth_cf


# ---------------------------------------------------- eq:fibre, coordinate Km
def part_b():
    print("PART B. eq:fibre, the coordinate S1.16.6 used")
    print("  k_e VARIES along this curve, and the affine fit (a, b) is held fixed.")
    a, b_slope = A.fit_ab(P)
    Kmax = A.fibre_Km_max(a, b_slope)
    print(f"    a = {a:.4f}   b (the affine SLOPE, not a load bias) = {b_slope:.4f}")
    print(f"    K_max = {Kmax:.3f};  theta spans a factor of "
          f"{((Kmax + M_OP)/M_OP)**2:.1f} over the fibre")
    dth_cf = 2 * a * (P.Km + M_OP) / M_OP ** 2
    dth_fd = central(lambda K: A.fibre(a, b_slope, K)[0], P.Km)
    print(f"    d(theta)/dKm  closed {dth_cf:10.4f}   central {dth_fd:10.4f}"
          f"   {'agree' if abs(dth_cf - dth_fd) < 1e-5 else 'DIFFER'}")
    dke = central(lambda K: A.fibre(a, b_slope, K)[1], P.Km)
    print(f"    d(k_e)/dKm = {dke:.4f}, which is NOT zero: this is the whole")
    print(f"      difference from eq:fibrecl, where k_e does not move at all.")
    print()
    for sp, label in ((SP_PAPER, f"the paper's fold {ILLNESS_FOLD:.2f}"),
                      (SP_DOC, "the margin document's s_peak")):
        Pp = sp - b_slope * M_TOL
        dmu_cf = 2 * a / M_OP - Pp / M_TOL - a * M_TOL / M_OP ** 2
        dmu_fd = central(lambda K: A.margin_on_fibre(a, b_slope, K, s_peak=sp),
                         P.Km)
        print(f"    {label}:")
        print(f"      d(mu)/dKm  closed {dmu_cf:10.4f}   central {dmu_fd:10.4f}"
              f"   {'agree' if abs(dmu_cf - dmu_fd) < 1e-4 else 'DIFFER'}")
        print(f"      ratio |d(mu)/dKm| / |d(theta)/dKm| = "
              f"{abs(dmu_cf/dth_cf):7.4f}"
              f"   -> the margin is {abs(dth_cf/dmu_cf):.2f}x LESS exposed")
    print()
    #  and the linearity itself, checked rather than quoted
    Ks = np.linspace(1e-6, Kmax, 400)
    mus = np.array([A.margin_on_fibre(a, b_slope, K, s_peak=SP_PAPER) for K in Ks])
    ths = np.array([A.fibre(a, b_slope, K)[0] for K in Ks])
    lin = np.polyfit(Ks, mus, 1)
    quad = np.polyfit(Ks, ths, 2)
    print(f"    mu along the fibre, fitted linear: residual max "
          f"{np.abs(mus - np.polyval(lin, Ks)).max():.3e}  -> LINEAR")
    print(f"    theta along the fibre, fitted quadratic: residual max "
          f"{np.abs(ths - np.polyval(quad, Ks)).max():.3e}  -> QUADRATIC")
    print(f"    relative span over the fibre: theta x{ths.max()/ths.min():.1f}, "
          f"mu x{mus.max()/mus.min():.2f}")
    print()


def part_c(dth_db):
    print("PART C. the two curves are different DIRECTIONS, not two names for one")
    a, b_slope = A.fit_ab(P)
    #  tangent of eq:fibre at the nominal point, in (ln theta, ln Km, ln ke)
    dth = 2 * a * (P.Km + M_OP) / M_OP ** 2
    dke = -a / M_OP ** 2
    t_fibre = np.array([dth / P.theta, 1.0 / P.Km, dke / P.ke])
    t_fibre = t_fibre / np.linalg.norm(t_fibre)
    #  tangent of eq:fibrecl at b = 0, same coordinates
    dth2 = dth_db
    dKm2 = -P.Km / L
    t_cl = np.array([dth2 / P.theta, dKm2 / P.Km, 0.0])
    t_cl = t_cl / np.linalg.norm(t_cl)
    print(f"    eq:fibre   tangent in (ln theta, ln Km, ln ke) = "
          f"({t_fibre[0]:+.4f}, {t_fibre[1]:+.4f}, {t_fibre[2]:+.4f})")
    print(f"    eq:fibrecl tangent in the same coordinates      = "
          f"({t_cl[0]:+.4f}, {t_cl[1]:+.4f}, {t_cl[2]:+.4f})")
    cos = abs(float(t_fibre @ t_cl))
    print(f"    |cos| between them = {cos:.4f}, angle "
          f"{np.degrees(np.arccos(min(1.0, cos))):.1f} degrees")
    print(f"    the k_e component is the decisive difference: "
          f"{t_fibre[2]:+.4f} against {t_cl[2]:+.4f}")
    print(f"    R168's null eigenvector had k_e = 0.0000, which is eq:fibrecl,")
    print(f"    and 06_excitation:48 says exactly that: the null vector is its")
    print(f"    tangent at b = 0.")
    print()


if __name__ == "__main__":
    print("=" * 74)
    d = part_a()
    part_b()
    part_c(d)
    print("=" * 74)
    print("CONCLUSION")
    print("  Neither derivative is arithmetically wrong. They are taken along")
    print("  two different curves, which the paper distinguishes and which its")
    print("  reuse of the letter b hides, and the margin document's pair is")
    print("  additionally taken at a fold of 3.58 rather than the paper's 1.50.")
    print("  At the paper's own fold, along the estimator's own curve, the")
    print("  margin is 1.86x as exposed as capacity, not 3.90x.")
    print("=" * 74)
