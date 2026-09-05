"""
Consistency checks on the model definitions.

These are not unit tests of numerical accuracy; they assert the structural
claims made in the paper, so that a change to the model that breaks a claim
fails loudly instead of silently changing a figure.

Run:  pytest tests/ -q
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dataclasses import replace

import numpy as np
import pytest

from common.models import (ParamsA, ParamsB, ParamsC, g, h,
                           clearance_A, clearance_B, clearance_C_quasistatic,
                           alt_flux_C, rhs_C,
                           steady_state, equilibria_nonmonotone,
                           sensitivity_A, condition_number, theta_crit)


# ---------------------------------------------------------------- kinetics
def test_g_saturates():
    p = ParamsA()
    assert g(0.0, p.Km) == 0.0
    assert g(1e9, p.Km) == pytest.approx(1.0, abs=1e-6)
    assert g(p.Km, p.Km) == pytest.approx(0.5)


# ---------------------------------------------------------------- Model A
def test_A_clearance_unbounded():
    """Model A's known defect: clearance grows without limit, so no crisis."""
    p = ParamsA()
    assert clearance_A(1e7, p) > 10 * p.theta


def test_A_equilibrium_always_exists():
    """Follows from unbounded clearance."""
    p = ParamsA()
    for s in (10.0, 1e3, 1e5, 1e7):
        assert np.isfinite(steady_state(clearance_A, p, s))


def test_A_rank_needs_three_levels():
    """Proposition: three distinct metabolite levels are needed."""
    p = ParamsA()
    for n, levels in [(1, [10.0]), (2, [10.0, 12.0]), (3, [10.0, 12.0, 14.0])]:
        S = np.array([sensitivity_A(m, p) for m in levels])
        assert np.linalg.matrix_rank(S, tol=1e-10) == min(n, 3)


def test_A_conditioning_scales_inverse_square():
    """Condition number should grow roughly as 1/excursion^2."""
    p = ParamsA()
    m0 = 10.0
    c = []
    for sp in (0.05, 0.5):
        c.append(condition_number(sensitivity_A, p,
                                  m0 * np.array([1.0, 1 + sp, 1 + 2 * sp])))
    ratio = c[0] / c[1]
    assert 50 < ratio < 200, f"expected ~100x, got {ratio:.1f}x"


def test_A_degeneracy_two_patients_indistinguishable():
    """Two very different patients with the same steady state."""
    from scipy.optimize import brentq
    Km = 100.0
    pA = ParamsA(3000.0, Km, 5.0)
    s = clearance_A(10.0, pA)
    keB = brentq(lambda k: clearance_A(10.0, ParamsA(1000.0, Km, k)) - s, 1e-6, 1e4)
    pB = ParamsA(1000.0, Km, keB)
    assert steady_state(clearance_A, pA, s) == pytest.approx(
           steady_state(clearance_A, pB, s), rel=1e-9)
    assert pA.theta / pB.theta == pytest.approx(3.0)


def test_A_anchor_recovers_parameters_exactly():
    """The anchor proposition, in the noise-free case."""
    from scipy.optimize import least_squares
    Km, ke = 100.0, 5.0
    thetas = [0.0, 2600.0, 1300.0, 450.0]      # patient 0 is null
    loads = [420.0, 360.0, 480.0, 445.0]
    ms = [steady_state(clearance_A, ParamsA(t, Km, ke), s)
          for t, s in zip(thetas, loads)]

    ke_hat = loads[0] / ms[0]                   # anchor (i)
    assert ke_hat == pytest.approx(ke, rel=1e-9)

    i, s2 = 1, loads[1] * 1.6                   # anchor (ii)
    m2 = steady_state(clearance_A, ParamsA(thetas[i], Km, ke), s2)
    th, Km_hat = least_squares(
        lambda v: [loads[i] - ke_hat * ms[i] - v[0] * ms[i] / (v[1] + ms[i]),
                   s2 - ke_hat * m2 - v[0] * m2 / (v[1] + m2)],
        [1000.0, 50.0], bounds=([0, 1e-3], [1e6, 1e5])).x
    assert Km_hat == pytest.approx(Km, rel=1e-6)
    assert th == pytest.approx(thetas[i], rel=1e-6)


# ---------------------------------------------------------------- Model B
def test_B_clearance_bounded_by_ceiling():
    p = ParamsB()
    assert clearance_B(1e9, p) == pytest.approx(p.theta + p.ke, rel=1e-6)


def test_B_no_equilibrium_above_ceiling():
    p = ParamsB()
    ceiling = p.theta + p.ke
    assert np.isfinite(steady_state(clearance_B, p, ceiling * 0.9))
    assert not np.isfinite(steady_state(clearance_B, p, ceiling * 1.1))


def test_B_is_harder_to_identify_than_A():
    """The biologically faithful correction costs identifiability."""
    from common.models import sensitivity_B
    m0, sp = 10.0, 0.1
    lv = m0 * np.array([1.0, 1 + sp, 1 + 2 * sp, 1 + 3 * sp])
    cA = condition_number(sensitivity_A, ParamsA(), lv[:3])
    cB = condition_number(sensitivity_B, ParamsB(), lv)
    assert cB > 10 * cA


# ---------------------------------------------------------------- Model C
def test_C_clearance_is_non_monotone():
    """Depletion makes clearance rise, peak, then fall. This is what makes
    bistability possible, and distinguishes C from B."""
    p = ParamsC()
    m = np.linspace(1.0, 2000.0, 4000)
    C = np.array([clearance_C_quasistatic(x, p) for x in m])
    peak = int(np.argmax(C))
    assert 0 < peak < len(m) - 1, "clearance should peak in the interior"
    assert C[-1] < C[peak] * 0.9, "clearance should decline after the peak"


def test_C_has_two_equilibria_at_moderate_load():
    p = ParamsC()
    grid = np.linspace(1e-3, 2500.0, 8000)
    roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, p, 600.0, grid)
    assert len(roots) == 2, f"expected 2 equilibria, got {len(roots)}"
    assert stab[0] and not stab[1], "lower should be stable, upper unstable"


def test_C_loses_equilibria_above_the_fold():
    p = ParamsC()
    grid = np.linspace(1e-3, 2500.0, 8000)
    roots, _ = equilibria_nonmonotone(clearance_C_quasistatic, p, 900.0, grid)
    assert len(roots) == 0


def test_C_intrinsic_theta_crit_is_where_two_equilibria_become_none():
    """
    Proposition (intrinsic threshold): theta_crit_intrinsic(p, s) is the value
    of theta at which the attractor and the separatrix merge at a saddle-node.
    Just above it there are two equilibria; just below it there are none.
    """
    from model_C.run import theta_crit_intrinsic

    p = ParamsC()
    s = 600.0
    thc = theta_crit_intrinsic(p, s)
    assert np.isfinite(thc) and thc > 0, f"expected a positive threshold, got {thc}"

    grid = np.linspace(1e-3, 30000.0, 90000)

    def n_equilibria(theta):
        q = replace(p, theta=float(theta))
        roots, _ = equilibria_nonmonotone(clearance_C_quasistatic, q, s, grid)
        return len(roots)

    assert n_equilibria(thc * 1.01) == 2, "above the threshold: attractor and separatrix"
    assert n_equilibria(thc * 0.99) == 0, "below the threshold: no equilibrium at all"


def test_C_intrinsic_theta_crit_does_not_use_a_tolerable_level():
    """
    The point of the construction: it takes the load and the parameters, and
    nothing else. Changing the clinical convention m_bar cannot change it,
    because m_tol is not one of its arguments. Contrast theta_crit() above.
    """
    from model_C.run import theta_crit_intrinsic
    import inspect

    names = list(inspect.signature(theta_crit_intrinsic).parameters)
    assert "m_bar" not in names
    assert names[:2] == ["p", "s_base"]


def test_C_margin_closes_monotonically_towards_the_threshold():
    """The separatrix-to-attractor distance is the safety margin, and it
    shrinks to zero as theta falls to the intrinsic threshold."""
    from model_C.run import theta_crit_intrinsic

    p = ParamsC()
    s = 600.0
    thc = theta_crit_intrinsic(p, s)
    grid = np.linspace(1e-3, 30000.0, 90000)

    margins = []
    for f in (1.02, 1.2, 1.6):
        q = replace(p, theta=thc * f)
        roots, st = equilibria_nonmonotone(clearance_C_quasistatic, q, s, grid)
        assert len(roots) == 2
        margins.append(roots[1] - roots[0])

    assert margins[0] < margins[1] < margins[2], f"margin should widen: {margins}"
    assert margins[0] < 0.15 * margins[2], "margin should close near the fold"


# ---------------------------------------------------------------- threshold
def test_theta_crit_depends_on_unmeasured_quantities():
    """A sanity check on the definition, and a reminder of its weakness."""
    p = ParamsA()
    t1 = theta_crit(s_peak=500.0, m_tol=40.0, p=p)
    t2 = theta_crit(s_peak=500.0, m_tol=60.0, p=p)
    assert t1 != t2, "theta_crit must depend on the assumed tolerable level"


# ------------------------------------------------- the two-state system
def test_rhs_C_reduces_to_alt_flux_C():
    """
    The paper states Model C as two states and then works with the
    quasi-steady reduction alpha(m). The two must agree, or Section 5's
    propositions are about a curve the stated model does not produce, which
    is exactly the defect this test exists to prevent reopening.

    Under the acylcarnitine-competition form the reduction also preserves the
    constants: the same R0 and the same Kc appear on both sides, with no
    rescaling, so the agreement is exact rather than up to a reparametrisation.
    """
    from scipy.optimize import brentq

    p = ParamsC()
    k = 1.0
    for m in (0.5, 5.0, 50.0, 173.0, 1036.0, 5000.0, 20000.0):
        c = brentq(lambda c: rhs_C(0.0, [m, c], p, lambda t: 0.0, k)[1], 1e-12, 1e12)
        flux = k * c * m / (p.Kc + m)
        assert flux == pytest.approx(alt_flux_C(m, p), rel=1e-10), m


def test_rhs_C_is_callable_with_ParamsC():
    """
    It once read p.k and p.R, which ParamsC does not define, so it raised on
    any call. Nothing called it and nothing caught it.
    """
    out = rhs_C(0.0, [100.0, 50.0], ParamsC(), lambda t: 600.0)
    assert len(out) == 2 and all(np.isfinite(out))


def test_conjugation_alone_gives_a_ceiling_not_a_decline():
    """
    Section 5.1's first claim. If conjugation were the only route out of the
    carnitine pool, dc/dt = 0 would fix the flux at R0 for every m: no vanishing
    at the origin and no decline. The second exit is what buys both.
    """
    p = ParamsC()
    m = np.array([1.0, 10.0, 100.0, 1000.0, 1e5])
    # rho0 -> 0 removes the second exit
    none_ = replace(p, rho0=1e-9)
    flux = alt_flux_C(m, none_)
    assert np.all(np.abs(flux - p.R0) / p.R0 < 0.02), \
        "with no second exit the flux should sit at R0 everywhere"
    # with it, the flux starts at zero and comes back down
    real = alt_flux_C(m, p)
    assert real[0] < 0.1 * real.max(), "must vanish as m does"
    assert real[-1] < 0.1 * real.max(), "and must come back down"
    assert alt_flux_C(np.array([0.0]), p)[0] == 0.0


def test_alt_flux_ceiling_fixes_R0():
    """
    R0 is not free. It is set so the alternative route alone cannot carry the
    baseline load, which is what makes theta_crit positive (Corollary 1). If
    Kc, Mi or rho0 change, R0 has to be re-derived and this fails.
    """
    p = ParamsC()
    m = np.linspace(1e-6, 2e5, 400000)
    assert alt_flux_C(m, p).max() == pytest.approx(500.0, rel=1e-4)


def test_theta_crit_positive_iff_alt_route_cannot_carry_the_load():
    """Corollary 1, both directions."""
    from model_C.run import theta_crit_intrinsic
    p = ParamsC()
    m = np.linspace(1e-6, 2e5, 400000)
    s_base = 600.0
    assert alt_flux_C(m, p).max() < s_base
    assert theta_crit_intrinsic(p, s_base) > 0
    # lift the alternative route above the load and the threshold goes vacuous
    strong = replace(p, R0=p.R0 * 2.0)
    assert alt_flux_C(m, strong).max() > s_base
    assert np.isnan(theta_crit_intrinsic(strong, s_base))


# ================================================================= Model A
#  The fibre of observationally identical patients, and what survives on it.
#
#  A band fixes the affine fit (a, b) of Proposition 1 and nothing more. Both
#  coefficients carry Km, so a known (a, b) leaves a one-parameter family.
#  These tests pin the claims the paper makes about that family.
# =========================================================================
from model_A.run import (fit_ab, fibre, fibre_Km_max, margin,
                         margin_on_fibre, margin_window, s_peak_default,
                         M_OP, M_TOL)


def test_fibre_inverts_the_affine_fit():
    """Every point of the fibre reproduces the same (a, b) it came from."""
    p = ParamsA()
    a, b = fit_ab(p)
    for Km in (5.0, 25.0, 60.0, 100.0, 115.0):
        th, ke = fibre(a, b, Km)
        a2, b2 = fit_ab(ParamsA(th, Km, ke))
        assert a2 == pytest.approx(a, rel=1e-12), Km
        assert b2 == pytest.approx(b, rel=1e-12), Km


def test_fibre_extent_is_wide_with_noiseless_data():
    """
    Proposition 1 says the intercept is free of ke. It is not free of Km, and
    the span that leaves is the point of Remark 4.
    """
    p = ParamsA()
    a, b = fit_ab(p)
    Kmax = fibre_Km_max(a, b)
    assert Kmax == pytest.approx(b * M_OP ** 2 / a, rel=1e-12)
    assert fibre(a, b, Kmax)[1] == pytest.approx(0.0, abs=1e-9), "ke vanishes there"
    lo, hi = fibre(a, b, 1e-9)[0], fibre(a, b, Kmax)[0]
    assert hi / lo > 100.0, f"span was only {hi / lo:.1f}x"
    assert lo < p.theta < hi, "the true patient must lie on its own fibre"


def test_margin_is_linear_in_Km_where_theta_is_quadratic():
    """
    The correction to Section 3.6 item 2. theta is quadratic along the fibre
    and the margin is linear: the quadratic terms cancel identically, so the
    second difference of mu is zero to machine precision while that of theta
    is not.
    """
    p = ParamsA()
    a, b = fit_ab(p)
    Ks = np.linspace(1.0, fibre_Km_max(a, b) * 0.99, 9)
    h = Ks[1] - Ks[0]
    mu = np.array([margin_on_fibre(a, b, K) for K in Ks])
    th = np.array([fibre(a, b, K)[0] for K in Ks])
    d2mu = np.diff(mu, 2)
    d2th = np.diff(th, 2)
    scale = max(abs(mu).max(), 1.0)
    assert np.max(np.abs(d2mu)) / scale < 1e-10, "mu must be linear in Km"
    assert np.max(np.abs(d2th)) > 1e-6 * abs(th).max(), "theta must not be"
    # and the closed form of (18) agrees with going round through (theta, ke)
    P = s_peak_default() - b * M_TOL
    slope = 2 * a / M_OP - P / M_TOL - a * M_TOL / M_OP ** 2
    for K in Ks:
        assert margin_on_fibre(a, b, K) == pytest.approx(slope * K + (a - P),
                                                         rel=1e-10), K


def test_matched_patients_separate_as_dke_times_mtol_minus_mbar():
    """
    Equation (19), exactly. Two patients matched in value at the operating
    point and sharing Km differ in margin by (ke_A - ke_B)(m_tol - m_bar),
    and coincide when the tolerable level is the operating point.
    """
    from scipy.optimize import brentq
    p = ParamsA()
    thB = 1000.0
    keB = brentq(lambda k: clearance_A(M_OP, ParamsA(thB, p.Km, k))
                 - clearance_A(M_OP, p), 1e-9, 1e5)
    for m_tol in (M_OP, 21.0, 40.0, 75.0):
        lhs = (margin(p.theta, p.ke, m_tol=m_tol)
               - margin(thB, keB, m_tol=m_tol))
        rhs = (p.ke - keB) * (m_tol - M_OP)
        assert lhs == pytest.approx(rhs, abs=1e-9), m_tol


def test_null_anchor_closes_the_fibre_without_a_second_steady_state():
    """
    The strengthened form of Proposition 5. One band plus a known ke recovers
    Km and theta exactly. No patient is moved to a second load.
    """
    p = ParamsA()
    a, b = fit_ab(p)
    Km_rec = (b - p.ke) * M_OP ** 2 / a
    th_rec, ke_rec = fibre(a, b, Km_rec)
    assert Km_rec == pytest.approx(p.Km, rel=1e-12)
    assert th_rec == pytest.approx(p.theta, rel=1e-12)
    assert ke_rec == pytest.approx(p.ke, rel=1e-12)


def test_two_patients_with_bands_close_the_system():
    """
    Figure 8's counting assumes one observation per patient. Under the band
    model two patients supply four coefficients for three unknowns, and they
    close it without either of them moving.
    """
    p = ParamsA()
    pA = ParamsA(3000.0, p.Km, p.ke)
    pB = ParamsA(1000.0, p.Km, p.ke)
    a1, b1 = fit_ab(pA)
    a2, b2 = fit_ab(pB)
    Km = (b2 - b1) / (a2 / M_OP ** 2 - a1 / M_OP ** 2)
    ke = b1 - a1 * Km / M_OP ** 2
    assert Km == pytest.approx(p.Km, rel=1e-10)
    assert ke == pytest.approx(p.ke, rel=1e-10)
    assert fibre(a1, b1, Km)[0] == pytest.approx(pA.theta, rel=1e-10)
    assert fibre(a2, b2, Km)[0] == pytest.approx(pB.theta, rel=1e-10)


def test_det_M_factorises_exactly():
    """
    Equation (14) with no affine hypothesis. The narrow-band form g(m_bar)^2
    is its midpoint approximation, and it is not free over a wide band.
    """
    Km = ParamsA().Km
    for m, mp in [(8.0, 14.0), (10.0, 10.5), (5.0, 200.0), (40.0, 90.0)]:
        lhs = g(m, Km) * mp - g(mp, Km) * m
        assert lhs == pytest.approx(g(m, Km) * g(mp, Km) * (mp - m), rel=1e-13)
    # the approximation error grows with the excursion
    approx8 = g(11.0, Km) ** 2 * 6.0
    exact8 = g(8.0, Km) * g(14.0, Km) * 6.0
    assert 0.05 < (approx8 - exact8) / exact8 < 0.11, "about 8% at m=8, m'=14"


def test_load_error_amplification_is_a_flux_ratio_and_diverges():
    """
    The 1.4 of Figure 6 is C(m)/(theta g(m)), not a constant of the model.
    """
    Km, ke, m = 100.0, 5.0, 20.0
    amps = []
    for th in (3000.0, 1500.0, 400.0, 100.0):
        q = ParamsA(th, Km, ke)
        amps.append(clearance_A(m, q) / (th * g(m, Km)))
    assert amps[1] == pytest.approx(1.4, rel=1e-9), "the figure's value"
    assert amps == sorted(amps), "must grow as theta falls"
    assert amps[-1] > 5.0, "and diverge at the severe end"


def test_margin_window_endpoints_match_a_numerical_search():
    """
    The closed form of the undetermined window against a bisection on the
    sign pattern at the two fibre ends. Inside it the two ends disagree, so
    the decision is not determined by the data however wide the band.
    """
    from model_A.run import margin_window, fibre_Km_max, margin_on_fibre
    p = ParamsA()
    a, b = fit_ab(p)
    Kmax = fibre_Km_max(a, b)

    def disagree(sp, m_tol):
        e0 = margin_on_fibre(a, b, 1e-9, m_tol=m_tol, s_peak=sp)
        e1 = margin_on_fibre(a, b, Kmax, m_tol=m_tol, s_peak=sp)
        return np.sign(e0) != np.sign(e1)

    for m_tol in (20.0, 40.0, 80.0):
        lo, hi = margin_window(a, b, m_tol)
        assert not disagree(lo * 0.999, m_tol), m_tol
        assert disagree(0.5 * (lo + hi), m_tol), m_tol
        assert not disagree(hi * 1.001, m_tol), m_tol
        # and the closed form is where mu vanishes at each end
        assert margin_on_fibre(a, b, Kmax, m_tol=m_tol, s_peak=lo) == \
            pytest.approx(0.0, abs=1e-8), m_tol
        assert margin_on_fibre(a, b, 1e-12, m_tol=m_tol, s_peak=hi) == \
            pytest.approx(0.0, abs=1e-6), m_tol


def test_window_narrows_as_m_tol_approaches_the_operating_point():
    """
    The window is what replaces a tolerance on theta. Its width is set by the
    convention m_tol, and it closes as that convention approaches the band the
    patient is actually observed in.
    """
    from model_A.run import margin_window, M_OP
    p = ParamsA()
    a, b = fit_ab(p)
    widths = []
    for m_tol in (15.0, 20.0, 40.0, 80.0, 160.0):
        lo, hi = margin_window(a, b, m_tol)
        widths.append(100.0 * (hi - lo) / (0.5 * (lo + hi)))
    assert widths == sorted(widths), f"must widen with m_tol: {widths}"
    assert widths[0] < 2.0, "narrow when m_tol is close to m_bar"
    assert widths[-1] > 50.0, "wide when it is far"


def test_illness_amplitude_is_a_multiple_of_this_patients_baseline():
    """
    s_peak stated absolutely disagreed across the note. As a fold of baseline
    it is one physiological claim applied to each model's own parameters, and
    at Model C's amplitude the Model A decision is determined and safe.
    """
    from model_A.run import (baseline_load, s_peak_default, margin_window,
                             margin_on_fibre, fibre_Km_max, ILLNESS_FOLD)
    p = ParamsA()
    a, b = fit_ab(p)
    base = baseline_load()
    assert base == pytest.approx(clearance_A(M_OP, p), rel=1e-12)
    assert s_peak_default() == pytest.approx(ILLNESS_FOLD * base, rel=1e-12)
    lo, hi = margin_window(a, b)
    assert s_peak_default() < lo, "Model C's illness sits below the window"
    ends = [margin_on_fibre(a, b, K) for K in (1e-9, fibre_Km_max(a, b))]
    assert all(e > 0 for e in ends), "and the decision there is determined, safe"


# ================================================================= design
def test_relative_precision_of_theta_is_identified_along_the_fibre():
    """
    Pass 2, item 1. theta_hat = a_hat/g(m_bar)^2 is a_hat times a constant,
    so sd(theta_hat)/theta = sd(a_hat)/a and the fibre coordinate cancels.
    Equation (7)'s g(m_bar)^-2 is not a separate unknown: it is absorbed into
    the fitted intercept. The RELATIVE precision is identified even though
    theta is not.
    """
    from model_A.run import rel_err_theta, fibre_Km_max
    p = ParamsA()
    a, b = fit_ab(p)
    sigma, n, Delta = 5.0, 6.0, 0.10 * M_OP
    direct = []
    for Km in np.linspace(1.0, fibre_Km_max(a, b) * 0.99, 12):
        th, _ = fibre(a, b, Km)
        direct.append((sigma / th) * np.sqrt(12.0 / n) * (M_OP / Delta)
                      * g(M_OP, Km) ** -2)
    assert np.std(direct) / np.mean(direct) < 1e-12, "must be constant on the fibre"
    assert direct[0] == pytest.approx(rel_err_theta(sigma, n, Delta), rel=1e-12)


def test_the_band_the_design_needs_exceeds_what_the_instruments_reach():
    """
    Pass 2, items 1 and 4. Both criteria land near 1.5x m_bar, and the
    quasi-exogenous instruments of Section 3.5 reach 70%. This is the answer
    to Table 2's question and it should fail loudly if a parameter moves.
    """
    from model_A.run import (delta_for_precision, delta_for_curvature,
                             CURVATURE_CONST, ASSAY_CV)
    assert CURVATURE_CONST == pytest.approx(np.sqrt(180.0), rel=1e-12)
    d_prec = delta_for_precision(0.20, ASSAY_CV, 6.0)
    d_curv = delta_for_curvature(ASSAY_CV, 6.0)
    assert 1.3 < d_prec / M_OP < 1.6, d_prec
    assert 1.4 < d_curv / M_OP < 1.7, d_curv
    assert min(d_prec, d_curv) > 0.70 * M_OP, "instruments must fall short"


def test_abscissa_noise_biases_the_margin_downward():
    """
    Pass 2, item 3. Classical attenuation inflates theta_hat, which alone
    would make a patient look safer. It does not: ke_hat collapses, so
    theta_crit rises and the margin is UNDERestimated, flipping sign around
    sigma_m/Delta ~ 0.2. The direction is measured, not assumed.
    """
    from model_A.run import eiv_margin_bias
    D = 0.5 * M_OP
    prev = None
    for r in (0.05, 0.10, 0.20, 0.40):
        mu_hat, mu_true, th_ratio = eiv_margin_bias(r * D, D, reps=800, seed=3)
        assert th_ratio > 1.0, f"theta_hat inflated at r={r}"
        assert mu_hat < mu_true, f"but the margin is biased DOWN at r={r}"
        if prev is not None:
            assert mu_hat < prev, "and the bias grows with the noise"
        prev = mu_hat
    mu_hat, mu_true, _ = eiv_margin_bias(0.30 * D, D, reps=800, seed=3)
    assert mu_true > 0 > mu_hat, "the sign flips: safe is estimated as unsafe"


def test_null_anchor_attenuates_load_error_into_theta():
    """
    Pass 2, item 5. The anchored construction with a variance. A relative
    error in the null patient's load enters ke 1:1 and reaches theta
    attenuated by about a third, so the anchor does not relocate the
    conditioning problem; it is the band that remains binding.
    """
    from model_A.run import anchor_propagation
    for rel_s in (0.05, 0.10, 0.20):
        d_Km, d_th = anchor_propagation(rel_s, reps=8000)
        assert d_th < rel_s, f"must attenuate at {rel_s}"
        assert 0.25 < d_th / rel_s < 0.55, f"about a third, got {d_th / rel_s}"
        assert d_Km < d_th, "and Km moves less than theta"
