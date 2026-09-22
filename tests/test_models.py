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


def test_sigma_is_stated_relative_to_the_recorded_load():
    """
    The criteria were once quoted at sigma = 5 "metabolite units", which is a
    load known to 1.55% of baseline and not defensible. sigma belongs on the
    ordinate's scale, and the conclusion is reported as a function of it.
    """
    from model_A.run import (sigma_from_load, baseline_load, rel_err_theta,
                             delta_for_curvature, SIGMA_FRACS)
    assert sigma_from_load(0.05) == pytest.approx(0.05 * baseline_load())
    prev_r = prev_d = 0.0
    for sf in SIGMA_FRACS:
        sg = sigma_from_load(sf)
        r = rel_err_theta(sg, 6.0, 0.10 * M_OP)
        d = delta_for_curvature(sg, 6.0)
        assert r > prev_r and d > prev_d, "both must worsen with sigma"
        assert d > 0.70 * M_OP, "and stay outside what the instruments reach"
        prev_r, prev_d = r, d
    assert rel_err_theta(sigma_from_load(0.05), 6.0, 0.10 * M_OP) == \
        pytest.approx(9.20, rel=1e-2), "920% at a 5% load error"


def test_anchor_factor_inverts_below_a_computable_theta():
    """
    The anchor cost is a ratio scaling as 1/theta, like the load-error
    amplification. It attenuates for a mild patient and amplifies for a severe
    one, and Model C's own theta is on the amplifying side.
    """
    from model_A.run import anchor_factor, anchor_crossover
    p = ParamsA()
    x_km, x_th = anchor_crossover()
    assert x_km == pytest.approx(p.ke * (p.Km + M_OP) ** 2 / p.Km, rel=1e-12)
    assert x_th == pytest.approx(2 * p.ke * (p.Km + M_OP), rel=1e-12)
    assert anchor_factor(theta=p.theta)[1] < 1.0, "attenuates at theta = 3000"
    assert anchor_factor(theta=300.0)[1] > 1.0, "amplifies at Model C's theta"
    for th in (x_km, x_th):
        i = 0 if th == x_km else 1
        assert anchor_factor(theta=th)[i] == pytest.approx(1.0, rel=1e-9)


def test_abscissa_noise_sign_flip_is_an_inadmissible_fit():
    """
    The unconstrained margin appears to change sign. It does so by returning
    ke < 0, which is not a parameter set. Under ke >= 0 the bias is small, and
    the inadmissible fraction is what makes the failure detectable.
    """
    from model_A.run import eiv_fit
    D = 0.5 * M_OP
    free = eiv_fit(0.20 * D, D, reps=500, seed=7, constrain=False)
    cons = eiv_fit(0.20 * D, D, reps=500, seed=7, constrain=True)
    assert free["mu"] < 0 < cons["mu"], "the flip is the unconstrained fit only"
    assert cons["mu"] < cons["mu_true"], "a small downward bias survives"
    assert cons["mu"] > 0.85 * cons["mu_true"], "and it is small"
    assert free["frac_inadmissible"] > 0.9, "and it is detectable"


# =========================================================================
#  PASS 3 -- the five blockers.
#
#  Each block below pins a number that was wrong, missing, or quoted under a
#  convention the model forbids. The failure mode they are written against is
#  the one this repository keeps producing: a quantity that is correct in its
#  own terms and is then compared with a quantity computed on a different
#  footing.
# =========================================================================

# ------------------------------------------------- 0.1  the shape constants
def test_params_table_is_the_model_and_not_a_transcription_of_it():
    """
    The parameter table is generated from the dataclasses, not typed beside
    them. M_i was missing from the paper entirely, which is why Figure 16
    could not be reproduced; a table that can drift would reintroduce that.
    """
    from common.params import (LOCAL, PARENT, M_OP, M_TOL, ILLNESS_FOLD,
                               KE_PRIOR_LOG, LOAD_FRAC)
    import model_A.run as A
    import design.run as DR
    import inspect
    assert (M_OP, M_TOL, ILLNESS_FOLD) == (A.M_OP, A.M_TOL, A.ILLNESS_FOLD)
    # Phase 7 X2: two constants reached the paper's numbers without a row.
    # They are registered now, so they must be checked against the module that
    # actually uses them, exactly as M_OP and its neighbours are.
    assert LOAD_FRAC == DR.LOAD_FRAC
    assert KE_PRIOR_LOG == inspect.signature(
        DR._filter).parameters["ke_prior"].default
    assert all(r.kind for r in LOCAL + PARENT), "a row with no kind stated" 
    by_symbol = {r.symbol: r for r in PARENT}
    p = ParamsC()
    # The Value column carries units since September 2026, so the check is that
    # the dataclass number APPEARS in the rendered value, not that it is the
    # whole of it. The drift this guards against is unchanged: move Mi in the
    # dataclass without regenerating and the substring stops matching.
    assert f"{p.Mi:.0f}" in by_symbol[r"$M_i$"].value
    assert f"{p.Kc:.0f}" in by_symbol[r"$K_c$"].value
    assert f"{p.rho0:.2f}" in by_symbol[r"$\rho_0$"].value
    # and a value that is only a number, with no unit, is now the exception
    dimensionless = {r"$\rho_0$", r"$\kappa$", r"$k$", r"$\mathrm{CV}$",
                     r"$\spk/s_{\mathrm{base}}$", r"$s_{\ln k_e}$",
                     r"$\sigma_L/L$"}
    for r in LOCAL + PARENT:
        if r.symbol in dimensionless or r.value in ("not assigned", "cancels"):
            continue
        assert ("mu$M" in r.value or "d$^{-1}$" in r.value
                or r.value.endswith("~d")), (r.symbol, r.value)
    assert by_symbol[r"$R_0$"].source == "derived"
    # every other constant of the parent model is assumed, and says so
    assert {r.source for r in PARENT} <= {"derived", "assumed", "literature"}
    assert all(r.note for r in LOCAL + PARENT), "a row with no cost stated"


def test_Mi_is_the_value_the_text_arithmetic_already_uses():
    """
    M_i was never stated, but it is recoverable: section 5.1 says absolute
    free-carnitine excretion exceeds normal 'above m ~ 1650', and that
    threshold is (37/5.7 - 1) * M_i. Stating M_i = 300 reproduces it.
    """
    p = ParamsC()
    depletion = 37.0 / 5.7
    assert (depletion - 1.0) * p.Mi == pytest.approx(1650.0, rel=0.01)
    # and the two other figures in the same paragraph
    from model_C.run import PC, S_BASE, margin_at
    roots, st = equilibria_nonmonotone(
        clearance_C_quasistatic, PC, S_BASE, np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    attractor, separatrix = float(roots[st][0]), float(roots[~st][0])
    assert 1.0 + attractor / p.Mi == pytest.approx(1.3, abs=0.03)
    assert (depletion - 1.0) * p.Mi / separatrix == pytest.approx(1.6, abs=0.05)


def test_holding_R0_while_moving_a_shape_constant_breaks_its_own_constraint():
    """
    R0 is defined by a constraint: the alternative route alone must fall short
    of the baseline load, so a null genotype is not maintainable. Moving Kc,
    Mi or rho0 with R0 held moves the peak too, so the swept parameter sets
    are not admissible and the sensitivity reported is mostly the calibration.
    """
    from model_C.run import alt_peak, R0_for_peak, S_BASE
    p = ParamsC()
    assert alt_peak(p)[0] == pytest.approx(500.0, rel=1e-4)
    for field, v in (("Mi", 900.0), ("rho0", 0.35), ("Kc", 30.0)):
        q = replace(p, **{field: v})
        assert abs(alt_peak(q)[0] - 500.0) > 20.0, f"{field} left the peak alone"
        q_cal = replace(q, R0=R0_for_peak(q))
        assert alt_peak(q_cal)[0] == pytest.approx(500.0, rel=1e-4)
        assert alt_peak(q_cal)[0] < S_BASE, "the constraint must still hold"


def test_on_one_footing_M_i_is_cheaper_than_K_c_not_nine_times_dearer():
    """
    The claim under review was a ~9x swing in theta_crit over an 8x range in
    M_i against 1.56x for K_c. It does not survive: 1.56 is an R0-re-derived
    number and the 9 is an R0-held one. On the admissible footing M_i costs
    LESS than K_c, and rho_0 costs almost nothing.
    """
    from model_C.run import shape_sensitivity
    d = shape_sensitivity()
    fac = lambda name: (lambda xs: max(xs) / min(xs))(
        [r[1] for r in d[name] if np.isfinite(r[1])])
    assert fac("K_c") == pytest.approx(1.56, abs=0.03)
    assert fac("M_i") == pytest.approx(1.37, abs=0.05)
    assert fac("rho_0") == pytest.approx(1.02, abs=0.02)
    assert fac("M_i") < fac("K_c"), "M_i is the cheaper of the two"


def test_theta_crit_is_nearly_a_function_of_the_calibration_alone():
    """
    Three unmeasured shape constants, swept over their whole ranges with the
    calibration held, move theta_crit by well under a factor of two. The
    exposure is concentrated in the constraint that fixes R0, which is a
    clinical claim rather than a shape constant.
    """
    from model_C.run import collapse_check
    lo, hi = collapse_check(n=8, seed=1)
    assert hi / lo < 1.6


def test_the_margin_is_far_more_exposed_to_the_shape_than_the_threshold():
    """
    Figure 16 quotes a margin of 954 flat. The separatrix sits where alpha has
    decayed, alpha decays on the scale of M_i, so the margin inherits M_i
    almost linearly where theta_crit does not. This is the number that needed
    the sensitivity analysis, and it is not the one the old remark bounded.
    """
    from model_C.run import shape_sensitivity
    d = shape_sensitivity()
    for name in ("K_c", "M_i"):
        tc = [r[1] for r in d[name] if np.isfinite(r[1])]
        mg = [r[3] for r in d[name] if np.isfinite(r[3])]
        assert max(mg) / min(mg) > 2.0 * (max(tc) / min(tc)), name


# ------------------------------------------------- 0.2  the band, and its floor
def test_the_intercept_variance_has_a_floor_no_band_width_beats():
    """
    Var(a_hat) = sigma^2 [1/n + m_bar^2/S]. Equation (7) keeps only the second
    term. The first one binds in exactly the wide-band regime the note
    recommends buying, and it does not move with the width.
    """
    from model_A.run import rel_err_exact, rel_err_floor, M_OP, sigma_total
    sig = sigma_total(0.05)
    n = 6
    for frac in (0.1, 0.5, 1.0, 1.9, 1.999):
        lo, hi = M_OP * (1 - frac / 2), M_OP * (1 + frac / 2)
        got = rel_err_exact(np.linspace(lo, hi, n), sig)
        assert got >= rel_err_floor(sig, n, "optimal") * (1 - 1e-9), frac
    # and the floor is approached as the band's foot reaches zero
    got = rel_err_exact(np.linspace(1e-9, 2 * M_OP, 400), sig)
    assert got == pytest.approx(rel_err_floor(sig, 400), rel=0.02)


def test_the_two_point_closed_form_is_the_best_use_of_an_excursion():
    """sd(a_hat) = sigma (H+L)/(sqrt(n)(H-L)), checked against brute force."""
    from model_A.run import sd_a_two_point, rel_err_exact, fit_ab, P
    a = fit_ab(P)[0]
    for L, H, n in ((0.5, 10.0, 6), (3.0, 10.0, 6), (2.0, 50.0, 20)):
        brute = min(rel_err_exact(np.array([L] * k + [H] * (n - k)), 1.0) * a
                    for k in range(1, n))
        closed = sd_a_two_point(L, H, n, 1.0)
        assert closed <= brute * (1 + 1e-9)
        assert closed == pytest.approx(brute, rel=0.06)


def test_the_old_design_table_quoted_bands_that_cannot_exist():
    """
    A band is [m_bar - Delta/2, m_bar + Delta/2] and m cannot be negative, so
    Delta/m_bar <= 2. The published table quoted 460%, 920% and 1841%.
    """
    from model_A.run import delta_for_precision, sigma_from_load, M_OP
    impossible = [f for f in (0.05, 0.10, 0.20)
                  if delta_for_precision(0.20, sigma_from_load(f), 6.0) / M_OP > 2.0]
    assert impossible == [0.05, 0.10, 0.20]
    assert delta_for_precision(0.20, sigma_from_load(0.02), 6.0) / M_OP < 2.0


def test_the_curvature_criterion_is_a_top_of_excursion_not_a_band_width():
    """
    The bending of C is at high m, so the band has to run up. Restated that
    way the criterion is a number a clinician can be asked for, and it is
    monotone in the noise.
    """
    from model_A.run import top_for_curvature, sigma_total, M_OP
    tops = [top_for_curvature(sigma_total(f), 6) for f in (0.02, 0.05, 0.10, 0.20)]
    assert all(t is not None for t in tops)
    assert all(b > a for a, b in zip(tops, tops[1:])), tops
    assert 2.0 < tops[0] / M_OP < 3.5 and 4.5 < tops[-1] / M_OP < 6.5


def test_the_tangent_and_least_squares_intercepts_part_company_on_the_convention():
    """
    Footnote 2's numbers, both readings. Delta is the FULL width: at
    Delta/m_bar = 100% the gap is 8.8%. Read as a half-width the same phrase
    means [0, 2 m_bar] and the gap is 35.4%, four times larger.
    """
    from model_A.run import fit_ab, P, M_OP
    a_tan = fit_ab(P)[0]

    def ls(lo, hi):
        m = np.linspace(lo, hi, 20001)
        C = clearance_A(m, P)
        return np.linalg.lstsq(np.vstack([np.ones_like(m), m]).T, C, rcond=None)[0][0]

    assert 100 * (a_tan - ls(0.95 * M_OP, 1.05 * M_OP)) / a_tan == pytest.approx(0.09, abs=0.02)
    assert 100 * (a_tan - ls(0.5 * M_OP, 1.5 * M_OP)) / a_tan == pytest.approx(8.8, abs=0.2)
    assert 100 * (a_tan - ls(1e-9, 2.0 * M_OP)) / a_tan == pytest.approx(35.4, abs=0.3)


# ------------------------------------------------- 0.3  the Gramian
def _matched(amp, T=40.0):
    """A carnitine step, and the therapy schedule reproducing its trajectory."""
    from scipy.interpolate import interp1d
    from common.dynamics import simulate, step
    p = ParamsA()
    L = lambda t: clearance_A(10.0, p)
    tt, mc, Qc = simulate(p.theta, p.Km, p.ke, L, kappa_u=step(amp), T=T)
    mf = interp1d(tt, mc, kind="cubic", fill_value="extrapolate")
    d = lambda t: step(amp)(t) * (p.Km + float(mf(np.clip(t, tt[0], tt[-1]))))
    tt2, mt, Qt = simulate(p.theta, p.Km, p.ke, L, delta=d, T=T)
    return (tt, mc, Qc), (tt2, mt, Qt)


def test_channel_blindness_survives_the_gramian_given_Km():
    """
    Proposition 2 generalises from the determinant to the Gramian exactly as
    far as it is stated. At a matched trajectory the (theta, ke) sub-Gramian
    differs between the channels by a few per cent, always in the therapy's
    favour, so posing the question on carnitine is conservative.
    """
    from common.dynamics import gramian, sd_log_theta
    for amp in (3.0, 30.0, 100.0):
        (tc, mc, Qc), (tt, mt, Qt) = _matched(amp)
        assert np.max(np.abs(mc - mt)) < 0.01 * (mc.max() - mc.min())
        Gc = gramian(tc, Qc, idx=(0, 2))
        Gt = gramian(tt, Qt, idx=(0, 2))
        r = np.linalg.det(Gt) / np.linalg.det(Gc)
        assert 1.0 < r < 1.10, (amp, r)
        assert 1.0 < sd_log_theta(Gc) / sd_log_theta(Gt) < 1.05


def test_the_gap_between_the_channels_is_the_decay_rate_and_is_one_signed():
    """
    J_therapy - J_carnitine = -kappa u g(m) at a matched trajectory. The
    therapy decays the sensitivities more slowly, so it always returns the
    larger Gramian, and the size of the gap is bounded by the input.
    """
    from scipy.interpolate import interp1d
    from common.dynamics import simulate, step
    p = ParamsA()
    amp = 30.0
    L = lambda t: clearance_A(10.0, p)
    tt, mc, _ = simulate(p.theta, p.Km, p.ke, L, kappa_u=step(amp))
    mf = interp1d(tt, mc, kind="cubic", fill_value="extrapolate")
    d = lambda t: step(amp)(t) * (p.Km + float(mf(np.clip(t, tt[0], tt[-1]))))
    gp = p.Km / (p.Km + mc) ** 2
    u = np.array([step(amp)(t) for t in tt])
    dd = np.array([d(t) for t in tt])
    J_c = p.theta * gp + p.ke + u
    J_t = (p.theta + dd) * gp + p.ke
    predicted = -u * (mc / (p.Km + mc))
    assert np.max(np.abs((J_t - J_c) - predicted)) < 1e-3
    assert np.all(J_t <= J_c + 1e-9)


def test_channel_blindness_fails_once_Km_joins_the_parameter_vector():
    """
    THE GATE. The K_m sensitivity is forced by -(theta + Delta(t)) Km dg/dKm,
    and that factor is modulated by the therapeutic input and constant under
    carnitine. So the therapy separates theta from K_m and carnitine does not,
    by orders of magnitude, at the same excursion.
    """
    from common.dynamics import gramian, collinearity, sd_log_theta
    for amp in (10.0, 30.0, 100.0):
        (tc, mc, Qc), (tt, mt, Qt) = _matched(amp)
        Gc, Gt = gramian(tc, Qc), gramian(tt, Qt)
        assert np.linalg.det(Gt) / np.linalg.det(Gc) > 1e4, amp
        assert sd_log_theta(Gc) / sd_log_theta(Gt) > 20.0, amp
        assert collinearity(Gc, 0, 1) > 0.999, "carnitine cannot separate theta from Km"
        assert collinearity(Gt, 0, 1) < 0.99, "the therapy can"


def test_curvature_is_what_restores_the_K_m_direction_to_carnitine():
    """
    Carnitine reaches K_m only through the state. Carrying the patient from
    m/K_m ~ 0.08 to m/K_m ~ 1.5 breaks the collinearity and raises the Gramian
    determinant by many orders of magnitude. This is the excitation condition
    in the form the design has to be stated in.
    """
    from common.dynamics import simulate, step, gramian, collinearity
    from scipy.optimize import brentq
    p = ParamsA()
    out = []
    for m0 in (10.0, 100.0, 200.0):
        L0 = clearance_A(m0, p)
        L = lambda t, L0=L0: L0
        f = lambda A: (lambda m: (m.max() - m.min()) / m.max())(
            simulate(p.theta, p.Km, p.ke, L, kappa_u=step(A), m0=m0)[1]) - 0.5
        A = brentq(f, 1e-4, 1e4, xtol=1e-8)
        tt, m, Q = simulate(p.theta, p.Km, p.ke, L, kappa_u=step(A), m0=m0)
        G = gramian(tt, Q)
        out.append((m0 / p.Km, collinearity(G, 0, 1), np.linalg.det(G)))
    assert out[0][1] > out[1][1] > out[2][1], "collinearity must fall"
    assert out[-1][2] / out[0][2] > 1e10, "and the information must rise"


# ------------------------------------------------- 0.4  which threshold
def test_the_parent_model_has_no_fibre_so_the_cancellation_has_no_premise():
    """
    The K_m^2 cancellation is a statement about a degeneracy. In the parent
    model the alternative channel is curved, so it contributes to both
    coefficients and, given alpha, (a, b) determines (theta, K_m) outright.
    There is no fibre for the margin to be linear along.
    """
    from model_C.run import PC, S_BASE
    p = PC
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, p, S_BASE,
                                       np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    mbar = float(roots[st][0])
    hh = 1e-4 * mbar
    C = lambda m: p.theta * g(m, p.Km) + alt_flux_C(m, p)
    b = (C(mbar + hh) - C(mbar - hh)) / (2 * hh)
    a = C(mbar) - mbar * b
    ap = (alt_flux_C(mbar + hh, p) - alt_flux_C(mbar - hh, p)) / (2 * hh)
    a_e, b_e = a - (alt_flux_C(mbar, p) - mbar * ap), b - ap
    Km_rec = mbar ** 2 * b_e / a_e
    th_rec = a_e * (Km_rec + mbar) ** 2 / mbar ** 2
    assert Km_rec == pytest.approx(p.Km, rel=1e-3)
    assert th_rec == pytest.approx(p.theta, rel=1e-3)


def test_a_equals_theta_g_squared_is_the_local_model_and_not_the_parent():
    """
    Model A's intercept identity holds only when the alternative channel is
    linear through the origin. alpha is not, and its own intercept term
    dominates: 374 against 61 at the nominal constants.
    """
    from model_C.run import PC, S_BASE
    p = PC
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, p, S_BASE,
                                       np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    mbar = float(roots[st][0])
    hh = 1e-4 * mbar
    C = lambda m: p.theta * g(m, p.Km) + alt_flux_C(m, p)
    b = (C(mbar + hh) - C(mbar - hh)) / (2 * hh)
    a = C(mbar) - mbar * b
    a_enz = p.theta * g(mbar, p.Km) ** 2
    ap = (alt_flux_C(mbar + hh, p) - alt_flux_C(mbar - hh, p)) / (2 * hh)
    a_alt = alt_flux_C(mbar, p) - mbar * ap
    assert a == pytest.approx(a_enz + a_alt, rel=1e-4)
    assert a_alt > 5.0 * a_enz, "the alternative channel owns the intercept here"


# ------------------------------------------------- 0.5 / 2.2  the noise budget
def test_the_assay_term_dominates_at_a_well_recorded_load():
    """
    Section 3.7 asserted assay variability is the smaller term. At steady
    state an abscissa error acts like an ordinate error b*sigma_m, and
    b*m_bar is within 8% of the baseline load, so the crossover sits at an
    assay CV of about 1.08 times the relative error in the load.
    """
    from model_A.run import assay_load_crossover, sigma_total, P, M_OP, fit_ab, baseline_load
    k = assay_load_crossover()
    assert k == pytest.approx(1.08, abs=0.03)
    b = fit_ab(P)[1]
    assert b * 0.05 * M_OP > 0.02 * baseline_load(), "assay wins at a 2% load"
    assert b * 0.05 * M_OP < 0.10 * baseline_load(), "load wins at a 10% load"
    assert sigma_total(0.02) > 2.0 * 0.02 * baseline_load()


def test_a_marker_panel_averages_the_assay_term_and_not_the_load_term():
    """
    One latent pool, channel-specific gains. Independent assay noise averages
    as sqrt(J); the load error is common to every channel and does not average
    at all. So a panel pays only while the load term is small, which is the
    answer to whether it is affordable and it is conditional, not general.
    """
    from model_A.run import sigma_total
    gain_good = sigma_total(0.02, channels=1) / sigma_total(0.02, channels=4)
    gain_bad = sigma_total(0.20, channels=1) / sigma_total(0.20, channels=4)
    assert gain_good > 1.5, gain_good
    assert gain_bad < 1.1, gain_bad
    # and it never beats what a perfectly recorded load would give
    assert sigma_total(0.02, channels=64) > 0.02 * 322.0


# =========================================================================
#  PASS 4 -- the restructure. The excitation condition, the reduction
#  hierarchy, the design, and the decision.
# =========================================================================

# ------------------------------------------------- 1.1  the hierarchy
def test_the_two_state_system_reduces_to_the_quasi_steady_flux_symbolically():
    """
    Solving dc/dt = 0 in the two-state system returns alpha exactly, and the
    conjugation rate constant k survives only inside rho_0 = delta_0/k.
    """
    import sympy as sp
    m, c, k, R0, Kc, Mi, d0 = sp.symbols("m c k R_0 K_c M_i delta_0", positive=True)
    conj = k * c * m / (Kc + m)
    c_star = sp.solve(sp.Eq(R0 - conj - d0 * (1 + m / Mi) * c, 0), c)[0]
    alpha = sp.simplify(conj.subs(c, c_star))
    u = m / (Kc + m)
    alpha_paper = R0 * u / (u + (d0 / k) * (1 + m / Mi))
    assert sp.simplify(alpha - alpha_paper) == 0
    # k appears only as the ratio: substituting rho_0 removes it entirely
    rho0 = sp.Symbol("rho_0", positive=True)
    assert k not in sp.simplify(alpha.subs(d0, rho0 * k)).free_symbols


def test_the_two_state_system_settles_on_the_reduced_equilibria():
    """
    Integrating the full two-state system from below the separatrix lands on
    the attractor of the reduced clearance; from above it, it escapes. This is
    the equilibrium reduction verified numerically rather than asserted.
    """
    from scipy.integrate import solve_ivp
    from model_C.run import PC, S_BASE
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE,
                                       np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    attractor, separatrix = float(roots[st][0]), float(roots[~st][0])
    for m0, expect_escape in ((0.5 * attractor, False), (1.2 * separatrix, True)):
        sol = solve_ivp(lambda t, y: rhs_C(t, y, PC, lambda _: S_BASE),
                        (0, 4000.0), [m0, PC.R0], rtol=1e-10, atol=1e-12,
                        method="LSODA")
        end = float(sol.y[0, -1])
        if expect_escape:
            assert end > 10 * separatrix, end
        else:
            assert end == pytest.approx(attractor, rel=1e-6), end


def test_the_ode_equilibrium_matches_the_closed_form_reduction():
    """Pass 2.1: an integrated steady state equals the algebraic one."""
    from scipy.integrate import solve_ivp
    from model_C.run import PC, S_BASE
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE,
                                       np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    attractor = float(roots[st][0])
    assert clearance_C_quasistatic(attractor, PC) == pytest.approx(S_BASE, rel=1e-9)
    sol = solve_ivp(lambda t, y: [S_BASE - clearance_C_quasistatic(y[0], PC)],
                    (0, 2000.0), [0.3 * attractor], rtol=1e-11, atol=1e-13,
                    method="LSODA")
    assert float(sol.y[0, -1]) == pytest.approx(attractor, rel=1e-8)


def test_halving_the_step_does_not_move_the_quasi_steady_level():
    """
    Pass 2.1: the integrator's step must not be carrying the answer. Euler is
    for the write-up; computation uses an adaptive stiff-capable solver, and
    this is the check that says so.
    """
    from scipy.integrate import solve_ivp
    import design.run as D
    Lb = D.load(0.0)
    rhs = lambda t, y: [Lb - (D.P.theta + float(D.delta(t))) * y[0]
                        / (D.P.Km + y[0]) - D.P.ke * y[0]]
    ends = []
    for ms in (0.05, 0.025):
        sol = solve_ivp(rhs, (0.0, 10.0), [D.quasi_steady_level(0.0)],
                        rtol=1e-9, atol=1e-12, method="Radau", max_step=ms)
        ends.append(float(sol.y[0, -1]))
    assert abs(ends[0] - ends[1]) / ends[0] < 1e-6, ends
    assert ends[0] == pytest.approx(float(D.quasi_steady_level(10.0)), rel=0.01)


# ------------------------------------------------- 1.4  the fibre span
def test_the_fibre_span_is_the_recorded_load_over_the_fitted_intercept_squared():
    """
    (1 + b m_bar/a)^2 = (C(m_bar)/a)^2, both observable. This is what turns
    the illustrative 169 into a per-patient diagnostic.
    """
    from model_A.run import (fit_ab, fibre, fibre_Km_max, fibre_span,
                             baseline_load, P, M_OP)
    a, b = fit_ab(P)
    Kmax = fibre_Km_max(a, b)
    brute = fibre(a, b, Kmax)[0] / fibre(a, b, 1e-12)[0]
    assert fibre_span(a, b) == pytest.approx(brute, rel=1e-9)
    assert fibre_span(a, b) == pytest.approx((baseline_load() / a) ** 2, rel=1e-12)
    assert fibre_span(a, b) == pytest.approx(169.4, rel=1e-3)


# ------------------------------------------------- 1.6  the excitation
def test_both_asymptotes_are_degenerate_and_they_fail_differently():
    """
    Proposition: both ends of the saturation curve are degenerate, and the
    near-null direction is different at each. The degeneracy is asymptotic, so
    the statement is about conditioning and about WHICH direction goes flat,
    not about exact rank at a finite level.
    """
    from common.excitation import sensitivities
    th, Km, ke = 3000.0, 100.0, 5.0

    def svd_of(levels):
        Q = sensitivities(np.asarray(levels, float), th, Km, ke)
        Q = Q / np.linalg.norm(Q, axis=0)
        _, sv, Vt = np.linalg.svd(Q, full_matrices=False)
        return sv, Vt

    sv_lo, V_lo = svd_of(np.geomspace(1e-4, 1e-3, 12) * Km)
    sv_mid, _ = svd_of(np.geomspace(0.5, 2.0, 12) * Km)
    sv_hi, V_hi = svd_of(np.geomspace(1e4, 1e5, 12) * Km)
    assert sv_lo[0] / sv_lo[-1] > 1e5, sv_lo
    assert sv_hi[0] / sv_hi[-1] > 1e5, sv_hi
    assert sv_mid[0] / sv_mid[-1] < 1e3, sv_mid

    # LOW end: TWO directions go flat, and neither involves the load bias --
    # theta, K_m and k_e all collapse onto m and only their lumped ratio
    # survives. The clearance stays constrained; its decomposition does not.
    assert sv_lo[0] / sv_lo[-2] > 1e3, sv_lo
    assert abs(V_lo[-1][3]) < 0.1 and abs(V_lo[-2][3]) < 0.1

    # HIGH end: ONE direction goes flat and it is theta against the load bias.
    # The target parameter is what is lost, which is the worse failure.
    assert sv_hi[0] / sv_hi[-2] < 1e3, sv_hi
    v = V_hi[-1]
    assert abs(v[0]) > 0.3 and abs(v[3]) > 0.3, v
    assert np.sign(v[0]) != np.sign(v[3]), v
    assert abs(v[1]) < 0.1, v


def test_the_excitation_is_maximised_straddling_K_m_geometrically():
    """
    E_max(r) = (sqrt(r)-1)/(sqrt(r)+1), attained at sqrt(m_lo m_hi) = K_m.
    Checked against a brute-force search over placement.
    """
    from common.excitation import excitation, excitation_max, optimal_band, fold_range_for
    Km = 100.0
    for r in (1.5, 2.0, 10.0, 100.0):
        lo_star, hi_star = optimal_band(r, Km)
        assert np.sqrt(lo_star * hi_star) == pytest.approx(Km, rel=1e-12)
        best = max(excitation(c / np.sqrt(r), c * np.sqrt(r), Km)
                   for c in np.geomspace(1e-3 * Km, 1e3 * Km, 4001))
        assert excitation_max(r) == pytest.approx(best, rel=1e-3)
        assert excitation(lo_star, hi_star, Km) == pytest.approx(excitation_max(r), rel=1e-12)
        assert fold_range_for(excitation_max(r)) == pytest.approx(r, rel=1e-9)


def test_placement_beats_width_by_orders_of_magnitude():
    """
    The band the illustrative patient is held in throws away a factor of three
    in E by placement alone, and far more than that in sd(ln theta).
    """
    from common.excitation import excitation, excitation_max, sd_log_theta, band_levels
    th, Km, ke, m_bar = 3000.0, 100.0, 5.0, 10.0
    here = excitation(m_bar / np.sqrt(2), m_bar * np.sqrt(2), Km)
    assert excitation_max(2.0) / here == pytest.approx(3.0, abs=0.15)
    best = min(sd_log_theta(band_levels(c, 2.0, 9), th, Km, ke)
               for c in np.geomspace(0.02, 100.0, 400) * Km)
    assert sd_log_theta(band_levels(m_bar, 2.0, 9), th, Km, ke) / best > 1000


def test_the_optimum_is_interior_only_because_the_load_is_unrecorded():
    """
    With the load known exactly, sd(ln theta) falls monotonically into
    saturation and there is no interior optimum. It is the unrecorded load
    that creates the upper obstruction.
    """
    from common.excitation import sd_log_theta, band_levels
    th, Km, ke = 3000.0, 100.0, 5.0
    cen = np.geomspace(0.05, 50.0, 40) * Km
    known = [sd_log_theta(band_levels(c, 2.0, 9), th, Km, ke, load_bias=False)
             for c in cen]
    biased = [sd_log_theta(band_levels(c, 2.0, 9), th, Km, ke, load_bias=True)
              for c in cen]
    assert all(b <= a * (1 + 1e-9) for a, b in zip(known, known[1:])), "must be monotone"
    i = int(np.argmin(biased))
    assert 0 < i < len(cen) - 1, "must be interior"


# ------------------------------------------------- 1.8 / 1.9  the design
def test_the_Km_sensitivity_carries_the_therapeutic_increment():
    """
    The regression guard that would have caught the 1.8 error.

    The ordinate is (theta + Delta) g(m) + k_e m + b, so its derivative in
    ln K_m carries (theta + Delta). Finite-difference the clearance itself and
    check which analytic column it agrees with. The uncorrected column is wrong
    by tens of per cent wherever the dose has just landed, which is most of a
    dose cycle.
    """
    import design.run as D
    th, Km, ke = D.P.theta, D.P.Km, D.P.ke
    h = 1e-6
    worst_pub = 0.0
    for t in (0.0, 1.0, 3.0, 7.0, 13.0, 40.0, 729.0):
        dl = float(D.delta(t))
        m = float(D.quasi_steady_level(t))
        cl = lambda k: (th + dl) * m / (k + m) + ke * m
        fd = (cl(Km * np.exp(h)) - cl(Km * np.exp(-h))) / (2 * h)
        corrected = -(th + dl) * Km * m / (Km + m) ** 2
        published = -th * Km * m / (Km + m) ** 2
        assert abs(fd - corrected) < 1e-6 * abs(fd), (t, fd, corrected)
        worst_pub = max(worst_pub, abs(fd - published) / abs(fd))
    assert worst_pub > 0.5, "the uncorrected column must be visibly wrong, else this guards nothing"


def test_the_dense_design_collapses_once_the_increment_is_carried():
    """
    THE 1.8 GATE, as it stands after the phase 2C audit.

    With Delta carried, two-daily sampling reaches the 20% target and so does a
    single seven-sample dose cycle: the excitation is bought by the dose train,
    not by the marker band. Placement still helps, but by single figures rather
    than by the factor of 122 the uncorrected column reported.
    """
    import design.run as D
    import model_A.run as LOC
    sigma = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
    t_sp = np.arange(180.0, 731.0, 180.0)
    t_dn = np.arange(0.0, 730.0, 2.0)
    t_cy = np.arange(0.0, D.DOSE_INTERVAL, 2.0)
    sd = lambda t, lv=None: float(np.sqrt(D._filter(
        D.quasi_steady_level(t) if lv is None else lv, sigma,
        np.asarray(D.delta(t), dtype=float))[-1][0, 0]))
    s_sparse, s_dense, s_cycle = sd(t_sp), sd(t_dn), sd(t_cy)
    placed = np.geomspace(D.P.Km / 2.0, D.P.Km * 2.0, len(t_dn))
    s_placed = sd(t_dn, placed)

    assert s_dense < 0.20, "the gate passes on the dense design"
    assert s_cycle > 0.20, (
        "but NOT on one dose cycle of seven samples. That result was reported in "
        "the phase 2A audit and is retracted: it held only with k_e known, which "
        "is not the footing this paper claims. See REVIEW_LOG.md.")
    assert s_sparse > 0.20, "nor on six-monthly visits"
    assert s_sparse / s_dense > 2.0, (
        "and on this footing density DOES buy something, against the published claim")
    assert 3.0 < s_dense / s_placed < 20.0, (
        "placement still helps, but by single figures, not by 122")

    # and the uncorrected column is what produced the opposite conclusion
    stale_dense = float(np.sqrt(D._filter(D.quasi_steady_level(t_dn), sigma,
                                          ke_prior=None)[-1][0, 0]))
    assert stale_dense > 0.20 and stale_dense / s_dense > 3.0


def test_the_corrected_design_matrix_is_rank_deficient_at_constant_load():
    """
    The price of the correction, proved so it cannot be forgotten.

    At a CONSTANT recorded load the load-bias column lies exactly in the span of
    the other two, at c_theta = (1 + K_m k_e/L)/theta and c_Km = -1/L. The
    covariance in that direction is then filled in by the load-record prior
    rather than by the data. Any real variation in recorded intake breaks it.
    """
    import design.run as D
    th, Km, ke = D.P.theta, D.P.Km, D.P.ke
    L = D.load(0.0)
    t = np.arange(0.0, 730.0, 2.0)
    lv = D.quasi_steady_level(t)
    dd = np.asarray(D.delta(t), dtype=float)
    Q = np.stack([th * lv / (Km + lv),
                  -(th + dd) * Km * lv / (Km + lv) ** 2,
                  np.ones_like(lv)], axis=1)
    c = np.array([(1.0 + Km * ke / L) / th, -1.0 / L, -1.0])
    assert np.linalg.norm(Q @ c) < 1e-9 * np.linalg.norm(Q), "the identity must be exact"
    assert np.linalg.matrix_rank(Q, tol=1e-10) == 2

    # a varying recorded load breaks it
    rng = np.random.default_rng(0)
    Ls = L * (1.0 + 0.05 * rng.standard_normal(t.size))
    B = ke * Km + th + dd - Ls
    lv2 = (-B + np.sqrt(B ** 2 + 4.0 * ke * Ls * Km)) / (2.0 * ke)
    Q2 = np.stack([th * lv2 / (Km + lv2),
                   -(th + dd) * Km * lv2 / (Km + lv2) ** 2,
                   np.ones_like(lv2)], axis=1)
    assert np.linalg.matrix_rank(Q2, tol=1e-10) == 3


def test_the_constant_load_deficiency_is_an_exact_fibre_not_a_curvature():
    """
    What the rank deficiency above actually is, asserted so that the
    first-order reading cannot come back.

    Solving the balance for the increment that produces a level m gives, with
    no remainder,

        Delta(m) = A/m + B + C m,   A = K_m (L + b),
                                    B = L + b - theta - k_e K_m,
                                    C = -k_e,

    and {1/m, 1, m} are independent as functions of m, so (A, B, C) is
    everything a constant-load record of (m_k, Delta_k) pairs contains. Four
    unknowns against three invariants leave a ONE-PARAMETER family of exact
    observational twins,

        K_m(b) = K_0 L/(L + b),   theta(b) = theta_0 + b + k_e (K_0 - K_m(b)),

    with k_e unchanged. The family is EXACT, not second order: the noiseless
    misfit stays at machine zero along it at any displacement. The straight
    line tangent to it at b = 0 is the null vector of the test above, and
    moving along THAT line gives a quartic misfit, which is the curvature of
    the fibre relative to its own tangent rather than information the record
    holds. audit/f5_fibre.py, audit/PHASE4_PI.md, and REVIEW_LOG.md R14.
    """
    import design.run as D
    th0, K0, ke = D.P.theta, D.P.Km, D.P.ke
    L = D.load(0.0)
    t = np.arange(0.0, 730.0, 30.4)
    dd = np.asarray(D.delta(t), dtype=float)

    def fibre(b):
        Km = K0 * L / (L + b)
        return th0 + b + ke * (K0 - Km), Km

    def levels(th, Km, b):
        B = ke * Km + th + dd - (L + b)
        return (-B + np.sqrt(B ** 2 + 4.0 * ke * (L + b) * Km)) / (2.0 * ke)

    def misfit(th, Km, b, m):
        return (th + dd) * m / (Km + m) + ke * m - (L + b)

    m0 = levels(th0, K0, 0.0)
    scale = np.linalg.norm(np.full_like(m0, L))

    # along the exact fibre the misfit is machine zero at every displacement,
    # and the twin generates the same levels, which is the observational claim
    for b in (1.0, 10.0, 30.0, 100.0, 200.0):
        th, Km = fibre(b)
        assert np.linalg.norm(misfit(th, Km, b, m0)) < 1e-9 * scale, b
        assert np.allclose(levels(th, Km, b), m0, rtol=1e-11, atol=1e-11), b

    # along the straight tangent the misfit is QUARTIC in the displacement,
    # which is the quantity the first-order reading measured
    n = np.array([(L + ke * K0) / L, -K0 / L, 1.0])
    ss = {}
    for eps in (1.0, 10.0):
        ss[eps] = float(np.sum(misfit(th0 + eps * n[0], K0 + eps * n[1],
                                      eps * n[2], m0) ** 2))
    assert ss[1.0] > 0.0, "the line must NOT be flat, else this guards nothing"
    power = np.log(ss[10.0] / ss[1.0]) / np.log(10.0)
    assert 3.6 < power < 4.4, (power, ss)

    # and that line is the tangent to the fibre at b = 0, which is the null
    # vector of prop:rank
    d = 1e-6
    thp, Kmp = fibre(d)
    thm, Kmm = fibre(-d)
    tang = np.array([(thp - thm) / (2 * d), (Kmp - Kmm) / (2 * d), 1.0])
    assert np.allclose(tang, n, rtol=1e-6), (tang, n)

    # a varying recorded load closes the fibre: the same twin no longer
    # reproduces the record
    rng = np.random.default_rng(0)
    Ls = L * (1.0 + 0.05 * rng.standard_normal(t.size))
    B = ke * K0 + th0 + dd - Ls
    mv = (-B + np.sqrt(B ** 2 + 4.0 * ke * Ls * K0)) / (2.0 * ke)
    th, Km = fibre(30.0)
    r = (th + dd) * mv / (Km + mv) + ke * mv - (Ls + 30.0)
    assert np.linalg.norm(r) > 1e-3 * np.linalg.norm(Ls), np.linalg.norm(r)


def test_the_marker_is_quasi_steady_except_just_after_a_dose():
    """
    The separation of time scales that makes the information filter exact,
    and the sampling exclusion window that goes with it.
    """
    import design.run as D
    settle, worst = D.quasi_steady_check()
    assert 0.0 < settle * 24 < 6.0, settle * 24
    assert worst < 0.06, worst


def test_the_two_decision_errors_are_asymmetric_and_the_costly_one_is_larger():
    """
    1.9. At equal distance from the boundary the error of failing to treat is
    far larger than the error of treating unnecessarily, because the abscissa
    noise biases theta_hat upward and mu_hat inherits it.
    """
    import design.run as D
    from common.models import clearance_A
    base = clearance_A(D.M_BAR, D.P)
    lv = np.geomspace(D.M_BAR, 5 * D.M_BAR, 9)
    safe = D.sign_error(lv, D.P.theta, D.P.Km, D.P.ke, reps=200, seed=7,
                        s_peak=2.98 * base)
    unsafe = D.sign_error(lv, D.P.theta, D.P.Km, D.P.ke, reps=200, seed=7,
                          s_peak=3.58 * base)
    assert safe["mu"] > 0 > unsafe["mu"]
    assert abs(safe["mu"] + unsafe["mu"]) < 0.05 * abs(safe["mu"]), \
        "the two points must be near-matched in distance from the boundary"
    assert unsafe["p"] > 3 * safe["p"], (safe["p"], unsafe["p"])


def test_a_gaussian_tail_on_mu_misrepresents_the_decision():
    """
    The delta method is nearly flat at about 40% while the simulated error
    runs from under 1% to over 50%. Keeping both is what makes the gap
    visible, and the gap is the margin theorem.
    """
    import design.run as D
    import model_A.run as LOC
    from common.models import clearance_A
    base = clearance_A(D.M_BAR, D.P)
    sigma = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
    lv = np.geomspace(D.M_BAR, 5 * D.M_BAR, 9)
    sim = D.sign_error(lv, D.P.theta, D.P.Km, D.P.ke, reps=200, seed=8,
                       s_peak=1.5 * base)["p"]
    lin = D.sign_error_linear(lv, D.P.theta, D.P.Km, D.P.ke, sigma,
                              s_peak=1.5 * base)["p"]
    assert sim < 0.05 and lin > 0.30, (sim, lin)


def test_the_cohort_size_diverges_as_the_sign_error_approaches_a_half():
    """A record that cannot decide the sign cannot be rescued by more patients."""
    import design.run as D
    assert D.cohort_size(0.02) < D.cohort_size(0.2) < D.cohort_size(0.45)
    assert not np.isfinite(D.cohort_size(0.5))
    assert not np.isfinite(D.cohort_size(0.6))


def test_only_a_decompensation_approaches_the_required_excitation():
    """
    The excursion table's conclusion, asserted rather than left to a caption.
    """
    import design.run as D
    from common.excitation import excitation, excitation_max
    rows = {name: (lo, hi) for name, lo, hi, _ in D.excursion_sources()}
    Es = {n: excitation(lo, hi, D.P.Km) for n, (lo, hi) in rows.items()}
    worst = max(v for k, v in Es.items() if "decompensation" not in k)
    best = Es[[k for k in Es if "decompensation" in k][0]]
    assert best > 4 * worst, (best, worst)
    # and the routine sources are in the wrong PLACE, not merely too narrow
    for name, (lo, hi) in rows.items():
        if "decompensation" in name:
            continue
        frac = excitation(lo, hi, D.P.Km) / excitation_max(hi / lo)
        assert frac < 0.5, (name, frac)


# ------------------------------- the reduction away from the illustrative centre
#
#  audit/e10_qss_corners.py runs the full sweep across the corners of the
#  cohort's sampled region.  These three lock in what it found, because the
#  finding is easy to lose: the equilibrium agreement above looks like a check
#  on the quasi-steady assumption and is not one.

def _c_star(m, p, k):
    """The slow manifold: the carnitine level the reduction assumes."""
    u = m / (p.Kc + m)
    return p.R0 / (k * (u + p.rho0 * (1.0 + m / p.Mi)))


@pytest.mark.parametrize("ft,fk,fl", [(1.0, 1.0, 1.0), (1 / 3, 1 / 3, 0.5),
                                      (3.0, 3.0, 0.5), (3.0, 1 / 3, 2.0),
                                      (1.0, 3.0, 1.0)])
def test_the_reduction_is_exact_at_equilibrium_across_the_cohort_corners(ft, fk, fl):
    """
    The cohort of audit/e8_cohort_unblocked.py spans theta and K_m nine-fold and
    the load four-fold, and the reduction had only ever been checked at the
    centre.  At the corners it is not merely within tolerance, it is exact:
    alpha is defined by solving dc/dt = 0, and dc/dt = 0 holds at any fixed
    point, so the parent's equilibrium condition and s = theta g + alpha are one
    equation.  Asserted here far tighter than the 1e-6 used elsewhere, because
    1e-6 would pass on an answer that had genuinely drifted.
    """
    from scipy.integrate import solve_ivp
    from model_C.run import PC, S_BASE
    p, s = replace(PC, theta=PC.theta * ft, Km=PC.Km * fk), S_BASE * fl
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, p, s,
                                       np.linspace(1e-3, 3e4, 120000))
    st = np.asarray(st, dtype=bool)
    if not st.any():
        pytest.skip("no equilibrium at this corner: the patient decompensates")
    attr = float(roots[st][0])
    sol = solve_ivp(lambda t, y: rhs_C(t, y, p, lambda _: s),
                    (0, 20000.0), [0.5 * attr, _c_star(0.5 * attr, p, 1.0)],
                    rtol=1e-10, atol=1e-12, method="LSODA")
    assert float(sol.y[0, -1]) == pytest.approx(attr, rel=1e-12)


@pytest.mark.parametrize("k", [0.1, 1.0, 100.0])
def test_equilibrium_agreement_does_not_depend_on_the_timescale_separation(k):
    """
    Why the equilibrium checks cannot report degradation, stated as a test.

    k is the conjugation constant.  It cancels out of alpha, so the reduced
    curve is the same for every k -- but it alone sets how fast the carnitine
    state is, so the reduction's PREMISE is false at small k and true at large
    k.  The equilibrium agreement is unmoved either way, which is exactly why it
    is not evidence that carnitine is fast.
    """
    from scipy.integrate import solve_ivp
    from model_C.run import PC, S_BASE
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE,
                                       np.linspace(1e-3, 3e4, 120000))
    attr = float(roots[np.asarray(st, dtype=bool)][0])
    sol = solve_ivp(lambda t, y: rhs_C(t, y, PC, lambda _: S_BASE, k),
                    (0, 40000.0 / k), [0.5 * attr, _c_star(0.5 * attr, PC, k)],
                    rtol=1e-10, atol=1e-12, method="LSODA")
    assert float(sol.y[0, -1]) == pytest.approx(attr, rel=1e-12)


def test_the_reduction_is_not_exact_off_equilibrium_and_improves_with_k():
    """
    Where the quasi-steady assumption actually costs something.  Along the
    transient the reduced and parent trajectories differ by a large fraction at
    the k = 1.0 rhs_C defaults to, and the error falls roughly in proportion to
    1/k.  If someone later makes this agree at k = 1 they have changed the
    model, not fixed a tolerance, and this should fail.
    """
    from scipy.integrate import solve_ivp
    from model_C.run import PC, S_BASE
    roots, st = equilibria_nonmonotone(clearance_C_quasistatic, PC, S_BASE,
                                       np.linspace(1e-3, 3e4, 120000))
    attr = float(roots[np.asarray(st, dtype=bool)][0])
    m0, T = 0.5 * attr, 30.0
    tt = np.linspace(0.0, T, 2001)
    red = solve_ivp(lambda t, y: [S_BASE - clearance_C_quasistatic(y[0], PC)],
                    (0, T), [m0], dense_output=True, rtol=1e-10, atol=1e-12,
                    method="LSODA").sol(tt)[0]
    errs = []
    for k in (1.0, 10.0, 100.0):
        par = solve_ivp(lambda t, y: rhs_C(t, y, PC, lambda _: S_BASE, k),
                        (0, T), [m0, _c_star(m0, PC, k)], dense_output=True,
                        rtol=1e-10, atol=1e-12, method="LSODA").sol(tt)[0]
        errs.append(float(np.max(np.abs(par - red) / red)))
    assert errs[0] > 0.05, errs          # not a rounding-level disagreement
    assert errs[0] > errs[1] > errs[2]   # and it is the timescale ratio that drives it
    assert errs[2] < 0.02, errs          # separation bought by k restores it


def test_the_params_table_asks_for_no_font_shape_latin_modern_lacks():
    """
    CI fails the paper build if paper/main.log contains the string "undefined",
    and that grep catches font warnings as well as references. Latin Modern has
    no bold small-caps, so a \\textsc nested inside a \\textbf emits

        LaTeX Font Warning: Font shape `T1/lmr/bx/sc' undefined

    which builds a correct PDF and fails the build anyway. The retraction note
    on the row for k did exactly that, having set the withdrawn source class in
    small caps inside the bold lead-in. Checked on the rendered row rather than
    on the source string, because the note is assembled from several literals
    and only the rendered form is what LaTeX sees.
    """
    import re
    from common.params import LOCAL, PARENT, _tex

    for line in _tex(LOCAL + PARENT):
        for bold in re.findall(r"\\textbf\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", line):
            assert "\\textsc" not in bold, (
                f"bold small-caps asks for a font shape Latin Modern lacks: {bold}")


# ------------------------------------------------------- physical admissibility
#
# A fit can balance the data and describe no patient. These assert the
# constraints that say so, one test per constraint, so that a change to a
# constant that makes the model unphysical fails here rather than producing a
# figure. audit/h1_admissibility.py reports the same checks with their numbers.


def _design():
    """The design model's constants, in one place so the tests cannot drift."""
    import design.run as D
    peak = D.D_AMP / (1.0 - np.exp(-D.LAMBDA * D.DOSE_INTERVAL))
    return D, D.P, float(D.load(0.0)), peak


def test_every_rate_in_every_parameter_set_is_positive():
    """
    Rates are capacities and clearances. A negative one is not a worse fit, it
    is a different physics, and nothing downstream checks for it.
    """
    for p in (ParamsA(), ParamsB(), ParamsC()):
        for name in ("theta", "Km"):
            assert getattr(p, name) > 0.0, (type(p).__name__, name)
    for p in (ParamsA(), ParamsB()):
        assert p.ke > 0.0, type(p).__name__
    for name in ("Kc", "Mi", "R0"):
        assert getattr(ParamsC(), name) > 0.0, name

    D, P, Lb, _ = _design()
    for name, v in (("D_AMP", D.D_AMP), ("LAMBDA", D.LAMBDA),
                    ("DOSE_INTERVAL", D.DOSE_INTERVAL), ("load", Lb),
                    ("M_BAR", D.M_BAR)):
        assert v > 0.0, (name, v)


def test_the_increment_is_non_negative_and_never_exceeds_its_steady_peak():
    """
    Delta is a sum of decaying doses, so it cannot be negative and cannot
    exceed the geometric limit D/(1 - e^{-lambda T}). A schedule that broke
    either would be adding capacity the therapy never delivers.
    """
    D, _, _, peak = _design()
    t = np.linspace(0.0, 60 * D.DOSE_INTERVAL, 20001)
    d = np.asarray(D.delta(t), dtype=float)
    assert d.min() >= 0.0, d.min()
    assert d.max() <= peak + 1e-9, (d.max(), peak)
    # and it does reach the limit, so the bound is tight rather than vacuous
    assert d.max() > 0.999 * peak, (d.max(), peak)


def test_the_balance_has_exactly_one_positive_root_at_every_increment():
    """
    k_e m^2 + (k_e K_m + theta + Delta - L) m - L K_m = 0. The product of the
    roots is -L K_m / k_e, which is negative for any positive load, so the
    roots straddle zero whatever Delta does and the quasi-steady level is
    well defined and unique. quasi_steady_level takes the + branch; this
    asserts that branch is the positive one.
    """
    D, P, Lb, peak = _design()
    assert -Lb * P.Km / P.ke < 0.0
    for d in (0.0, 1.0, peak, 10.0 * peak):
        B = P.ke * P.Km + P.theta + d - Lb
        disc = B ** 2 + 4.0 * P.ke * Lb * P.Km
        assert disc > 0.0, d
        plus = (-B + np.sqrt(disc)) / (2.0 * P.ke)
        minus = (-B - np.sqrt(disc)) / (2.0 * P.ke)
        assert plus > 0.0 > minus, (d, plus, minus)
        res = P.ke * plus ** 2 + B * plus - Lb * P.Km
        assert abs(res) < 1e-6 * Lb * P.Km, (d, res)
    t = np.linspace(0.0, 4 * D.DOSE_INTERVAL, 2001)
    assert np.asarray(D.quasi_steady_level(t), dtype=float).min() > 0.0


def test_the_total_load_stays_positive_across_the_load_prior():
    """
    L + b is an intake and cannot be negative. The load record's own 95%
    interval must sit well inside that, or the prior is placing mass on
    patients who eat a negative amount.
    """
    D, _, Lb, _ = _design()
    half = 1.96 * D.LOAD_FRAC * Lb
    assert Lb - half > 0.0, (Lb, half)
    # and with a lot of room: positivity is nowhere near binding
    assert Lb / half > 5.0, Lb / half


@pytest.mark.parametrize("frac", [0.05, 0.10])
def test_the_capacity_ceiling_holds_at_the_adopted_reference_and_the_check(frac):
    """
    theta + Delta_peak <= normal capacity.

    `normal` is not a constant of this repository. It is implied by a stated
    fraction. The author settled the reference at 5% in September 2026
    (audit/h1_admissibility.py REFERENCE_FRAC), and the 10% reading is kept
    here as a check that the verdict does not turn on the choice. It cannot:
    the used fraction is frac (1 + Delta_peak/theta), so theta cancels and the
    answer is the same in the local and parent blocks. BOTH MUST STILL HOLD,
    because a reference that changed the verdict would mean the ceiling was
    being carried by the choice rather than by the model.
    """
    _, P, _, peak = _design()
    normal = P.theta / frac
    assert P.theta + peak <= normal, (frac, P.theta + peak, normal)
    used = (P.theta + peak) / normal
    assert used == pytest.approx(frac * (1.0 + peak / P.theta), rel=1e-12)


def test_admissibility_bounds_the_constant_load_fibre_on_both_sides():
    """
    The exact fibre of the constant-load case is unidentifiable from the data,
    so admissibility is the only thing that bounds it. It does, on both sides:
    positivity of theta from below, the capacity ceiling from above.

    The bound is real and it is loose. The interval it leaves is orders of
    magnitude wider than the load record's own 95% interval, which is why the
    record and not physics carries the identification at a constant load.
    """
    from scipy.optimize import brentq
    D, P, Lb, peak = _design()
    K0, TH0, KE = P.Km, P.theta, P.ke
    Km_b = lambda b: K0 * Lb / (Lb + b)
    th_b = lambda b: TH0 + b + KE * (K0 - Km_b(b))

    lo = brentq(th_b, -Lb + 1e-9, 0.0, xtol=1e-12)
    assert -Lb < lo < 0.0, lo
    assert Km_b(lo) > 0.0 and th_b(lo) == pytest.approx(0.0, abs=1e-6)

    half = 1.96 * D.LOAD_FRAC * Lb
    for frac in (0.05, 0.10):
        hi = brentq(lambda b: th_b(b) + peak - TH0 / frac, 0.0, 1e7, xtol=1e-9)
        assert hi > 0.0, (frac, hi)
        assert th_b(hi) + peak == pytest.approx(TH0 / frac, rel=1e-9)
        assert hi - lo > 100.0 * (2.0 * half), (frac, hi - lo, 2.0 * half)

    # the marker is bit-identical along the fibre, so every admissible b
    # describes the same record: that is why only physics can bound it
    for b in (lo / 2.0, 0.0, 10.0, 1000.0):
        km, th = Km_b(b), th_b(b)
        B0 = KE * K0 + TH0 - Lb
        m0 = (-B0 + np.sqrt(B0 ** 2 + 4.0 * KE * Lb * K0)) / (2.0 * KE)
        B = KE * km + th - (Lb + b)
        m = (-B + np.sqrt(B ** 2 + 4.0 * KE * (Lb + b) * km)) / (2.0 * KE)
        assert m == pytest.approx(m0, rel=1e-12), (b, m, m0)
