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
    because m_bar is not one of its arguments. Contrast theta_crit() above.
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
    t1 = theta_crit(s_peak=500.0, m_bar=40.0, p=p)
    t2 = theta_crit(s_peak=500.0, m_bar=60.0, p=p)
    assert t1 != t2, "theta_crit must depend on the assumed tolerable level"
