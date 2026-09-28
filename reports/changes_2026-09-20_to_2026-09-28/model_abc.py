"""
The invariants A, B, C of the local model at the Table 2 constants, shared by
make_numbers.py and fig_abc.py so that the report's check values and its
figures come from one computation.

Everything is taken from the repository's model code rather than restated:

    theta, K_m, k_e   common.models.ParamsA            (3000, 100, 5)
    m_bar             design.run.M_BAR                 (10)
    D, lambda, T      design.run.D_AMP, LAMBDA, DOSE_INTERVAL
    Delta(t)          design.run.delta                 (the biweekly train)
    load fraction     design.run.LOAD_FRAC             (the 5% diary record)
    L                 common.models.clearance_A(m_bar) = theta g(m_bar) + k_e m_bar

The steady-state balance L + b = (theta + Delta) m/(K_m + m) + k_e m solved for
the increment that holds the marker at m is

    Delta(m) = A/m + B + C m,
    A = K_m (L + b),   B = L + b - theta - k_e K_m,   C = -k_e,

which is equation (5) of the paper (Section 3.1). The twins of equation (6) keep
A, B, C fixed while b moves.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from common.models import ParamsA, g, clearance_A          # noqa: E402
from design import run as DR                                # noqa: E402

P = ParamsA()
THETA, KM, KE = P.theta, P.Km, P.ke
M_BAR = DR.M_BAR
D_AMP, LAMBDA, T_DOSE = DR.D_AMP, DR.LAMBDA, DR.DOSE_INTERVAL
LOAD_FRAC = DR.LOAD_FRAC
L = float(clearance_A(M_BAR, P))          # theta g(m_bar) + k_e m_bar
THETA_DISEASE = 0.05 * THETA              # the 5% row of Table 14
M_WATERFALL = 6.0                         # the level the waterfall is drawn at
TWIN_BS = (-200.0, 0.0, 200.0, 400.0)     # the twins drawn in Figure 5b
NONTWIN_DTHETA = 400.0                    # theta moved alone, K_m and b held

#: the periodic steady state of the train, for a patient already on therapy
DELTA_PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * T_DOSE))
DELTA_TROUGH = DELTA_PEAK * np.exp(-LAMBDA * T_DOSE)
#: 60 earlier doses, so draws are read on the filled train (design.run.delta
#: starts the train at t = 0 and sums up to 200 doses)
ON_THERAPY_SHIFT = 60 * T_DOSE


def invariants(theta=THETA, Km=KM, ke=KE, b=0.0, load=L):
    A = Km * (load + b)
    B = load + b - theta - ke * Km
    C = -ke
    return A, B, C


def b_pieces(theta=THETA, Km=KM, ke=KE, b=0.0, load=L):
    """B split into its three physical pieces: load, capacity, backup saving."""
    return {"load": load + b, "theta": -theta, "keKm": -ke * Km}


def delta_of_m(m, theta=THETA, Km=KM, ke=KE, b=0.0, load=L):
    A, B, C = invariants(theta, Km, ke, b, load)
    return A / m + B + C * m


def level_of_delta(delta, theta=THETA, Km=KM, ke=KE, b=0.0, load=L):
    """Invert Delta(m): the positive root of C m^2 + (B - Delta) m + A = 0."""
    A, B, C = invariants(theta, Km, ke, b, load)
    q = B - np.asarray(delta, dtype=float)
    return (-q - np.sqrt(q * q - 4.0 * C * A)) / (2.0 * C)


def untreated_level(theta=THETA, **kw):
    return float(level_of_delta(0.0, theta=theta, **kw))


def twin(b, theta0=THETA, Km0=KM, ke=KE, load=L):
    """The observational twin at load bias b, equation (6): same A, B, C."""
    Km = Km0 * load / (load + b)
    theta = theta0 + b - ke * (Km - Km0)
    return theta, Km


def estimate_under_bias(b, theta=THETA, Km=KM, ke=KE, load=L):
    """What a fit returns for a record at load L + b fitted with the load at L
    (Section 3.1 of the 27 Sept draft): affine in b."""
    return theta - b * (1.0 + Km * ke / load), Km * (1.0 + b / load)


def relaxation_rate(m, theta, Km, ke=KE, delta=0.0):
    """(theta + Delta) g'(m) + k_e, the marker's own relaxation rate."""
    return (theta + delta) * Km / (Km + m) ** 2 + ke


def train_delta(times):
    """The repository's biweekly train, read on a patient already on therapy."""
    return np.asarray(DR.delta(np.asarray(times, float) + ON_THERAPY_SHIFT))


def schedule(step, n=None, span=730.0):
    t = np.arange(0.0, span, step)
    return t if n is None else t[:n]
