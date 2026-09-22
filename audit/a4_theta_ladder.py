"""Candidate generator for the theta_hat ladder of Section 8's boundary paragraph."""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import lsq_linear
import model_C.run as C
from common.models import clearance_C_quasistatic, equilibria_nonmonotone, ParamsA

def misfit_ladder(tops=(1.6, 5.0, 8.0, 12.6, 20.0), n=9, p=None, s=None, PA=None):
    """
    theta_hat the DESIGN's estimator returns from a noiseless parent patient.

    Band: n geometrically spaced levels from the parent's own attractor m_bar
    to top*m_bar -- the same reference level and the same band as the 42%
    figure in the same paragraph. Ordinate: the parent's noiseless quasi-static
    clearance. Model: the design's admissible form theta*g(m,K_m) + k_e*m with
    K_m at the population value and the unknown ADDITIVE CONSTANT in the
    recorded load that fit_ensemble always carries, fitted by least squares
    under the design's own admissibility bounds theta >= 0, k_e >= 0.
    """
    p = C.PC if p is None else p
    s = C.S_BASE if s is None else s
    PA = ParamsA() if PA is None else PA
    m_bar = float(equilibria_nonmonotone(clearance_C_quasistatic, p, s,
                                         np.linspace(1e-3, 3000.0, 120000))[0][0])
    out = []
    for top in tops:
        lv = np.geomspace(m_bar, top * m_bar, int(n))
        ys = clearance_C_quasistatic(lv, p)
        X = np.column_stack([lv / (PA.Km + lv), lv, np.ones_like(lv)])
        r = lsq_linear(X, ys, bounds=([0.0, 0.0, -np.inf], [np.inf] * 3))
        out.append((top, float(r.x[0]), float(r.x[1]), float(r.x[2])))
    return m_bar, out

if __name__ == "__main__":
    tc = C.theta_crit_intrinsic(C.PC, C.S_BASE)
    for n in (9, 200):
        m_bar, rows = misfit_ladder(n=n)
        print(f"n={n}  m_bar={m_bar:.5f}  theta_crit_intrinsic={tc:.4f}")
        for top, th, ke, c0 in rows:
            print(f"   top={top:5} m_bar  theta_hat={th:8.3f}  k_e_hat={ke:.3f}  "
                  f"load const={c0:8.3f}  mu_hat={th-tc:+8.3f}")
        print()
