import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D
import model_A.run as LOC
from common.excitation import sensitivities
P = D.P
sigma = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
Lb = D.load(0.0)

def filt(levels, deltas, sigma, use_delta):
    Finv = np.linalg.inv(np.diag([0.8**2, 0.8**2, (D.LOAD_FRAC*Lb)**2]))
    for m, dd in zip(np.atleast_1d(levels), np.atleast_1d(deltas)):
        eff = P.theta + (dd if use_delta else 0.0)
        q = np.array([P.theta*m/(P.Km+m),
                      -eff*P.Km*m/(P.Km+m)**2,
                      1.0])
        Finv = Finv + np.outer(q,q)/sigma**2
    return np.linalg.inv(Finv)

T = 730.0
designs = {
 "six-monthly (n=4)": np.arange(180.0, T+1.0, 180.0),
 "one dose cycle (n=7)": np.arange(0.0, D.DOSE_INTERVAL, 2.0),
 "every 2 days x 2y (n=365)": np.arange(0.0, T, 2.0),
}
print(f"sigma_total = {sigma:.4f}   theta={P.theta} Km={P.Km} ke={P.ke}  Lb={Lb:.2f}")
for name, t in designs.items():
    lv = D.quasi_steady_level(t)
    dd = np.asarray(D.delta(t), dtype=float)
    Pa = filt(lv, dd, sigma, False); Pb = filt(lv, dd, sigma, True)
    r=lambda M: M[0,1]/np.sqrt(M[0,0]*M[1,1])
    print(f"{name:>28}: levels {lv.min():.3f}-{lv.max():.3f}  Delta/theta {dd.min()/P.theta:.3f}-{dd.max()/P.theta:.3f}")
    print(f"{'':>28}  sd(ln th) published {np.sqrt(Pa[0,0]):8.4f} (corr {r(Pa):+.4f})   Delta-corrected {np.sqrt(Pb[0,0]):8.4f} (corr {r(Pb):+.4f})")

lv3 = np.geomspace(P.Km/2, P.Km*2, 365)
d3 = np.zeros_like(lv3)
print(f"{'straddling K_m (n=365)':>28}: published {np.sqrt(filt(lv3,d3,sigma,False)[0,0]):.5f}  (load-driven, Delta=0 by construction)")
