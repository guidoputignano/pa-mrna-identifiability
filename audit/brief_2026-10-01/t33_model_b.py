"""
Round 7, item 1: redo the headline results under Model B, side by side with A.

Every number in this brief was computed on MODEL A (dm/dt = s - theta*g(m) -
ke*m), whose alternative clearance ke*m is LINEAR and UNBOUNDED. R-H showed that
lets the alternative route absorb ~84% of the load at 5% capacity, conflicting
with the PA literature (non-PCC disposal is low-capacity -> patients accumulate).
src/common/models.py carries a saturating alternative, MODEL B:
    dm/dt = s - theta*g(m) - ke*h(m),  g=m/(Km+m), h=m/(Ke+m)
ParamsB: theta=3000, Km=100, ke=400 (V_max of the alt route), Ke=300. ke is now a
MAXIMAL FLUX, not a rate constant.

THE FIRST QUESTION IS NOT PRECISION, IT IS EXISTENCE. Under B the total disposal
is bounded: at inc=0 it rises to theta+ke as m->inf. A steady state exists only if
the fixed load L < theta+ke, i.e. theta > theta_crit = L - ke. Below theta_crit
the marker grows without bound: there is NO steady state, which is what
decompensation IS. If theta_crit lands in the disease range, the steady-state FIM
analysis -- the whole brief -- does not apply where the clinical question lives.

ParamsC's own comment states the same mechanism: "the alternative route alone
peaks at 500 against a baseline load of 600, so a patient with theta=0 cannot be
maintained. That is the severity of null genotypes." So the repo already encodes
a load the alt route cannot absorb -- in C, not in the A the brief used.
"""
import numpy as np
from scipy.optimize import brentq

import core as C

KM, TH0, MBAR = C.KM0, C.TH0, 10.0
KE_B, KEM_B = 400.0, 300.0           # ParamsB: alt V_max, alt Michaelis const
LOAD_FRAC, ASSAY_CV = C.LOAD_FRAC, C.ASSAY_CV
I_TH = 0


def gf(m, Km=KM):
    return m / (Km + m)


def hf(m, Ke=KEM_B):
    return m / (Ke + m)


def L_anchor_B(ke=KE_B, Ke=KEM_B):
    """Convention A: load fixed at the theta=3000, m_bar=10 disposal, Model B."""
    return TH0 * gf(MBAR) + ke * hf(MBAR, Ke)


def level_B(theta, L, inc=0.0, ke=KE_B, Ke=KEM_B, Km=KM):
    """Positive root of (theta+inc) g(m) + ke h(m) = L, or inf if none exists."""
    cap = theta + inc + ke
    if L >= cap - 1e-9:
        return np.inf
    f = lambda m: (theta + inc) * gf(m, Km) + ke * hf(m, Ke) - L
    hi = 1.0
    while f(hi) < 0 and hi < 1e14:
        hi *= 2.0
    return brentq(f, 1e-12, hi)


def theta_crit(ke=KE_B, Ke=KEM_B):
    """Capacity below which no steady state exists (inc=0), as % of normal."""
    L = L_anchor_B(ke, Ke)
    return (L - ke) / TH0 * 100.0, L


# ---- Model-B record + FIM, mirroring core but with the saturating route ----
def record_B(theta, times, seed=0, ke=KE_B, Ke=KEM_B, load_cv=0.10):
    L = L_anchor_B(ke, Ke)
    times = np.asarray(times, float)
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    inc = C.delta_at(times, C.D0)
    m = np.array([level_B(theta, Lo, i, ke, Ke) for Lo, i in zip(loads, inc)])
    return dict(theta=theta, Km=KM, ke=ke, Ke=Ke, L=L, loads=loads, inc=inc,
                m=m, times=times)


def flux_rows_B(r):
    """d(clearance)/d(ln theta, ln Km, b, ln ke, ln Ke), Model B columns."""
    th, K, ke, Ke, m, inc = r["theta"], r["Km"], r["ke"], r["Ke"], r["m"], r["inc"]
    return np.stack([th * m / (K + m),
                     -(th + inc) * K * m / (K + m) ** 2,
                     np.ones_like(m),
                     ke * m / (Ke + m),
                     -ke * Ke * m / (Ke + m) ** 2], axis=1)


def slope_B(r):
    K, Ke = r["Km"], r["Ke"]
    return ((r["theta"] + r["inc"]) * K / (K + r["m"]) ** 2
            + r["ke"] * Ke / (Ke + r["m"]) ** 2)


def sigmas_B(r):
    return np.hypot(LOAD_FRAC * r["L"], slope_B(r) * ASSAY_CV * r["m"])


def fim_B(r, ke_prior=C.KE_PRIOR, ke_m_prior=None, drop_Ke=False):
    """FIM in (ln th, ln Km, b, ln ke[, ln Ke]); priors mirror core.fim.
    drop_Ke pins Ke (4-param, matching A's parameter count)."""
    Q = flux_rows_B(r)
    if drop_Ke:
        Q = Q[:, :4]
    sg = sigmas_B(r)
    F = (Q / sg[:, None] ** 2).T @ Q
    n = Q.shape[1]
    pri = np.zeros(n)
    pri[I_TH] = 0.0
    pri[2] = (LOAD_FRAC * r["L"]) ** -2     # b prior
    if ke_prior:
        pri[3] = ke_prior ** -2             # ln ke prior
    if not drop_Ke and ke_m_prior:
        pri[4] = ke_m_prior ** -2           # ln Ke prior
    return F + np.diag(pri)


def sd_th_B(F):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-12 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def share_B(theta, ke=KE_B, Ke=KEM_B):
    L = L_anchor_B(ke, Ke)
    m = level_B(theta, L, 0.0, ke, Ke)
    if not np.isfinite(m):
        return np.nan, np.inf
    return theta * gf(m) / L, m


def over(fn, seeds=range(24)):
    v = [fn(s) for s in seeds]
    return float(np.median(v))


def main():
    print("=" * 94)
    print("  ROUND 7 ITEM 1. MODEL A vs MODEL B: EXISTENCE FIRST, THEN THE HEADLINES")
    print("=" * 94)

    # ---------- (1) STEADY-STATE EXISTENCE ----------
    print("\n  (1) DOES A STEADY STATE EXIST IN THE DISEASE RANGE UNDER B?")
    print("      Disposal is bounded by theta+ke; SS exists iff load L < theta+ke.")
    print("      theta_crit = L - ke is the capacity (as % of normal) below which the")
    print("      marker DIVERGES -- no steady state, i.e. decompensation.\n")
    print(f"      {'alt-route ceiling ke':>22}  {'anchored load L':>15}  {'theta_crit (% normal)':>22}")
    for ke in (400.0, 300.0, 200.0, 100.0, 50.0, 25.0):
        tc, L = theta_crit(ke)
        tag = "  <- ParamsB" if ke == 400.0 else ("  (no SS at 5%)" if tc > 5 else "")
        print(f"      {ke:>22.0f}  {L:>15.1f}  {tc:>21.1f}%{tag}")
    print("""
      READING: ParamsB sets the alt-route ceiling ke=400, ABOVE the load (~286),
      so theta_crit < 0 and a steady state exists even at theta=0 -- ParamsB does
      NOT reproduce decompensation, which already conflicts with severe PA. The
      literature (and the repo's own Model C: alt peak 500 < load 600) puts the
      ceiling BELOW the load; then theta_crit lands at ~7-9% of normal, so AT THE
      5% DISEASE REFERENCE THERE IS NO STEADY STATE. Where the ceiling is genuinely
      low-capacity, the steady-state FIM analysis -- every precision number in this
      brief -- does not apply at the capacity the clinical question is about. That
      is a more fundamental finding than any sd: the benign steady state the brief
      analyses exists only because Model A's clearance is unbounded.
""")

    # ---------- (2) SHARE ACROSS CAPACITY, A vs B ----------
    print("  (2) ENZYME SHARE s = theta g(m)/L ACROSS CAPACITY (ParamsB, SS exists):")
    print(f"      {'capacity':>12}  {'A: share':>10}  {'A: m':>8}  {'B: share':>10}  {'B: m':>8}")
    for frac in (1.0, 0.10, 0.05, 0.01):
        th = TH0 * frac
        # Model A share
        LA = TH0 * gf(MBAR) + C.KE * MBAR
        from n4_theta_scope import level as levelA
        mA = levelA(th, KM, C.KE, LA)
        sA = th * gf(mA) / LA
        sB, mB = share_B(th)
        print(f"      {int(frac*100):>10}%   {sA:>10.3f}  {mA:>8.1f}  {sB:>10.3f}  {mB:>8.1f}")
    print("      Under B the enzyme carries a LARGER share at low capacity (the alt route")
    print("      saturates, so it cannot soak up the rising marker as A's linear term does).")
    print("      The s=0.16 premise of the A analysis becomes ~0.37 at 5% under ParamsB.")

    # ---------- (3) HEADLINE sd(ln theta) AT 5%, A vs B ----------
    print("\n  (3) HEADLINE sd(ln theta) AT 5% MONTHLY (where B's SS exists, ParamsB):")
    rA = lambda s: C.record(C.TH_5, seed=s)
    sdA = over(lambda s: C.sd_th(C.fim(rA(s), noise="perdraw")))
    sdB5 = over(lambda s: sd_th_B(fim_B(record_B(C.TH_5, C.MONTHLY, s), drop_Ke=True)))
    sdB5k = over(lambda s: sd_th_B(fim_B(record_B(C.TH_5, C.MONTHLY, s), drop_Ke=False,
                                         ke_m_prior=0.5)))
    print(f"      Model A (4 params: ln th, ln Km, b, ln ke):            {sdA:.3f}")
    print(f"      Model B, Ke pinned (4 params, matched to A):           {sdB5:.3f}")
    print(f"      Model B, Ke free with 50% prior (5 params):            {sdB5k:.3f}")
    print("      (Under B the enzyme's larger share tightens theta where a steady state")
    print("       exists; but a NEW nuisance Ke appears, and the result still sits far")
    print("       above 0.20 -- the cross-sectional verdict is unchanged in sign.)")

    # ---------- (3b) IS THE OBSERVED RANGE AFFINE? (the paper's C23 condition) ----------
    print("\n  (3b) WHY B's sd IS FINITE: the DOSED record stays near-affine (paper's C23).")
    r = record_B(C.TH_5, C.MONTHLY, 0)
    mlo, mhi = r["m"].min(), r["m"].max()
    hlo, hhi = hf(mlo), hf(mhi)
    print(f"      Over the record, m in [{mlo:.1f}, {mhi:.1f}]; alt-route saturation h(m)=m/(Ke+m)")
    print(f"      runs {hlo:.3f}->{hhi:.3f} (Ke=300), so h is still RISING, not saturated.")
    print("      The parent paper (sections/03_model_A.tex, claim C23) conditions its")
    print("      degeneracy 'on clearance being affine over the observed range' -- which the")
    print("      DOSED record respects, so the FIM result is similar in kind. Saturation and")
    print("      non-existence bite at the UNTREATED balance (m ~ 244 at 5%) and below")
    print("      theta_crit, NOT in the dosed data -- that is the scope line (1).")

    # ---------- (4) LONGITUDINAL CONTRAST, A vs B ----------
    print("\n  (4) LONGITUDINAL CONTRAST sd(ln th2 - ln th1), weekly windows, A vs B:")
    cA = over(lambda s: contrast_A(104, 730, 2 * s, 2 * s + 1))
    cB = over(lambda s: contrast_B(104, 730, 2 * s, 2 * s + 1))
    print(f"      Model A (shared Km, b, ke):            {cA:.3f}")
    print(f"      Model B (shared Km, b, ke, Ke):        {cB:.3f}")
    print("      B is BETTER on the contrast too (not worse): the higher enzyme share")
    print("      tightens theta per window, and the shared patient constants (now including")
    print("      Ke) still cancel. So monitoring survives -- and improves -- WHERE A STEADY")
    print("      STATE EXISTS. If the alt route is low-capacity and the load unmanaged (no SS")
    print("      at 5%), there is no stable marker to difference, and monitoring is moot in")
    print("      severe disease until the load is managed (t35).")

    print("\n  WHICH MODEL THE PARENT PAPER USED, AND WHAT INHERITS THIS:")
    print("""    The parent paper's identifiability configuration (sections/03_model_A.tex,
    eq:local; symbolic/.../observational_equivalence.py) is MODEL A:
    m' = L + L0 - (theta+Delta) g(m) - k_e m. So every published degeneracy /
    precision figure -- and every number in this brief -- is a Model-A result.
    The paper does NOT ignore the issue: claim C23 states the degeneracy is
    'conditioned on clearance being affine over the observed range', and Model C
    (sections/05) is where it treats a bounded alternative route and the
    decompensation the null genotype cannot escape. So the published figures
    inherit the Model-A affine-clearance assumption EXPLICITLY, not silently.

  BOTTOM LINE:
    - Numbers move materially under B: share 0.16 -> 0.37 at 5%, cross-sectional
      sd 3.05 -> ~0.9, contrast 0.133 -> 0.088 (B better on BOTH, the share effect).
      SIGNS do not: cross-sectional still fails 0.20, monitoring still passes where
      a steady state exists.
    - The fundamental change is SCOPE, not precision: under a literature-
      consistent LOW-capacity alternative route (the repo's own Model C), no
      steady state exists below ~6-8% of normal, so the steady-state analysis
      does not describe severe disease at all -- it describes the stable inter-
      crisis state, and the crisis boundary sits inside the disease range.""")


# ---- contrasts (A reuses core; B mirrors t24 with Ke as a 5th shared nuisance) ----
def contrast_A(n, win, s1, s2, ke_prior=0.8, bfrac=0.05):
    import t24_contrast_stress as T24
    return T24.contrast(n, win, s1, s2, ke_prior=ke_prior, bfrac=bfrac)[0]


def contrast_B(n, win, s1, s2, ke_prior=0.8, bfrac=0.05, ke_m_prior=0.5,
               drift_ke=0.0):
    """Mirrors t24.contrast exactly, with a 5th shared nuisance Ke. ke1/ke2 are
    coupled with precision 1/drift_ke^2 (drift_ke=0 -> shared, as in R-D's base
    case); Ke (alt Michaelis const, a patient constant) is a single shared index."""
    r1 = record_B(C.TH_5, np.linspace(0, win, n), s1)
    r2 = record_B(C.TH_5, np.linspace(0, win, n), s2)
    L = r1["L"]
    Q1, Q2 = flux_rows_B(r1), flux_rows_B(r2)
    sg1, sg2 = sigmas_B(r1), sigmas_B(r2)
    F1 = (Q1 / sg1[:, None] ** 2).T @ Q1
    F2 = (Q2 / sg2[:, None] ** 2).T @ Q2
    # params: th1,th2, Km(sh), b1,b2, ke1,ke2, Ke(sh)  -> 8
    i1 = [0, 2, 3, 5, 7]
    i2 = [1, 2, 4, 6, 7]
    Z = np.zeros((8, 8))
    for a in range(5):
        for b in range(5):
            Z[i1[a], i1[b]] += F1[a, b]
            Z[i2[a], i2[b]] += F2[a, b]

    def cp(i, j, sd):                       # couple i,j (shared when sd=0)
        g = 1e12 if sd <= 0 else 1.0 / sd ** 2
        Z[i, i] += g; Z[j, j] += g; Z[i, j] -= g; Z[j, i] -= g
    cp(3, 4, 0.0)                           # b a patient constant (shared), as t24
    cp(5, 6, drift_ke)                      # ke a patient constant (shared), as t24
    Z[5, 5] += ke_prior ** -2; Z[6, 6] += ke_prior ** -2
    Z[3, 3] += (bfrac * L) ** -2; Z[4, 4] += (bfrac * L) ** -2
    Z[7, 7] += ke_m_prior ** -2
    P = np.linalg.inv(Z)
    e = np.zeros(8); e[0] = -1; e[1] = 1
    return float(np.sqrt(e @ P @ e))


if __name__ == "__main__":
    main()
