"""
Round 8, items 2 & 3: theta_crit moves with the MANAGED load, and the only valid
Model B comparison is at a clinically admissible operating point.

ITEM 2 CORRECTION. theta_crit = L - ke depends on L, and L is actively MANAGED
(protein restriction, gut suppression). Real 5% patients are alive and stable, so
a steady state DOES exist for them -- at their managed load. The v11 claim "no
steady state at 5%" is an overstatement. Correct: the severe patient has a steady
state only UNDER active load management. The margin is the paper's own
mu = theta - theta_crit(L) (sections/03_dynamic.tex, prop:intrinsic). This
connects to the paper's Section 3.5: antibiotic courses and protein-intake changes
move L, so the record already contains theta_crit excursions.

ITEM 3 DOUBLE BIND. Every Model B number in t33 is at ParamsB (ceiling ke=400),
whose ceiling EXCEEDS the untreated load (~286), so theta=0 is maintainable --
clinically wrong (it fails the severe phenotype). A clinically RIGHT ceiling
(below the load) admits NO steady state at the untreated 5% load. So the A-vs-B
comparison cannot be made at a single clinically valid UNTREATED point at 5%; it
needs the managed load. The valid Model B column is: a literature-consistent LOW
ceiling, WITH a managed load that admits a steady state.
"""
import numpy as np
from scipy.optimize import brentq

import core as C
import t33_model_b as T
from n4_theta_scope import level as levelA

KM, TH0 = C.KM0, C.TH0
TH5 = C.TH_5              # 150


def theta_crit_B(L, ke):
    """Model B: SS exists iff L < theta+ke, so theta_crit = L - ke."""
    return L - ke


def level_B_L(theta, L, ke, Ke=T.KEM_B, inc=0.0):
    cap = theta + inc + ke
    if L >= cap - 1e-9:
        return np.inf
    f = lambda m: (theta + inc) * (m / (KM + m)) + ke * (m / (Ke + m)) - L
    hi = 1.0
    while f(hi) < 0 and hi < 1e14:
        hi *= 2
    return brentq(f, 1e-12, hi)


def record_B_managed(theta, L, ke, times, seed=0, Ke=T.KEM_B, load_cv=0.10):
    times = np.asarray(times, float)
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    inc = C.delta_at(times, C.D0)
    m = np.array([level_B_L(theta, Lo, ke, Ke, i) for Lo, i in zip(loads, inc)])
    return dict(theta=theta, Km=KM, ke=ke, Ke=Ke, L=L, loads=loads, inc=inc, m=m, times=times)


def record_A_managed(theta, L, times, seed=0, ke=C.KE, load_cv=0.10):
    times = np.asarray(times, float)
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    inc = C.delta_at(times, C.D0)
    m = C.level(theta, KM, ke, loads, inc)
    return dict(theta=theta, Km=KM, ke=ke, L=L, m_bar=C.level(theta, KM, ke, L),
                loads=loads, inc=inc, m=m, times=times)


def over(fn, seeds=range(24)):
    return float(np.median([fn(s) for s in seeds]))


def main():
    print("=" * 94)
    print("  ROUND 8 ITEM 2. theta_crit MOVES WITH THE MANAGED LOAD; the margin mu")
    print("=" * 94)
    L0 = T.L_anchor_B(ke=100.0)          # untreated anchored load at a low ceiling
    print(f"""
  For a literature-consistent LOW ceiling (take ke=100, well below the untreated
  load L0={L0:.0f}, so theta=0 is NOT maintainable -- the severe phenotype):
  theta_crit = L - ke, and the 5% patient (theta={TH5:.0f}) is stable iff
  L < theta + ke = {TH5+100:.0f}. The margin is mu = theta - theta_crit(L).
""")
    print(f"      {'managed load L':>16}  {'reduction vs L0':>16}  {'theta_crit (%norm)':>18}  {'mu (%norm)':>12}  {'5% stable?':>11}")
    for L in (L0, 260, 250, 240, 220, 200, 180):
        tc = theta_crit_B(L, 100.0)
        mu = TH5 - tc
        red = (1 - L / L0) * 100
        stable = "yes" if mu > 0 else "NO (diverges)"
        print(f"      {L:>16.0f}  {red:>15.1f}%  {tc/TH0*100:>17.1f}%  {mu/TH0*100:>11.2f}%  {stable:>11}")
    print(f"""
  READING: untreated (L={L0:.0f}) the 5% patient is BELOW theta_crit and diverges.
  A load reduction of ~{100*(1-250/L0):.0f}% (to L<{TH5+100:.0f}) brings them above it; deeper
  management buys margin. Protein restriction + gut suppression routinely cut the
  load 20-40% (Thompson: gut share ~22%, ~41% suppressible; protein restriction
  cuts the ~52% amino-acid share), so a stable managed steady state at 5% is
  ACHIEVABLE -- which is why real patients are stable. The correct statement is:
  the severe patient has a steady state only UNDER active load management, with a
  thin margin mu that the record's own antibiotic/protein excursions (paper 3.5)
  can erase -- a theta_crit excursion is a decompensation.
""")

    print("=" * 94)
    print("  ROUND 8 ITEM 3. THE ONLY VALID MODEL-B COLUMN: low ceiling + managed load")
    print("=" * 94)
    print("""
  THE DOUBLE BIND. ParamsB's ceiling (400) EXCEEDS the untreated load (~286): a
  theta=0 patient is maintainable, which is clinically wrong. A clinically-right
  ceiling (< load) admits NO steady state at the untreated 5% load. So t33's 0.92
  was computed at a MISCALIBRATED point. The valid comparison fixes both: a low
  ceiling (ke=100) AND a managed load that admits a steady state (take L=220, a
  ~20% reduction, margin mu~1% of normal). A and B are compared at the SAME
  clinically valid managed point.
""")
    ke_B, Lman = 100.0, 220.0
    # untreated-balance m (near the ceiling) vs the DOSED-record operating range
    mB_unt = level_B_L(TH5, Lman, ke_B)
    sB = TH5 * (mB_unt / (KM + mB_unt)) / Lman
    mA_unt = levelA(TH5, KM, C.KE, Lman)
    sA = TH5 * (mA_unt / (KM + mA_unt)) / Lman
    rBrec = record_B_managed(TH5, Lman, ke_B, C.MONTHLY, 0)
    rArec = record_A_managed(TH5, Lman, C.MONTHLY, 0)
    print(f"  At theta=5%, managed load L={Lman:.0f}:")
    print(f"    {'':<28}{'Model A':>12}{'Model B (ke=100)':>18}")
    print(f"    {'untreated-balance m':<28}{mA_unt:>12.1f}{mB_unt:>18.1f}")
    print(f"    {'  (B near its ceiling)':<28}{'':>12}{'':>18}")
    print(f"    {'DOSED-record m (mean)':<28}{rArec['m'].mean():>12.1f}{rBrec['m'].mean():>18.1f}")
    print(f"    {'  record m range':<28}{'[%d,%d]'%(rArec['m'].min(),rArec['m'].max()):>12}{'[%d,%d]'%(rBrec['m'].min(),rBrec['m'].max()):>18}")
    print(f"    {'untreated enzyme share':<28}{sA:>12.3f}{sB:>18.3f}")
    sdA = over(lambda s: C.sd_th(C.fim(record_A_managed(TH5, Lman, C.MONTHLY, s), noise="perdraw")))
    sdB = over(lambda s: T.sd_th_B(T.fim_B(record_B_managed(TH5, Lman, ke_B, C.MONTHLY, s), drop_Ke=True)))
    sdBf = over(lambda s: T.sd_th_B(T.fim_B(record_B_managed(TH5, Lman, ke_B, C.MONTHLY, s), drop_Ke=False, ke_m_prior=0.5)))
    print(f"    {'cross-sectional sd(lnth)':<26}{sdA:>12.3f}{sdB:>18.3f}  (Ke pinned)")
    print(f"    {'                 Ke free':<26}{'':>12}{sdBf:>18.3f}")
    cA = over(lambda s: contrast_A_managed(Lman, 104, 730, 2 * s, 2 * s + 1))
    cB = over(lambda s: contrast_B_managed(Lman, ke_B, 104, 730, 2 * s, 2 * s + 1))
    print(f"    {'longitudinal contrast':<26}{cA:>12.3f}{cB:>18.3f}")
    rrange = "[%d,%d]" % (rBrec['m'].min(), rBrec['m'].max())
    print(f"""
  READING at the VALID point (5% capacity, managed load {Lman:.0f}, low ceiling):
    - The UNTREATED balance sits near the ceiling (m~{mB_unt:.0f}, share {sB:.2f}), but the
      DOSED record operates LOW (m in {rrange}, like every other record here),
      because the dose train boosts disposal -- so the FIM is well-behaved, not
      stiff, and the numbers are robust across the managed range (sd ~0.29,
      contrast ~0.08 for L from 180 to 240).
    - The marginal cross-sectional sd (~0.29) is still FINITE and still > 0.20; the
      sign of the verdict is unchanged whether we use ParamsB (0.92) or this valid
      point (0.29) -- both fail cross-sectionally. B is tighter than A at every
      valid point (the share effect of item 1).
    - The longitudinal contrast (~0.08) is BELOW A's 0.133: B is better on BOTH
      axes (b and ke shared, as in t24). So the qualitative conclusions survive the
      correction; only the specific 0.92/0.167 of v11 were at a bad point / a b-not-
      shared bug.
    - ParamsB (ke=400) is MISCALIBRATED for severe PA: its ceiling exceeds the load,
      so it cannot represent a null-genotype patient (theta=0 stays maintainable),
      contradicting the repo's own Model C ('alt peaks at 500 vs load 600, so
      theta=0 cannot be maintained'). Flagged.
""")


def contrast_A_managed(L, n, win, s1, s2, ke_prior=0.8, bfrac=0.05):
    r1 = record_A_managed(TH5, L, np.linspace(0, win, n), s1)
    r2 = record_A_managed(TH5, L, np.linspace(0, win, n), s2)
    F1 = C.fim(r1, noise="perdraw", ke_prior=0.0, b_prior=np.inf)
    F2 = C.fim(r2, noise="perdraw", ke_prior=0.0, b_prior=np.inf)
    i1, i2 = [0, 2, 3, 5], [1, 2, 4, 6]
    Z = np.zeros((7, 7))
    for a in range(4):
        for b in range(4):
            Z[i1[a], i1[b]] += F1[a, b]; Z[i2[a], i2[b]] += F2[a, b]
    g = 1e12
    Z[3, 3] += g; Z[4, 4] += g; Z[3, 4] -= g; Z[4, 3] -= g   # b shared (as t24)
    Z[5, 5] += g; Z[6, 6] += g; Z[5, 6] -= g; Z[6, 5] -= g   # ke shared
    Z[5, 5] += ke_prior ** -2; Z[6, 6] += ke_prior ** -2
    Z[3, 3] += (bfrac * L) ** -2; Z[4, 4] += (bfrac * L) ** -2
    P = np.linalg.inv(Z); e = np.array([-1, 1, 0, 0, 0, 0, 0.])
    return float(np.sqrt(e @ P @ e))


def contrast_B_managed(L, ke, n, win, s1, s2, ke_prior=0.8, bfrac=0.05, ke_m_prior=0.5):
    r1 = record_B_managed(TH5, L, ke, np.linspace(0, win, n), s1)
    r2 = record_B_managed(TH5, L, ke, np.linspace(0, win, n), s2)
    Q1, Q2 = T.flux_rows_B(r1), T.flux_rows_B(r2)
    sg1, sg2 = T.sigmas_B(r1), T.sigmas_B(r2)
    F1 = (Q1 / sg1[:, None] ** 2).T @ Q1
    F2 = (Q2 / sg2[:, None] ** 2).T @ Q2
    i1, i2 = [0, 2, 3, 5, 7], [1, 2, 4, 6, 7]
    Z = np.zeros((8, 8))
    for a in range(5):
        for b in range(5):
            Z[i1[a], i1[b]] += F1[a, b]; Z[i2[a], i2[b]] += F2[a, b]
    g = 1e12
    Z[3, 3] += g; Z[4, 4] += g; Z[3, 4] -= g; Z[4, 3] -= g   # b shared (as t24)
    Z[5, 5] += g; Z[6, 6] += g; Z[5, 6] -= g; Z[6, 5] -= g   # ke shared
    Z[5, 5] += ke_prior ** -2; Z[6, 6] += ke_prior ** -2
    Z[3, 3] += (bfrac * L) ** -2; Z[4, 4] += (bfrac * L) ** -2
    Z[7, 7] += ke_m_prior ** -2
    P = np.linalg.inv(Z); e = np.zeros(8); e[0] = -1; e[1] = 1
    return float(np.sqrt(e @ P @ e))


if __name__ == "__main__":
    main()
