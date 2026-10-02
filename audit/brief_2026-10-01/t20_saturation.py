"""
Round-3 close-out A. The saturation regime, the integral-estimator collapse,
and whether 'needs an excursion' and 'needs saturation' are one constraint.

The enzyme is NOT saturated in the disease range: at 5% m/K_m ~ 0.54, g ~ 0.35,
so lambda = (theta)g'(m) + k_e is ~1.13 k_e with a theta-dependence of order the
precision we want -- not '~k_e with no theta'. This script makes that explicit,
explains the TB-A integral collapse as a g-structure (near-collinearity) problem,
and tests the excursion = saturation hypothesis.
"""
import numpy as np

import core as C
import dyn as D
import t6_estimators as T6


def g_of(m, Km=C.KM0):
    return m / (Km + m)


def design_cond(m, t, ex, Km):
    """condition number of the integral design [tau, -int g, -int m] (b')."""
    X, _ = T6._vp_design(m, t, ex, Km, cumulative=False)
    return float(np.linalg.cond(X))


def main():
    print("=" * 96)
    print("  CLOSE-OUT A. SATURATION REGIME, THE INTEGRAL COLLAPSE, EXCURSION = SATURATION")
    print("=" * 96)

    # ---- (1) the m/K_m regime across capacity -----------------------------
    print("\n  (1) THE SATURATION REGIME ACROSS CAPACITY (K_m = 100).")
    print("      The enzyme is sub-saturated through the disease range; it only")
    print("      approaches g -> 1 at very low capacity or under a big load excursion.")
    print(f"      {'capacity':>9}  {'untreated m':>11}  {'m/K_m':>6}  {'g=m/(Km+m)':>10}"
          f"  {'eps=Km/(Km+m)':>13}  {'on-therapy g (24 recs)':>24}")
    for pct in (1, 5, 10, 20, 50, 100):
        th = pct / 100 * C.TH0
        L, m0 = C.N4.anchor(th, C.KM0, C.KE, "A")
        gs = np.concatenate([g_of(C.record(th, seed=s)["m"]) for s in C.SEEDS])
        print(f"      {pct:>7}%   {m0:>11.2f}  {m0/C.KM0:>6.2f}  {g_of(m0):>10.3f}"
              f"  {C.KM0/(C.KM0+m0):>13.3f}  {gs.min():.3f}-{gs.max():.3f} (med {np.median(gs):.3f})")
    lam5 = (C.TH_5 * C.KM0 / (C.KM0 + 54.02) ** 2 + C.KE)
    print(f"      At 5%: lambda = theta g'(m) + k_e = {lam5:.3f} = {lam5/C.KE:.3f} k_e;")
    print(f"      the theta term is {100*(lam5-C.KE)/lam5:.0f}% of lambda -- NOT negligible.")

    # ---- (2) the integral-estimator collapse on TB-A ----------------------
    print("\n  (2) WHY THE INTEGRAL ESTIMATORS COLLAPSE ON TB-A.")
    print("      (b) 'saturated' assumes g = 1, so theta enters only as L0 - theta and")
    print("      is structurally unidentifiable. (b')/(b'') use the proper")
    print("      g = m/(K_m+m), but theta's regressor is INT g: it is identifiable only")
    print("      when g VARIES across the record. When g is small and nearly flat, the")
    print("      INT g column is collinear with the constant (tau) and INT m columns.")
    for nm in ("A", "B"):
        b = T6.bench(nm)
        g = g_of(b["m"])
        cond = design_cond(b["m"], b["t"], b["ex"], C.KM0)
        est = T6.est_int_gen(b["m"], b["t"], b["ex"])          # noiseless
        print(f"      TB-{nm} (theta={b['th']:.0f}): g spans {g.min():.3f}-{g.max():.3f} "
              f"(range {g.max()-g.min():.3f}); design cond[tau,∫g,∫m] = {cond:.1e};")
        print(f"             noiseless (b') theta-hat = {est['th']:.1f} against truth "
              f"{b['th']:.0f}  ({'RECOVERED' if abs(est['th']-b['th'])<0.1*b['th'] else 'COLLAPSED'})")
    print("      NUANCE (correcting the brief's framing): (b')/(b'') are NOT g=1-")
    print("      misspecified -- they use the proper g and RECOVER theta noiselessly even")
    print("      on TB-A. Their TB-A collapse is ILL-CONDITIONING: cond ~2e4, so the 5%")
    print("      assay noise is amplified ~2e4-fold and theta blows up. g is nearly flat")
    print("      (0.21-0.35), so the ∫g column barely separates. On TB-B the excursion")
    print("      sweeps g to 0.75, cond falls ~7x, and theta survives the noise. Only (b)")
    print("      assumes g=1 (a ~1/0.35 ~ 3 misspecification at 5%), and it cannot")
    print("      identify theta regardless. So: ill-conditioning from sub-saturation, not")
    print("      a misspecified flux, is what sinks the integral estimators on TB-A.")

    # ---- (3) excursion = saturation ---------------------------------------
    print("\n  (3) IS 'NEEDS AN EXCURSION' THE SAME AS 'NEEDS SATURATION'?")
    print("      Load multiple f -> peak marker, peak g, design cond, transient sd(ln theta)")
    print("      (theta=5%, dense 5-min 48 h, k_e prior 0.8 + K_m from record + L0 to 1%):")
    print(f"      {'f':>5}  {'peak m':>7}  {'peak g':>6}  {'cond[τ,∫g,∫m]':>14}  {'sd(ln θ)':>9}  admissible")
    import t7_design as T7
    L = D.LTOT
    tt = np.arange(0, 2.0 + 1e-9, 5 / 60 / 24)
    for f in (0.5, 1.0, 1.177, 1.5, 2.0, 2.3, 3.0, 6.0):
        if f <= 1.0:
            ex = D.Experiment([(-1, 0.5 * L), (0, (f - 0.5) * L)], L0=0.5 * L)
        else:
            ex = D.Experiment([(-1, 0.5 * L), (0, 0.5 * L + (f - 1) * L)], L0=0.5 * L)
        m, _ = D.simulate(D.params(ex=ex), ex, tt, sens=False)
        gpk = g_of(m.max())
        cond = design_cond(m, tt, ex, C.KM0)
        Fd = D.fim(D.params(ex=ex), ex, tt)
        sdth = T7.joint_record(Fd, {"ke": 0.8, "Km": T7.KM_REC, "L0": 0.01 * L})
        adm = "yes" if f <= 1.177 + 1e-9 else "NO"
        print(f"      {f:>5.2f}  {m.max():>7.1f}  {gpk:>6.3f}  {cond:>14.1e}  {sdth:>9.3f}  {adm}")
    print("      theta becomes identifiable (sd and cond both fall) as the excursion")
    print("      drives the peak g UP from its sub-saturated baseline ~0.35. Clearing")
    print("      0.20 needs f ~ 2.3 (peak g ~ 0.57, marker ~131, well above K_m=100) --")
    print("      the enzyme pushed toward saturation, though not all the way to g=1. That")
    print("      f ~ 2.3 is ~2x past the 1.177x admissible edge. So 'needs an excursion'")
    print("      and 'needs to drive the enzyme toward saturation' are ONE constraint:")
    print("      theta's integral regressor ∫g only separates when g is swept upward, and")
    print("      sweeping it means adding substrate past the crisis fold. The dynamic")
    print("      route fails for the same reason it would be unsafe.")


if __name__ == "__main__":
    main()
