"""
3.0.5 and 3.0.6 -- does the model earn its place, and where does the rate
estimand stop working?

3.0.5 asks for the baseline: read the rate straight off the observed marker
through the published rate ratio, with no model, no theta, no identification.
The model-based route only earns its place if it beats that.

Three things have to be in the comparison for it to be fair:

  (i)   at the patient's OWN operating point, where the direct read-off has an
        observation and the model merely reconstructs one;
  (ii)  under a COUNTERFACTUAL load or dose, where the direct read-off has no
        observation at all and cannot answer;
  (iii) the uncertainty in the MAP itself, which both routes inherit equally
        and which neither can reduce.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D

P = D.P
L0 = D.load(0.0)

# --- the map, and its own uncertainty -------------------------------------
# RR 0.70 (95% CI 0.60-0.90) for a 50% LOWER 3-HP  =>  rate ~ m^beta
RR, RR_LO, RR_HI = 0.70, 0.60, 0.90
BETA = -np.log2(RR)
BETA_LO, BETA_HI = -np.log2(RR_HI), -np.log2(RR_LO)     # note the swap
SD_BETA = (BETA_HI - BETA_LO) / (2 * 1.96)


def level_at(th, Km, ke, L):
    B = ke * Km + th - L
    d = B ** 2 + 4.0 * ke * L * Km
    if not np.isfinite(d) or d < 0 or ke <= 0:
        return np.nan
    return (-B + np.sqrt(d)) / (2.0 * ke)


def main():
    print(__doc__)
    print("=" * 84)
    print("THE MAP'S OWN UNCERTAINTY -- inherited by both routes")
    print("=" * 84)
    print(f"  RR = {RR} (95% CI {RR_LO}-{RR_HI}) per halving")
    print(f"  beta = -log2(RR) = {BETA:.4f}, 95% CI {BETA_LO:.4f} to {BETA_HI:.4f}")
    print(f"  so the EXPONENT is known to a factor of {BETA_HI/BETA_LO:.1f}")
    print(f"  sd(beta) ~ {SD_BETA:.4f}")
    print()
    print("  A rate prediction at a level m carries ln(m/m_ref)*sd(beta) from this")
    print("  alone, and NEITHER route can reduce it:")
    print(f"{'m / m_ref':>12} {'from the map alone':>22}")
    for r in (1.0, 1.25, 1.5, 2.0, 3.0, 5.0):
        print(f"{r:12.2f} {abs(np.log(r))*SD_BETA:22.4f}")

    print()
    print("=" * 84)
    print("3.0.5  (i) AT THE PATIENT'S OWN OPERATING POINT")
    print("=" * 84)
    t = np.arange(0.0, 730.0, 30.4)
    lv, inc = D.design_record(t)
    n = len(lv)
    m_ref = float(level_at(P.theta, P.Km, P.ke, L0))

    ens = D.fit_ensemble(lv, P.theta, P.Km, P.ke, reps=400, seed=11, increments=inc)
    th, Km, ke = ens[:, 0], ens[:, 1], ens[:, 2]
    mh = np.array([level_at(a, b, c, L0) for a, b, c in zip(th, Km, ke)])
    ok = np.isfinite(mh) & (mh > 0)
    sd_lnm_model = float(np.std(np.log(mh[ok]), ddof=1))

    # direct read-off: average the observed levels, assay noise only
    sd_lnm_direct = D.ASSAY_CV / np.sqrt(n)

    print(f"  design: monthly, n = {n}, assay CV = {D.ASSAY_CV:.0%}")
    print()
    print(f"{'route':>34} {'sd(ln m)':>10} {'x beta':>9} {'+ map':>9} {'95% width':>12}")
    for lab, s in (("direct read-off, no model", sd_lnm_direct),
                   ("model-based, via theta-hat", sd_lnm_model)):
        core = BETA * s
        tot = float(np.hypot(core, 0.0))     # at m = m_ref the map term vanishes
        print(f"{lab:>34} {s:10.4f} {core:9.4f} {tot:9.4f} "
              f"{np.exp(1.96*tot)/np.exp(-1.96*tot):12.2f}x")
    print()
    print(f"  ratio of interval widths, model / direct: "
          f"{(BETA*sd_lnm_model)/(BETA*sd_lnm_direct):.1f}x WIDER")
    print("  At the operating point the model reconstructs a level that was")
    print("  measured directly, and pays parameter uncertainty to do it.")

    print()
    print("=" * 84)
    print("3.0.5  (ii) UNDER A COUNTERFACTUAL -- what the direct route cannot do")
    print("=" * 84)
    print("  'What would this patient's rate be if the dose were halved, or if")
    print("   their load rose 50%?' There is no observation at that condition.")
    print()
    print(f"{'counterfactual':>34} {'model sd(ln rate)':>18} {'direct':>22}")
    for lab, Lf, dsc in (("load 1.5x baseline", 1.5, None),
                         ("load 3.0x baseline", 3.0, None)):
        mcf = np.array([level_at(a, b, c, Lf * L0) for a, b, c in zip(th, Km, ke)])
        o = np.isfinite(mcf) & (mcf > 0)
        s = BETA * float(np.std(np.log(mcf[o]), ddof=1))
        print(f"{lab:>34} {s:18.4f} {'no answer possible':>22}")
    print()
    print("  This is the only column where the model has anything the data do not.")

    print()
    print("=" * 84)
    print("3.0.5  (iii) THE ANCHOR -- and this blocks BOTH routes")
    print("=" * 84)
    print("  rate(m) = rate_ref * (m/m_ref)^beta needs rate_ref AND m_ref.")
    print("  The natural history study reports a cohort mean annualised MDE")
    print("  frequency of 1.0 (SD 1.6), so rate_ref exists. It reports NO")
    print("  concentration range and no reference 3-HP value, so m_ref does not.")
    print("  Without m_ref neither route yields a rate: both yield only a RATIO")
    print("  between two levels, and the model adds nothing to that.")


if __name__ == "__main__":
    main()


def scope_of_the_rate_estimand():
    """
    3.0.6 -- where the theta -> marker map stops carrying information.

    Stated in the units a clinician reads: the MARKER LEVEL, not theta. As
    capacity falls the marker rises towards the ceiling L0/k_e set by the
    alternative channel alone, and the elasticity d ln m / d ln theta goes to
    zero there. Above some observable level, two patients with very different
    residual capacity are indistinguishable on the marker.
    """
    from scipy.optimize import brentq
    print()
    print("=" * 84)
    print("3.0.6  WHERE THE RATE ESTIMAND STOPS WORKING, in observable units")
    print("=" * 84)
    ceil = L0 / P.ke
    m_bar = float(level_at(P.theta, P.Km, P.ke, L0))
    print(f"  marker ceiling as theta -> 0:  m -> L0/k_e = {ceil:.2f} "
          f"= {ceil/m_bar:.2f} x the managed level")
    print()
    h = 1e-5

    def elast(th):
        return abs((np.log(level_at(th * np.exp(h), P.Km, P.ke, L0))
                    - np.log(level_at(th * np.exp(-h), P.Km, P.ke, L0))) / (2 * h))

    print(f"{'theta':>11} {'% of nominal':>13} {'marker m':>10} {'m / m_bar':>10} "
          f"{'|elasticity|':>13}")
    for th in (3000., 1500., 900., 600., 300., 150., 60., 30., 10., 3., 1.):
        print(f"{th:11.0f} {100*th/P.theta:12.1f}% {level_at(th,P.Km,P.ke,L0):10.2f} "
              f"{level_at(th,P.Km,P.ke,L0)/m_bar:10.2f} {elast(th):13.4f}")
    print()
    for cut in (0.50, 0.25, 0.10, 0.05):
        try:
            th_c = brentq(lambda x: elast(x) - cut, 1e-6, P.theta * 10, xtol=1e-9)
            m_c = level_at(th_c, P.Km, P.ke, L0)
            print(f"  elasticity falls below {cut:.2f} at theta = {th_c:8.1f} "
                  f"({100*th_c/P.theta:5.1f}% of nominal), where the marker reads "
                  f"{m_c/m_bar:.2f} x managed")
        except ValueError:
            print(f"  elasticity never falls below {cut:.2f} in range")
    print()
    print("  THE CLINICIAN-FACING STATEMENT. The rate estimand is informative about")
    print("  capacity only while the marker sits well below its ceiling. Once a")
    print("  patient's managed marker is a few times the reference level, a further")
    print("  fall in capacity barely moves it, and the estimate of capacity behind")
    print("  it is worthless -- while the RATE prediction itself stays fine, because")
    print("  the rate is read off the marker and the marker is what is measured.")
    print()
    print("  So the failure is asymmetric and must be stated that way: the rate")
    print("  survives where the capacity does not. That is an argument for the rate")
    print("  estimand and against reading theta out of it for severe patients.")


scope_of_the_rate_estimand()
