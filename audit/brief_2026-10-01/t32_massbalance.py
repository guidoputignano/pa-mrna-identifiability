"""
Item 4: reconcile the model's clearance split with the PA literature.

The model (Model A, src/common/models.py) is dm/dt = s - theta*g(m) - ke*m, with
ke*m the LINEAR, UNBOUNDED alternative (non-PCC) clearance. At 5% capacity the
enzyme's share s = theta g(m)/L ~ 0.16 (R-A), so ke*m absorbs ~84% of the load.

The literature (user) says measured urinary propionate-metabolite output in PA is
1.5-3 umol/kg/h (Thompson & Chalmers 1990), one to TWO orders of magnitude below
estimated production (~1-2 mmol/kg/day = 42-83 umol/kg/h). i.e. in real PA the
non-PCC routes clear only a few percent of the load, which is WHY patients
accumulate and decompensate. The model's 84% non-PCC clearance at low capacity
looks far too high.

Caveat (user, correct): urinary output is not all non-PCC disposal -- odd-chain
lipid incorporation, protein propionylation, the methylcitrate shunt and acyl-
carnitine/CoA buffering are not excreted (research doc 4). So the urinary 1.5-3
is a floor on non-PCC disposal, not the whole of it. But the research doc is
explicit that the remainder is "buffered transiently ... or accumulated", NOT
cleared at high capacity. The conflict is real at the level of CAPACITY.
"""
import numpy as np

import core as C
import n4_theta_scope as N4

KM, KE = C.KM0, C.KE


def partition(theta):
    """At capacity theta, convention A (same diet): L, m, enzymatic flux theta*g,
    alternative flux ke*m, and the non-PCC fraction."""
    L, _ = N4.anchor(theta, KM, KE, "A")
    # re-solve m at this theta with the SAME load L (convention A)
    from numpy import roots
    # L = theta*m/(Km+m) + ke*m  ->  ke*m^2 + (ke*Km + theta - L)*m - L*Km = 0
    a = KE
    b = KE * KM + theta - L
    c = -L * KM
    m = max((-b + np.sqrt(b * b - 4 * a * c)) / (2 * a), 1e-9)
    enz = theta * m / (KM + m)
    alt = KE * m
    return L, m, enz, alt, alt / (enz + alt)


def main():
    print("=" * 90)
    print("  ITEM 4. THE CLEARANCE SPLIT vs THE PA LITERATURE")
    print("=" * 90)
    print("\n  Model A disposal partition, convention A (same diet, L fixed at nominal):")
    print(f"    {'capacity':>22}  {'L':>7}  {'m':>8}  {'enzymatic':>10}  {'ke*m (non-PCC)':>15}  {'non-PCC %':>9}")
    for frac, lab in ((1.00, "100% (normal)"), (0.10, "10%"), (0.05, "5% (disease)"),
                      (0.01, "1%")):
        L, m, enz, alt, fr = partition(C.TH0 * frac)
        print(f"    {lab:>22}  {L:7.1f}  {m:8.1f}  {enz:10.1f}  {alt:15.1f}  {100*fr:8.1f}%")

    print("""
  READING: at NORMAL capacity the enzyme carries ~85% and the alternative route
  ~15%. As capacity falls the split INVERTS: at 5% the alternative route ke*m
  carries ~84% and rising, because Model A's ke*m is LINEAR and UNBOUNDED -- as m
  climbs, ke*m climbs without limit and silently disposes of the whole load. That
  is the structural reason a 5% patient still reaches a benign steady state in the
  model.

  THE CONFLICT. The PA literature says the opposite: non-PCC disposal is LOW-
  CAPACITY (urinary metabolite output 1.5-3 umol/kg/h, 1-2 orders below
  production), so a real PCC-deficient patient does NOT clear 84% by alternative
  routes -- they accumulate and decompensate. Model A's unbounded ke*m overstates
  residual disposal at low capacity, exactly where the disease lives.

  THE REPO ALREADY CONTAINS THE FIX. models.py Model B replaces ke*m with
  ke*h(m), h(m)=m/(Ke+m): a SATURATING alternative route with its own V_max = ke.
  Capping the alternative route is precisely what the literature's low-capacity
  non-PCC disposal implies. The FIM/identifiability work in this brief used Model
  A throughout (the linear, analytically tractable form); the clearance-structure
  question (A vs B vs C) is the repo's own PHASE5/hypothesis question, and the
  literature here comes down on the side of a SATURATING (or otherwise bounded)
  alternative route.
""")
    # What a bounded alternative route would change: ke*m at 5% in Model A
    L5, m5, enz5, alt5, _ = partition(C.TH_5)
    print(f"  Model A at 5%: ke*m = {alt5:.0f} (vs ke*m_bar at nominal = {KE*10:.0f}),"
          f" a {alt5/(KE*10):.1f}x rise in the alternative flux as m goes 10 -> {m5:.0f}.")
    print("  A saturating route (Model B) would cap that rise, so the marker would")
    print("  climb far higher (or the steady state would not exist) at 5% -- which is")
    print("  the clinically realistic picture of accumulation. The s=0.16 fraction")
    print("  oxidised is specific to Model A's unbounded clearance and would MOVE under")
    print("  a bounded route; much of the cross-sectional analysis is conditioned on it.")

    print("\n  PHYSICAL UNITS. The model is in abstract units (theta=3000 'V_max',")
    print("  L~323/time); the repository supplies NO bridge to umol/kg/h (it is")
    print("  calibrated to a mouse 3-HP marker baseline, not an absolute propionate")
    print("  flux). So the model's ke*m CANNOT be put in umol/kg/h to compare numerically")
    print("  with the 1.5-3 umol/kg/h urinary figure -- only the FRACTION (84% non-PCC)")
    print("  is comparable, and that is what conflicts. An absolute comparison needs a")
    print("  unit anchor the repo does not contain (flagged, not fabricated).")

    print("\n  PROVENANCE CONFIRMED: no human measurement of k_e exists in the literature")
    print("  (research doc 4: 'no direct human measurement of the total fraction cleared")
    print("  by non-PCC routes'; 'a precise human plasma propionate half-life ... is not")
    print("  well established'). This SUPPORTS the parent paper's claim that no published")
    print("  source supplies k_e, and is why k_e is carried with a wide prior here.")


if __name__ == "__main__":
    main()
