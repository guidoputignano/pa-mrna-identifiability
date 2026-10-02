"""
Item 1, redone: the ratio route is right about kappa and WRONG about the estimand.

t28 showed that if the estimand is an ex vivo activity RATIO (A_patient / mean
A_healthy), the population mapping kappa cancels. True. But the repository's
theta is NOT in vitro activity: src/common/models.py defines theta as "V_max of
PCC, the maximal enzymatic velocity" in the mass balance dm/dt = s - theta*g(m)
- ke*m. That is WHOLE-BODY enzymatic capacity. The ex vivo assay measures IN
VITRO residual activity. The literature (provided by the user) says these two
are only weakly related:

  - Lehnert/Suormala/Baumgartner 1994: residual PCC activity is "neither a
    predictive marker for severity nor for outcome".
  - The ONLY same-subject in vitro / in vivo study is Manoli 2021 -- in MMA, not
    PA, using 14C-propionate incorporation (not a direct PCC assay): Spearman
    r = 0.682, R^2 = 0.327, n = 11. A correlation, not a calibration.
  - Field view: the relationship is "difficult or impossible to reduce to a
    single scaling factor" (tissue-specific expression, hepatic vs fibroblast
    flux, alternative disposal).

CONSEQUENCE. Write ln(in-vitro activity) = a + b*ln(theta_invivo) + eta, with
eta the PER-PATIENT departure of the subject's in-vitro reading from what their
whole-body capacity predicts. Reading the in-vitro ratio to estimate the in-vivo
ratio leaves a residual sd(ln theta_invivo | in-vitro) = sqrt(1 - R^2) *
sd(ln theta_between). This eta is INDEPENDENT across patients, so -- unlike kappa
-- it does NOT cancel in a ratio. My item-3 (a ~15-20% systematic kappa offset)
understated the problem: the real obstruction is per-patient SCATTER.

ALL FIGURES BELOW ARE LITERATURE, FROM DIFFERENT POPULATIONS/ASSAYS THAN THE
REPO MODEL (MMA not PA; 14C-incorporation not a direct PCC assay). They are used
to TEST the ratio route's assumption and are swept, not substituted for any repo
value. The repo itself contains no in-vitro/in-vivo mapping.
"""
import numpy as np

import core as C


def sd_ratio(Y, sig_healthy, N, sig_eta, nread=1, seeds=C.SEEDS):
    """sd(ln theta_invivo) from the monthly record + an ex vivo ratio reading.
    Reading precision = hypot(patient repeatability Y/sqrt(nread),
    healthy-mean SE sig_healthy/sqrt(N), mapping scatter sig_eta). kappa cancels
    (t28); sig_eta does NOT (per-patient, independent)."""
    P = np.hypot(np.hypot(Y / np.sqrt(nread), sig_healthy / np.sqrt(N)), sig_eta)

    def one(s):
        r = C.record(C.TH_5, seed=s)
        row = np.array([1.0, 0, 0, 0])
        return C.sd_th(C.fim(r, noise="perdraw") + np.outer(row, row) / P ** 2)
    return C.over_seeds(one, seeds)[0]


def main():
    print("=" * 94)
    print("  ITEM 1 REDONE. TWO ESTIMANDS: in-vitro ACTIVITY vs WHOLE-BODY CAPACITY")
    print("=" * 94)
    print("""
  The repo's theta is WHOLE-BODY V_max (models.py: 'V_max of PCC, the maximal
  enzymatic velocity'; mass balance dm/dt = s - theta*g(m) - ke*m). The ex vivo
  assay reads IN VITRO activity. The ratio route of t28 cancels kappa but NOT the
  per-patient in-vitro/in-vivo mapping scatter eta.
""")
    # Between-patient spread of ln(in-vivo capacity), from Thompson 1990 in-vivo
    # oxidation (user: mean 36.3%, range 3.0-91.1%). ln(91.1/3.0)=3.41; as +/-2SD
    # that is sd~0.85, so 0.6-1.0 brackets it. SWEPT, flagged, not a repo value.
    spreads = (0.6, 0.8, 1.0)
    # R^2 of the in-vitro -> in-vivo mapping. Manoli 2021 (MMA, 14C): R^2=0.327.
    # 0.5 and 0.8 are optimistic what-ifs (a better, PA-specific, direct-PCC map).
    r2s = (0.327, 0.50, 0.80)

    print("  (0) BASELINE -- estimand = in-vitro activity (t28, sig_eta = 0):")
    print(f"      patient read 10%, 30 controls, spread 50%: sd = {sd_ratio(0.10,0.50,30,0.0):.3f}"
          "  <= GATE (t28 stands for THIS estimand)")

    print("\n  (1) estimand = WHOLE-BODY CAPACITY: the mapping scatter sig_eta enters.")
    print("      sig_eta = sqrt(1-R^2) * sd(ln theta_between); it does NOT cancel.\n")
    print("      patient read 10%, N=30 controls, healthy ex-vivo spread 50%:")
    h1 = "R^2 \\ between-patient sd"
    print(f"      {h1:>26}  " + "".join(f"sd={s:>4}" + "  " for s in spreads))
    for r2 in r2s:
        cells = []
        for sd_bw in spreads:
            eta = np.sqrt(1 - r2) * sd_bw
            v = sd_ratio(0.10, 0.50, 30, eta)
            cells.append(f"{v:7.3f}" + ("* " if v <= 0.20 else "  "))
        print(f"      R^2={r2:>5}  (sig_eta={np.sqrt(1-r2)*0.8:.2f} at sd=0.8)   " + "".join(cells))

    print("""
  READING: at the only measured mapping (Manoli R^2=0.327, MMA) and a between-
  patient spread of ~0.8 log units (Thompson in-vivo oxidation range), sig_eta
  ~= 0.66 -- about FIVE TIMES the 0.135 t28 reported, and far above the 0.20
  gate. Even an optimistic R^2=0.8 with a tight 0.6 spread gives sig_eta ~= 0.27,
  still over the gate. For WHOLE-BODY CAPACITY the ratio route does NOT reach
  target: the floor is set by how well an in-vitro assay predicts in-vivo flux,
  which no study has shown to be better than ~1/3 of variance, and none at all
  in PA.
""")

    print("  (2) DOES A TIGHTER HEALTHY CONTROL SPREAD RESCUE IT? (the other axis)")
    print("      Gravel 1977 fibroblast PCC normal 450-850 pmol/min/mg -> ~16% if read")
    print("      as +/-2SD (may be min-max; treat as a BOUND). Sweeping it against sig_eta:")
    eta_fixed = np.sqrt(1 - 0.327) * 0.8   # Manoli R^2, spread 0.8
    print(f"      (mapping scatter held at sig_eta = {eta_fixed:.2f}: Manoli R^2, spread 0.8)\n")
    h2 = "healthy spread \\ N"
    print(f"      {h2:>22}  {'N=10':>8}{'N=30':>8}{'N=100':>8}")
    for sig in (0.16, 0.30, 0.50):
        cells = "".join(f"{sd_ratio(0.10, sig, N, eta_fixed):8.3f}" for N in (10, 30, 100))
        print(f"      spread {int(sig*100):>3}%           {cells}")
    print("""
      A tighter control spread (0.16) barely moves the answer: once sig_eta
      dominates, neither more controls nor better controls help. The control
      panel was never the binding term for the whole-body estimand -- the mapping
      is.
""")

    print("=" * 94)
    print("  THE DECIDING QUESTION: which estimand is the clinical one?")
    print("=" * 94)
    print("""
  - If the question is IN VITRO RESIDUAL ENZYME ACTIVITY (% of control), the
    ratio route answers it: kappa cancels, sd ~= 0.11-0.19, t28/v9 stands.
  - If the question is WHOLE-BODY CAPACITY (the parent model's theta, a flux
    parameter in a mass balance), the ratio route INHERITS the weak in-vitro/
    in-vivo mapping: floor ~= 0.5-0.8, it does NOT reach 0.20.

  The repo's theta is whole-body V_max. So for the estimand the model actually
  carries, the ex vivo ratio route does not reach target, and the cross-sectional
  conclusion reverts closer to the original 'no route reaches 0.20 by inference'
  -- the one exception being a DIRECT whole-body flux measure (the 1-13C breath
  oxidation test), which reads in-vivo capacity itself rather than an in-vitro
  proxy. The scatter here is a LITERATURE estimate from MMA; the honest statement
  is that the mapping is unmeasured in PA and assumed tight by the ratio route,
  and the only same-disease-adjacent datum puts it far too loose for 0.20.
""")

    print("=" * 94)
    print("  ITEM 2. THE FLOOR IS BORROWED, NOT MEASURED -- and the study that would settle it")
    print("=" * 94)
    print("""
  The 0.5-0.8 floor rests entirely on Manoli R^2 = 0.327 -- MMA, 11 patients,
  14C-incorporation, not a direct PCC assay. There is NO PA/PCC same-subject
  in-vitro/in-vivo study at all. So the honest conclusion is NOT 'the floor is
  0.66'; it is: THERE IS NO EVIDENCE THAT AN IN-VITRO ASSAY REACHES TARGET FOR
  WHOLE-BODY CAPACITY, AND THE STUDY THAT WOULD SETTLE IT HAS NEVER BEEN DONE IN
  PA. The 0.66 is an indicative calculation from the nearest datum.

  WHAT R^2 THE RATIO ROUTE WOULD NEED. To bring the mapping scatter sig_eta =
  sqrt(1-R^2)*sd_between under the 0.20 gate, the mapping must explain""")
    for sd_bw in (0.6, 0.8, 1.0):
        r2_need = 1 - (0.20 / sd_bw) ** 2
        print(f"      between-patient sd {sd_bw}: R^2 >= {r2_need:.3f}  "
              f"(Manoli's 0.327 is far below)")
    print("""
  i.e. R^2 ~= 0.89-0.96 -- the in-vitro assay would have to predict in-vivo
  capacity almost perfectly, versus the ~1/3 of variance actually measured.

  THE STUDY: paired DIRECT ex vivo PCC activity AND 1-13C-propionate oxidation
  in the SAME PA patients (what Manoli did in MMA, with a direct PCC assay added).
  N to decide, from the Fisher-z SE of a correlation (1/sqrt(N-3)):""")
    for N in (12, 20, 30, 60):
        se_z = 1 / np.sqrt(N - 3)
        # if true r ~ 0.57 (R^2=0.33), upper 95% bound on R^2:
        import math
        z = math.atanh(0.57)
        r_hi = math.tanh(z + 1.96 * se_z)
        print(f"      N={N:>2}: if true R^2~0.33, 95% upper bound on R^2 ~ {r_hi**2:.2f} "
              f"-> {'REJECTS the route' if r_hi**2 < 0.89 else 'inconclusive'}")
    print("""
  So ~20-30 PA patients with both measurements would REJECT the ratio route if
  the PA mapping is as loose as MMA's; CONFIRMING it (a tight CI around R^2 >=
  0.9) would take more, but the point estimate alone would almost certainly
  settle it. This is the single experiment that decides whether the ex vivo
  route is viable, and it is small.
""")


if __name__ == "__main__":
    main()
