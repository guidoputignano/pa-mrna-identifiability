"""
Item 1: the kappa question for ex vivo PCC. And item 2: scoping the propionyl-CoA
model extension. Repo values only; gaps flagged, not filled.

EX VIVO PCC reads A = kappa * theta * exp(eps), eps ~ N(0, tau_map)
(PHASE10_SENSORS.md 1.2). The repo's 20% is "an ex vivo assay's REPEATABILITY"
(q15_sensors.py step 3) = the eps/tau_map scatter. kappa, the cell-to-whole-body
mapping, is a NUISANCE "neither of which any human study reports (R4)" -- i.e.
UNCALIBRATED. kappa is structurally the breath test's recovery factor c: a
multiplicative factor that a single reading cannot separate from theta.
"""
import numpy as np

import core as C


def sd_exvivo(Y, Z, nread=1, seeds=C.SEEDS):
    """sd(ln theta): record + ex vivo reading (repeatability Y, nread reads),
    kappa carried free with prior sd Z. Z<=0 means kappa KNOWN (column dropped)."""
    def one(s):
        r = C.record(C.TH_5, seed=s)
        F4 = C.fim(r, noise="perdraw")
        sd = Y / np.sqrt(nread)
        if Z <= 0:                      # kappa known: direct theta reading, 4-param
            row = np.array([1.0, 0, 0, 0])
            F = F4 + np.outer(row, row) / sd ** 2
        else:                           # kappa free + prior Z, 5-param
            F = np.zeros((5, 5)); F[:4, :4] = F4
            row = np.array([1, 0, 0, 0, 1.0])
            F += np.outer(row, row) / sd ** 2
            if np.isfinite(Z):
                F[4, 4] += 1.0 / Z ** 2
        w = np.linalg.eigvalsh(F)
        return np.inf if w[0] <= 1e-12 * w[-1] else float(np.sqrt(np.linalg.inv(F)[0, 0]))
    return C.over_seeds(one, seeds)[0]


def main():
    print("=" * 92)
    print("  ITEM 1. THE kappa QUESTION FOR EX VIVO PCC")
    print("=" * 92)
    print("\n  PROVENANCE: the 20% is REPEATABILITY only (q15 step 3: 'an ex vivo assay's")
    print("  repeatability'). kappa (cell-to-whole-body mapping) is unreported in any human")
    print("  study (PHASE10 1.2, R4), so it is UNCALIBRATED -- the same object as the breath")
    print("  test's recovery c. q15's 0.200 adds the 20% as a direct theta-prior, which")
    print("  silently assumes kappa = 1 known.")

    print("\n  sd(ln theta) at 5% monthly, over (repeatability Y, calibration Z = sd(ln kappa)):")
    Ys = [0.10, 0.15, 0.20, 0.25, 0.30]
    print("    Z = sd(ln kappa) \\ Y   " + "".join(f"{y*100:7.0f}%" for y in Ys))
    for Z in (0.0, 0.05, 0.10, 0.20, np.inf):
        lab = "exact (kappa known)" if Z <= 0 else ("free (uncalibrated)" if not np.isfinite(Z) else f"{Z:.2f}")
        print(f"    {lab:>21}   " + "".join(f"{sd_exvivo(y, Z):8.3f}" for y in Ys))

    print("\n  READING:")
    print("  - kappa KNOWN: sd(ln theta) ~ the repeatability Y (0.20 at 20%, 0.30 at 30%).")
    print("    So the 0.200 is knife-edge: it sits exactly on the gate at the stated 20%,")
    print("    and 25%/30% repeatability give 0.25/0.30 -- a result on the gate should never")
    print("    have been reported without this.")
    print("  - kappa UNCALIBRATED (free): sd(ln theta) = 3.045, the record's own value. A")
    print("    single ex vivo reading is then worth NOTHING for absolute theta, exactly as")
    print("    a single breath reading is (R11). The ex vivo route has the SAME unresolved")
    print("    calibration problem, and at face value the conclusion REVERTS.")
    print("  - The 0.20 contour is sqrt(Y^2 + Z^2) <= 0.20 (hypot, as in R11). Maximum kappa")
    print("    uncertainty compatible with target, by repeatability:")
    for y in Ys:
        z2 = 0.20 ** 2 - y ** 2
        print(f"      Y = {int(y*100)}%: sd(ln kappa) <= "
              f"{('%.3f' % np.sqrt(z2)) if z2 > 0 else 'NONE (repeatability alone >= 0.20)'}")
    print("    At the repo's 20% repeatability, kappa must be known essentially EXACTLY;")
    print("    any tolerance for kappa needs the repeatability pushed below ~15% first.")
    print("  - MONITORING escapes this: kappa is a patient constant, so it CANCELS in the")
    print("    longitudinal contrast (like k_e, b). Two ex vivo reads give a contrast of")
    print(f"    ~sqrt(2)*0.20 = {np.sqrt(2)*0.20:.2f} (worse than the marker's 0.13), so ex vivo")
    print("    does not improve monitoring, but the kappa problem is a CLASSIFICATION")
    print("    problem only -- the same split as everything else in this study.")

    print("\n" + "=" * 92)
    print("  ITEM 2. SCOPING THE PROPIONYL-CoA MODEL EXTENSION (PHASE10 1.1 route b)")
    print("=" * 92)
    print("""
  THE REPO DOES NOT PRICE C3. C3 is dropped at step 1 (inexpressible: no
  propionyl-CoA state). So 'compute sd(ln theta) with C3' cannot be done on repo
  values -- it requires building route (b), which the repo states it 'cannot
  itself run'. Scoping it, repo-grounded where possible and flagged where not:

  NEW STRUCTURE (PHASE10 1.1b): a propionyl-CoA pool P, the enzyme acting on P,
  the marker m formed from P 'by overflow'. New parameters beyond (theta, K_m,
  k_e, L, b):
    - k_P : non-enzymatic disposal of P              (NEW, not in any repo file)
    - the overflow kinetics P -> m (form + rate)     (NEW, functional form unstated)
    - kappa_c : the C3 assay's mapping constant        (NEW, a THIRD uncalibrated
                multiplier, exactly like kappa and c)

  IDENTIFIABILITY OF THE NEW PARAMETERS FROM EXISTING DATA: NO. The marker record
  observes only m, which is downstream of P; P and its disposal split (theta g(P)
  vs k_P P vs overflow) are unobserved. Adding k_P and the overflow kinetics as
  free parameters, with only m observed, adds unidentified parameters to an
  already ill-posed model -- it can only make theta worse, not better. So the
  extension is worth running ONLY together with a P-observable (C3).

  WOULD C3 THEN BREAK THE DEGENERACY, OR MOVE IT? Structural answer:
    - C3 ∝ kappa_c * P reads the pool directly. With P pinned, the enzyme flux
      theta g(P) is read more directly than through the overflow marker -- that is
      the genuine promise, and why the enzyme's SUBSTRATE pool is the right place
      to look.
    - BUT two obstructions reappear on P: (i) the theta-vs-k_P confound on the P
      balance (L = theta g(P) + k_P P + overflow) is the SAME structure as the
      theta-vs-k_e confound on m, so it breaks theta out only if non-enzymatic P
      disposal k_P is small or known; (ii) kappa_c is a new uncalibrated
      multiplier, so an absolute C3 reading inherits the kappa/c calibration
      problem above.
    - VERDICT: the extension MOVES the degeneracy from (theta vs k_e, on m) to
      (theta vs k_P, on P) and adds a third calibration unknown. It BREAKS the
      degeneracy only if (a) non-enzymatic propionyl-CoA disposal is negligible or
      measured and (b) C3's mapping kappa_c is calibrated -- neither of which the
      repository provides. A definitive number needs the structural model built
      with those priors; that is out of the current model and is the repo's own
      stated main open item (PHASE10 'route (b) is what this phase should
      recommend and cannot itself run'). Not fabricated here.
""")


def calibration():
    print("\n" + "=" * 92)
    print("  ITEM 1b. CALIBRATING kappa (a POPULATION constant) -- the shortest path")
    print("=" * 92)
    print("""
  kappa is a population constant (cell-to-whole-body mapping), so it is pinned
  ONCE, not per subject. Calibration subgroup: N patients whose theta is
  independently known to sd(ln theta) = sig_th, each giving an ex vivo reading
  ln A = ln kappa + ln theta + eps (repeatability Y_cal, n_read reads). Each
  patient estimates ln kappa with variance Y_cal^2/n_read + sig_th^2; over N,
      sd(ln kappa) = sqrt(Y_cal^2/n_read + sig_th^2) / sqrt(N).
  Target for the 10%-repeatability deployment: sd(ln kappa) <= 0.173 (so that
  hypot(0.10, 0.173) = 0.20).

  COMPLETE-DEFICIENCY (theta ~ 0) DOES NOT SERVE. At theta = 0 the reading
  A = kappa*theta = 0 for every kappa: zero information. And 'theta near 0' with
  an unknown value leaves ln theta unbounded, so kappa is not separable. The
  anchor must be patients with KNOWN, NONZERO theta.""")
    Z = 0.173
    print("  N patients needed (Y_cal = 10%, 1 read each):")
    for sig in (0.10, 0.20, 0.30, 0.50):
        per = 0.10 ** 2 + sig ** 2
        print(f"    theta independently known to {int(sig*100):>2}%:  N = {int(np.ceil(per/Z**2))}")
    print("""
  DELIVERABLE: kappa calibrated to sd(ln kappa) <= 0.173 from N ~ 2 patients with
  theta independently known to 20% (or ~9 at 50%) makes the ex vivo route reach
  0.20 at <=10% repeatability -- PROVIDED a source of independently-known nonzero
  theta exists. The repository identifies none (ex vivo itself is circular;
  complete deficiency fails), so that source is the open requirement, not the N.

  DOES THE SAME SERVE THE BREATH TEST's c? NO. u3 carries c as a SUBJECT-SPECIFIC
  recovery factor, so a population calibration study does not transfer -- each
  patient would need their own c. Only kappa (population) is fixable once. That is
  the one structural advantage ex vivo has over the breath test.""")


if __name__ == "__main__":
    main()
    calibration()
