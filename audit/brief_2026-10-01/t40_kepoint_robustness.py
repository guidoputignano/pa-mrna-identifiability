"""
Round 8 close: is the Model-B 0.19 a finding or a single constructed point?

(1) Redo the ceiling x managed-load sweep under the conditions that produced 0.19
    (tightened k_e prior, K_e known), and report sd against the steady-state margin
    mu = theta - theta_crit, so the ~1%-margin chosen point is not the only evidence.
(2) K_e tolerance: sd(ln theta) as a function of sd(ln K_e), with the provenance of
    K_e stated (it is a bare default, unmeasured in humans; the exchange degeneracy
    is exactly (theta,K_m)<->(k_e,K_e), so "K_e known" collapses it by fiat).
"""
import numpy as np

import core as C
import t33_model_b as T
import t35_managed_load as M

TH5, TH0 = C.TH_5, C.TH0
SEEDS = range(24)


def sd_B(L, ke, ke_prior, Ke_sd=None):
    """sd(ln theta) at 5%, managed load L, ceiling ke. Ke_sd None -> K_e pinned
    (known); else K_e free with prior sd(ln K_e) = Ke_sd."""
    def one(s):
        r = M.record_B_managed(TH5, L, ke, C.MONTHLY, s)
        if Ke_sd is None:
            return T.sd_th_B(T.fim_B(r, ke_prior=ke_prior, drop_Ke=True))
        return T.sd_th_B(T.fim_B(r, ke_prior=ke_prior, ke_m_prior=Ke_sd, drop_Ke=False))
    return float(np.median([one(s) for s in SEEDS]))


def main():
    print("=" * 92)
    print("  (1) THE 0.19 CONDITIONS ACROSS THE CEILING x LOAD GRID, AND vs THE MARGIN")
    print("=" * 92)
    for kep, lab in ((0.2, "k_e prior 0.2 (the value required)"),
                     (0.054, "k_e prior 0.054 (the IV tracer's precision)")):
        print(f"\n  {lab}, K_e known.  cell = sd(ln theta)  [margin mu, %% of normal]")
        header = "      " + f"{'ceiling':>8}"
        for L in (180, 200, 220, 240):
            header += f"{'L=' + str(L):>17}"
        print(header)
        allv = []
        for ke in (50.0, 100.0, 150.0, 200.0):
            row = f"      {ke:>8.0f}"
            for L in (180, 200, 220, 240):
                mu = TH5 - (L - ke)
                if mu <= 0:
                    row += f"{'(no SS)':>17}"
                    continue
                v = sd_B(L, ke, kep)
                allv.append((v, mu, ke))
                flag = "*" if v <= 0.20 else " "
                row += f"{v:>9.3f}{flag}[{100 * mu / TH0:4.1f}%]"
            print(row)
        sds = [v for v, _, _ in allv]
        npass = sum(v <= 0.20 for v in sds)
        print(f"      RANGE {min(sds):.3f} to {max(sds):.3f};  "
              f"{npass}/{len(sds)} valid cells reach 0.20")
        bymu = sorted(allv, key=lambda x: x[1])
        print("      sd vs margin (ascending mu):  " +
              ", ".join(f"{100*mu/TH0:.1f}%->{v:.3f}" for v, mu, _ in bymu))

    print("""
  READING (1): the 0.19 is a KNIFE-EDGE AT THE GATE, and how robust it is to the
  ceiling depends entirely on the (unmeasured) k_e prior:
    - at the tracer prior (0.054), the whole ceiling sweep 50-200 sits at the gate
      (0.18-0.20, 10/13 cells <= 0.20); the result is marginal everywhere, not a
      clear pass and not a clear miss.
    - at the looser 0.2 prior it degrades off the gate: only the low ceilings
      (k_e <= 100) pass, k_e = 150-200 miss, though only mildly (up to 0.228) --
      NOT into any 0.3+ range (the whole item-1 grid spans 0.18-0.23).
  Within a ceiling, sd is nearly flat in the managed load, so it does NOT track the
  margin mu -- it tracks the enzyme SHARE, which the ceiling sets; and the cells
  nearest the gate are the low-ceiling = smallest-mu = most-fragile ones. So "0.19"
  is a marginal, gate-hugging result conditional on BOTH a tight k_e prior and a
  low ceiling (and, from item 2, K_e known) -- a point held up by three unmeasured
  assumptions, not a finding that stands on its own.
""")

    print("=" * 92)
    print("  (2) K_e TOLERANCE AND PROVENANCE  (k_e=100, L=220, k_e prior = tracer 0.054)")
    print("=" * 92)
    print("""
  PROVENANCE: K_e is the half-saturation of the alternative (non-PCC) route. In
  the repo it is a bare default (ParamsB.Ke = 300 in src/common/models.py) with NO
  entry in the parameter-provenance registry (src/common/params.py) and NO human
  measurement behind it -- the same gap as k_e itself (round-6: no human non-PCC
  disposal constant exists). And the paper's exchange degeneracy is EXACTLY
  (theta, K_m) <-> (k_e, K_e): the pair is unidentified without a selective,
  time-varying input (A1_supplementary). So "K_e known" does not tighten a
  nuisance -- it collapses the degeneracy by fiat, assuming the very constant whose
  unknownness is the degeneracy.
""")
    print("  sd(ln theta) as a function of sd(ln K_e):")
    print(f"      {'sd(ln K_e)':>14}  {'sd(ln theta)':>14}")
    for Ke_sd in (1e-4, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, np.inf):
        lab = "known (pinned)" if Ke_sd <= 1e-4 else ("free (no prior)" if not np.isfinite(Ke_sd) else f"{Ke_sd:.2f}")
        v = sd_B(220, 100.0, 0.054, Ke_sd=(None if Ke_sd <= 1e-4 else Ke_sd))
        flag = "  <= 0.20" if v <= 0.20 else ""
        print(f"      {lab:>14}  {v:>14.3f}{flag}")

    # tolerance: largest sd(ln K_e) keeping sd(ln theta) <= 0.20
    grid = np.linspace(0.0, 0.6, 61)
    tol = None
    prev = None
    for z in grid:
        v = sd_B(220, 100.0, 0.054, Ke_sd=(None if z < 1e-6 else z))
        if v > 0.20 and prev is not None:
            tol = 0.5 * (z + prev)
            break
        prev = z
    print(f"""
  READING (2): sd(ln theta) crosses 0.20 at sd(ln K_e) ~ {('%.2f' % tol) if tol else '<0.05'}.
  So the gate is reached only while K_e is known to within ~{('%.0f%%' % (100*tol)) if tol else 'a few %'};
  beyond that the exchange degeneracy reopens and sd climbs back toward the
  record-only value. Since no human K_e exists, this tolerance is not met by any
  data the repo or the literature supplies -- the 0.19 rests on assuming it.

  THE PATTERN (belongs in the decision): every route that reaches target
  cross-sectionally does so CONDITIONAL ON AN UNMEASURED MULTIPLICATIVE/SHAPE
  CONSTANT -- kappa (ex vivo mapping), c (breath recovery), and now K_e (the
  alternative route's half-saturation). Three different routes, three different
  unmeasured constants, each of which the repository and the human literature
  leave unfixed. That recurrence is itself the finding: the cross-sectional
  target is reachable only by importing a constant nobody has measured.
""")


if __name__ == "__main__":
    main()
