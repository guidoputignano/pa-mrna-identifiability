"""
Round 15: lambda is no longer ASSUMED. Two literature searches supply a measured
HUMAN value and a mouse-fitted one, so the repo's provenance ("rodent PK scaled
allometrically") is wrong twice over -- a human measurement exists, and allometry
is the wrong operation for an intensive tissue-level rate -- and the decision must
be re-run against the real numbers and their uncertainty band.

Literature (see REPORT.md R-N for full provenance):
  - Chandler & Ballard 1985, HUMAN fibroblast PCC, k_deg = 0.19/d  -> t_half ~ 3.6 d
    (single-exponential [3H]biotin pulse-chase; the ONLY human PCC turnover number).
  - Attarwala 2023, MOUSE-fitted PK/PD: k_deg = 0.0075/h -> ln2/k_deg = 3.8 d intrinsic;
    a 2-compartment protein layer (k_q = 0.00474/h) gives an effective terminal ~7 d.
  - Jiang 2020, mouse-liver exogenous human PCC: 5.3-6.7 d (corroboration).
  - Recommended sensitivity: 3-fold base case (t_half 1.2-11.5 d), 10-fold worst case.

Round-11/13 thresholds (half-life at which a headline flips) are recomputed here and
the new values placed against them. Same-rate dosing throughout (the repo's sec-7
convention): Delta_trough(tau) = D14*(tau/14)*e^{-lam*tau}/(1-e^{-lam*tau}).
"""
import sys
sys.path.insert(0, "../../src")
import numpy as np
from scipy.optimize import brentq, least_squares

import model_C.run as MC
from common.models import ParamsC
from model_A.run import ParamsA, baseline_load, M_TOL
import t46_timing as T46

D14, TAU0 = 6000.0, 14.0
KM, KE = 100.0, 5.0
THC_INT = MC.theta_crit_intrinsic(ParamsC(), MC.S_BASE)           # 150.8 (Model C, no m_tol)
THC_ROUT = (baseline_load(ParamsA()) - KE * M_TOL) * (KM + M_TOL) / M_TOL   # 429.5 (fold-1.0)
REPO_HL = 3.5                                                      # ln2/LAMBDA, the old assumption


def trough(tau, thalf):
    lam = np.log(2.0) / thalf
    x = np.exp(-lam * tau)
    return D14 * (tau / TAU0) * x / (1.0 - x)


def main():
    print("=" * 92)
    print("  ROUND 15. lambda RE-GROUNDED ON LITERATURE: a measured human 3.6 d, a mouse 7 d")
    print("=" * 92)

    # ---------- ITEM 1: which lambda does the trough see? ----------
    print("""
  ITEM 1. WHICH lambda DOES THE TROUGH SEE? The human value (Chandler 1985) is a single
  exponential (3.6 d) -- there is no two-compartment subtlety in it. The mouse Attarwala
  model is 2-compartment; the trough sits at 14 d, deep in terminal decline, so it sees the
  SLOW (terminal) rate, not the intrinsic k_deg. The exact 2-comp equations are not public
  (Suppl Table S1 unavailable), so two bracketing symmetric reconstructions are shown; the
  paper's stated effective terminal is ~7 d.""")
    kdeg, kq = 0.0075, 0.00474     # /h
    print(f"\n      {'2-comp reconstruction':>34}  {'terminal t_half':>15}")
    for lab, A in [("elimination from central only", np.array([[-(kdeg+kq), kq], [kq, -kq]])),
                   ("elimination from both comps", np.array([[-(kdeg+kq), kq], [kq, -(kdeg+kq)]]))]:
        w = np.linalg.eigvals(A); slow = -np.max(w)
        print(f"      {lab:>34}  {np.log(2)/slow/24:>13.1f} d")
    print(f"      {'paper stated effective terminal':>34}  {'~7.0 d':>15}")
    print(f"      {'intrinsic ln2/k_deg (1-comp view)':>34}  {'3.8 d':>15}")
    # effective single-exp half-life that reproduces the STEADY-STATE trough accumulation
    # (the decision-relevant match): single-exp trough = D*x/(1-x), x=e^{-lam*tau}; 2-comp
    # steady-state trough = sum_{n>=1} P(n*tau), P = central-comp impulse response.
    from scipy.linalg import expm
    def eff_halflife(A):
        imp = lambda t: (expm(A*t) @ np.array([1.0, 0.0]))[0]
        ss = sum(imp(n*14*24) for n in range(1, 60))     # steady-state trough (shape)
        # solve x/(1-x) = ss for x, then lam = -ln(x)/tau
        x = ss / (1.0 + ss)
        lam = -np.log(x) / (14*24)
        return np.log(2) / lam / 24
    for lab, A in [("central-only", np.array([[-(kdeg+kq), kq], [kq, -kq]])),
                   ("both-comp", np.array([[-(kdeg+kq), kq], [kq, -(kdeg+kq)]]))]:
        print(f"      effective single-exp t_half matching the STEADY-STATE trough ({lab} recon): "
              f"{eff_halflife(A):.1f} d")
    print("""      READING, and a correction to the "trough sees ~7 d" intuition. The effective single
      exponential that reproduces the mouse model's STEADY-STATE trough is ~3-5 d under both
      reconstructions -- BELOW the paper's 7-d terminal, not at it. The reason: the 14-day
      interval samples the decay before terminal dominance sets in, so the accumulated trough
      falls faster than the terminal rate alone implies. Mapping the mouse model into the
      repo's single-exponential increment via the decision-relevant quantity (the trough)
      therefore gives ~3-5 d, which CONVERGES with the human single-exp value (3.6 d) rather
      than the 7-d terminal. So reporting the decision "at 7 d" (below) uses the naive terminal
      mapping and OVERSTATES the trough; the proper mapping of either source lands at ~3.5-5 d
      -- right at the thresholds. The exact 2-comp value is structure-dependent (equations
      unpublished), so the three stated points (3.6 / 3.8 / 7 d) are carried forward with 7 d
      flagged as the optimistic end.""")

    # ---------- ITEM 2: decision at the candidate half-lives, and the fail half-lives ----------
    print("\n" + "=" * 92)
    print("  ITEM 2. THE DECISION AT 3.6 / 3.8 / 7 d, AND WHERE EACH HEADLINE FAILS")
    print("=" * 92)
    print(f"\n  theta_crit^int = {THC_INT:.1f} (Model C, convention-free); theta_crit^rout = {THC_ROUT:.1f} (fold 1.0).")
    print(f"  Delta_trough and Model-C safe-for-all status (Delta_trough(tau) > theta_crit^int):\n")
    print(f"      {'t_half':>8}  {'Dtr14':>7}  {'Dtr15':>7}  {'Dtr21':>7}   {'14d':>4} {'15d':>4} {'21d':>4}  note")
    for th, note in [(REPO_HL, "repo's old ASSUMED value"), (3.6, "human fibroblast (primary)"),
                     (3.8, "ln2/k_deg (mouse intrinsic)"), (7.0, "mouse effective terminal")]:
        t14, t15, t21 = trough(14, th), trough(15, th), trough(21, th)
        mk = lambda v: " ok " if v > THC_INT else "FAIL"
        print(f"      {th:>7.1f}d  {t14:>7.0f}  {t15:>7.0f}  {t21:>7.0f}   {mk(t14)} {mk(t15)} {mk(t21)}  {note}")

    print("\n  Half-life at which each headline FAILS (flips), and the margin at the human 3.6 d:")
    # each headline with a bracket that avoids the coarse-grid singularity at trough14=THC_ROUT (~3.57 d)
    heads = [("Model C: 14-day safe for all", lambda h: trough(14, h) - THC_INT, (0.5, 49.0)),
             ("Model C: 15-day safe for all", lambda h: trough(15, h) - THC_INT, (0.5, 49.0)),
             ("Model C: 21-day safe for all (floor theta=0)", lambda h: trough(21, h) - THC_INT, (0.5, 49.0)),
             ("coarse grid 14-vs-10 resolvable at record 3.05 (fold 1.0)",
              lambda h: (trough(10, h) - trough(14, h)) / (THC_ROUT - trough(14, h)) - 3.05, (2.5, 3.55))]
    print(f"      {'headline':>56}  {'fails below':>11}  {'margin at 3.6 d':>15}")
    for lab, f, br in heads:
        try:
            hl = brentq(f, *br)
        except Exception:
            hl = float("nan")
        marg = 100 * (3.6 - hl) / hl if np.isfinite(hl) else float("nan")
        print(f"      {lab:>56}  {hl:>9.2f} d  {marg:>+13.0f}%")
    print(f"""
  READING: the correction runs the WRONG WAY. The repo's old 3.5 d sat just BELOW the
  21-day floor (3.55 d), so round 13 found 21-day needs theta (and the record cannot supply
  it); the human-measured 3.6 d sits just ABOVE it, flipping 21-day to "free for all" -- but
  by only ~1%, a knife-edge. The Model-C 14-/15-day safe-for-all headlines clear 3.6 d by a
  comfortable ~30-37%; the coarse-grid resolvability clears by only ~8%. And human HEPATIC
  PCC very likely turns over FASTER than fibroblast (liver is the fastest-turnover tissue,
  Rolfs 2021), so the true hepatic half-life may sit at or below 3.6 d -- toward every
  threshold, not away. The honest statement is a margin, not a pass: at 3.6 d the 14-/15-day
  steps hold, the coarse grid barely holds, and the 21-day step is on its boundary.""")

    # ---------- ITEM 3: sensitivity over the 3-fold band ----------
    print("\n" + "=" * 92)
    print("  ITEM 3. SENSITIVITY OVER THE 3-FOLD BAND (t_half 1.2-11.5 d; 10-fold worst case)")
    print("=" * 92)
    print(f"\n      {'t_half':>8}  {'Dtr(14)':>8}   {'14d':>4} {'15d':>4} {'21d':>4} {'coarse':>6}")
    for th in (1.2, 2.0, 2.62, 2.75, 3.32, 3.55, 3.6, 5.0, 7.0, 11.5):
        t14 = trough(14, th)
        c14 = "ok" if trough(14, th) > THC_INT else "FAIL"
        c15 = "ok" if trough(15, th) > THC_INT else "FAIL"
        c21 = "ok" if trough(21, th) > THC_INT else "FAIL"
        req = (trough(10, th) - trough(14, th)) / (THC_ROUT - trough(14, th)) if THC_ROUT - trough(14, th) > 0 else 1e9
        cc = "ok" if req >= 3.05 else "FAIL"
        print(f"      {th:>7.1f}d  {t14:>8.0f}   {c14:>4} {c15:>4} {c21:>4} {cc:>6}")
    print("""
  READING: EVERY headline flips inside the 3-fold base-case band (1.2-11.5 d): 14-day at
  2.62, 15-day at 2.75, coarse-grid at 3.32, 21-day at 3.55 d. The within-species-and-method
  spread alone is ~3-fold (2.7-8.5 d), so this is measurement spread, not cross-species
  conservatism. Per the source's own rule, a conclusion that flips inside the 3-fold band
  makes lambda the DOMINANT model uncertainty -- which it now is for every interval headline.
  The 15-day step (round 13's "free for all") is NOT exempt: it fails below 2.75 d, inside
  the band, so even it is lambda-conditional, though it holds at all three literature points.""")

    # ---------- ITEM 5: the record discriminates the hypotheses (a test, not self-calibration) ----------
    print("\n" + "=" * 92)
    print("  ITEM 5. THE RECORD'S lambda IS NOW A TEST: can it tell 3.6 d from 7 d?")
    print("=" * 92)
    print(f"\n  Hypotheses are ln-separated by: 3.6 vs 7 d -> {np.log(7/3.6):.3f}; 3.6 vs 3.8 d -> {np.log(3.8/3.6):.3f}.")
    print(f"  Separation in sigma (d' = |ln ratio| / sd(ln lambda)) at the on-therapy record's precision:\n")
    h1, h2, h3, h4 = "record", "sd(ln lam)", "d'(3.6 vs 7)", "d'(3.6 vs 3.8)"
    print(f"      {h1:>16}  {h2:>11}  {h3:>12}  {h4:>14}")
    for lab, mo, cad in [("6 mo monthly", 6, 30.4), ("12 mo monthly", 12, 30.4),
                         ("12 mo weekly", 12, 7.0), ("24 mo weekly", 24, 7.0)]:
        sl = T46.sd_lnlam(mo, cad)
        print(f"      {lab:>16}  {sl:>11.3f}  {np.log(7/3.6)/sl:>12.1f}  {np.log(3.8/3.6)/sl:>14.1f}")
    print("""
  READING: with the human (3.6 d) and mouse-terminal (7 d) values nearly 2x apart and the
  on-therapy record pinning lambda to ~5-9% within the first year (round 12), the record
  separates the two hypotheses at ~7-15 sigma -- it can decisively tell which regime a real
  patient is in. So round 12's self-calibration is better stated as a PREDICTION the trial's
  own marker record could settle: measure lambda from the dose-timing response and it will
  say whether human hepatic PCC turns over near 3.6 d (headlines on the edge) or near 7 d
  (headlines comfortable). The 3.6-vs-3.8 d distinction is NOT resolvable (<1 sigma), and
  need not be -- both sit on the same side of the knife-edge pattern. The decision does not
  need lambda assumed; it needs lambda READ, and the record reads it.""")


if __name__ == "__main__":
    main()
