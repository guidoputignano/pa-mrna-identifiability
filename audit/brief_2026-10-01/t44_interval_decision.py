"""
Round 11: the DECISION the paper was written for -- the dosing INTERVAL -- is a
trough condition on capacity, not a measurement of theta.

The therapeutic increment Delta(t) adds to capacity in (theta+Delta) g(m). On a
train of decaying doses (rate lambda, interval tau) the steady-state trough
increment is lowest just before the next dose, so the safety constraint binds at
the trough (paper sec 7: "A dosing interval is chosen to protect the trough, the
lowest enzymatic increment the patient sees before the next dose"):

    theta_base + Delta_trough(tau) > theta_crit.

TWO dosing regimes, both reported, because they change the trough:
  (a) fixed per-dose amplitude D (the user's formula):
        Delta_trough = D e^-lt/(1-e^-lt)
  (b) same delivered RATE, D proportional to tau (the paper's sec-7 convention:
      "the same delivered rate", "the same mean dose" -- mean increment = D/(lambda tau),
      so holding it fixed forces D ~ tau):
        Delta_trough = D_14 (tau/14) e^-lt/(1-e^-lt)
Regime (b) is the repository's own: it reproduces the paper's "87% cut in the
trough at a 28-day interval" exactly, and its "doubles at 7 days" approximately.
Regime (a) gives wider interval-to-interval gaps, so (b) is the conservative one.

All constants are the repository's own:
    lambda = ln2/3.5 /day, D_14 = 6000 uM/d, grid {7,10,14,15} d  (design/run.py)
    theta_crit^A = (s_pk - k_e m_tol)(K_m+m_tol)/m_tol            (model_A/run.py)
      s_pk = fold * clearance_A(M_OP) at theta=3000; fold=1.0 (routine load) and
      the repo default 1.5 (illness amplitude) BOTH reported -- the convention is
      exactly what the paper flags as unmeasured.
    theta_crit^int(s_base) = theta s.t. max_m C(m;theta) = s_base  (model_C, no m_tol)
"""
import sys; sys.path.insert(0, "../../src")
import numpy as np

import core as C
from model_A.run import ParamsA, M_TOL, baseline_load
from common.models import ParamsC
import model_C.run as MC

LAM = np.log(2.0) / 3.5          # /day, repo LAMBDA
D14 = 6000.0                      # uM/d, repo D_AMP at the 14-day interval
KE, KM = C.KE, 100.0             # 5, 100
TAU0 = 14.0
GRID = [7, 10, 14, 15]
Z = 1.6449                       # one-sided 95% for "determined" headroom

SD_RECORD = 3.05                 # Model A, monthly record (the physical reading)
SD_MANAGED = 0.29                # Model B, managed load
SD_CONTRAST = (0.13, 0.48)       # within-patient monitoring contrast


def trough(tau, regime="rate", lam=LAM, d14=D14):
    x = np.exp(-lam * tau)
    d = d14 * (tau / TAU0) if regime == "rate" else d14
    return d * x / (1.0 - x)


def theta_crit_A(fold):
    return (fold * baseline_load(ParamsA()) - KE * M_TOL) * (KM + M_TOL) / M_TOL


def max_safe_interval(theta_base, thc, regime, grid=GRID):
    safe = [t for t in grid if theta_base + trough(t, regime) >= thc]
    return max(safe) if safe else None


# ---------------------------------------------------------------------------
def item0_validate():
    print("=" * 92)
    print("  VALIDATION: the trough formula reproduces the repository's own numbers")
    print("=" * 92)
    # train_trough_frac = 0.134, train_peak_frac = 2.13 at theta=3000 (results/numbers/design.tex)
    print(f"  repo train_trough_frac 0.134*theta=3000 = 402 uM/d;  formula Delta_trough(14) = {trough(14,'fixed'):.1f}")
    print(f"  repo train_peak_frac   2.13 *theta=3000 = 6390;      formula Delta_peak(14)   = {D14/(1-np.exp(-LAM*14)):.1f}")
    cut28 = 1 - trough(28, "rate") / trough(14, "rate")
    cut7 = trough(7, "rate") / trough(14, "rate")
    print(f"  paper sec 7 'a 28-day interval ... 87% cut in the trough' (same rate): formula = {100*cut28:.0f}% cut")
    print(f"  paper sec 7 'shortening to seven days DOUBLES the trough'  (same rate): formula = {cut7:.2f}x")
    print()


def item1_tolerance():
    print("=" * 92)
    print("  ITEM 1. TOLERANCE ON theta TO CHOOSE BETWEEN ADJACENT INTERVALS")
    print("=" * 92)
    print("\n  Trough increment and the convention-FREE capacity GAP between adjacent")
    print("  intervals, under both dosing regimes (the gap is the ABSOLUTE tolerance;")
    print("  it does not use theta_crit):\n")
    print(f"      {'tau (d)':>8}  {'Delta_trough rate(b)':>20}  {'Delta_trough fixed(a)':>22}")
    for t in GRID:
        print(f"      {t:>8}  {trough(t,'rate'):>20.1f}  {trough(t,'fixed'):>22.1f}")
    pairs = [(14, 15), (14, 10), (14, 7)]
    gaps = {}
    print(f"\n      {'adjacent pair':>16}  {'gap rate(b) uM/d':>18}  {'gap fixed(a) uM/d':>18}")
    for a, b in pairs:
        gr = abs(trough(a, "rate") - trough(b, "rate"))
        gf = abs(trough(a, "fixed") - trough(b, "fixed"))
        gaps[(a, b)] = {"rate": gr, "fixed": gf}
        print(f"      {f'{a} vs {b}':>16}  {gr:>18.1f}  {gf:>18.1f}")

    print("""
  Relative tolerance sd(ln theta) = gap / theta*, theta* = theta_crit - Delta_trough
  of the LONGER interval of the pair (the boundary patient's residual capacity).
  theta_crit rests on the unmeasured s_pk/m_tol convention, so it is reported at the
  routine-load reading (fold 1.0, the "drug dominates" regime the user names) and the
  conservative illness reading (fold 1.5, the repo default). Primary regime: same-rate (b).
""")
    for fold in (1.0, 1.5):
        thc = theta_crit_A(fold)
        th_book = thc - trough(14, "rate")
        ratio = (th_book + trough(14, "rate")) / th_book if th_book > 0 else float("inf")
        print(f"  --- Model A, fold={fold}: theta_crit = {thc:.1f}; booked-boundary theta* = {th_book:.1f}; "
              f"(theta+Delta_trough)/theta* = {ratio:.1f}x")
        print(f"      {'pair':>10}  {'gap(b)':>7}  {'req sd(b)':>10}  {'deriv-formula':>14}  {'req sd(a)':>10}")
        for a, b in pairs:
            hi = max(a, b)
            thb = thc - trough(hi, "rate")
            thb_a = thc - trough(hi, "fixed")
            req_b = gaps[(a, b)]["rate"] / thb if thb > 0 else float("inf")
            req_a = gaps[(a, b)]["fixed"] / thb_a if thb_a > 0 else float("inf")
            lin = LAM * trough(hi, "rate") * (hi - min(a, b)) / thb if thb > 0 else float("inf")
            fb = f"{req_b:.2f}" if np.isfinite(req_b) else "inf"
            fl = f"{lin:.2f}" if np.isfinite(lin) else "inf"
            fa = f"{req_a:.2f}" if np.isfinite(req_a) else "inf"
            print(f"      {f'{a} vs {b}':>10}  {gaps[(a,b)]['rate']:>7.0f}  {fb:>10}  {fl:>14}  {fa:>10}")
        print()
    print("""  READING: the user's linearisation dtau = dtheta/(lambda Delta_trough) is a tangent
  approximation, good to tens of per cent but not exact: against the same-rate column it
  OVERSTATES the tolerance at the fine and mid grids (deriv/req ~1.31 at 14v15, ~1.11 at
  14v10) and UNDERSTATES at the coarse 14v7 (~0.92), because Delta_trough is convex in tau
  and the secant (the true capacity gap) departs from the tangent in both directions as the
  step widens. The qualitative message is unchanged -- the coarse tolerances are large --
  but the exact gaps, not the tangent, are what the table reports. Under the routine-load
  reading the drug dominates theta_crit at the trough and every tolerance is large; under
  the conservative fold-1.5 reading the drug no longer dominates and the fine grid tightens
  to ~0.08-0.13. The fixed-amplitude regime (a) gives wider gaps, so same-rate (b) is
  conservative.""")
    return gaps


# ---------------------------------------------------------------------------
def item2_compare(gaps):
    print("\n" + "=" * 92)
    print("  ITEM 2. REQUIRED TOLERANCE (same-rate b) vs WHAT THE RECORD DELIVERS")
    print("=" * 92)
    print(f"\n  Achievable: record {SD_RECORD}, managed {SD_MANAGED}, "
          f"contrast {SD_CONTRAST[0]}-{SD_CONTRAST[1]}.\n")
    for fold in (1.0, 1.5):
        thc = theta_crit_A(fold)
        print(f"  --- Model A, fold={fold} (theta_crit={thc:.0f}):")
        print(f"      {'pair':>10}  {'req sd':>8}   determined by...")
        for a, b in [(14, 15), (14, 10), (14, 7)]:
            thb = thc - trough(max(a, b), "rate")
            req = gaps[(a, b)]["rate"] / thb if thb > 0 else float("inf")
            if not np.isfinite(req):
                verdict = "trivially (longer interval safe for every PA patient)"
            else:
                who = [n for n, s in [("record 3.05", SD_RECORD), ("managed 0.29", SD_MANAGED),
                                      ("contrast", SD_CONTRAST[1])] if s <= req]
                verdict = (", ".join(who) + " -- and anything tighter") if who else \
                          "NOTHING achievable (needs tighter than contrast)"
            rq = f"{req:.2f}" if np.isfinite(req) else "inf"
            print(f"      {f'{a} vs {b}':>10}  {rq:>8}   {verdict}")
        print()
    print("""  HEADLINE, with its honest limit. This is a GRID-RESOLUTION tolerance: how far apart
  two adjacent intervals sit in theta-space. In that sense the coarse grid (14 vs 10,
  14 vs 7) is determined at the record 3.05 under the routine-load reading and needs no
  new measurement -- the intervals are far apart. That is NOT the same as every patient's
  choice being determined: a PA population straddles the one safety boundary
  theta_crit - Delta_trough(14), and a patient sitting on it is undecidable at any sd
  (item 3). The two reconcile under Model C, where the convention-free theta_crit is low
  enough that Delta_trough(14) alone clears it and the boundary leaves the admissible
  range. So "the therapeutic decision needs no new measurement" holds cleanly where the
  drug trough dominates theta_crit (Model C, routine load) and softens to a boundary
  problem under the conservative fold-1.5 Model-A reading.""")


# ---------------------------------------------------------------------------
def determined_fraction(thc, sd, regime="rate", grid=(7, 10, 14),
                        th_lo=5.0, th_hi=150.0, npat=4000, ndraw=400, seed=0):
    rng = np.random.default_rng(seed)
    th = np.exp(rng.uniform(np.log(th_lo), np.log(th_hi), npat))
    true = [max_safe_interval(t, thc, regime, grid) for t in th]
    hit = tot = 0
    for t, tr in zip(th, true):
        draws = t * np.exp(rng.normal(0.0, sd, ndraw))
        hit += sum(max_safe_interval(d, thc, regime, grid) == tr for d in draws)
        tot += ndraw
    bnds = [thc - trough(t, regime) for t in grid if thc - trough(t, regime) > 0]
    if not bnds:
        unamb = 1.0
    else:
        logd = np.array([min(abs(np.log(t) - np.log(b)) for b in bnds) for t in th])
        unamb = float(np.mean(logd > Z * sd))
    return hit / tot, unamb


def item3_full_uncertainty():
    print("\n" + "=" * 92)
    print("  ITEM 3. THE WHOLE POSTERIOR, NOT JUST THE FIBRE (fraction determined, regime b)")
    print("=" * 92)
    print("""
  The sd's are the brief's FULL-uncertainty numbers: 3.05 already bundles the fibre
  (K_m), the k_e prior, the load-bias b and the assay noise (prior-free NLS record sd);
  0.29 is the managed posterior. So this IS the proper version, not just along the fibre.
  Population: theta_base log-uniform over 0.17%-5% of normal; grid {7,10,14} d.
""")
    thc_int = MC.theta_crit_intrinsic(ParamsC(), MC.S_BASE)
    rows = [("Model A, routine load (fold 1.0)", theta_crit_A(1.0)),
            ("Model A, illness load (fold 1.5)", theta_crit_A(1.5)),
            ("Model C, intrinsic margin (no m_tol)", thc_int)]
    print(f"      {'model / convention':>38}  {'sd':>6}  {'correct-choice':>14}  {'unambiguous':>12}")
    for lab, thc in rows:
        for sd in (SD_RECORD, SD_MANAGED):
            acc, unamb = determined_fraction(thc, sd)
            print(f"      {lab:>38}  {sd:>6.2f}  {acc:>13.1%}  {unamb:>11.1%}")
    print(f"""
  READING: under Model C's intrinsic margin (theta_crit^int(s_base=600) = {thc_int:.1f} uM/d,
  carrying NO m_tol convention) the trough increment Delta_trough(14)=400 alone EXCEEDS
  theta_crit, so every PA patient clears the bar at the booked interval regardless of
  residual capacity: 100% determined at any precision -- the convention-free headline.
  Under Model A a PA population straddles the safety boundary theta_crit-Delta_trough(14),
  so at the record sd 3.05 classification is weak (61% at fold 1.0, 78% at fold 1.5; 0%
  unambiguous either way) and only at the managed sd 0.29 reaches 93-100%. The
  decision survives the full uncertainty where the PARAMETER never does -- but under
  Model A that survival still needs the managed precision for boundary patients; it is
  measurement-free only under Model C's convention-free margin.""")


# ---------------------------------------------------------------------------
def item4_lambda():
    print("\n" + "=" * 92)
    print("  ITEM 4. SENSITIVITY TO lambda (which the repo ASSUMES, and nothing measures)")
    print("=" * 92)
    print("""
  PROVENANCE: design/run.py sets LAMBDA = ln2/3.5 /day, labelled ASSUMED -- "Nothing
  measures it." common/params.py calls it "Partly corroborated, and only early": the
  sole support is a translational PK model's k_deg = 0.0075/h for PCC protein, fitted in
  mouse/rat/monkey and allometrically scaled, inside a two-compartment model whose own
  protein half-life is ~7 d. lambda is not a human measurement, and it enters the trough
  -- and so every tolerance above -- directly.
""")
    print(f"  Delta_trough(14), same-rate, and the fold-1.0 14-vs-15 required sd, across half-life:\n")
    print(f"      {'half-life (d)':>14}  {'lambda (/d)':>11}  {'Delta_trough(14)':>16}  {'req sd 14v15':>13}")
    thc = theta_crit_A(1.0)
    for hl in (2.0, 3.0, 3.5, 5.0, 7.0):
        lam = np.log(2.0) / hl
        tr14 = trough(14, "rate", lam); tr15 = trough(15, "rate", lam)
        thb = thc - tr15
        req = abs(tr14 - tr15) / thb if thb > 0 else float("inf")
        rq = f"{req:.2f}" if np.isfinite(req) else "inf"
        print(f"      {hl:>14.1f}  {lam:>11.4f}  {tr14:>16.1f}  {rq:>13}")
    print("""
  READING: a shorter half-life (faster decay) drains the trough, lowering the whole
  increment toward the safety boundary and tightening the fine grid; a longer half-life
  raises the trough so high the 15-day interval is safe for everyone (req sd infinite,
  trivially determined). The coarse-grid headline is robust across the plausible range
  (the gaps still dwarf a deep-PA theta_base), but the exact boundary interval and the
  fine-grid precision move sharply with an UNMEASURED human parameter. The single
  measurement the interval decision turns on is not sd(ln theta) at all -- it is the human
  decay rate lambda of the therapeutic increment.""")


def main():
    item0_validate()
    gaps = item1_tolerance()
    item2_compare(gaps)
    item3_full_uncertainty()
    item4_lambda()


if __name__ == "__main__":
    main()
