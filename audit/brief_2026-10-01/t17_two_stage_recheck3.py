"""
Adversarial independent recheck of ITEM 15 (out/t17_two_stage.txt).

Imports ONLY core + dyn primitives (NOT t17_two_stage, NOT t7_design, NOT
t14_ke_sampling). Everything the claim relies on -- the monthly-record FIM,
the record_n builder, the joint transient+record FIM, the two-stage sweep --
is rebuilt here from first principles. Constants that live in source files
(KM_REC = 0.062 sd(ln K_m) from the record; FOLD = 1.177 from
results/numbers/boundary.tex) are both re-derived / re-read below rather than
trusted.

Claims under test (from the brief):
  (a) safe withdrawal transient, K_m prior 0.062 (dyn idx 2), L0 to 5%/1% (idx 3):
      marginal sd(ln k_e)=sqrt(inv F[1,1]) ~0.20 @5%, ~0.06 @1%; sd(ln theta)
      from the SAME transient stays ~1.2.
  (b) two-stage: k_e fixed @0.20 then n monthly theta-draws gives sd(ln theta)
      ~0.96(n=24)..0.81(n=365); @k_e 0.06 ~0.49..0.33. Does NOT reach 0.20.
  (c) fold value 1.177; every load-up f that drives sd(ln theta) below 0.20
      needs f well above 1.177; f<=1.177 never reaches 0.20 in a joined
      transient+record FIM.
"""
import numpy as np

import core as C
import dyn as D

np.set_printoptions(linewidth=140)

L = D.LTOT
TH = C.TH_5
KM_REC = 0.062               # sd(ln K_m) the record delivers (Task 1); re-derived below
FOLD = 1.177                 # results/numbers/boundary.tex fold_load_ratio; re-read below
HORIZON = 730.0              # two-year window, as in the monthly record


def sdmarg(F, idx):
    """Marginal sd of coord idx from an information matrix, with an eig guard."""
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-12 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[idx, idx]))


# ============================================================ rebuilt helpers
#  The monthly record FIM in core coords (ln theta, ln K_m, b, ln k_e), WITHOUT
#  a k_e prior (it is added at join time), b prior = 5% of L (core default).
FR = C.fim(C.record(TH, seed=0), ke_prior=None)


def joint_fim(Fd, ke_sd=None, Km_sd=None, L0_sd=None, kap_sd=None):
    """
    Rebuild the joined transient+record information matrix from scratch.

    Joint coords  z = (ln theta, ln k_e, ln K_m, L0, b, [ln kappa]).
      transient (dyn) coords (ln theta, ln k_e, ln K_m, L0[, ln kap])
                       -> joint indices [0, 1, 2, 3, (5)]
      record   (core) coords (ln theta, ln K_m, b, ln k_e)
                       -> joint indices [0, 2, 4, 1]
    L0 (dynamic unrecorded load) and b (record load-bias) are DISTINCT nuisance
    parameters, so they occupy separate joint slots. The record already carries
    its own b prior inside FR. A K_m prior is NOT added (the record supplies K_m).
    """
    has_kap = Fd.shape[0] > 4
    n = 6 if has_kap else 5
    Z = np.zeros((n, n))
    dyn_map = [0, 1, 2, 3] + ([5] if has_kap else [])
    Z[np.ix_(dyn_map, dyn_map)] += Fd
    rec_map = [0, 2, 4, 1]                      # th, Km, b, ke
    Z[np.ix_(rec_map, rec_map)] += FR
    if ke_sd:
        Z[1, 1] += ke_sd ** -2
    if L0_sd:
        Z[3, 3] += L0_sd ** -2
    if kap_sd and has_kap:
        Z[5, 5] += kap_sd ** -2
    # Km_sd deliberately ignored: record supplies K_m (matches the claim's design)
    return sdmarg(Z, 0)                          # sd(ln theta)


def record_n(n, seed):
    """n draws spread uniformly over the two-year window, dose train on."""
    times = np.linspace(0.0, HORIZON, int(n))
    return C.record(TH, seed=seed, times=times)


def two_stage_theta(ske, n):
    """median-over-seeds sd(ln theta): n monthly draws, k_e pinned to sd ske."""
    def one(s):
        return C.sd_th(C.fim(record_n(n, s), noise="perdraw", ke_prior=ske))
    return C.over_seeds(one)[0]


def main():
    print("=" * 96)
    print("  INDEPENDENT RECHECK OF ITEM 15  (t17_two_stage.txt)")
    print("=" * 96)
    print(f"  L = dyn.LTOT = {L:.5f}   (expect 322.73)")
    print(f"  TH = 5% theta = {TH:.4f}   K_m = {C.KM0}   k_e = {C.KE}")

    # ---- side check: is KM_REC = 0.062 really what the record delivers? ----
    sd_km_record = sdmarg(C.fim(C.record(TH, seed=0)), 1)   # core idx 1 = ln K_m
    print(f"\n  SIDE CHECK  sd(ln K_m) from the monthly record (core idx 1, with")
    print(f"              the full default priors incl. k_e 0.8) = {sd_km_record:.4f}")
    print(f"              KM_REC used by the claim = {KM_REC}")

    # 5-min sampling over 48 h, as in the claim
    tt = np.arange(0, 2.0 + 1e-9, 5 / 60 / 24)
    print(f"  transient grid: {tt.size} samples, {tt[0]:.4f}..{tt[-1]:.4f} d (5 min, 48 h)")

    # ================================================================ PART (a)
    print("\n" + "-" * 96)
    print("  (a) SAFE FULL WITHDRAWAL: dyn.Experiment([(-1,0.5L),(0,0.0)], L0=0.5L)")
    print("      K_m prior 0.062 (dyn idx 2), L0 prior at 5% / 1% / exact (dyn idx 3)")
    print("-" * 96)
    exW = D.Experiment([(-1, 0.5 * L), (0, 0.0)], L0=0.5 * L)
    mW, _ = D.simulate(D.params(ex=exW), exW, tt, sens=False)
    print(f"      marker over the transient: {mW.min():.2f} .. {mW.max():.2f}"
          f"  (drops from {mW.max():.1f}; untreated 5% marker ~54 -> DROPS, safe)")
    a_tab = []
    for l0, lab in ((0.05 * L, "L0 5%"), (0.01 * L, "L0 1%"), (1e-6, "L0 exact")):
        F = D.fim(D.params(ex=exW), exW, tt, priors={2: KM_REC, 3: l0})
        sdke = sdmarg(F, 1)
        sdth = sdmarg(F, 0)
        a_tab.append((lab, sdke, sdth))
        print(f"        {lab:>9}:  sd(ln k_e) = {sdke:.4f}   sd(ln theta) = {sdth:.4f}")
    # mild 25%-of-diet withdrawal
    mild = D.Experiment([(-1, 0.5 * L), (0, 0.375 * L)], L0=0.5 * L)
    mm, _ = D.simulate(D.params(ex=mild), mild, tt, sens=False)
    Fm = D.fim(D.params(ex=mild), mild, tt, priors={2: KM_REC, 3: 0.05 * L})
    print(f"      mild 25% withdrawal (marker {mm.min():.1f}-{mm.max():.1f}): "
          f"sd(ln k_e) @L0 5% = {sdmarg(Fm, 1):.4f}")

    # ================================================================ PART (c)
    print("\n" + "-" * 96)
    print("  (c) FOLD + ADMISSIBILITY SWEEP (joined transient+record FIM)")
    print("      fold check and whether any f<=1.177 reaches sd(ln theta)=0.20")
    print("-" * 96)
    # re-read the fold from the repository's generated number file
    import re
    from pathlib import Path
    btex = (Path(__file__).resolve().parents[2] / "results" / "numbers" / "boundary.tex").read_text()
    fold_src = float(re.search(r"fold_load_ratio\\endcsname\{([\d.]+)\}", btex).group(1))
    foldm_src = float(re.search(r"fold_m\\endcsname\{([\d.]+)\}", btex).group(1))
    print(f"      boundary.tex fold_load_ratio = {fold_src}  (claim uses {FOLD})")
    print(f"      boundary.tex fold_m          = {foldm_src}  (marker at the saddle)")
    print(f"      {'f x load':>9}  {'marker':>7}  {'admissible':>10}  {'sd(ln theta)':>12}")
    sweep = []
    for f in (0.5, 0.75, 1.0, 1.177, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        m_f = float(C.m_ss(TH, C.KM0, C.KE, f * L))
        adm = f <= FOLD + 1e-9
        if f <= 1.0:
            ex = D.Experiment([(-1, 0.5 * L), (0, (f - 0.5) * L)], L0=0.5 * L)
        else:
            ex = D.Experiment([(-1, 0.5 * L), (0, 0.5 * L + (f - 1) * L)], L0=0.5 * L)
        Fd = D.fim(D.params(ex=ex), ex, tt)                 # transient FIM, no priors
        sdth = joint_fim(Fd, ke_sd=0.8, L0_sd=0.01 * L)     # join to record
        sweep.append((f, m_f, adm, sdth))
        print(f"      {f:>9.3f}  {m_f:>7.1f}  {'yes' if adm else 'NO':>10}  {sdth:>12.4f}")
    reach = [f for f, m, a, s in sweep if s <= 0.20]
    adm_reach = [f for f, m, a, s in sweep if s <= 0.20 and a]
    adm_min = min((s for f, m, a, s in sweep if a), default=np.inf)
    print(f"      f reaching sd<=0.20: {reach}; of these ADMISSIBLE (f<=1.177): "
          f"{adm_reach if adm_reach else 'NONE'}")
    print(f"      best (lowest) sd(ln theta) over the ADMISSIBLE band f<=1.177: {adm_min:.4f}")

    # ================================================================ PART (b)
    print("\n" + "-" * 96)
    print("  (b) TWO-STAGE COST: k_e pinned at stage 1, then n monthly theta-draws")
    print("      median over 24 intake seeds, perdraw noise")
    print("-" * 96)
    print(f"      {'draws':>7}  {'k_e~0.20':>10}  {'k_e~0.06':>10}  {'k_e exact':>10}")
    b_tab = []
    for n in (24, 52, 104, 365):
        a = two_stage_theta(0.20, n)
        b = two_stage_theta(0.06, n)
        c = two_stage_theta(1e-4, n)
        b_tab.append((n, a, b, c))
        print(f"      {n:>7}  {a:>10.4f}  {b:>10.4f}  {c:>10.4f}")
    reaches20 = any(a <= 0.20 or b <= 0.20 for n, a, b, c in b_tab)
    print(f"      any safe two-stage (k_e 0.20 or 0.06) reaching sd(ln theta)<=0.20? "
          f"{reaches20}")

    print("\n" + "=" * 96)
    print("  DONE")
    print("=" * 96)


if __name__ == "__main__":
    main()
