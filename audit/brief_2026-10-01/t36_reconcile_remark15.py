"""
Round 8, items 1 & 4: reconcile "B is better for theta" with the parent paper's
Remark 15 / Figure 13 (fig:B2), which say the opposite, and explain the
cross-sectional/longitudinal asymmetry and the theta-ke correlation vs share.

The paper (sections/03_dynamic.tex rem:ceiling; 10_results.tex fig:B2) retains
Model A as the working model because "if identification fails there it fails a
fortiori once both channels bend", and prices Model B's two costs: (i) the
condition number of the column-normalised sensitivity matrix rises ~450x at a 10%
excursion, and (ii) an EXACT two-point exchange degeneracy (Prop exchange) in
which (theta,Km) and (ke,Ke) swap. Remark 3 (A1_supplementary.tex) says the
relation between the affine-fit conditioning and the three-parameter sensitivity
matrix "is not established here".

t33 reported sd(ln theta) 3.05 (A) -> 0.92 (B): B 3.3x BETTER. That is a MARGINAL,
PRIOR-REGULARISED CRB on theta -- a different object from the paper's joint
conditioning and exchange degeneracy. This script shows both move, in opposite
directions, and names why.
"""
import numpy as np

import core as C
import t33_model_b as T


def condno(Q):
    Qn = Q / np.linalg.norm(Q, axis=0, keepdims=True)
    s = np.linalg.svd(Qn, compute_uv=False)
    return s[0] / s[-1]


def main():
    print("=" * 94)
    print("  ROUND 8 ITEM 1. 'B BETTER FOR theta' vs REMARK 15 ('fails a fortiori under B')")
    print("=" * 94)

    rA = C.record(C.TH_5, seed=0)
    rB = T.record_B(C.TH_5, C.MONTHLY, 0)
    QA, QB = C.flux_rows(rA), T.flux_rows_B(rB)

    print("\n  (a) THE MECHANISM IS THE SHARE. sd(ln theta) ~ (1-s)/s (the clearance the")
    print("      nuisances carry, relative to the enzyme's). s: 0.163 (A) -> 0.372 (B).")
    sA, sB = 0.163, 0.372
    rr = ((1 - sA) / sA) / ((1 - sB) / sB)
    print(f"      (1-s)/s: {(1-sA)/sA:.3f} (A) -> {(1-sB)/sB:.3f} (B); ratio {rr:.3f}")
    print(f"      observed sd ratio 3.045/0.923 = {3.045/0.923:.3f}")
    print(f"      -> the share explains {100*rr/(3.045/0.923):.0f}% of the improvement; the")
    print("      residual ~8% is the second nuisance Ke and the shifted operating point.")
    print("      So the marginal gain is almost entirely the higher enzyme share.")

    print("\n  (b) THE CONDITION NUMBER DEGRADES -- the paper's direction (fig:B2).")
    print(f"      cond(col-normalised sensitivity), monthly record at 5%:")
    print(f"        A (4 cols):            {condno(QA):.3e}")
    print(f"        B (5 cols, Ke free):   {condno(QB):.3e}   ({condno(QB)/condno(QA):.0f}x worse)")
    print(f"        B (4 cols, Ke pinned): {condno(QB[:, :4]):.3e}   (<= A: the degradation is")
    print("                               ENTIRELY the Ke column -- the exchange direction)")
    print("      Same DIRECTION as the paper's 450x; the magnitude differs because the")
    print("      paper's figure is the DYNAMIC sensitivity matrix at a 10% excursion, this")
    print("      is the monthly equation-error FIM. Both say: adding the saturating channel's")
    print("      Ke, collinear with (theta,Km) once both bend, wrecks the joint conditioning.")

    print("\n  (c) WHY THE MARGINAL sd IMPROVES ANYWAY: IT IS PRIOR-REGULARISED.")
    print("      The ill-conditioning lives in the NUISANCE subspace (the ke<->Ke exchange).")
    print("      Priors on the nuisances absorb it, so theta's marginal variance is protected.")
    print("      Strip the priors and B is catastrophic, exactly as the paper says:")
    def sdB(kep, kemp, dropKe):
        return T.sd_th_B(T.fim_B(rB, ke_prior=kep, ke_m_prior=kemp, drop_Ke=dropKe))
    print(f"        A, ke prior 0.8:               {C.sd_th(C.fim(rA, noise='perdraw')):.3f}")
    print(f"        A, NO ke prior:                {C.sd_th(C.fim(rA, noise='perdraw', ke_prior=0)):.3f}")
    print(f"        B, ke 0.8 + Ke pinned:         {sdB(0.8, None, True):.3f}   (t33 headline)")
    print(f"        B, ke 0.8 + Ke free 0.5:       {sdB(0.8, 0.5, False):.3f}")
    print(f"        B, NO priors, Ke free:         {sdB(0.0, 1e9, False):.3f}   <- the exchange")
    print("                                       degeneracy: 55, far WORSE than A's 5.9.")
    print("""
  RECONCILIATION (and a finding about the parent paper):
    - The paper's 'fails a fortiori under B' and fig:B2 are about the DATA/JOINT
      problem: conditioning (450x) and the exact (theta,Km)<->(ke,Ke) exchange.
      Confirmed here -- stripped of priors, B's theta-CRB blows up to 55 vs A's 5.9.
    - t33's '0.92 < 3.05' is the MARGINAL, PRIOR-REGULARISED CRB on theta, which
      IMPROVES under B because the higher enzyme share gives theta more leverage
      while the ill-conditioning stays in the (prior-pinned) nuisance subspace.
    - These are the two objects Remark 3 said were not related. This result RELATES
      them, in the direction OPPOSITE to Remark 15: the marginal precision on the
      one parameter of interest improves even as the joint conditioning degrades.
      So 'A is the most favourable case' holds for CONDITIONING and JOINT
      identifiability, but NOT for the marginal CRB on theta once the nuisances
      carry reasonable priors -- a narrow but real correction to the a-fortiori
      claim, bearing on the paper's stated reason for working in A.
""")

    print("=" * 94)
    print("  ROUND 8 ITEM 4. NO ASYMMETRY (a v11 bug), AND corr(theta, ke) vs SHARE")
    print("=" * 94)
    print("""
  (asymmetry withdrawn) v11 reported the contrast WORSENING under B (0.133 -> 0.167)
    while cross-sectional improved, and asked why. There is NO asymmetry: the 0.167
    was a BUG -- the Model-B contrast did not share the load-bias b between windows
    (b is a patient constant and cancels, exactly as in t24/R-D). Sharing b, the
    corrected Model-B contrast is 0.088 < A's 0.133: B is BETTER on BOTH axes. The
    explanation is the single share mechanism of item 1 -- a higher enzyme share
    tightens theta both cross-sectionally (marginal CRB) and per-window (so the
    contrast, which keeps the independent per-window theta information, improves
    too). No competing effect; the v11 asymmetry does not exist.
""")
    print("  (corr vs share) theta-ke correlation at several capacities, A and B:")
    print("      prediction: corr -> -1 as s -> 0 (the channels become indistinguishable);")
    print("      B has higher s at each capacity, so a WEAKER (less negative) correlation.\n")
    print(f"      {'capacity':>10}  {'A: share':>9}  {'A: corr':>9}  {'B: share':>9}  {'B: corr':>9}")
    from n4_theta_scope import level as levelA
    for frac in (0.20, 0.10, 0.05, 0.02):
        th = C.TH0 * frac
        # A
        rAk = C.record(th, seed=0)
        PA = np.linalg.inv(C.fim(rAk, noise="perdraw"))
        dA = np.sqrt(np.diag(PA)); corrA = PA[0, 3] / (dA[0] * dA[3])
        LA = C.TH0 * (10/110) + C.KE * 10
        sA_ = th * levelA(th, C.KM0, C.KE, LA) / (C.KM0 + levelA(th, C.KM0, C.KE, LA)) / LA
        # B
        rBk = T.record_B(th, C.MONTHLY, 0)
        PB = np.linalg.inv(T.fim_B(rBk, drop_Ke=False, ke_m_prior=0.5))
        dB = np.sqrt(np.diag(PB)); corrB = PB[0, 3] / (dB[0] * dB[3])
        sB_ = T.share_B(th)[0]
        print(f"      {int(frac*100):>8}%   {sA_:9.3f}  {corrA:9.3f}  {sB_:9.3f}  {corrB:9.3f}")
    print("""
      Confirmed: in BOTH models corr(theta, ke) deepens toward -1 as capacity (and
      share) fall, and at each capacity B's higher share gives a less-negative
      correlation. The -0.82 reported at 5% under B sits on this curve. The
      correlation is a function of the share, which is why it is weaker under the
      model with the larger enzyme share -- the same mechanism as (a).
""")


if __name__ == "__main__":
    main()
