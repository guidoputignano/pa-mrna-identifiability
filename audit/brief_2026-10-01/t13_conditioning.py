"""
Item 9 (ten-minute check). Does the likelihood condition on the per-draw
recorded load L(t_k) and on the per-draw drug increment Delta(t_k), or does it
collapse the load to one unknown constant per patient?

Read straight off core.record / core.flux_rows / core.fim, with the numbers
that prove each clause. Nothing new is estimated.
"""
import numpy as np

import core as C
from n4_theta_scope import delta_at


def main():
    print("=" * 92)
    print("  ITEM 9.  WHAT THE LIKELIHOOD CONDITIONS ON")
    print("=" * 92)

    r = C.record(C.TH_5, seed=0)
    L = r["L"]

    print(f"  one 5%-capacity record, 25 monthly draws; anchored load L = {L:.2f}")
    print()
    print("  (1) THE LOAD IS PER-DRAW, NOT ONE CONSTANT.")
    print(f"      recorded loads L(t_k): min {r['loads'].min():.1f}, "
          f"median {np.median(r['loads']):.1f}, max {r['loads'].max():.1f} "
          f"(sd {r['loads'].std():.1f} = {100*r['loads'].std()/L:.1f}% of L)")
    print(f"      each marker m(t_k) solves the balance at ITS OWN load and increment:")
    print(f"        m(t_k): min {r['m'].min():.2f}, median {np.median(r['m']):.2f}, "
          f"max {r['m'].max():.2f}")
    #  m is a deterministic function of (load_k, inc_k); show the coupling
    mm = C.m_ss(r["theta"], r["Km"], r["ke"], r["loads"], r["inc"])
    print(f"      m recomputed from (load_k, inc_k) matches r['m'] to "
          f"{np.abs(mm - r['m']).max():.1e}: the record IS this coupling")
    print()

    print("  (2) THE DRUG DOSE IS A PER-DRAW KNOWN INCREMENT Delta(t_k).")
    print(f"      increment on the enzyme channel: min {r['inc'].min():.1f}, "
          f"median {np.median(r['inc']):.1f}, max {r['inc'].max():.1f}")
    print(f"      it enters the mean as (theta + Delta_k) g(m_k): see core.flux_rows "
          f"col 0 = theta g, and _resid 'pred = (th + inc)*g + ke*m + b'")
    #  dose linearity check
    inc2 = delta_at(r["times"], 12000.0)
    print(f"      doubling the booked dose doubles the increment: "
          f"max ratio {np.nanmax(inc2 / np.where(r['inc'] > 0, r['inc'], np.nan)):.3f}")
    print()

    print("  (3) WHAT *IS* ONE UNKNOWN CONSTANT: THE LOAD-RECORD BIAS b.")
    Q = C.flux_rows(r)
    print(f"      flux_rows has four columns {C.COORDS}.")
    print(f"      the b-column is constant across draws: "
          f"unique values {np.unique(np.round(Q[:, C.I_B], 9))}")
    print(f"      b is a single systematic offset between recorded and true load,")
    print(f"      with a prior sd = 5% of L = {C.LOAD_FRAC*L:.2f} (flux units).")
    print(f"      The per-draw scatter (CV 10%) is RECORDED, so it is conditioned on,")
    print(f"      not estimated: it moves m_k around and is already in the design.")
    print()

    print("  (4) CONSEQUENCE FOR THE 3.07.")
    base = C.sd_th(C.fim(r, noise="perdraw"))
    #  the covariation of L and m is already used: kill it (constant load) and
    #  the record gets no better -- it gets a touch worse, the spread was signal.
    r_flat = C.record(C.TH_5, seed=0, load_cv=0.0)
    flat = C.sd_th(C.fim(r_flat, noise="perdraw"))
    #  the only load knowledge NOT in the 3.07 is a tighter bias b
    tight = C.sd_th(C.fim(r, noise="perdraw", b_prior=1e-6 * L))
    print(f"      sd(ln theta) on the per-draw record (the 3.07):        {base:.3f}")
    print(f"      same record with the load held flat (no draw spread):  {flat:.3f}")
    print(f"      same record with the bias b pinned (L known exactly):  {tight:.3f}")
    print()
    print("  VERDICT: the likelihood conditions on L(t_k) and Delta(t_k) per draw.")
    print("  The L-m covariation is already inside the 3.07. The ONLY load")
    print("  information not yet used is a reduction of the systematic bias b")
    print("  below its 5% prior -- that is item 13's 'measure L to X%'.")


if __name__ == "__main__":
    main()
