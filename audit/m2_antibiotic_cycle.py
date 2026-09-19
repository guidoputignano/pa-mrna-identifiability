"""
m2 -- the periodic input the design does not know about.

Two facts verified against primary sources in this pass, neither of which
appeared in the September 2026 review pack:

  Baumgartner et al. 2014 (PMC4180313, verbatim): "Metronidazole, 10-20 mg/kg/day
  divided in 2-3 times... 1-2 weeks of therapy alternating with 2-3 weeks off or
  alternating every month is advisable."

  Thompson GN, Chalmers RA, Walter JH et al., Eur J Pediatr 1990;149(11):792-6
  (PMID 2226555, verified; 4 PA and 5 MMA children): on metronidazole "urinary
  excretion of propionate metabolites fell significantly... the mean decrease
  being 41% (range 12-76, P less than 0.01), while mean plasma propionate was
  reduced from 45.0 mumol/l to 25.1 mumol/l".

So in a patient on the standard regimen the non-dietary load is not constant.
It is close to a square wave of order 40% amplitude on a 21 to 35 day period,
sitting alongside a 14 or 21 day dose train and a 30.4 day sampling grid. The
paper treats cyclical antibiotics as a pre-therapy instrument (A1, 07 "three
readings") and never as a concurrent periodic input.

Two questions, neither with an answer written in advance:
  (a) if the antibiotic cycle IS in the load record, what does it buy?
  (b) if it is NOT (a diet prescription does not record antibiotic weeks), what
      does an unrecorded 40% periodic swing cost?
"""
import sys
import numpy as np
sys.path.insert(0, 'src')
import matplotlib
matplotlib.use('Agg')
import design.run as D, model_A.run as LOC

P = D.P
L0 = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
T = np.arange(0.0, 730.0, 30.4)          # the paper's monthly schedule, n = 24
AMP = 0.41                                # Thompson 1990, Eur J Pediatr


def abx(t, period, amp=AMP, duty=0.4, phase=0.0):
    """On/off square wave: load suppressed by `amp` while on."""
    frac = ((np.asarray(t, float) / period) + phase) % 1.0
    return 1.0 - amp * (frac < duty)


def steady(Ltrue, dd):
    B = P.ke * P.Km + P.theta + dd - Ltrue
    return (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ltrue * P.Km)) / (2.0 * P.ke)


def columns(lv, dd):
    return np.array([[P.theta * m / (P.Km + m),
                      -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                      1.0, P.ke * m] for m, d in zip(lv, dd)])


def sd_from(Q, sigma):
    F = np.linalg.inv(np.diag([0.8 ** 2, 0.8 ** 2, (D.LOAD_FRAC * L0) ** 2, 0.8 ** 2]))
    F = F + Q.T @ Q / sigma ** 2
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def run(period, recorded, seeds=24, load_cv=0.10):
    """recorded=True: the antibiotic weeks are in the load record, so they are
    excitation. recorded=False: the record is a prescription and misses them, so
    the swing lands in the residual."""
    out = []
    for s in range(seeds):
        rng = np.random.default_rng(200 + s)
        dd = np.asarray(D.delta(T), float)
        base = L0 * (1.0 + load_cv * rng.standard_normal(T.size))
        Ltrue = base * abx(T, period, phase=rng.random())
        lv = steady(Ltrue, dd)
        if recorded:
            sigma = SIG
        else:
            # what the record says is `base`; the abx swing is unmodelled, and
            # its size relative to the load sets the extra ordinate error
            extra = np.std(Ltrue - base, ddof=1)
            sigma = float(np.hypot(SIG, extra))
        out.append(sd_from(columns(lv, dd), sigma))
    return np.array(out)


def main():
    print(__doc__)
    dd = np.asarray(D.delta(T), float)
    ref = []
    for s in range(24):
        rng = np.random.default_rng(200 + s)
        base = L0 * (1.0 + 0.10 * rng.standard_normal(T.size))
        ref.append(sd_from(columns(steady(base, dd), dd), SIG))
    ref = np.array(ref)
    print(f"no antibiotic cycle at all: median {np.median(ref):.4f} "
          f"[{ref.min():.4f}, {ref.max():.4f}]")
    print()
    print(f"{'abx period (d)':>15} {'in the record':>26} {'NOT in the record':>26}")
    for per in (21.0, 28.0, 30.4, 35.0, 60.0):
        a = run(per, True)
        b = run(per, False)
        print(f"{per:15.1f} "
              f"{format(np.median(a), '.4f') + ' [' + format(a.min(), '.4f') + ', ' + format(a.max(), '.4f') + ']':>26} "
              f"{format(np.median(b), '.4f') + ' [' + format(b.min(), '.4f') + ', ' + format(b.max(), '.4f') + ']':>26}")
    print()
    print("30.4 d is the sampling interval itself, which is the aliasing case:")
    print("every draw lands at the same phase of the antibiotic cycle.")


if __name__ == "__main__":
    main()
