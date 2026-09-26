"""
S1, the last piece: Gap B on two baseline footings, and where ke0 comes from.

The reconciliation (`q1b_reconcile.py`) found that the t = 0 draw of the monthly
record carries much of Gap B, and the author has set out what that means. Read
post-dose, the draw sits inside the marker's settling window, which is wrong
for any kernel. Read pre-dose, it becomes an UNTREATED BASELINE draw, and that
one draw is what halves the bias. That is a change of DESIGN, not an error
correction: a patient already on therapy has no such draw. So Gap B is
reported here on two footings, under all three amplitude conventions:

    baseline     therapy starts within the record; the t = 0 draw is taken
                 before the first infusion, an untreated baseline. n = 25.
    on therapy   no untreated draw. The probe's operational version: the
                 t = 0 draw removed, n = 24.
    steady       a sensitivity row for "on therapy": the same 25 monthly draws
                 starting 100 dose intervals into the train, so every draw is
                 on therapy at steady state and none is first-cycle.

Truths: single exponentials at 2.5, 3.5, 5 and 7 d, and the two-compartment
protein model under BOTH readings of where the dose enters (S1.R): `rising`,
through the effect compartment, which is q1's; and `jump`, straight into PCC,
which is the probe's. Analyst: the committed 3.5 d shape, amplitude fixed.
Noisy records on q1's own noise model, 24 seeds, with the seed range and the
coverage of the 95% profile interval.

Part B records the provenance of ke0 and what the reading it rests on omits.

Nothing here records an expected result in advance.

Run:  python audit/q1c_baseline_footings.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import q1_increment_shape as Q

SEEDS = Q.SEEDS
T = Q.T_DOSE
STEADY = 100 * T
CONVS = ("dose", "peak", "mean")


def kernel_jump(kq=None, kdeg=None):
    """The probe's two-compartment kernel: a unit impulse straight into PCC."""
    kq = Q.K1.KQ if kq is None else kq
    kdeg = Q.K1.KDEG if kdeg is None else kdeg
    A = np.array([[-(kq + kdeg), kq], [kq, -kq]]) * 24.0
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.array([1.0, 0.0]))
    return Q.ExpSum(np.real(V[0] * c), np.real(w))


def kernel_with_mrna(t_half_mrna_h, kq=None, kdeg=None, ke0=None):
    """
    The published structure with a plasma-mRNA stage restored: a first-order
    plasma mRNA C13 decaying with half-life `t_half_mrna_h`, feeding
    d(Ce)/dt = ke0 (C13 - Ce), feeding the protein block. Four states, an
    impulse into C13. Neither q1 nor the probe carries this stage.
    """
    kq = Q.K1.KQ if kq is None else kq
    kdeg = Q.K1.KDEG if kdeg is None else kdeg
    ke0 = Q.K1.KE0 if ke0 is None else ke0
    km = np.log(2.0) / t_half_mrna_h
    A = np.array([[-km, 0.0, 0.0, 0.0],
                  [ke0, -ke0, 0.0, 0.0],
                  [0.0, Q.K1.SLOPE, -(kq + kdeg), kq],
                  [0.0, 0.0, kq, -kq]], float)
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.array([1.0, 0.0, 0.0, 0.0]))
    return Q.ExpSum(np.real(V[2] * c), np.real(w) * 24.0)


TRUTHS = {
    "exp 2.5 d": Q.kernel_exp(2.5),
    "exp 3.5 d": Q.kernel_exp(3.5),
    "exp 5.0 d": Q.kernel_exp(5.0),
    "exp 7.0 d": Q.kernel_exp(7.0),
    "2-comp, rising": Q.kernel_attarwala(),
    "2-comp, jump": kernel_jump(),
}


def banner(t):
    print("\n" + "=" * 96)
    print("  " + t)
    print("=" * 96)


def footing_times(name):
    """(times, dose train) for each footing. The long train covers `steady`."""
    tn = np.arange(0.0, STEADY + 800.0, T)
    wn = np.ones(tn.size)
    if name == "baseline":
        return Q.MONTHLY, (tn, wn)
    if name == "on therapy":
        return Q.MONTHLY[1:], (tn, wn)
    return STEADY + Q.MONTHLY, (tn, wn)


def truth_inc(kernel, conv, train, times, ref):
    A = Q.amplitude(kernel, conv, ref)
    return A * Q.assemble(kernel, train[0], train[1], Q.at_draw(times))


# --------------------------------------------------------------------------
def part_a(seeds=SEEDS):
    banner("A. GAP B ON TWO BASELINE FOOTINGS, ALL THREE CONVENTIONS")
    print("  Fixed 3.5 d analyst, amplitude fixed, theta = 3000, noisy, "
          f"{len(seeds)} seeds. Each cell")
    print("  is median theta_hat/theta [min, max], the RMS residual in sigma "
          "and the coverage")
    print("  of the 95% profile interval. The analyst matches the "
          "convention's amplitude.")
    print(f"  Dose-day draws are read {Q.DOSE_DAY_DRAW}; the `baseline` t = 0 "
          "draw is therefore")
    print("  untreated, and it is the only draw that differs between the "
          "first two footings.")
    ref = Q.committed_ref()
    out = {}
    for foot in ("baseline", "on therapy", "steady"):
        times, train = footing_times(foot)
        print(f"\n  FOOTING: {foot}   (n = {times.size}, first draw at "
              f"t = {times[0]:.1f} d)")
        print(f"    {'truth':16s}" + "".join(f"{c:>34s}" for c in CONVS))
        for tname, k in TRUTHS.items():
            cells = []
            for conv in CONVS:
                A = Q.make_analysts(train, ref, conv,
                                    which=("fixed 3.5 d",))["fixed 3.5 d"]
                inc = truth_inc(k, conv, train, times, ref)
                c = Q.cell(A, inc, seeds=seeds, times=times)
                out[(foot, tname, conv)] = c
                r = c["ratio"]
                cov = [x for x in c["cov"] if x is not None]
                cells.append(f"{np.median(r):6.3f} [{r.min():5.3f},"
                             f"{r.max():6.3f}] {np.median(c['resid']):4.2f}s "
                             f"{100 * sum(cov) / len(cov):3.0f}%")
            print(f"    {tname:16s}" + "".join(f"{x:>34s}" for x in cells))
    print("\n  What the baseline draw is worth, median theta_hat/theta:")
    print(f"    {'truth':16s} {'conv':5s} {'baseline':>9s} {'on therapy':>11s} "
          f"{'steady':>8s}")
    for tname in TRUTHS:
        for conv in CONVS:
            v = [np.median(out[(f, tname, conv)]["ratio"])
                 for f in ("baseline", "on therapy", "steady")]
            print(f"    {tname:16s} {conv:5s} {v[0]:9.3f} {v[1]:11.3f} "
                  f"{v[2]:8.3f}")
    return out


# --------------------------------------------------------------------------
def part_b():
    banner("B. WHERE ke0 COMES FROM, AND WHAT THE RISING READING OMITS")
    print(f"""  ke0 = {Q.K1.KE0} /h is Attarwala et al. 2023, Nucleic Acid Ther
  33(2):141-147, Table 1, RSE 3.35% (transcribed in audit/k1_sources.py, KE0).
  SPECIES: MOUSE. The paper scales only PK parameters to human by allometry;
  ke0, Slope, kdeg and kq are fixed from the Pcca-/- (A138T) mouse fit, which
  k1_sources.py finding 2 and brief R2 both state.

  ITS ROLE is not the decay of a bolus. Suppl. Table S1(B) writes

      d(Ce)/dt = ke0 (C13 - Ce),   Ksyn = Slope Ce,

  an effect compartment EQUILIBRATING WITH PLASMA mRNA, C13, which has its own
  two-compartment PK. q1's rising kernel starts Ce at 1 and lets it decay at
  ke0, which is that model with the plasma mRNA stage removed: the protein
  rise then lasts hours because ke0 is fast (half-life {np.log(2) / Q.K1.KE0:.2f} h).
  The probe's jump kernel removes the effect compartment as well. So BOTH
  readings drop a stage the published model has, and the S1.8 timing
  conclusion rests on how fast the protein rises, which is exactly what the
  dropped stage would set.

  What restoring it does, as a SWEEP, because no human plasma-mRNA half-life
  is established in this repository. Brief R2 quotes a human population-PK
  PCCA mRNA clearance of 0.266 L/h and volume 13.9 L (Gao/Zhang, ACoP 2025),
  which as one compartment is a half-life of
  {np.log(2) * 13.9 / 0.266:.1f} h; it is NOT verified here and is carried
  as one point of the sweep.""")
    ref = Q.committed_ref()
    rows = [("jump (probe)", kernel_jump()),
            ("rising, no mRNA stage (q1)", Q.kernel_attarwala())]
    for th in (2.0, 12.0, 36.2, 72.0):
        rows.append((f"rising, mRNA t1/2 {th:g} h", kernel_with_mrna(th)))
    u = np.linspace(0.0, 30.0, 30001)
    print(f"\n    {'kernel':30s} {'time to single-dose peak':>25s} "
          f"{'family sep/CV, functional, window from 1 h / 6 h / 1 d':>56s}")
    for name, k in rows:
        y = k(u)
        tpk = u[int(np.argmax(y))]
        seps = []
        for lo in (0.04, 0.25, 1.0):
            #  the same phase grid as q1b R6 and q1 Part 5, so the no-mRNA row
            #  reproduces S1.8's separations
            g = np.round(np.arange(lo, T - lo + 1e-9, 0.02), 3)

            def ly(kk):
                return np.log(Q.TH + Q.steady_profile(kk, ref, "peak", g))
            d = ly(k) - ly(Q.kernel_exp(Q.TH_HALF0))
            seps.append(abs(d.max() - d.min()) / np.sqrt(2.0))
        print(f"    {name:30s} {24 * tpk:22.1f} h "
              f"{seps[0]:18.4f} {seps[1]:10.4f} {seps[2]:10.4f}")
    print("""
  Time to peak is 24 x the argmax of the single-dose kernel on a 0.001 d grid.
  The separation columns are S1.8's statistic, the best opposed pair of phases
  after each shape has its best free amplitude, in units of the measurement
  CV, peak-matched, readout theta + Delta.""")


def main():
    print("=" * 96)
    print("  S1: GAP B ON TWO BASELINE FOOTINGS, AND THE PROVENANCE OF ke0")
    print("=" * 96)
    part_b()
    part_a()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
