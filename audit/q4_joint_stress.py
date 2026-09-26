"""
S4. One joint stress scenario at normal capacity.

S1 to S3 varied one thing at a time. A reviewer's objection is that the design
survives each insult alone and would not survive them together, so this script
builds ONE scenario in which every stress the audit has found is present at
once, and measures bias, RMSE, coverage and residual on it.

THE SCENARIO, at normal capacity, theta = 3000:

    increment   the two-compartment protein model, under BOTH readings of where
                the dose enters (S1.R): `rising`, through the effect
                compartment, and `jump`, straight into PCC. The analyst carries
                the committed 3.5 d exponential with a fixed amplitude, so the
                shape is misspecified in the way S1 measured.
    marker CV   15%, with 10% beside it. S3 places the crossover with a 10%
                load record at about CV 11%, so 15% is just past the point
                where the marker, not the diet record, carries the noise.
    calibration the TROUGH reading, which `k4_admissibility.py` finds is the
                most internally coherent combination with a 7 d interval and
                which sets a much larger amplitude for the same target.
    train       the 7 d variant, and the trial's real irregular C5 train with
                its delayed dose, skipped dose and de-escalation.
    k_e         drifting within the patient, at S2b.5's swept 0, 5 and 10%.

CARRIED THROUGHOUT, because S1.R and S2b showed they matter:

    footing     `baseline`, with an untreated pre-dose draw at t = 0, against
                `on therapy`, with no untreated draw at all.
    convention  `dose` is PRIMARY, the author's decision: the single-dose peak
                is what a per-dose calibration pins. `peak`, `mean` and the
                scenario's own `trough` are the sensitivity column.

TWO QUESTIONS THE AUTHOR ADDED, both possible specification items:

    (1) does the baseline draw's advantage survive the noise and the 14.5 d
        schedule, where the draws walk in phase against the 14 d train;
    (2) do two or three pre-treatment draws do better than one.

THE t = 0 CONVENTION. Dose-day draws are read PRE-DOSE, so the `baseline`
footing's t = 0 draw is untreated and the two footings differ in that draw
alone.

THE DRIFT CONVENTION. The record is generated with the patient's own k_e drawn
from the prior and the analyst carries the reference centre, which is S2's
`truth` convention. At these widths, 0.15 and below, S2b found the two honest
conventions coincide; `q9` check 15 measures that rather than assuming it.

Nothing here records an expected result in advance.

Run:  python audit/q4_joint_stress.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import q1_increment_shape as Q
import q1c_baseline_footings as Q1C

SEEDS = tuple(range(24))
T = Q.T_DOSE                           # 14 d
MONTHLY = Q.MONTHLY                    # arange(0, 730, 30.4), n = 25, t = 0 in
SCHED_145 = np.arange(0.5, 730.0, 14.5)
CONVS = ("dose", "peak", "mean", "trough")
rebind = Q.rebind


def banner(t):
    print("\n" + "=" * 100)
    print("  " + t)
    print("=" * 100)


# ==========================================================================
#  The trough reading, which q1's three conventions do not carry.
# ==========================================================================
def steady_trough(kernel, period=T):
    """The steady-state increment just BEFORE the next dose, phase tau = period.
    k4_admissibility's trough reading."""
    return float(np.atleast_1d(kernel.steady(np.array(float(period)),
                                              period))[0])


def amplitude(kernel, conv, ref, ref_trough):
    """q1's amplitude, extended with `trough`."""
    if conv != "trough":
        return Q.amplitude(kernel, conv, ref)
    return ref_trough / steady_trough(kernel)


def committed_trough():
    return steady_trough(Q.kernel_exp(Q.TH_HALF0)) * Q.D_AMP


# ==========================================================================
#  Trains and footings.
# ==========================================================================
def uniform_train(period, t_end=900.0):
    t = np.arange(0.0, t_end, period)
    return t, np.ones(t.size)


TRAINS = {"nominal 14 d": uniform_train(T),
          "7 d": uniform_train(7.0),
          "irregular C5": Q1C.Q.irregular_train("C5 de-escalated")}


def footing_times(name, sched=MONTHLY, n_pre=1):
    """
    The draw times on each footing.

    `baseline` keeps the untreated t = 0 draw. `on therapy` drops it. `pre n`
    adds n - 1 further untreated draws BEFORE t = 0, one dose interval apart,
    which is what a run-in visit schedule would look like; the treated draws
    are unchanged, so the comparison is what the extra untreated draws buy and
    not what more draws buy.
    """
    if name == "on therapy":
        return sched[sched > 0.0]
    if name == "baseline":
        return sched if sched[0] == 0.0 else np.concatenate([[0.0], sched])
    pre = -T * np.arange(n_pre - 1, 0, -1.0)
    base = sched if sched[0] == 0.0 else np.concatenate([[0.0], sched])
    return np.concatenate([pre, base])


def truth_inc(kernel, conv, train, times, ref, ref_trough):
    A = amplitude(kernel, conv, ref, ref_trough)
    return A * Q.assemble(kernel, train[0], train[1], Q.at_draw(times))


def analyst_for(train, ref, conv, ref_trough):
    """The committed 3.5 d analyst, amplitude fixed, on this train and
    convention. For `trough` the reference triple is rescaled so the analyst's
    own shape carries the same trough, which is what a trough calibration
    means."""
    if conv == "trough":
        k = Q.kernel_exp(Q.TH_HALF0)
        scale = ref_trough / steady_trough(k)
        r = tuple(np.array(Q.norm_stats_uniform(k)) * scale)
        a = Q.make_analysts(train, r, "dose",
                            which=("fixed 3.5 d",))["fixed 3.5 d"]
        return a
    return Q.make_analysts(train, ref, conv,
                           which=("fixed 3.5 d",))["fixed 3.5 d"]


# ==========================================================================
#  A cell, with the marker CV and the k_e drift applied.
# ==========================================================================
def cell(kernel, conv, train, times, cv, drift, ref, ref_trough,
         seeds=SEEDS):
    """
    One (kernel, convention, train, schedule, CV, drift) cell over seeds.

    The drift is applied by drawing the patient's own k_e per record and
    building the record's footing from it, while the analyst keeps the
    reference footing. Its prior width is widened to the drift so the analyst
    is not told more than the study would tell it.
    """
    with rebind(Q, CV=cv):
        A = analyst_for(train, ref, conv, ref_trough)
        F_ref = Q.Footing()
        if drift:
            F_ref.s_ke = float(np.hypot(0.02, drift) * Q.KE)
        ratio, lerr, resid, cov = [], [], [], []
        for sd in seeds:
            if drift:
                z = np.random.default_rng(90_001 * sd + 11).normal(0.0, drift)
                F_t = Q.Footing(ke=Q.KE * float(np.exp(z)))
            else:
                F_t = F_ref
            inc = truth_inc(kernel, conv, train, times, ref, ref_trough)
            m, s, _ = Q.make_record(inc, sd, F=F_t, times=times)
            f = A.fit(m, s, times=times, F=F_ref)
            ratio.append(f["theta"] / Q.TH)
            lerr.append(np.log(max(f["theta"], 1e-300)) - np.log(Q.TH))
            resid.append(Q.rms_residual(f["theta"], f["Km"], f["ke"], f["b"],
                                        f["kappa"], m, s, f["phi"], F=F_ref))
            lo, hi, _, _ = A.interval(m, s, f, times=times, F=F_ref)
            cov.append(bool(lo <= Q.TH <= hi))
    r = np.array(ratio)
    e = np.array(lerr)
    return dict(ratio=r, med=float(np.median(r)), lo=float(r.min()),
                hi=float(r.max()),
                rmse=float(np.sqrt(np.mean(e ** 2))),
                bias=float(np.median(e)),
                resid=float(np.median(resid)),
                cov=100.0 * sum(cov) / len(cov), n=len(seeds))


def fmt(c):
    return (f"{c['med']:6.3f} [{c['lo']:5.3f},{c['hi']:6.3f}] "
            f"{c['rmse']:5.3f} {c['resid']:4.2f}s {c['cov']:3.0f}%")


HEAD = f"{'median ratio [min, max]':>22s} {'RMSE':>5s} {'res':>5s} {'cov':>4s}"


# ==========================================================================
def part_a(ref, ref_trough):
    banner("A. THE STRESSES ONE AT A TIME, THEN TOGETHER")
    print(f"  theta = 3000, dose convention, {len(SEEDS)} seeds, monthly draws. "
          "Each row adds one stress")
    print("  to the row above it, so the last row IS the joint scenario. "
          "Analyst: committed 3.5 d,")
    print("  amplitude fixed. Coverage is of the 95% profile interval.")
    ladder = [("committed truth, CV 5%", "exp 3.5 d", 0.05, "nominal 14 d", 0.0),
              ("+ 2-comp truth", None, 0.05, "nominal 14 d", 0.0),
              ("+ CV 10%", None, 0.10, "nominal 14 d", 0.0),
              ("+ CV 15%", None, 0.15, "nominal 14 d", 0.0),
              ("+ 7 d train", None, 0.15, "7 d", 0.0),
              ("+ irregular C5", None, 0.15, "irregular C5", 0.0),
              ("+ k_e drift 10%", None, 0.15, "irregular C5", 0.10)]
    out = {}
    for foot in ("baseline", "on therapy"):
        times = footing_times(foot)
        print(f"\n  FOOTING: {foot}  (n = {times.size})")
        print(f"    {'row':24s}" +
              "".join(f"{f'2-comp {k}':>44s}" for k in ("rising", "jump")))
        print(f"    {'':24s}" + "".join(f"{HEAD:>44s}" for _ in range(2)))
        for label, override, cv, tr, dr in ladder:
            cells = []
            for kname in ("2-comp, rising", "2-comp, jump"):
                k = Q1C.TRUTHS[override if override else kname]
                c = cell(k, "dose", TRAINS[tr], times, cv, dr, ref, ref_trough)
                out[(foot, label, kname)] = c
                cells.append(fmt(c).rjust(44))
            print(f"    {label:24s}" + "".join(cells))
    return out


def part_b(ref, ref_trough):
    banner("B. THE AMPLITUDE CONVENTION, INCLUDING THE TROUGH READING")
    print("  The joint scenario: 2-comp truth, CV 15%, irregular C5 train, "
          "k_e drift 10%. `dose` is")
    print("  primary and the other three are the sensitivity column. The "
          "trough reading is k4's.")
    out = {}
    for foot in ("baseline", "on therapy"):
        times = footing_times(foot)
        print(f"\n  FOOTING: {foot}")
        for kname in ("2-comp, rising", "2-comp, jump"):
            print(f"    {kname:16s}" +
                  "".join(f"{c:>44s}" for c in CONVS))
            print(f"    {'':16s}" + "".join(f"{HEAD:>44s}" for _ in CONVS))
            cells = []
            for conv in CONVS:
                c = cell(Q1C.TRUTHS[kname], conv, TRAINS["irregular C5"],
                         times, 0.15, 0.10, ref, ref_trough)
                out[(foot, kname, conv)] = c
                cells.append(fmt(c).rjust(44))
            print(f"    {'':16s}" + "".join(cells))
    print("\n  And what each convention sets the amplitude to, per kernel:")
    print(f"    {'kernel':16s}" + "".join(f"{c:>14s}" for c in CONVS))
    for kname in ("exp 3.5 d", "2-comp, rising", "2-comp, jump"):
        k = Q1C.TRUTHS[kname]
        print(f"    {kname:16s}" +
              "".join(f"{amplitude(k, c, ref, ref_trough):14.1f}"
                      for c in CONVS))
    return out


def part_c(ref, ref_trough):
    banner("C. k_e DRIFT AT NORMAL CAPACITY")
    print("  S2b.5 found drift binding below 10% of normal. The same prior "
          "carries the headline")
    print("  precision, so it is swept here too: the joint scenario at drift "
          "0, 5 and 10%.")
    for foot in ("baseline", "on therapy"):
        times = footing_times(foot)
        print(f"\n  FOOTING: {foot}")
        print(f"    {'kernel':16s} {'CV':>4s}" +
              "".join(f"{f'drift {int(100 * d)}%':>44s}"
                      for d in (0.0, 0.05, 0.10)))
        print(f"    {'':16s} {'':4s}" + "".join(f"{HEAD:>44s}" for _ in range(3)))
        for kname in ("2-comp, rising", "2-comp, jump"):
            for cv in (0.10, 0.15):
                cells = []
                for dr in (0.0, 0.05, 0.10):
                    c = cell(Q1C.TRUTHS[kname], "dose",
                             TRAINS["irregular C5"], times, cv, dr, ref,
                             ref_trough)
                    cells.append(fmt(c).rjust(44))
                print(f"    {kname:16s} {100 * cv:3.0f}%" + "".join(cells))


def part_d(ref, ref_trough):
    banner("D. QUESTION 1: DOES THE BASELINE DRAW SURVIVE NOISE AND THE 14.5 d "
           "SCHEDULE")
    print("""  The monthly schedule's draws all fall at the same cycle phase; the 14.5 d
  schedule walks in phase against the 14 d train, which is the richer record and
  the one where a single untreated draw has the least to add. Both schedules,
  both kernels, CV 5, 10 and 15%, joint train and drift held at the scenario's
  irregular C5 and 10%. The cell is the median ratio and, beside it, what the
  baseline draw moves it by.""")
    for sname, sched in (("monthly", MONTHLY), ("14.5 d", SCHED_145)):
        print(f"\n  SCHEDULE: {sname}  (n = {sched[sched > 0].size} on therapy)")
        print(f"    {'kernel':16s} {'CV':>4s} {'baseline':>22s} "
              f"{'on therapy':>22s} {'gain':>8s} {'cov base':>9s} "
              f"{'cov ther':>9s}")
        for kname in ("2-comp, rising", "2-comp, jump"):
            for cv in (0.05, 0.10, 0.15):
                cb = cell(Q1C.TRUTHS[kname], "dose", TRAINS["irregular C5"],
                          footing_times("baseline", sched), cv, 0.10, ref,
                          ref_trough)
                ct = cell(Q1C.TRUTHS[kname], "dose", TRAINS["irregular C5"],
                          footing_times("on therapy", sched), cv, 0.10, ref,
                          ref_trough)
                g = abs(ct["med"] - 1.0) - abs(cb["med"] - 1.0)
                print(f"    {kname:16s} {100 * cv:3.0f}% "
                      f"{cb['med']:6.3f} [{cb['lo']:5.3f},{cb['hi']:6.3f}] "
                      f"{ct['med']:6.3f} [{ct['lo']:5.3f},{ct['hi']:6.3f}] "
                      f"{g:+8.3f} {cb['cov']:8.0f}% {ct['cov']:8.0f}%")
    print("""
  `gain` is how much closer to 1 the baseline footing's median sits, in units of
  the ratio itself. A positive number means the untreated draw helps.""")


def part_e(ref, ref_trough):
    banner("E. QUESTION 2: DO TWO OR THREE PRE-TREATMENT DRAWS BEAT ONE")
    print("""  The extra untreated draws are placed one dose interval apart BEFORE t = 0, so
  the treated record is identical across the three columns and what is measured
  is the run-in, not the record length. Joint scenario, both kernels, both
  schedules, CV 15% and 10%, drift 10%.""")
    for sname, sched in (("monthly", MONTHLY), ("14.5 d", SCHED_145)):
        print(f"\n  SCHEDULE: {sname}")
        print(f"    {'kernel':16s} {'CV':>4s} {'none (on therapy)':>22s}" +
              "".join(f"{f'{n} pre-draw(s)':>22s}" for n in (1, 2, 3)))
        for kname in ("2-comp, rising", "2-comp, jump"):
            for cv in (0.10, 0.15):
                row = []
                c0 = cell(Q1C.TRUTHS[kname], "dose", TRAINS["irregular C5"],
                          footing_times("on therapy", sched), cv, 0.10, ref,
                          ref_trough)
                row.append(c0)
                for n in (1, 2, 3):
                    t = (footing_times("baseline", sched) if n == 1
                         else footing_times("pre", sched, n_pre=n))
                    row.append(cell(Q1C.TRUTHS[kname], "dose",
                                    TRAINS["irregular C5"], t, cv, 0.10, ref,
                                    ref_trough))
                print(f"    {kname:16s} {100 * cv:3.0f}%" +
                      "".join(f"{c['med']:6.3f} [{c['lo']:5.3f},"
                              f"{c['hi']:6.3f}]".rjust(22) for c in row))
                print(f"    {'':16s} {'RMSE':>4s}" +
                      "".join(f"{c['rmse']:.3f}".rjust(22) for c in row))


def main():
    print("=" * 100)
    print("  S4. ONE JOINT STRESS SCENARIO AT NORMAL CAPACITY")
    print("=" * 100)
    ref = Q.committed_ref()
    ref_trough = committed_trough()
    print(f"\n  committed reference: single {ref[0]:.1f}, peak {ref[1]:.1f}, "
          f"mean {ref[2]:.1f}, trough {ref_trough:.1f}")
    part_a(ref, ref_trough)
    part_b(ref, ref_trough)
    part_c(ref, ref_trough)
    part_d(ref, ref_trough)
    part_e(ref, ref_trough)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
