"""
Phase 3, what the E2 curves are conditional on.

g2's curves are not unconditional statements about a schedule. Four things are
being assumed, and each one is doing work:

    (a) a 5% prior on the load bias `b`. At a constant recorded load that
        prior is the ONLY thing closing P5's exact fibre, so any comparison
        across load CV with the prior on is a statement about the prior.
    (b) the dose scale is known exactly. The design columns have no dose-scale
        entry, which asserts that `Delta(t_i)` is right to infinite precision.
        design/run.py:943 carries that column; g1 to g4 did not.
    (c) therapy starts at t = 0. A fortnightly schedule is only nearly
        singular because the dose train is still accumulating. A patient
        ALREADY ON THERAPY has no such transient.
    (d) doses land exactly every 14 days at a fixed time of day. The whole
        phase-to-weekday argument of g3 rests on it.

This script prices all four. Writes results/figures/fig_E4_conditioning.{pdf,png}.

Footing as g2's, and stated again beside anything that departs from it:
information filter in (ln theta, ln K_m, b, ln k_e), load record 5%,
sigma = sigma_total(0.05, 0.05), k_e log-prior 0.8, half-day visit grid,
median of 8 seeds with the min-max range quoted.

Nothing here records an expected result in advance. The one number taken from
elsewhere is the repository's own required dose-scale CV, recomputed from
design/run.py rather than typed in.

Run:  python audit/g5_conditions.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from common.style import save, figsize, INK, MUTED, RED, GREEN, BLUE, AMBER, GREY
from design.run import D_AMP, LAMBDA, DOSE_INTERVAL, LOAD_FRAC
import g_design_common as G
import g3_calendar as g3

NS = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
               25, 30, 36, 44, 52, 60])
CVS = (0.0, 0.02, 0.10)
#: the schedules that are not already dead; the two fortnightly ones are
#: carried separately in part (c), where they are the subject.
LIVE = ("14.5 d", "monthly (30.4 d)", "phase-spread")
COLOUR = {"14.5 d": AMBER, "monthly (30.4 d)": BLUE, "phase-spread": GREEN}


def banner(t):
    print("\n" + "=" * 96); print("  " + t); print("=" * 96)


def sd_at(name, n, cv, **kw):
    return G.sd_over_seeds(G.SCHEDULES[name](int(n)), cv, **kw)


def visits_to_gate(name, cv, ns=NS, **kw):
    """The smallest n tried whose MEDIAN sd is at or below the 0.20 target."""
    for n in ns:
        md = sd_at(name, n, cv, **kw)[1]
        if np.isfinite(md) and md <= G.GATE:
            return int(n)
    return None


def fmt(v):
    return "singular" if not np.isfinite(v) else f"{v:.4f}"


def fmt_n(n):
    return f"{n:>5d}" if n else f"{'>' + str(NS[-1]):>5}"


# ------------------------------------------------------------------ part (a)
def part_a_load_prior():
    banner("(a) THE LOAD PRIOR. g2's curves with it, and without it")
    print("  With the prior ON, the three load-CV rows of fig_E2 agree to the")
    print("  fourth decimal, and g2 read that as `load variability buys no")
    print("  precision`. That reading is withdrawn: it is conditional on the")
    print("  prior. Below, the same designs with the prior on `b` DROPPED (the")
    print("  column stays, only the prior goes), which is design/run.py:965's")
    print("  own `free_` construction.\n")
    print(f"    {'schedule':>20} {'load CV':>8} | {'prior ON, n=12':>15} "
          f"{'n to 0.20':>10} | {'prior OFF, n=12':>16} {'n to 0.20':>10}")
    print("  " + "-" * 94)
    for name in LIVE:
        for cv in CVS:
            on = sd_at(name, 12, cv)[1]
            off = sd_at(name, 12, cv, load_prior=False)[1]
            print(f"    {name:>20} {cv:8.0%} | {fmt(on):>15} "
                  f"{fmt_n(visits_to_gate(name, cv)):>10} | {fmt(off):>16} "
                  f"{fmt_n(visits_to_gate(name, cv, load_prior=False)):>10}")
        print()
    print("  At a CONSTANT load and no prior the design is exactly singular, at")
    print("  every schedule and every n. That is P5's fibre, and it is what the")
    print("  prior was filling. Once the load really does vary, the data fill the")
    print("  same direction and the prior is worth a few per cent.")
    print("\n  So the honest statement is the paper's own rank-versus-precision")
    print("  one and not g2's: variation restores RANK, and the 5% prior is an")
    print("  assumption standing in for variation that the record may not have.")


# ------------------------------------------------------------------ part (b)
def required_kappa(name="monthly (30.4 d)", n=25, cv=0.10, seeds=range(8)):
    """
    The dose-scale log-prior at which this design just reaches 0.20, per seed.

    design/run.py:971 solves the same equation on one seed, which is the
    repository's printed requirement. Here it is solved on eight, so the
    requirement can be quoted with its spread.
    """
    t = G.SCHEDULES[name](n)
    out = []
    for s in seeds:
        m, inc, _ = G.record(t, cv, seed=s)
        f = lambda x: G.sd_ln_theta(m, inc, kappa_prior=x) - G.GATE
        try:
            out.append(brentq(f, 1e-3, 3.0, xtol=1e-6))
        except ValueError:
            pass
    return np.array(out)


def part_b_dose_scale():
    banner("(b) THE DOSE SCALE. What happens when Delta is not known exactly")
    kap = required_kappa()
    print("  A dose-scale prior adds the column Delta g(m), which is exactly")
    print("  design/run.py:943's. Leaving it out asserts the infusion delivers")
    print("  its nominal increment to infinite precision.\n")
    print(f"  The repository's own requirement, recomputed: the monthly n=25")
    print(f"  design reaches 0.20 at a dose-scale CV of {np.median(kap):.4f}")
    print(f"  (8 seeds {kap.min():.4f} to {kap.max():.4f}). The paper's printed")
    print(f"  Tier 1 row carries 0.1167, which is the same solve on its own")
    print(f"  single seed; it is not changed here.\n")
    KAPS = [(None, "exact (g1-g4)"), (float(np.median(kap)), None), (0.20, None)]
    for i, (k, lab) in enumerate(KAPS):
        if lab is None:
            KAPS[i] = (k, f"dose-scale CV {k:.1%}")
    print(f"    {'schedule':>20} {'load CV':>8} | " +
          " | ".join(f"{lab:>26}" for _, lab in KAPS))
    print(f"    {'':>20} {'':>8} | " +
          " | ".join(f"{'n=12':>12} {'n to 0.20':>13}" for _ in KAPS))
    print("  " + "-" * 118)
    for name in LIVE:
        for cv in CVS:
            cells = []
            for k, _ in KAPS:
                v = sd_at(name, 12, cv, kappa_prior=k)[1]
                cells.append(f"{fmt(v):>12} {fmt_n(visits_to_gate(name, cv, kappa_prior=k)):>13}")
            print(f"    {name:>20} {cv:8.0%} | " + " | ".join(cells))
        print()
    print("  The `n to 0.20` column steps on the sweep grid, and at the")
    print("  required CV the monthly design is sitting exactly ON the target")
    print("  rather than well past it, because that is how the required CV was")
    print("  defined. At the repository's own n = 25:")
    for k, lab in KAPS:
        v = sd_at("monthly (30.4 d)", 25, 0.10, kappa_prior=k)[1]
        print(f"      monthly n=25, {lab:>22}: {fmt(v)}")
    print("  So `30` above means `25 misses by a thousandth`, not `25 is far")
    print("  off`. Read it as: the repository's requirement buys 0.20 at n = 25")
    print("  and has no margin, and anything worse than it does not.\n")
    print("  A 20% dose-scale CV is not an exotic case: the repository's own")
    print("  between-visit CV row is `mean 14.0%, range 0.4 to 50.8%`, and that")
    print("  is an upper bound on the same-day figure the design actually needs.")
    return KAPS


def part_b_on_the_recommendation(KAPS):
    banner("(b, continued) THE 12-VISIT RECOMMENDATION, RE-PRICED")
    settle, _ = G.settling_days(tol=0.01)
    r = g3.reachable(g3.AFTERNOON, settle=settle)
    t = np.asarray(G.phase_spread(12, gap=30.4, phases=g3.spread_over(r, 12)),
                   float)
    print("  g3's booked design, 12 visits, afternoon dose, every visit on a")
    print("  clinic day. Its headline was 0.1929, which is BELOW 0.20 by 3.5%")
    print("  and inside the seed spread of the target. Here is what it becomes")
    print("  under each conditioning, at load CV 10%:\n")
    print(f"    {'conditioning':>36} {'median':>9} {'8 seeds':>21} {'vs 0.20':>9}")
    print("  " + "-" * 82)
    rows = [("as reported (prior on, dose exact)", dict()),
            ("no load prior", dict(load_prior=False))]
    rows += [(lab, dict(kappa_prior=k)) for k, lab in KAPS if k is not None]
    rows += [("no load prior AND dose-scale CV 20%",
              dict(load_prior=False, kappa_prior=0.20))]
    for lab, kw in rows:
        lo, md, hi = G.sd_over_seeds(t, 0.10, **kw)
        v = "passes" if np.isfinite(md) and md <= G.GATE else "FAILS"
        print(f"    {lab:>36} {fmt(md):>9} "
              f"{'(' + fmt(lo) + ' to ' + fmt(hi) + ')':>21} {v:>9}")
    print("\n  The recommendation is not robust to the dose scale. It should be")
    print("  quoted as `12 visits reach 0.20 IF the dose scale is known to")
    print("  better than about 12%`, with the visit count that follows if it")
    print("  is not.")
    for k, lab in KAPS:
        if k is None:
            continue
        need = None
        for n in NS:
            tt = np.asarray(G.phase_spread(
                int(n), gap=30.4, phases=g3.spread_over(r, int(n))), float)
            if G.sd_over_seeds(tt, 0.10, kappa_prior=k)[1] <= G.GATE:
                need = int(n); break
        print(f"    at {lab}, the booked design needs "
              f"{need if need else '>' + str(NS[-1])} visits")


# ------------------------------------------------------------------ part (c)
def part_c_enrolment():
    banner("(c) WHEN THERAPY STARTS. t = 0 against a patient already on it")
    print("  g1 and g2 ran every schedule from t = 0, so the dose train was")
    print("  still accumulating through the record. For a fortnightly schedule")
    print("  that transient is the ENTIRE excitation: it is why sd came out")
    print("  around 5 rather than infinite. A patient already on therapy at")
    print("  enrolment has none of it.\n")
    print(f"    {'schedule':>28} {'Delta/theta range':>18} {'sd, CV 0%':>11} "
          f"{'sd, CV 10% (8 seeds)':>32}")
    print("  " + "-" * 96)
    for name in G.SCHEDULES:
        t = G.SCHEDULES[name](12)
        lo, md, hi = G.sd_over_seeds(t, 0.10)
        print(f"    {name:>28} {G.delta_range(t):18.4e} "
              f"{fmt(G.sd_over_seeds(t, 0.0)[1]):>11} "
              f"{fmt(md) + ' (' + fmt(lo) + ' to ' + fmt(hi) + ')':>32}")
    print("\n  The already-on-therapy fortnightly design has a Delta/theta range")
    print("  of EXACTLY zero and is exactly singular at a constant load, at")
    print("  every n tried. The t = 0 version's finite sd is an artefact of")
    print("  enrolling at the start of therapy, and no trial should count on it.")
    print("\n  Nothing else changes: the other three schedules cover the cycle")
    print("  by construction and are the same design at either enrolment point.")
    for name in LIVE:
        a = G.sd_over_seeds(G.SCHEDULES[name](12), 0.0)[1]
        b = G.sd_over_seeds(G.fixed(12, {"14.5 d": 14.5,
                                         "monthly (30.4 d)": 30.4}[name],
                                    t0=G.STEADY), 0.0)[1] \
            if name != "phase-spread" else \
            G.sd_over_seeds(G.phase_spread(12, gap=30.4, t0=G.STEADY), 0.0)[1]
        print(f"    {name:>20}: from t=0 {fmt(a)}, already on therapy {fmt(b)}")


# ------------------------------------------------------------------ part (d)
def delta_slipped(times, doses):
    """The increment from a dose train whose doses land at `doses`."""
    times = np.asarray(times, float)[:, None]
    doses = np.asarray(doses, float)[None, :]
    on = times >= doses
    return (D_AMP * np.exp(-LAMBDA * np.where(on, times - doses, 0.0))
            * on).sum(axis=1)


def sd_with_slips(times, slip, cv=0.10, mode="iid", seeds=range(8)):
    """
    sd(ln theta) when each dose lands within +/- `slip` days of schedule.

    mode `iid`   every dose is independently late or early, the schedule does
                 not drift. This is a dose given at a variable time of day.
    mode `walk`  each interval's error accumulates, so the train drifts away
                 from the calendar. This is a missed-and-rebooked infusion.
    """
    times = np.asarray(times, float)
    nd = int(times[-1] // G.T) + 3
    out = []
    for s in seeds:
        rng = np.random.default_rng(1000 + s)
        e = rng.uniform(-slip, slip, nd)
        doses = np.arange(nd) * G.T + (np.cumsum(e) if mode == "walk" else e)
        doses = np.sort(doses)
        inc = delta_slipped(times, doses)
        L = G.LB
        if cv:
            s2 = np.log1p(cv ** 2)
            L = G.LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
        B = G.KE * G.KM + G.TH + inc - L
        m = (-B + np.sqrt(B ** 2 + 4.0 * G.KE * L * G.KM)) / (2.0 * G.KE)
        out.append(G.sd_ln_theta(m, inc))
    v = np.array(out, float)
    v = v[np.isfinite(v)]
    return (np.inf, np.inf, np.inf) if not v.size else \
        (float(v.min()), float(np.median(v)), float(v.max()))


def part_d_slips():
    banner("(d) DOSING REGULARITY. What g3's weekday argument assumes")
    settle, _ = G.settling_days(tol=0.01)
    r = g3.reachable(g3.AFTERNOON, settle=settle)
    t = np.asarray(G.phase_spread(12, gap=30.4, phases=g3.spread_over(r, 12)),
                   float)
    print("  g3's claim that a phase fixes a weekday is arithmetic, and the")
    print("  arithmetic is `14 days is two whole weeks`. It therefore holds")
    print("  EXACTLY as long as doses land exactly 14 days apart at a fixed")
    print("  time of day, and not otherwise. A slip does not make a booked")
    print("  visit unbookable, because the clinic owns the dates; it moves the")
    print("  visit's PHASE, which is what the design was buying.\n")
    SLIPS = [(0.0, "iid", "exactly on schedule"),
             (0.125, "iid", "+/- 3 h, same day"),
             (0.25, "iid", "+/- 6 h, same day"),
             (0.5, "iid", "+/- 12 h, independent"),
             (1.0, "iid", "+/- 1 day, independent"),
             (0.5, "walk", "+/- 12 h, accumulating"),
             (1.0, "walk", "+/- 1 day, accumulating")]
    print(f"    {'dose timing':>24} | {'booked 12-visit design':>38} | "
          f"{'14.5 d, n=12 (control)':>24}")
    print(f"    {'':>24} | {'median':>8} {'spread':>7} {'|dphase|':>9} {'range':>6} | "
          f"{'median':>8} {'spread':>7} {'x wider':>7}")
    print("  " + "-" * 100)
    t145 = G.SCHEDULES["14.5 d"](12)
    base = None
    for slip, mode, lab in SLIPS:
        lo, md, hi = sd_with_slips(t, slip, mode=mode)
        b, mdb, c = sd_with_slips(t145, slip, mode=mode)
        dph, rng = phase_drift(t, slip, mode, settle)
        w, wb = hi - lo, c - b
        if base is None:
            base = (w, wb)
        print(f"    {lab:>24} | {md:8.4f} {w:7.4f} {dph:9.2f} {rng:6.2f} | "
              f"{mdb:8.4f} {wb:7.4f} {wb/base[1]:7.1f}")
    print(f"\n  `spread` is the width of the 8-seed min-to-max interval, and")
    print(f"  `x wider` is that width against the control's own exactly-on-")
    print(f"  schedule width ({base[1]:.4f}); the booked design's widens by "
          f"{(w/base[0]):.1f}x on the same rows.")
    print("\n  |dphase| is the mean absolute drift, in days, between the phase a")
    print("  visit was BOOKED for and the phase it actually lands on; `range` is")
    print("  the Delta/theta range the design realises.")
    print("\n  What the table says, which is not what was expected. The phases")
    print("  DO drift: by 1.5 days on average, a tenth of a cycle, in the worst")
    print("  row. The precision does not follow them. That is f4 once more and")
    print("  it is the reason the spread design is worth having: it is not")
    print("  buying any PARTICULAR phase, it is buying the RANGE, and a bounded")
    print("  shift applied to a well-spread set leaves the range where it was,")
    print("  1.60 to 1.72 across every row against 1.66 on schedule.")
    print("\n  The control is the informative column, and not in the direction")
    print("  of a simple penalty. The 14.5 d median moves BOTH WAYS, sometimes")
    print("  improving, because a slip can hand a thinly-covering schedule")
    print("  phases it would not otherwise have visited. What it loses is")
    print("  predictability: its seed spread widens by the factor in the last")
    print("  column while the booked design's barely moves. A design whose")
    print("  precision is a lottery on the infusion clerk is not a design.")
    print("\n  So the honest conditional on g3 is narrower than `the design")
    print("  needs exact dosing`. It is this: the PHASE-TO-WEEKDAY MAP needs")
    print("  exact 14-day dosing at a fixed time of day, because that map is")
    print("  the arithmetic of 14 = 2 x 7 and nothing else. The precision the")
    print("  design delivers does not. A clinic that loses the map still has")
    print("  the design, provided it re-derives the phases from the actual")
    print("  infusion dates, which the dosing-record row the paper already asks")
    print("  for is exactly what makes possible.")


def phase_drift(times, slip, mode, settle, cv=0.0, seeds=range(8)):
    """Mean |booked phase - realised phase|, and the realised Delta/theta range."""
    times = np.asarray(times, float)
    nd = int(times[-1] // G.T) + 3
    booked = times % G.T
    d, r = [], []
    for s in seeds:
        rng = np.random.default_rng(1000 + s)
        e = rng.uniform(-slip, slip, nd)
        doses = np.sort(np.arange(nd) * G.T +
                        (np.cumsum(e) if mode == "walk" else e))
        got = np.array([x - doses[doses <= x][-1] if (doses <= x).any() else x
                        for x in times])
        gap = np.abs(got - booked)
        d.append(np.mean(np.minimum(gap, G.T - gap)))
        inc = delta_slipped(times, doses)
        r.append((inc.max() - inc.min()) / G.TH)
    return float(np.mean(d)), float(np.mean(r))


# ------------------------------------------------------------------- figure
def figure(KAPS):
    fig, axes = plt.subplots(1, 3, figsize=figsize(1.0, 0.36))
    fig.set_layout_engine("none")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.40,
                        wspace=0.26)

    def dress(ax, title):
        ax.axhline(G.GATE, color=GREY, ls=":", lw=1.0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(3.6, 66); ax.set_ylim(0.09, 3.0)
        ax.set_xticks([4, 6, 10, 20, 40, 60])
        ax.set_xticklabels(["4", "6", "10", "20", "40", "60"])
        ax.set_yticks([0.1, 0.2, 0.5, 1.0, 2.0])
        ax.set_yticklabels(["0.1", "0.2", "0.5", "1", "2"])
        ax.minorticks_off()
        ax.grid(True, which="major", alpha=0.45)
        ax.set_xlabel("number of visits", fontsize=8)
        ax.set_title(title, fontsize=8.5, pad=5)
        ax.tick_params(labelsize=6.5)

    # (a) the load prior
    ax = axes[0]
    for name in LIVE:
        for cv, ls in ((0.02, "--"), (0.10, "-")):
            y = [sd_at(name, n, cv, load_prior=False)[1] for n in NS]
            ax.plot(NS, y, color=COLOUR[name], lw=1.4, ls=ls)
    ax.set_ylabel(r"sd$(\ln\hat\theta)$, median of 8 seeds", fontsize=8)
    dress(ax, "(a) no prior on the load bias")

    # (b) the dose scale
    ax = axes[1]
    for name in LIVE:
        for (k, _), ls in zip(KAPS, ("-", "--", ":")):
            y = [sd_at(name, n, 0.10, kappa_prior=k)[1] for n in NS]
            ax.plot(NS, y, color=COLOUR[name], lw=1.4, ls=ls)
    dress(ax, "(b) dose scale not known exactly")

    # (c) the enrolment point
    ax = axes[2]
    for name, col in (("14 d (the dose interval)", RED),
                      ("14 d, already on therapy", MUTED),
                      ("monthly (30.4 d)", BLUE)):
        y = [sd_at(name, n, 0.10)[1] for n in NS]
        y = [np.nan if not np.isfinite(v) else v for v in y]
        ax.plot(NS, y, color=col, lw=1.5, label=name)
    dress(ax, "(c) when therapy starts")
    ax.set_ylim(0.09, 40)
    ax.set_yticks([0.1, 0.2, 1.0, 5.0, 20.0])
    ax.set_yticklabels(["0.1", "0.2", "1", "5", "20"])
    ax.legend(fontsize=6.0, frameon=False, loc="center right")

    handles = [plt.Line2D([], [], color=COLOUR[k], lw=1.6) for k in LIVE]
    fig.legend(handles, list(LIVE), loc="lower center", ncol=3, fontsize=7.5,
               bbox_to_anchor=(0.5, 0.150), frameon=False)
    fig.text(0.5, 0.010,
             "fig_E2 holds the load prior on, the dose scale exact and therapy "
             "starting at t = 0. Each panel here releases one of those, at "
             r"$k_e$ prior 0.8." "\n"
             "(a) solid load CV 10%, dashed CV 2%; CV 0% is exactly singular "
             "at every n and is not drawn.   (b) load CV 10%; solid exact, "
             "dashed the required CV, dotted 20%.\n"
             "(c) load CV 10%; the grey schedule has a "
             r"$\Delta/\theta$ range of exactly zero and is singular at a "
             "constant load, so what is plotted for it is load variation alone.",
             ha="center", va="bottom", fontsize=6.3, color=MUTED,
             linespacing=1.6)
    save(fig, "fig_E4_conditioning")


if __name__ == "__main__":
    print(__doc__)
    part_a_load_prior()
    KAPS = part_b_dose_scale()
    part_b_on_the_recommendation(KAPS)
    part_c_enrolment()
    part_d_slips()
    figure(KAPS)
