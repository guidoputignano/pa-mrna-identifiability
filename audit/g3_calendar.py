"""
Phase 3, the inverse map: from target cycle phases to dates a clinic can book.

g1 and g2 argue in cycle phase. A clinic books in dates, so the design is only
usable if the map runs the other way: given a target phase, on what day does
the patient come in. That map is the whole of this script, and it turns out to
be more constrained than it looks, for an arithmetic reason:

    the dose interval is 14 days, which is exactly two weeks.

So the phase of a visit determines its WEEKDAY, once and for all. A target
phase is either always bookable or never bookable; there is no rescheduling it
into the working week. That is what makes this worth computing rather than
assuming, and it is the difference between a schedule that reads well on a
dial and one a clinic can actually run.

    Part A  which phases a weekday clinic can reach at all, and what
            Delta/theta each of them buys
    Part B  the booking itself: target phase -> date, for one design
    Part C  the comparison the phase asked for, phase-spread against fixed
            intervals at equal n
    Part D  results/figures/fig_E3_calendar.{pdf,png}

Footing, stated once and true of every number below. Information filter in the
paper's coordinates (ln theta, ln K_m, b, ln k_e), load record 5%, sigma =
sigma_total(0.05, 0.05), k_e log-prior 0.8 unless a column says otherwise;
`k_e known` DROPS the column, which is _filter's convention and e9's. Visits
lie on a half-day grid (morning or afternoon). sd figures at a non-zero load CV
are the median over 8 seeds with the range shown. The first dose is taken to
fall on a Monday; nothing depends on which Monday.

Nothing here records an expected result in advance.

Run:  python audit/g3_calendar.py
"""
import sys
from pathlib import Path
from datetime import date, timedelta

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from common.style import save, figsize, INK, MUTED, RED, GREEN, BLUE, AMBER, GREY
from design.run import delta
import g_design_common as G

#: an arbitrary Monday, so the printed table is a real calendar. Only the
#: weekday matters, and the dose interval being two whole weeks means the
#: weekday never drifts.
DOSE_DAY = date(2026, 1, 5)
DAYNAME = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
CLINIC = (0, 1, 2, 3, 4)          # a weekday clinic
HALF = {0.0: "am", 0.5: "pm"}

#: when the dose is given within its own day. 0.0 is a morning dose, 0.5 an
#: afternoon one. It shifts which phases land on a weekday and so is a design
#: choice, not a detail -- Part A prices it.
MORNING, AFTERNOON = 0.0, 0.5


# --------------------------------------------------------------------- utils
#: every interpolated number this module draws, classified for
#: tools/check_artefacts.py:check_figure_numbers. See that function for what
#: the classifications mean and for what the check can and cannot verify.
DRAWN = {
    "n_cal": "input: figure()'s own n_cal argument, the visit count this "
             "calendar was drawn at",
}


def day_index(t, dose_offset):
    """The calendar day a visit at cycle time t falls on, day 0 = DOSE_DAY."""
    return int(np.floor(t + dose_offset + 1e-9))


def weekday(t, dose_offset):
    return day_index(t, dose_offset) % 7


def inc_of(t):
    """Delta/theta at a visit, from the repository's own dose train."""
    return float(np.atleast_1d(delta(np.atleast_1d(float(t))))[0]) / G.TH


#: far enough in that the dose train has reached steady state. Part A talks
#: about a PHASE rather than a date, and a phase only has a single Delta once
#: the train has settled; read on the first cycle it is the transient, which
#: understates every increment by about 7%.
STEADY = 100 * G.T


def inc_at_phase(p):
    """Delta/theta at cycle phase p, at steady state."""
    return inc_of(STEADY + float(p))


def slot(t, dose_offset):
    """(date, 'am'|'pm') for a visit at cycle time t."""
    d = day_index(t, dose_offset)
    frac = (t + dose_offset) - d
    return DOSE_DAY + timedelta(days=d), HALF[0.5 if frac > 0.25 else 0.0]


# ------------------------------------------------------------- the inverse map
def date_for_phase(phase, earliest, dose_offset=MORNING, period=G.T,
                   grid=G.GRID, clinic=CLINIC):
    """
    THE INVERSE MAP. The first bookable slot at or after cycle time `earliest`
    whose dose-cycle phase is `phase`.

    Returns (t, date, half, weekday index), or None if the phase never falls on
    a clinic day -- which, because 14 days is two whole weeks, is a property of
    the phase alone and not of how long one is willing to wait.
    """
    phase = float(phase) % period
    if weekday(phase, dose_offset) not in clinic:
        return None
    k = int(np.floor((earliest - phase) / period)) + 1
    for k in range(max(k - 2, 0), k + 3):
        t = k * period + phase
        if t >= earliest - 1e-9 and abs(t / grid - round(t / grid)) < 1e-6:
            d, h = slot(t, dose_offset)
            return t, d, h, weekday(t, dose_offset)
    return None


def reachable(dose_offset, period=G.T, grid=G.GRID, clinic=CLINIC, settle=None):
    """
    The phases a weekday clinic can book, on the half-day grid, excluding any
    slot inside the post-dose settling window.
    """
    ph = np.arange(0.0, period, grid)
    out = []
    for p in ph:
        if weekday(p, dose_offset) not in clinic:
            continue
        if settle is not None and p < settle - 1e-9:
            continue
        out.append(float(p))
    return np.array(out)


def spread_over(reach, n):
    """
    n targets spread as evenly as the reachable set allows, INCLUDING both
    ends of it.

    The ends are not an aesthetic choice. f4 found the binding quantity is the
    RANGE of Delta/theta a schedule covers, so a chooser that spaces targets
    evenly but drops the extreme phases gives away the thing that is being
    bought. With n larger than the reachable set the interior targets repeat,
    which is allowed -- two visits at the same phase in different cycles are
    two observations.
    """
    R = len(reach)
    if R == 0:
        return np.array([])
    if n == 1:
        return np.array([reach[0]])
    return np.array([reach[int(round(i * (R - 1) / (n - 1)))] for i in range(n)])


def snap_to_clinic(times, dose_offset, settle, grid=G.GRID, clinic=CLINIC,
                   period=G.T):
    """
    Move each visit to the nearest bookable slot: a half-day on a clinic day,
    outside the settling window. Ties go to the earlier slot.

    This is what a clinic actually does with a schedule that lands on a
    Saturday, and it is the fair comparator for the fixed intervals. It is not
    free: moving a visit moves its phase, which is the quantity the design is
    spending its visits on.
    """
    out = []
    for t in np.asarray(times, float):
        best, bd = None, np.inf
        for k in np.arange(-6, 7):
            c = np.round((t + k * grid) / grid) * grid
            if c < -1e-9:
                continue
            if weekday(c, dose_offset) not in clinic:
                continue
            if (c % period) < settle - 1e-9:
                continue
            d = abs(c - t)
            if d < bd - 1e-12:
                best, bd = float(c), d
        out.append(t if best is None else best)
    return np.array(out)


# ------------------------------------------------------------------- part A
def part_a():
    settle, worst = G.settling_days(tol=0.01)
    print("=" * 96)
    print("  (A) WHICH PHASES A WEEKDAY CLINIC CAN REACH")
    print("=" * 96)
    print(f"  The dose interval is {G.T:g} d = {G.T/7:g} weeks exactly, so a visit's")
    print("  phase fixes its weekday. The reachable set is therefore a property")
    print("  of the schedule, not of the patient's flexibility.")
    print(f"  Slots inside the settling window ({settle*24:.2f} h at a 1% tolerance)")
    print("  are dropped: a draw there is not a steady state.\n")

    full = np.arange(0.0, G.T, G.GRID)
    fullok = full[full >= settle - 1e-9]
    frange = (max(inc_at_phase(p) for p in fullok)
              - min(inc_at_phase(p) for p in fullok))
    print(f"    {'dose given':>14} {'slots':>6} {'of':>4} {'phases reachable':>52}")
    print("  " + "-" * 94)
    rows = {}
    for off, lab in ((MORNING, "morning"), (AFTERNOON, "afternoon")):
        r = reachable(off, settle=settle)
        rows[lab] = r
        txt = ", ".join(f"{p:g}" for p in r)
        if len(txt) > 52:
            txt = txt[:49] + "..."
        print(f"    {lab:>14} {len(r):6d} {len(fullok):4d} {txt:>52}")
    print()
    print(f"    {'dose given':>14} {'Delta/theta reachable':>26} {'range':>8} "
          f"{'of the full':>12} {'trough?':>9}")
    print("  " + "-" * 94)
    for lab, r in rows.items():
        v = np.array([inc_at_phase(p) for p in r])
        rng = v.max() - v.min()
        trough = ("yes" if abs(v.min() - min(inc_at_phase(p) for p in fullok))
                  < 1e-9 else "no")
        print(f"    {lab:>14} [{v.min():7.4f}, {v.max():7.4f}] {rng:8.4f} "
              f"{rng/frange:11.1%} {trough:>9}")
    print(f"\n  The full half-day dial outside the settling window spans "
          f"{frange:.4f}, at steady state. Read on the FIRST cycle instead it\n"
          f"  comes out about 7% smaller, because the dose train has not\n"
          f"  finished accumulating; that is a real feature of the early\n"
          f"  visits, not of the phase, and Part C measures it where it\n"
          f"  belongs, at the actual visit times.")
    print("  Both dose timings leave the same number of slots; they are not the")
    print("  same slots. Which of the two reaches the trough is the whole of the")
    print("  difference, and the trough is one end of the quantity f4 found binds.")

    print("\n  What is lost, phase by phase (blank = not bookable):")
    print(f"    {'phase':>6} {'Delta/theta':>12} {'morning dose':>14} {'afternoon dose':>16}")
    print("  " + "-" * 94)
    for p in full:
        if p < settle - 1e-9:
            tag_m = tag_a = "settling"
        else:
            tag_m = DAYNAME[weekday(p, MORNING)]
            tag_a = DAYNAME[weekday(p, AFTERNOON)]
            tag_m = tag_m if weekday(p, MORNING) in CLINIC else ""
            tag_a = tag_a if weekday(p, AFTERNOON) in CLINIC else ""
        print(f"    {p:6g} {inc_at_phase(p):12.4f} {tag_m:>14} {tag_a:>16}")
    return settle, rows


# ------------------------------------------------------------------- part B
def part_b(settle, n=12, gap=30.4, dose_offset=AFTERNOON):
    print("\n" + "=" * 96)
    print(f"  (B) THE BOOKING: {n} TARGET PHASES TURNED INTO DATES")
    print("=" * 96)
    r = reachable(dose_offset, settle=settle)
    tgt = spread_over(r, n)
    print(f"  Dose in the afternoon, first dose {DOSE_DAY.isoformat()} "
          f"({DAYNAME[DOSE_DAY.weekday()]}), minimum gap {gap:g} d.")
    print(f"  Targets, chosen from the {len(r)} reachable phases: "
          f"{', '.join(f'{p:g}' for p in tgt)}")
    print("  The targets are taken NEAREST-REACHABLE FIRST rather than in")
    print("  numerical order: booking them in order forces a wait of up to a")
    print("  whole cycle between visits and stretches the record for a reason")
    print("  that has nothing to do with phase.\n")
    times = np.asarray(G.phase_spread(n, gap=gap, phases=tgt), float)
    # the inverse map must reproduce the same slots, visit by visit
    prev, agree = None, True
    for t in times:
        hit = date_for_phase(t % G.T, 0.0 if prev is None else prev + gap,
                             dose_offset)
        agree &= hit is not None
        prev = t
    print(f"    {'#':>3} {'target phase':>13} {'date':>12} {'day':>5} {'when':>5} "
          f"{'cycle':>6} {'realised':>9} {'Delta/theta':>12} {'gap (d)':>8}")
    print("  " + "-" * 94)
    prev = None
    for i, t in enumerate(times):
        d, h = slot(t, dose_offset)
        wd = weekday(t, dose_offset)
        print(f"    {i+1:3d} {t % G.T:13g} {d.isoformat():>12} {DAYNAME[wd]:>5} {h:>5} "
              f"{int(t//G.T):6d} {t % G.T:9g} {inc_of(t):12.4f} "
              f"{'' if prev is None else f'{t-prev:8.1f}'}")
        prev = t
    unbookable = [p for p in tgt if date_for_phase(p, 0.0, dose_offset) is None]
    print(f"\n  Every visit falls on a clinic day outside the settling window: "
          f"{all(weekday(t, dose_offset) in CLINIC and (t % G.T) >= settle - 1e-9 for t in times)}.")
    print(f"  Targets the inverse map cannot book at all: "
          f"{len(unbookable)} of {len(tgt)}.")
    print(f"\n  Horizon {times[-1]:.1f} d ({times[-1]/30.4:.1f} months), "
          f"{len(G.cycle_phases(times))} distinct phases, "
          f"Delta/theta range {G.delta_range(times):.4f}.")
    for cv in (0.0, 0.02, 0.10):
        lo, md, hi = G.sd_over_seeds(times, cv)
        print(f"  sd(ln theta) at load CV {cv:4.0%}: median {md:.4f} "
              f"(8 seeds {lo:.4f} to {hi:.4f})")
    return times


# ------------------------------------------------------------------- part C
def variants(settle, gap=30.4):
    """The schedules compared at equal n, each as a function of n."""
    rm = reachable(MORNING, settle=settle)
    ra = reachable(AFTERNOON, settle=settle)
    return {
        "14 d (the dose interval)": (G.SCHEDULES["14 d (the dose interval)"], RED),
        "14.5 d": (G.SCHEDULES["14.5 d"], AMBER),
        "monthly (30.4 d)": (G.SCHEDULES["monthly (30.4 d)"], BLUE),
        "phase-spread, whole dial": (lambda n: G.phase_spread(n, gap=gap), GREY),
        "phase-spread, clinic (am dose)":
            (lambda n: G.phase_spread(n, gap=gap, phases=spread_over(rm, n)), MUTED),
        "phase-spread, clinic (pm dose)":
            (lambda n: G.phase_spread(n, gap=gap, phases=spread_over(ra, n)), GREEN),
    }


def offset_of(name):
    """Which dose timing a schedule's weekday count should be read against."""
    return AFTERNOON if "pm dose" in name else MORNING


def part_c(settle, ns=(8, 12, 16, 25)):
    print("\n" + "=" * 96)
    print("  (C) PHASE-SPREAD AGAINST FIXED INTERVALS, AT EQUAL n")
    print("=" * 96)
    print("  `off-clinic` counts visits the schedule puts on a weekend or inside")
    print("  the settling window, i.e. ones a clinic would have to move -- and")
    print("  moving one changes its phase, so it is a cost, not a formality.\n")
    V = variants(settle)
    for n in ns:
        print(f"  n = {n}")
        print(f"    {'schedule':>32} {'horizon':>9} {'phases':>7} {'range':>7} "
              f"{'off-clinic':>11} | {'sd, CV 0%':>10} "
              f"{'sd, CV 10% (median, 8 seeds)':>32}")
        print("  " + "-" * 110)
        for name, (fn, _) in V.items():
            t0 = np.asarray(fn(n), float)
            off = offset_of(name)
            rows = [(name, t0)]
            if not name.startswith("phase-spread, clinic"):
                rows.append(("  ... snapped to the clinic",
                             snap_to_clinic(t0, off, settle)))
            for lab, t in rows:
                bad = sum(1 for x in t
                          if weekday(x, off) not in CLINIC
                          or (x % G.T) < settle - 1e-9)
                s0 = G.sd_over_seeds(t, 0.0)[1]
                lo, md, hi = G.sd_over_seeds(t, 0.10)
                f = lambda v: "singular" if not np.isfinite(v) else f"{v:.4f}"
                print(f"    {lab:>32} {t[-1]:8.1f}d {len(G.cycle_phases(t)):7d} "
                      f"{G.delta_range(t):7.4f} {bad:5d} /{len(t):4d} | {f(s0):>10} "
                      f"{f(md):>10}  ({f(lo)} to {f(hi)})")
        print()
    print("  Read the off-clinic column against the sd column. A schedule that")
    print("  wins on sd and loses on bookability has not won.")


# ------------------------------------------------------------------- part D
def figure(settle, n_cal=12, gap=30.4, ns=np.array([4, 6, 8, 10, 12, 16, 20, 25, 30])):
    fig = plt.figure(figsize=figsize(1.0, 0.52))
    fig.set_layout_engine("none")
    axc = fig.add_axes([0.075, 0.10, 0.30, 0.78])
    axs = fig.add_axes([0.605, 0.155, 0.370, 0.735])

    # ---- the calendar
    r = reachable(AFTERNOON, settle=settle)
    times = np.array([t for t in
                      G.phase_spread(n_cal, gap=gap, phases=spread_over(r, n_cal))])
    weeks = int(day_index(times[-1], AFTERNOON) // 7) + 1
    for w in range(weeks):
        for dcol in range(7):
            axc.add_patch(Rectangle((dcol, w), 1, 1, facecolor=
                                    ("#eeece5" if dcol >= 5 else "white"),
                                    edgecolor="#d8d5cc", lw=0.4))
    doses = [d for d in range(0, day_index(times[-1], AFTERNOON) + 1, int(G.T))]
    for d in doses:
        axc.add_patch(Rectangle((d % 7 + 0.05, d // 7 + 0.22), 0.11, 0.56,
                                facecolor=RED, edgecolor="none", alpha=0.75))
    v = np.array([inc_of(t) for t in times])
    norm = plt.Normalize(v.min(), v.max())
    cmap = plt.get_cmap("viridis")
    for t, val in zip(times, v):
        d = day_index(t, AFTERNOON)
        frac = (t + AFTERNOON) - d
        axc.scatter([d % 7 + (0.38 if frac < 0.25 else 0.66)], [d // 7 + 0.5],
                    c=[val], cmap=cmap, norm=norm, s=34, zorder=4,
                    edgecolors="white", linewidths=0.5)
    axc.set_xlim(0, 7); axc.set_ylim(weeks, 0)
    axc.set_xticks(np.arange(7) + 0.5); axc.set_xticklabels(DAYNAME, fontsize=6.5)
    step = max(1, int(np.ceil(weeks / 14)))
    axc.set_yticks(np.arange(0, weeks, step) + 0.5)
    axc.set_yticklabels([(DOSE_DAY + timedelta(days=7 * w)).strftime("%d %b")
                         for w in range(0, weeks, step)], fontsize=6.0)
    axc.tick_params(length=0)
    for s in axc.spines.values():
        s.set_visible(False)
    axc.set_title(f"{n_cal} visits, phase-spread on a weekday clinic\n"
                  f"red bar = dose day, dot = draw (left am, right pm)",
                  fontsize=8, pad=6, linespacing=1.5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axc, fraction=0.045, pad=0.08)
    cb.ax.set_title(r"$\Delta/\theta$", fontsize=7, pad=3)
    cb.ax.tick_params(labelsize=6, length=2, pad=1)

    # ---- what the constraint costs
    curves = []
    for name, (fn, col) in variants(settle, gap).items():
        if name.startswith("14"):
            continue
        curves.append((name, fn, col, "-"))
        if name == "monthly (30.4 d)":
            curves.append(("monthly, snapped to the clinic",
                           lambda k, f=fn: snap_to_clinic(f(k), MORNING, settle),
                           col, "--"))
    for name, fn, col, ls in curves:
        y = [G.sd_over_seeds(fn(int(k)), 0.10)[1] for k in ns]
        y = np.array([np.nan if not np.isfinite(q) else q for q in y])
        axs.plot(ns, y, color=col, lw=1.5, ls=ls, marker="o", ms=2.6, label=name)
    axs.axhline(G.GATE, color=GREY, ls=":", lw=1.0)
    axs.set_xscale("log"); axs.set_yscale("log")
    axs.set_xticks([4, 6, 8, 12, 16, 20, 30])
    axs.set_xticklabels(["4", "6", "8", "12", "16", "20", "30"])
    axs.set_yticks([0.18, 0.20, 0.25, 0.30])
    axs.set_yticklabels(["0.18", "0.20", "0.25", "0.30"])
    axs.minorticks_off()
    axs.text(0.98, G.GATE, " 0.20 target", transform=axs.get_yaxis_transform(),
             fontsize=6.5, color=GREY, ha="right", va="bottom")
    axs.set_xlabel("number of visits", fontsize=8)
    axs.set_ylabel(r"sd$(\ln\hat\theta)$, median of 8 seeds", fontsize=8)
    axs.set_title("what the weekday clinic costs\nload CV 10%, $k_e$ prior 0.8",
                  fontsize=8, pad=6, linespacing=1.5)
    axs.grid(True, which="major", alpha=0.45)
    axs.tick_params(labelsize=6.5)
    axs.legend(fontsize=6.2, frameon=False, loc="upper right")
    save(fig, "fig_E3_calendar")


if __name__ == "__main__":
    print(__doc__)
    settle, _ = part_a()
    part_b(settle)
    part_c(settle)
    figure(settle)
