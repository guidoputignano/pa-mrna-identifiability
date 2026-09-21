"""
Phase 3, the phase overlay: where in the dose cycle a schedule actually lands.

The design rule this paper turns on is invisible on a calendar axis. Sampling
every 14 days and every 14.5 days look alike on a timeline and differ by a
factor of twenty-seven in the precision of the answer. Folding the visit times
onto the dose cycle, t -> t mod 14, makes the difference the first thing a
reader sees.

WHY THIS IS AN OVERLAY AND NOT THE POLAR DIAL IT REPLACED. The dial plotted
the same fold in polar coordinates with RADIUS AS VISIT ORDER, and visit order
carries nothing about identifiability. Worse, it separated on the page points
that lie on one curve: because Delta is a deterministic function of phase,
every draw of every schedule falls exactly on the single decay curve drawn
here, and the dial's radius pushed them off it into a spiral that reads as a
trajectory. A reader seeing it for the first time could not tell that the
points were one set. On the phase axis the collapse is immediate: the
fortnightly schedule is twenty-five draws stacked on ONE POINT of the curve.

The dial could show one thing this cannot, that phase is circular, so a draw
at 13.9 d is adjacent to one at 0.1 d. The curve is periodic and the caption
says so. That is a small loss against a large gain.

TWO ROWS, AND THE SECOND IS THE ONE THE PAPER HAS NEVER DRAWN.

    top     Delta/theta against phase. Every point lies EXACTLY on the curve,
            because the increment is a deterministic function of phase. The
            only thing a schedule chooses is WHERE on the curve its draws sit.
    bottom  the marker m against phase. Here the points do NOT lie on the
            curve: they scatter around it by the day-to-day variation in the
            recorded load. That scatter is sigma, drawn rather than stated,
            and sigma is what Section 5.3 identifies as the binding constraint
            at nominal capacity.

The fortnightly panels carry the whole aliasing result between them. At
enrolment (t0 = 0) the dose train is still filling, so its draws smear over a
Delta/theta range of 0.13; for a patient ALREADY ON THERAPY the train is at
steady state and the range is EXACTLY ZERO. The dial showed both as one spoke
and could not tell them apart.

Writes results/figures/fig_E1_phaseclock.{pdf,png}.

Nothing here records an expected result in advance.

Run:  python audit/g1_phase_clock.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import matplotlib.pyplot as plt

from common.style import save, figsize, INK, MUTED, RED, GREEN, BLUE, AMBER
from design.run import delta, D_AMP, DOSE_INTERVAL, LAMBDA, LOAD_CV
import g_design_common as G

N = 25          # the accessible design's visit count, held across panels

#: panel titles. The registry names are written for a table and are too long
#: for five columns at text width; these say the same thing in the space there
#: actually is.
SHORT = {
    "14 d (the dose interval)": "14 d, the dose interval",
    "14 d, already on therapy": "14 d, on therapy already",
    "14.5 d": "14.5 d",
    "monthly (30.4 d)": "monthly, 30.4 d",
    "phase-spread": "phase-spread",
}
#: the steady-state peak of the dose train, which is the top of the curve
#: every panel is drawn against
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))


def _steady_curve(n=400):
    """
    The increment and the marker across one steady-state dose cycle.

    Drawn from the STEADY-STATE train rather than from the first cycle: it is
    the curve every schedule's late draws sit on, and a curve taken from the
    filling train would put the fortnightly panel's points visibly off it for
    a reason that has nothing to do with the schedule.
    """
    ph = np.linspace(0.0, G.T, n)
    inc = PEAK * np.exp(-LAMBDA * ph)
    return ph, inc / G.TH, G.level(G.TH, inc, G.LB)


def panel_increment(ax, times, name, curve):
    """Delta/theta against phase. Every draw lies exactly on the curve."""
    ph_c, d_c, _ = curve
    t = np.asarray(times, float)
    inc = np.asarray(delta(t), float)
    ax.plot(ph_c, d_c, color=MUTED, lw=1.2, zorder=1)
    ax.scatter(t % G.T, inc / G.TH, s=24, color=GREEN, edgecolors="white",
               linewidths=0.4, zorder=3, clip_on=False)
    nph = len(G.cycle_phases(t))
    rng = inc.max() / G.TH - inc.min() / G.TH
    #  "of N draws" is not decoration. The count is capped by the number of
    #  visits drawn, not by the schedule's phase lattice, and a header that
    #  says only "25 phases" invites a reader to think the lattice has 25.
    ax.set_title(f"{SHORT.get(name, name)}\n"
                 f"{nph} phase{'s' if nph != 1 else ''}, range {rng:.2f}",
                 fontsize=7.2, pad=4, linespacing=1.5)


def panel_marker(ax, times, curve, seed=3):
    """
    m against phase. The draws do NOT lie on the curve, and the gap is sigma.

    One intake record is drawn rather than a median over the usual
    twenty-four, because the object of this panel is the SCATTER of a single
    patient's record and a median over records would average exactly that
    away. The seed is fixed so the figure is reproducible.
    """
    ph_c, _, m_c = curve
    t = np.asarray(times, float)
    inc = np.asarray(delta(t), float)
    L = G.LB * (1.0 + LOAD_CV *
                np.random.default_rng(seed).standard_normal(t.size))
    m_obs, m_ref = G.level(G.TH, inc, L), G.level(G.TH, inc, G.LB)
    ax.plot(ph_c, m_c, color=MUTED, lw=1.2, zorder=1)
    ax.vlines(t % G.T, np.minimum(m_obs, m_ref), np.maximum(m_obs, m_ref),
              color=BLUE, lw=0.7, alpha=0.5, zorder=2)
    ax.scatter(t % G.T, m_obs, s=24, color=BLUE, edgecolors="white",
               linewidths=0.4, zorder=3, clip_on=False)


def figure():
    settle, worst = G.settling_days(tol=0.01)
    settle5, worst5 = G.settling_days(tol=0.05)
    print(f"  settling window: {settle*24:.2f} h at a 1% tolerance "
          f"(worst departure after it {100*worst:.3f}%),")
    print(f"                   {settle5*24:.2f} h at the repository's own 5% "
          f"(worst after {100*worst5:.3f}%).")

    designs = [(name, fn(N)) for name, fn in G.SCHEDULES.items()]
    curve = _steady_curve()

    fig, axes = plt.subplots(2, len(designs), figsize=figsize(1.0, 0.52),
                             sharex=True, sharey="row")
    for j, (name, t) in enumerate(designs):
        panel_increment(axes[0, j], t, name, curve)
        panel_marker(axes[1, j], t, curve)
        for ax in (axes[0, j], axes[1, j]):
            ax.set_xlim(-0.5, G.T + 0.5)
            ax.set_xticks([0, 7, 14])
            ax.grid(color="#eeece6", lw=0.5)
            ax.tick_params(labelsize=7)
            for sp in ax.spines.values():
                sp.set_color("#d8d5cc")
        #  the settling window, shaded on the phase axis rather than as a
        #  wedge. It is 0.4% of the cycle and the point of drawing it is
        #  that it removes the PEAK, not that it removes much time.
        for ax in (axes[0, j], axes[1, j]):
            ax.axvspan(0.0, settle, color=RED, alpha=0.13, lw=0, zorder=0)
        print(f"  {name:28} n={len(t):3d}  horizon {t[-1]:7.1f} d  "
              f"phases {len(G.cycle_phases(t)):3d}  "
              f"Delta/theta range {G.delta_range(t):.4f}")
    axes[0, 0].set_ylabel(r"increment  $\Delta/\theta$", fontsize=8.5)
    axes[1, 0].set_ylabel(r"marker  $m$  [$\mu$M]", fontsize=8.5)
    #  ONE shared x-label. Five copies is five times the ink for one fact and
    #  at text width they collide. tight_layout reserves no room for a figure
    #  label, so the room is reserved with `rect` and the text placed into it
    #  afterwards; supxlabel alone lands off the bottom of the canvas.
    fig.tight_layout(pad=0.6, rect=(0, 0.085, 1, 1))
    fig.text(0.5, 0.028, "phase of the dose cycle (days after a dose)",
             ha="center", va="bottom", fontsize=8.5)
    print(f"    settling window {100*settle/G.T:.1f}% of the cycle, "
          f"{settle*24:.1f} h at a 1% tolerance")
    save(fig, "fig_E1_phaseclock")


def numbers():
    print("\n  what the overlay is showing, at load CV 10%, 8 seeds, "
          "k_e log-prior 0.8, load record 5%:")
    print(f"    {'schedule':>28} {'n':>4} {'phases':>7} {'range':>8} "
          f"{'sd(ln theta): min':>18} {'median':>9} {'max':>9}")
    print("  " + "-" * 92)
    for name, fn in G.SCHEDULES.items():
        t = fn(N)
        lo, md, hi = G.sd_over_seeds(t, 0.10)
        f = lambda v: "singular" if not np.isfinite(v) else f"{v:.4f}"
        print(f"    {name:>28} {len(t):4d} {len(G.cycle_phases(t)):7d} "
              f"{G.delta_range(t):8.4f} {f(lo):>18} {f(md):>9} {f(hi):>9}")


if __name__ == "__main__":
    print(__doc__)
    figure()
    numbers()
