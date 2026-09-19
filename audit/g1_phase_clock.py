"""
Phase 3, the phase clock: where in the dose cycle a schedule actually lands.

The design rule this paper turns on is invisible on a calendar axis. Sampling
every 14 days and every 14.5 days look alike on a timeline and differ by a
factor of twenty in the precision of the answer. Folding the visit times onto
the dose cycle, t -> t mod 14, makes the difference the first thing a reader
sees: one schedule collapses to a single spoke, the other fills the dial.

    angle   2 pi (t mod 14)/14, dose at 12 o'clock, time running clockwise
    radius  visit order, so the record reads outwards
    colour  Delta/theta at that visit, which f4 found is the binding quantity
    shaded  the post-dose settling window, inside which a draw is not a
            steady state and does not belong in a steady-state fit

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
from design.run import delta
import g_design_common as G

N = 25          # the accessible design's visit count, held across panels


def panel(ax, times, title, settle, cmap, norm):
    inc = np.asarray(delta(times), dtype=float) / G.TH
    ph = np.asarray(times, float) % G.T
    ang = 2 * np.pi * ph / G.T
    rad = np.arange(1, len(times) + 1, dtype=float)

    # the settling window, a wedge from the dose instant
    ax.bar(x=[0.5 * 2 * np.pi * settle / G.T], height=[len(times) + 1.4],
           width=[2 * np.pi * settle / G.T], bottom=0.0,
           color=RED, alpha=0.16, edgecolor="none", zorder=0)

    ax.plot(ang, rad, color=MUTED, lw=0.5, alpha=0.45, zorder=1)
    ax.scatter(ang, rad, c=inc, cmap=cmap, norm=norm, s=17,
               edgecolors="white", linewidths=0.35, zorder=3)

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rticks([])
    ax.set_rlim(0, len(times) + 1.4)
    ax.set_xticks(np.arange(0, 2 * np.pi, 2 * np.pi / 7))
    ax.set_xticklabels([f"{d:g}" for d in np.arange(0, G.T, G.T / 7)],
                       fontsize=6.5, color=MUTED)
    ax.grid(color="#e6e4dd", lw=0.5)
    nph = len(G.cycle_phases(times))
    ax.set_title(f"{title}\n{nph} phase{'s' if nph != 1 else ''}, "
                 f"$\\Delta/\\theta$ range {G.delta_range(times):.2f}",
                 fontsize=7, pad=6, linespacing=1.5)


def figure():
    settle, worst = G.settling_days(tol=0.01)
    settle5, worst5 = G.settling_days(tol=0.05)
    print(f"  settling window: {settle*24:.2f} h at a 1% tolerance "
          f"(worst departure after it {100*worst:.3f}%),")
    print(f"                   {settle5*24:.2f} h at the repository's own 5% "
          f"(worst after {100*worst5:.3f}%).")
    print(f"  the shaded wedge is the 1% window, "
          f"{360*settle/G.T:.2f} degrees of the dial.\n")

    designs = [(name, fn(N)) for name, fn in G.SCHEDULES.items()]
    allinc = np.concatenate([np.asarray(delta(t), float) / G.TH
                             for _, t in designs])
    norm = plt.Normalize(allinc.min(), allinc.max())
    cmap = plt.get_cmap("viridis")

    fig, axes = plt.subplots(1, len(designs), figsize=figsize(1.0, 0.27),
                             subplot_kw={"projection": "polar"})
    fig.set_layout_engine("none")
    fig.subplots_adjust(left=0.015, right=0.985, top=0.80, bottom=0.30,
                        wspace=0.50)
    for ax, (name, t) in zip(axes, designs):
        panel(ax, t, name, settle, cmap, norm)
        print(f"  {name:28} n={len(t):3d}  horizon {t[-1]:7.1f} d  "
              f"phases {len(G.cycle_phases(t)):3d}  "
              f"Delta/theta range {G.delta_range(t):.4f}")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax = fig.add_axes([0.34, 0.035, 0.32, 0.045])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.ax.set_title(r"$\Delta/\theta$ at the visit", fontsize=7.5, pad=3)
    cb.ax.tick_params(labelsize=6.5, length=2, pad=1)
    #  the reading of this figure is the paper's caption, not a block
    #  of small type inside the axes: a figure that has to be annotated
    #  to be read is a figure that will be read wrongly when the caption
    #  is what gets quoted.
    print(f"    settling window {360*settle/G.T:.1f} deg of the dial, "
          f"{settle*24:.1f} h at a 1% tolerance")
    save(fig, "fig_E1_phaseclock")


def numbers():
    print("\n  what the dial is showing, at load CV 10%, 8 seeds, "
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
