"""
Figures 5a and 5b of the change report: the invariants A, B, C of the local
model, drawn from the repository's model code at the Table 2 constants.

5a  Delta(m) = A/m + B + C m decomposed, at theta = 3000 and at theta = 150,
    with a waterfall at m = 6 that builds Delta from its physical pieces.
5b  left:  the record in (m, Delta m) form, Delta m = A + B m + C m^2, with the
           25 monthly visits, four twins along equation (6), one non-twin, and
           the 14 d and 14.01 d designs as insets;
    right: the (K_m, theta) plane, the curved fibre and the straight locus of
           estimates under a constant diary error, tangent to it at the truth.

Every constant comes through model_abc.py. Run:
    python reports/changes_2026-09-20_to_2026-09-28/fig_abc.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_abc as M                                                  # noqa: E402
from common.style import (TEXTWIDTH_IN, INK, MUTED, GREEN, RED, BLUE,   # noqa: E402
                          AMBER, GREY)
import matplotlib.pyplot as plt                                        # noqa: E402

OUT = HERE / "fig"
#: colour by invariant: A, the three pieces of B, and C
COL = {"A": BLUE, "load": "#7a4fb5", "theta": RED, "keKm": "#c77dba",
       "C": AMBER, "total": INK}


# ------------------------------------------------------------------ 5a
def fig5a():
    fig, axes = plt.subplots(2, 2, figsize=(TEXTWIDTH_IN, 0.60 * TEXTWIDTH_IN),
                             gridspec_kw=dict(height_ratios=[1.45, 1.0]))
    m = np.linspace(1.5, 70, 600)
    for j, th in enumerate((M.THETA, M.THETA_DISEASE)):
        A, B, C = M.invariants(theta=th)
        ax = axes[0, j]
        ax.axhline(0, color=GREY, lw=0.6)
        lo, hi = M.level_of_delta(M.DELTA_PEAK, theta=th), M.level_of_delta(M.DELTA_TROUGH, theta=th)
        ax.axvspan(lo, hi, color=GREEN, alpha=0.13, lw=0)
        ax.text((lo + hi) / 2, 0.93, "dose train\n$\\Delta$ %.0f to %.0f" % (M.DELTA_TROUGH, M.DELTA_PEAK),
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7, color=GREEN)
        ax.plot(m, A / m, color=COL["A"], lw=1.2, label="$A/m$")
        ax.plot(m, np.full_like(m, B), color=COL["theta"], lw=1.2, ls="--", label="$B$")
        ax.plot(m, C * m, color=COL["C"], lw=1.2, label="$Cm$")
        ax.plot(m, A / m + B + C * m, color=COL["total"], lw=1.8, label="$\\Delta(m)$")
        mu = M.untreated_level(theta=th)
        ax.plot([mu], [0], "o", color=INK, ms=4, zorder=5)
        ax.annotate("untreated\n$m=%.1f$" % mu, (mu, 0), (mu + 4, 1900), fontsize=7,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
        ax.set_xlim(0, 70)
        ax.set_ylim(-4000, 9000)
        ax.set_xlabel("marker $m$  [$\\mu$M]")
        ax.set_title("$\\theta = %.0f$  (%s)" % (th, "normal capacity" if j == 0 else "5% of normal"),
                     fontsize=8.5)
        if j == 0:
            ax.set_ylabel("increment $\\Delta$  [$\\mu$M/d]")
            ax.legend(loc="upper right", fontsize=7, ncol=2, frameon=False)
        # the waterfall at m = 6
        axw = axes[1, j]
        mm = M.M_WATERFALL
        pcs = M.b_pieces(theta=th)
        steps = [("$+A/m$ (saturation cost)", A / mm, COL["A"]),
                 ("$+(L+b)$ (load)", pcs["load"], COL["load"]),
                 ("$-\\theta$ (own enzyme)", pcs["theta"], COL["theta"]),
                 ("$-k_eK_m$ (backup, fixed)", pcs["keKm"], COL["keKm"]),
                 ("$-k_e m$ (backup, $C m$)", C * mm, COL["C"])]
        x = 0.0
        for i, (lab, v, c) in enumerate(steps):
            y = len(steps) - i
            axw.barh(y, v, left=x, color=c, height=0.62)
            axw.text(max(x, x + v) + 120, y, "%s  %+.0f" % (lab, v), va="center", fontsize=6.6)
            x += v
        axw.barh(0, x, color=INK, height=0.62)
        axw.text(x + 120, 0, "$\\Delta(%.0f) = %.1f$" % (mm, x), va="center", fontsize=7)
        axw.axvline(0, color=GREY, lw=0.6)
        axw.set_yticks([])
        axw.set_xlim(-500, 11800)
        axw.set_xlabel("$\\mu$M/d, at $m = %.0f$" % mm)
        for s in ("left", "right", "top"):
            axw.spines[s].set_visible(False)
    fig.tight_layout(h_pad=0.8, w_pad=1.2)
    fig.savefig(OUT / "fig5a_abc.pdf")
    plt.close(fig)


# ------------------------------------------------------------------ 5b
def record(times, theta=M.THETA, Km=M.KM, b=0.0):
    """Noiseless record: the level each visit's increment holds the marker at.
    A twin generated at load L + b is fitted at the recorded L, so its levels
    are those of its own parameters at load L + b."""
    d = M.train_delta(times)
    m = M.level_of_delta(d, theta=theta, Km=Km, b=b)
    return m, d * m


def fig5b():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN, 0.43 * TEXTWIDTH_IN),
                                 gridspec_kw=dict(width_ratios=[1.15, 1.0]))
    # left: the parabola
    mm = np.linspace(2.5, 9.6, 300)
    styles = list(zip(M.TWIN_BS, (":", "-", "--", "-.")))
    for b, ls in styles:
        th, km = M.twin(b)
        A, B, C = M.invariants(theta=th, Km=km, b=b)
        a1.plot(mm, A + B * mm + C * mm ** 2, ls=ls, color=BLUE, lw=1.3 if b == 0 else 1.0,
                label="twin $b=%+d$: $\\theta=%.0f$, $K_m=%.0f$" % (b, th, km))
    t = M.schedule(30.4)
    m0, y0 = record(t)
    a1.plot(m0, y0, "o", color=INK, ms=3.2, zorder=5, label="%d monthly visits" % t.size)
    A, B, C = M.invariants(theta=M.THETA + M.NONTWIN_DTHETA)
    a1.plot(mm, A + B * mm + C * mm ** 2, color=RED, lw=1.0, label="non-twin: $\\theta%+.0f$ alone" % M.NONTWIN_DTHETA)
    a1.set_xlabel("marker $m$  [$\\mu$M]")
    a1.set_ylabel("$\\Delta\\cdot m$  [$\\mu$M$^2$/d]")
    a1.legend(fontsize=6.3, loc="lower left", frameon=False)
    a1.set_xlim(2.5, 9.6)
    a1.set_ylim(-9500, 31000)
    A0, B0, C0 = M.invariants()
    for k, (step, lab) in enumerate(((14.0, "14 d"), (14.01, "14.01 d"))):
        ins = a1.inset_axes([0.38 + 0.31 * k, 0.70, 0.27, 0.22])
        zm = np.linspace(3.0, 4.2, 50)
        ins.plot(zm, A0 + B0 * zm + C0 * zm ** 2, color=BLUE, lw=0.8)
        ts = M.schedule(step)
        mi, yi = record(ts)
        ins.plot(mi, yi, "o", color=INK, ms=2.2, alpha=0.6)
        ins.set_xlim(3.0, 4.2)
        ins.set_ylim(A0 + B0 * 4.2 + C0 * 4.2 ** 2, A0 + B0 * 3.0 + C0 * 9.0)
        ins.tick_params(labelsize=5.5, length=2, pad=1)
        ins.set_yticks([])
        ins.set_title("%s: %d visits" % (lab, ts.size), fontsize=6.3, pad=1.5)
    # right: the fibre and the estimate locus
    bb = np.linspace(-200, 800, 400)
    fth, fkm = M.twin(bb)
    a2.plot(fkm, fth, color=BLUE, lw=1.5, label="fibre, eq. (6): exact twins")
    eth, ekm = M.estimate_under_bias(bb)
    a2.plot(ekm, eth, color=RED, lw=1.2, ls="--", label="estimate under a diary error $b$")
    for b, _ in styles:
        th, km = M.twin(b)
        a2.plot(km, th, "o", color=BLUE, ms=3.5)
    bd = M.LOAD_FRAC * M.L
    hi, _ = M.estimate_under_bias(-bd)
    lo, _ = M.estimate_under_bias(bd)
    a2.axhspan(lo, hi, color=AMBER, alpha=0.35, lw=0,
               label="$\\pm%.1f$ diary bound: factor %.3f" % (bd, hi / lo))
    a2.plot([M.KM], [M.THETA], "*", color=INK, ms=8, zorder=6)
    a2.annotate("truth", (M.KM, M.THETA), (M.KM - 10, M.THETA - 420), fontsize=7)
    slope = -(M.L + M.KE * M.KM) / M.KM
    a2.text(0.97, 0.52, "tangent slope %.3f\n= null direction (%.1f, %.0f)"
            % (slope, M.L + M.KE * M.KM, -M.KM), transform=a2.transAxes, ha="right", fontsize=6.8)
    a2.set_xlabel("$K_m$  [$\\mu$M]")
    a2.set_ylabel("$\\theta$  [$\\mu$M/d]")
    a2.set_xlim(25, 280)
    a2.set_ylim(1700, 5000)
    a2.legend(fontsize=6.3, loc="upper right", frameon=False)
    fig.tight_layout(w_pad=1.2)
    fig.savefig(OUT / "fig5b_fibre.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig5a()
    fig5b()
    print("fig5a_abc.pdf, fig5b_fibre.pdf")
