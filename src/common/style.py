"""Shared plotting style, so every figure in the paper looks the same."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# repository root, two levels up from src/common/
ROOT = Path(__file__).resolve().parents[2]
FIGDIR = ROOT / "results" / "figures"
TABDIR = ROOT / "results" / "tables"
FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

# Text block of the paper: A4 with 2.1cm margins, in inches.
TEXTWIDTH_IN = (21.0 - 2 * 2.1) / 2.54


def figsize(frac, aspect=0.85):
    """
    Figure size for a figure printed at `frac` of the paper's text block.

    Figures are generated at the size they are printed. Drawing one larger and
    letting includegraphics shrink it scales the type down with the axes: at
    frac = 0.49 an 8pt annotation reaches the page at under 5pt.
    """
    w = frac * TEXTWIDTH_IN
    return (w, w * aspect)


INK = "#1a1a1a"
MUTED = "#5a5a5a"
GREEN = "#1baf7a"
RED = "#e34948"
BLUE = "#2a78d6"
AMBER = "#e8a33d"
GREY = "#9a9890"

plt.rcParams.update({
    "figure.dpi": 130,
    # Figures are now drawn at printed size, so panels are narrow and titles
    # collide without automatic spacing.
    "figure.autolayout": True,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    # Embed TrueType rather than matplotlib's default Type 3. Type 3 fonts do
    # not scale, and most journals reject them outright; the paper's own body
    # text is Type 1 via lmodern, so without this the figures were the only
    # Type 3 left in the PDF.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#b8b6ae",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 1.8,
    "grid.color": "#e6e4dd",
    "grid.linewidth": 0.6,
})


def save(fig, name):
    """Save a figure into results/figures and report the path."""
    path = FIGDIR / f"{name}.pdf"
    fig.savefig(path)
    fig.savefig(FIGDIR / f"{name}.png")
    plt.close(fig)
    print(f"  wrote {path.relative_to(ROOT)}")
    return path


def write_table(name, header, rows, note=""):
    """Write a LaTeX booktabs table into results/tables."""
    path = TABDIR / f"{name}.tex"
    ncol = len(header)
    lines = [r"\begin{tabular}{@{}" + "l" * ncol + r"@{}}", r"\toprule",
             " & ".join(f"\\textbf{{{h}}}" for h in header) + r" \\",
             r"\midrule"]
    for r in rows:
        lines.append(" & ".join(str(c) for c in r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    if note:
        lines.append("% " + note)
    path.write_text("\n".join(lines))
    print(f"  wrote {path.relative_to(ROOT)}")
    return path


# ====================================================================
#  Annotation placement
#
#  Every label below used to be written in DATA coordinates. That is
#  fine until a parameter changes: refitting Model C moved the
#  separatrix from 470 to 1033 and six annotations, pinned to numbers
#  chosen against the old geometry, ended up on top of curves, legends
#  and each other. Nudging the numbers would only defer the problem to
#  the next refit.
#
#  The rule these helpers enforce: a label that describes a REGION goes
#  in axes fractions, which no parameter can move, and a label that
#  points at a DATA POINT is offset from that point in typographic
#  points, which no parameter can move either. Nothing is written in
#  data units unless the data itself is the thing being located.
# ====================================================================
def corner(ax, text, where="upper left", pad=0.035, **kw):
    """
    A label in one corner of the axes, positioned in axes fractions.

    `where` is the matplotlib legend vocabulary, so a caller can say
    where a label belongs without knowing what the axis limits are.
    """
    vy, hx = where.split()
    x = {"left": pad, "center": 0.5, "right": 1 - pad}[hx]
    y = {"upper": 1 - pad, "center": 0.5, "lower": pad}[vy]
    kw.setdefault("ha", {"left": "left", "center": "center", "right": "right"}[hx])
    kw.setdefault("va", {"upper": "top", "center": "center", "lower": "bottom"}[vy])
    kw.setdefault("fontsize", 8)
    return ax.text(x, y, text, transform=ax.transAxes, **kw)


def over_x(ax, x, text, frac=0.94, **kw):
    """
    A label at a given position on the x axis, at a fixed height in the
    axes box. Use it for shaded bands and vertical markers, whose x is
    data and whose y is only "near the top".
    """
    kw.setdefault("ha", "left")
    kw.setdefault("va", "top")
    kw.setdefault("fontsize", 8)
    return ax.text(x, frac, text, transform=ax.get_xaxis_transform(), **kw)


def beside_y(ax, y, text, frac=0.5, **kw):
    """The transpose of over_x: y is data, x is a fraction of the axes."""
    kw.setdefault("ha", "center")
    kw.setdefault("va", "bottom")
    kw.setdefault("fontsize", 8)
    return ax.text(frac, y, text, transform=ax.get_yaxis_transform(), **kw)


def at_point(ax, xy, text, dx=0, dy=0, arrow=False, color=None, **kw):
    """
    A label anchored to a data point and offset in POINTS.

    The offset is typographic, so it is the same distance on the page
    whatever the axis limits are, and it stays clear of the marker when
    the point moves.
    """
    kw.setdefault("fontsize", 8)
    if color is not None:
        kw.setdefault("color", color)
    if arrow:
        kw.setdefault("arrowprops",
                      dict(arrowstyle="->", lw=0.9,
                           color=kw.get("color", "0.35")))
    return ax.annotate(text, xy=xy, xytext=(dx, dy),
                       textcoords="offset points", **kw)
