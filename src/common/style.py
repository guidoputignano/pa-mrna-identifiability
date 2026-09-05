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
