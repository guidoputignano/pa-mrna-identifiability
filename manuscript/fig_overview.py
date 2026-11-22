"""Figure 1 (overview), drawn as an SVG at print size (180 mm wide) and rendered to PDF/PNG.

The figure is hand-built SVG so that strokes, fonts and spacing are controlled exactly; it is
rendered with the pre-installed Chromium through Playwright. Units in the SVG are millimetres.
Schematic only: every statement on it is a result stated in the text.

Run:  python3 fig_overview.py      -> fig/fig_overview.svg, .html, .pdf, .png
"""
from pathlib import Path
import math

HERE = Path(__file__).resolve().parent
OUT = HERE / "fig"
OUT.mkdir(exist_ok=True)

W, H = 180.0, 134.0
# palette (one meaning per colour, used consistently in every panel)
INK, MUTED, LINE, PANEL = "#1f2933", "#5f6b7a", "#aab4bf", "#f3f5f8"
ENZ, ENZ_L = "#16866f", "#d6efe9"          # enzyme exit (own + delivered)
DEL = "#2f6fdf"                            # delivered mRNA / infusions
OTH, OTH_L = "#e07a2f", "#fbe3d1"          # other (carnitine/urine) route
MET, MET_L = "#d9a21b", "#fbefc9"          # metabolite pool
RED, GOOD = "#cf3f3f", "#2e9e5b"
BLOOD, BLOOD_L = "#c0392b", "#f6d6d2"

FS = 2.25      # body text (mm), about 6.4 pt
FS_S = 1.95    # small text, about 5.5 pt
FS_T = 2.6     # panel titles, about 7.4 pt

parts = []


def add(s):
    parts.append(s)


def text(x, y, s, size=FS, color=INK, anchor="start", weight=400, italic=False, lh=1.22):
    """Multi-line text; lines separated by '\n'. y is the baseline of the first line."""
    style = f'font-size:{size}px;fill:{color};font-weight:{weight};' + ("font-style:italic;" if italic else "")
    lines = s.split("\n")
    tsp = "".join(
        f'<tspan x="{x:.2f}" dy="{0 if i == 0 else size * lh:.2f}">{ln}</tspan>' for i, ln in enumerate(lines))
    add(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" style="{style}">{tsp}</text>')


def line(x1, y1, x2, y2, color=INK, w=0.3, dash=None, cap="round"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" '
        f'stroke-width="{w}" stroke-linecap="{cap}"{d}/>')


def path(d, fill="none", stroke=INK, w=0.3, extra=""):
    add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{w}" stroke-linejoin="round" '
        f'stroke-linecap="round" {extra}/>')


def rect(x, y, w, h, fill, stroke="none", sw=0.3, rx=0.0):
    add(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"/>')


def circle(cx, cy, r, fill, stroke="none", sw=0.3):
    add(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def arrow(x1, y1, x2, y2, color=INK, w=0.35, head=1.4):
    """Straight arrow with a filled triangular head drawn as a polygon (prints reliably)."""
    ang = math.atan2(y2 - y1, x2 - x1)
    bx, by = x2 - head * math.cos(ang), y2 - head * math.sin(ang)
    line(x1, y1, bx, by, color, w, cap="butt")
    px, py = -math.sin(ang), math.cos(ang)
    hw = max(head * 0.45, w * 0.9)
    add(f'<polygon points="{x2:.2f},{y2:.2f} {bx + px * hw:.2f},{by + py * hw:.2f} '
        f'{bx - px * hw:.2f},{by - py * hw:.2f}" fill="{color}"/>')


def drop(cx, cy, s=1.0, color=DEL):
    path(f"M{cx:.2f},{cy - 1.9 * s:.2f} C{cx + 0.2 * s:.2f},{cy - 1.2 * s:.2f} {cx + 1.1 * s:.2f},{cy - 0.6 * s:.2f} "
         f"{cx + 1.1 * s:.2f},{cy + 0.25 * s:.2f} A{1.1 * s:.2f},{1.1 * s:.2f} 0 1 1 {cx - 1.1 * s:.2f},{cy + 0.25 * s:.2f} "
         f"C{cx - 1.1 * s:.2f},{cy - 0.6 * s:.2f} {cx - 0.2 * s:.2f},{cy - 1.2 * s:.2f} {cx:.2f},{cy - 1.9 * s:.2f} Z",
         fill=color, stroke="none")


def tube(cx, y, s=1.0):
    """Blood collection tube, top at y."""
    w, h = 1.5 * s, 5.0 * s
    rect(cx - w / 2 - 0.25 * s, y - 0.7 * s, w + 0.5 * s, 0.9 * s, "#8a96a3", rx=0.2)
    path(f"M{cx - w / 2:.2f},{y:.2f} V{y + h - w / 2:.2f} A{w / 2:.2f},{w / 2:.2f} 0 0 0 {cx + w / 2:.2f},{y + h - w / 2:.2f} "
         f"V{y:.2f} Z", fill="#ffffff", stroke="#8a96a3", w=0.22)
    path(f"M{cx - w / 2 + 0.15:.2f},{y + h * 0.42:.2f} V{y + h - w / 2:.2f} A{w / 2 - 0.15:.2f},{w / 2 - 0.15:.2f} 0 0 0 "
         f"{cx + w / 2 - 0.15:.2f},{y + h - w / 2:.2f} V{y + h * 0.42:.2f} Z", fill=BLOOD, stroke="none")


def badge(cx, cy, kind, r=1.9):
    """Result badge: ok (tick), no (cross), model ('M')."""
    if kind == "ok":
        circle(cx, cy, r, GOOD)
        path(f"M{cx - 0.85:.2f},{cy + 0.05:.2f} L{cx - 0.2:.2f},{cy + 0.7:.2f} L{cx + 0.95:.2f},{cy - 0.7:.2f}",
             stroke="#ffffff", w=0.45)
    elif kind == "no":
        circle(cx, cy, r, RED)
        line(cx - 0.7, cy - 0.7, cx + 0.7, cy + 0.7, "#ffffff", 0.45)
        line(cx - 0.7, cy + 0.7, cx + 0.7, cy - 0.7, "#ffffff", 0.45)
    else:
        circle(cx, cy, r, "#ffffff", MUTED, 0.3)
        text(cx, cy + 0.75, "M", size=2.1, color=MUTED, anchor="middle", weight=700)


def panel(x, y, w, h, letter, title):
    rect(x, y, w, h, PANEL, rx=1.6)
    text(x + 2.6, y + 5.0, letter, size=3.3, weight=700)
    text(x + 6.3, y + 4.9, title, size=FS_T, weight=600)


# ------------------------------------------------------------------------------------------
# layout
G = 2.6
R1H, R2H = 71.0, H - 71.0 - G
AX, AW = 0.0, 57.0
BX, BW = AW + G, W - AW - G
SW = (W - 3 * G) / 4
CX, DX, EX, FX = 0.0, SW + G, 2 * (SW + G), 3 * (SW + G)
R2Y = R1H + G

# ---------------------------------------------------------------------- (a) therapy and data
panel(AX, 0, AW, R1H, "a", "Therapy and routine data")
bx, by = AX + 7.0, 11.0
line(bx, by - 1.6, bx, by + 21.5, "#9aa5b1", 0.35)                     # pole
line(bx - 2.2, by - 1.6, bx + 2.2, by - 1.6, "#9aa5b1", 0.35)
rect(bx - 3.0, by, 6.0, 8.6, "#e8f0fd", rx=1.0)
rect(bx - 3.0, by + 3.6, 6.0, 5.0, "#bcd2f7", rx=0.0)
rect(bx - 3.0, by, 6.0, 8.6, "none", DEL, 0.3, rx=1.0)
text(bx, by + 2.6, "mRNA", size=1.75, color=DEL, anchor="middle", weight=700)
text(bx, by + 6.9, "LNP", size=1.6, color=DEL, anchor="middle", weight=600)
path(f"M{bx},{by + 8.6} C{bx},{by + 13} {bx + 6},{by + 12.5} {bx + 10.2},{by + 11.3}", stroke=DEL, w=0.3)
px, py = AX + 21.0, 10.2
circle(px, py + 2.6, 2.5, "#f1c9a5")
path(f"M{px - 4.0},{py + 15.5} V{py + 9.6} C{px - 4.0},{py + 6.7} {px - 2.3},{py + 5.6} {px},{py + 5.6} "
     f"C{px + 2.3},{py + 5.6} {px + 4.0},{py + 6.7} {px + 4.0},{py + 9.6} V{py + 15.5} Z", fill="#3b6fd1", stroke="none")
lx, ly = AX + 41.5, 15.0
path(f"M{lx - 9.0},{ly - 0.6} C{lx - 8.4},{ly - 4.2} {lx - 2.0},{ly - 4.6} {lx + 3.0},{ly - 3.8} "
     f"C{lx + 7.2},{ly - 3.1} {lx + 9.4},{ly - 2.2} {lx + 8.6},{ly - 0.4} C{lx + 7.6},{ly + 1.8} {lx + 3.2},{ly + 3.6} "
     f"{lx - 1.4},{ly + 4.4} C{lx - 5.6},{ly + 5.0} {lx - 9.6},{ly + 2.8} {lx - 9.0},{ly - 0.6} Z",
     fill="#b9573b", stroke="none")
path(f"M{lx - 1.0},{ly - 4.2} C{lx - 0.4},{ly - 1.0} {lx - 0.8},{ly + 1.6} {lx - 2.0},{ly + 4.3}",
     stroke="#8f3f28", w=0.25)
path(f"M{px + 2.4},{py + 10.0} C{px + 6},{py + 9.5} {lx - 12},{ly + 3.5} {lx - 9.4},{ly + 1.2}",
     stroke=MUTED, w=0.22, extra='stroke-dasharray="0.6 0.6"')
text(lx, ly + 8.6, "the liver makes PCC", size=FS, anchor="middle")
text(lx, ly + 12.0, "own capacity θ", size=FS, color=ENZ, anchor="middle", weight=600)
text(lx, ly + 14.9, "+ delivered enzyme Δ(t)", size=FS, color=DEL, anchor="middle", weight=600)
ty = 47.0
x0, x1 = AX + 4.0, AX + AW - 4.0
line(x0, ty, x1, ty, INK, 0.3)
days = 56.0
for d in range(0, 57, 7):
    xx = x0 + (x1 - x0) * d / days
    line(xx, ty, xx, ty + 0.9, INK, 0.25)
for d in (0, 14, 28, 42, 56):
    xx = x0 + (x1 - x0) * d / days
    drop(xx, ty - 2.6, 0.95)
text(x0, ty - 6.6, "infusion every 2 weeks, each one dated", size=FS, color=DEL, weight=600)
for d in (5, 19, 33, 47):
    xx = x0 + (x1 - x0) * d / days
    tube(xx, ty + 2.6, 0.95)
text(x0, ty + 12.6, "blood draws (metabolite level)", size=FS, color=BLOOD, weight=600)
text(x0, ty + 15.7, "and a record of intake as eaten", size=FS, color=MUTED)
text(x1, ty + 3.6, "days", size=FS_S, color=MUTED, anchor="end")

# ---------------------------------------------------------------------- (b) biology -> model
panel(BX, 0, BW, R1H, "b", "From biology to model: one balance, two exits")
cx0, cy0, cw, ch = BX + 4.0, 12.0, 58.0, 33.0
text(cx0 + 0.6, cy0 - 1.4, "liver cell", size=FS_S, color=MUTED, italic=True)
rect(cx0, cy0, cw, ch, "#ffffff", "#c6ced8", 0.3, rx=3.0)
text(cx0 + 2.5, cy0 + 5.3, "protein (amino acids),\nodd-chain fatty acids,\ngut bacteria", size=FS_S, color=INK)
pcx, pcy = cx0 + 21.0, cy0 + 18.0
add(f'<ellipse cx="{pcx}" cy="{pcy}" rx="8.6" ry="5.0" fill="{MET_L}" stroke="{MET}" stroke-width="0.35"/>')
text(pcx, pcy - 0.4, "propionyl-CoA", size=FS, anchor="middle", weight=600)
text(cx0 + 2.5, pcy + 8.4, "toxic when\nit builds up", size=FS_S, color=MUTED, italic=True)
arrow(cx0 + 9.0, cy0 + 11.5, pcx - 5.0, pcy - 4.2, MET, 0.4, 1.3)
mx, my = cx0 + 44.5, cy0 + 9.5
add(f'<ellipse cx="{mx}" cy="{my}" rx="11.2" ry="5.6" fill="{ENZ_L}" stroke="{ENZ}" stroke-width="0.35"/>')
text(mx, my + 2.4, "own θ + delivered Δ", size=1.8, anchor="middle", color=ENZ)
text(mx, my - 0.9, "PCC", size=FS, anchor="middle", weight=700, color=ENZ)
text(mx, my + 8.2, "mitochondrion", size=FS_S, anchor="middle", color=MUTED, italic=True)
arrow(pcx + 8.0, pcy - 2.0, mx - 10.6, my + 1.8, ENZ, 0.55, 1.5)
text(mx + 4.0, my - 6.1, "→ energy", size=FS_S, color=ENZ, anchor="start")
ccx, ccy = cx0 + 44.0, cy0 + 25.5
rect(ccx - 12.0, ccy - 3.4, 24.0, 6.8, OTH_L, OTH, 0.35, rx=3.2)
text(ccx, ccy - 0.3, "propionylcarnitine", size=FS, anchor="middle", weight=600, color="#9c4a12")
text(ccx, ccy + 2.4, "(+ carnitine)", size=FS_S, anchor="middle", color="#9c4a12")
arrow(pcx + 7.4, pcy + 2.8, ccx - 12.2, ccy - 0.6, OTH, 0.55, 1.5)
vy = cy0 + ch + 2.6
rect(cx0, vy, cw, 4.2, BLOOD_L, rx=2.1)
text(cx0 + 2.0, vy + 2.95, "blood: metabolite marker (3-hydroxypropionate)", size=FS_S, color=BLOOD)
arrow(pcx - 1.0, pcy + 5.0, pcx - 1.0, vy - 0.2, MET, 0.35, 1.1)
arrow(ccx + 4.0, ccy + 3.4, ccx + 4.0, vy - 0.2, OTH, 0.35, 1.1)
kx, ky = cx0 + cw - 5.0, vy + 9.6
path(f"M{kx},{ky - 3.3} C{kx + 3.1},{ky - 3.5} {kx + 3.5},{ky + 3.4} {kx},{ky + 3.3} C{kx - 1.8},{ky + 3.3} "
     f"{kx - 2.3},{ky + 1.7} {kx - 1.1},{ky + 0.8} C{kx - 0.5},{ky + 0.3} {kx - 0.5},{ky - 0.3} {kx - 1.1},{ky - 0.8} "
     f"C{kx - 2.3},{ky - 1.7} {kx - 1.8},{ky - 3.3} {kx},{ky - 3.3} Z", fill="#c97a6a", stroke="none")
path(f"M{kx - 0.6},{ky} C{kx - 2.0},{ky} {kx - 2.6},{ky + 1.2} {kx - 3.4},{ky + 3.6}", stroke="#c97a6a", w=0.4)
arrow(kx + 0.6, vy + 4.3, kx + 0.6, ky - 3.6, OTH, 0.35, 1.0)
text(kx - 4.2, ky + 0.9, "kidney → urine", size=FS_S, color="#9c4a12", anchor="end")
text(cx0 + 1.0, ky + 5.6, "in the kidney, acylcarnitines compete with carnitine for reabsorption", size=FS_S, color=MUTED)
gx = BX + 64.4
path(f"M{gx},{25.5} l1.8,2.2 l-1.8,2.2", stroke=MUTED, w=0.45)
text(gx + 0.9, 33.6, "model", size=FS_S, color=MUTED, anchor="middle", italic=True)
tx0 = BX + 82.0
tw, th_ = 19.0, 23.0
tyy = 14.5
path(f"M{tx0},{tyy} V{tyy + th_ - 2.5} Q{tx0},{tyy + th_} {tx0 + 2.5},{tyy + th_} H{tx0 + tw - 2.5} "
     f"Q{tx0 + tw},{tyy + th_} {tx0 + tw},{tyy + th_ - 2.5} V{tyy}", fill="#ffffff", stroke=INK, w=0.4)
path(f"M{tx0 + 0.2},{tyy + 9.0} H{tx0 + tw - 0.2} V{tyy + th_ - 2.5} Q{tx0 + tw - 0.2},{tyy + th_ - 0.2} "
     f"{tx0 + tw - 2.5},{tyy + th_ - 0.2} H{tx0 + 2.5} Q{tx0 + 0.2},{tyy + th_ - 0.2} {tx0 + 0.2},{tyy + th_ - 2.5} Z",
     fill=MET_L, stroke="none")
text(tx0 + tw / 2, tyy + 15.6, "metabolite", size=FS, anchor="middle", weight=600)
text(tx0 + tw / 2, tyy + 18.7, "level m", size=FS, anchor="middle", weight=600)
arrow(tx0 + tw / 2, tyy - 7.0, tx0 + tw / 2, tyy + 4.5, INK, 0.5, 1.5)
text(tx0 + tw / 2 + 1.5, tyy - 4.2, "load L(t) + L₀", size=FS, anchor="start")
arrow(tx0 + tw, tyy + th_ - 5.0, tx0 + tw + 8.0, tyy + th_ - 5.0, ENZ, 1.1, 2.2)
text(tx0 + tw + 1.0, tyy + th_ - 12.2, "enzyme", size=FS, color=ENZ, weight=700)
text(tx0 + tw + 1.0, tyy + th_ - 9.1, "(θ + Δ) g(m)", size=FS, color=ENZ, weight=600)
arrow(tx0, tyy + th_ - 5.0, tx0 - 8.0, tyy + th_ - 5.0, OTH, 1.1, 2.2)
text(tx0 - 1.0, tyy + th_ - 12.2, "other route", size=FS, color=OTH, weight=700, anchor="end")
text(tx0 - 1.0, tyy + th_ - 9.1, "α(m)", size=FS, color=OTH, weight=600, anchor="end")
hy = 48.0
text(BX + 91.5, hy - 1.8, "three hypotheses for the other route α(m)", size=FS_S, color=MUTED, anchor="middle")
shapes = [
    ("carnitine", "rises, then falls", lambda u: 2.4 * u * math.exp(-1.5 * u) / 0.65),
    ("linear", "straight line", lambda u: 0.95 * u),
    ("saturable", "levels off", lambda u: 1.25 * u / (0.28 + u)),
]
for i, (lab, sub, f) in enumerate(shapes):
    ox, oy, ow, oh = BX + 70.5 + 15.8 * i, hy, 12.6, 9.0
    rect(ox, oy, ow, oh, "#ffffff", "#c6ced8", 0.25, rx=0.8)
    pts = []
    for k in range(41):
        u = k / 40
        v = min(f(u), 1.0)
        pts.append(f"{ox + 1.0 + (ow - 2.0) * u:.2f},{oy + oh - 1.0 - (oh - 2.4) * v:.2f}")
    path("M" + " L".join(pts), stroke=OTH, w=0.5,
         extra='stroke-dasharray="1.0 0.7"' if i == 2 else ("" if i == 0 else 'stroke-dasharray="0.25 0.7"'))
    text(ox + ow / 2, oy + oh + 2.6, lab, size=FS_S, anchor="middle", weight=600)
    text(ox + ow / 2, oy + oh + 4.9, sub, size=FS_S, anchor="middle", color=MUTED)

# ---------------------------------------------------------------------- (c) the experiment
panel(CX, R2Y, SW, R2H, "c", "Infusions as the experiment")
text(CX + SW / 2, R2Y + 9.6, "where blood draws fall in the 14-day cycle", size=FS_S, color=MUTED, anchor="middle")
for k, (ox, spread) in enumerate(((CX + 11.0, False), (CX + 31.6, True))):
    oy, r = R2Y + 26.0, 8.2
    circle(ox, oy, r, "#ffffff", "#c6ced8", 0.3)
    for d in range(14):
        a = math.radians(90 - d * 360 / 14)
        line(ox + (r - 0.9) * math.cos(a), oy - (r - 0.9) * math.sin(a),
             ox + r * math.cos(a), oy - r * math.sin(a), LINE, 0.25)
    drop(ox, oy - r - 2.1, 0.85)
    if spread:
        for j in range(13):
            a = math.radians(90 - j * 360 / 13 - 8)
            circle(ox + (r - 2.6) * math.cos(a), oy - (r - 2.6) * math.sin(a), 0.75, BLOOD)
        badge(ox + r * 0.86, oy - r * 0.86, "ok", 1.6)
    else:
        circle(ox, oy - r + 2.6, 0.75, BLOOD)
        text(ox, oy - r + 6.4, "all 13", size=FS_S, color=BLOOD, weight=600, anchor="middle")
        badge(ox + r * 0.86, oy - r * 0.86, "no", 1.6)
text(CX + 11.0, R2Y + 38.6, "every 14 days:\nalways the same\npoint", size=FS_S, anchor="middle")
text(CX + 31.6, R2Y + 38.6, "every 14.5 days:\nthe whole cycle\nis covered", size=FS_S, anchor="middle")
text(CX + SW / 2, R2Y + 52.6, "Book draws off the dose day;\nrecord every dose and the intake.", size=FS_S,
     anchor="middle", weight=600)

# ---------------------------------------------------------------------- (d) capacity in severe PA
panel(DX, R2Y, SW, R2H, "d", "Capacity in severe PA")
badge(DX + SW - 3.4, R2Y + 3.8, "no", 1.6)
text(DX + 3.0, R2Y + 10.4, "who clears the load (5% of normal)", size=FS_S, color=MUTED)
bx0, bw0, byy = DX + 3.0, SW - 6.0, R2Y + 12.6
rect(bx0, byy, bw0 * 0.84, 4.6, OTH, rx=0.6)
rect(bx0 + bw0 * 0.84 + 0.25, byy, bw0 * 0.16 - 0.25, 4.6, ENZ, rx=0.6)
text(bx0 + 1.2, byy + 3.1, "other route 84%", size=FS, color="#ffffff", weight=700)
text(bx0 + bw0, byy + 7.6, "own enzyme 16%", size=FS_S, color=ENZ, anchor="end", weight=700)
text(DX + SW / 2, R2Y + 26.0, "a 1% error in the large share\nis a 5% error in the small one", size=FS_S,
     anchor="middle")
gx0, gy0, gw, gh = DX + 5.0, R2Y + 46.5, SW - 9.0, 12.0
line(gx0, gy0, gx0 + gw, gy0, MUTED, 0.25)
line(gx0, gy0, gx0, gy0 - gh, MUTED, 0.25)
pts = []
for k in range(61):
    u = k / 60
    v = 0.02 if u < 0.62 else min(((u - 0.62) / 0.3) ** 2, 1.0)
    pts.append(f"{gx0 + gw * u:.2f},{gy0 - 0.5 - (gh - 1.5) * v:.2f}")
path("M" + " L".join(pts), stroke=INK, w=0.5)
text(gx0 + 1.2, gy0 - 2.0, "no lower bound", size=FS_S, color=MUTED)
text(gx0 + gw * 0.62 - 1.0, gy0 - gh + 1.8, "upper bound", size=FS_S, color=MUTED, anchor="end")
add(f'<text x="{gx0 - 1.2:.2f}" y="{gy0 - gh / 2:.2f}" text-anchor="middle" transform="rotate(-90 {gx0 - 1.2:.2f} '
    f'{gy0 - gh / 2:.2f})" style="font-size:{FS_S}px;fill:{MUTED}">fit worse</text>')
text(gx0 + gw, gy0 + 2.6, "candidate capacity", size=FS_S, color=MUTED, anchor="end")
text(DX + SW / 2, R2Y + 55.4, "only an upper bound on θ", size=FS_S, anchor="middle", weight=600)

# ---------------------------------------------------------------------- (e) therapeutic window
panel(EX, R2Y, SW, R2H, "e", "The therapeutic window")
ex0, ew, ey_top, ey_bot = EX + 3.0, SW - 6.0, R2Y + 14.0, R2Y + 37.0
crit = 0.07


def yy(v):
    return ey_bot - (ey_bot - ey_top) * v


rect(ex0, yy(1.0), ew, yy(crit) - yy(1.0), ENZ_L)
line(ex0, yy(1.0), ex0 + ew, yy(1.0), RED, 0.35, dash="1.0 0.7")
line(ex0, yy(crit), ex0 + ew, yy(crit), INK, 0.35, dash="1.0 0.7")
text(ex0 + ew, yy(1.0) - 1.1, "ceiling: never above normal", size=FS_S, color=RED, anchor="end")
text(ex0, yy(crit) + 2.6, "critical capacity", size=FS_S, color=INK)
lam = math.log(2) / 4.3
ts = [i * 42 / 300 for i in range(301)]
vals = [sum(math.exp(-lam * (t - p)) for p in (0, 14, 28) if t >= p) for t in ts]
vmax = max(vals)
pts = [f"{ex0 + ew * t / 42:.2f},{yy(v / vmax * 0.985):.2f}" for t, v in zip(ts, vals)]
path("M" + " L".join(pts), stroke=ENZ, w=0.55)
tr = vals[-1] / vmax * 0.985
arrow(ex0 + ew - 1.0, yy(tr) + 4.6, ex0 + ew - 1.0, yy(tr) + 0.9, INK, 0.3, 0.9)
text(ex0 + ew - 2.0, yy(tr) + 5.0, "trough", size=FS_S, anchor="end")
text(ex0 + 1.0, ey_bot + 6.4, "delivered enzyme over three infusions", size=FS_S, color=ENZ)
text(EX + SW / 2, R2Y + 50.0, "safe for every patient if the trough\nclears the line: this needs the\nhalf-life and the dose, not θ",
     size=FS_S, anchor="middle", weight=600)

# ---------------------------------------------------------------------- (f) what the record identifies
panel(FX, R2Y, SW, R2H, "f", "What the record identifies")
rows = [("ok", "the enzyme half-life,", "from dated infusions"),
        ("ok", "change within one patient,", "if the marker is stable"),
        ("no", "absolute capacity", "in severe PA"),
        ("model", "how late an infusion can be:", "a model result")]
for i, (k, a1, a2) in enumerate(rows):
    ry = R2Y + 13.0 + 9.4 * i
    badge(FX + 5.0, ry, k, 1.8)
    text(FX + 8.4, ry - 0.5, a1, size=FS_S, weight=600)
    text(FX + 8.4, ry + 2.1, a2, size=FS_S, color=MUTED)
text(FX + SW / 2, R2Y + 55.4, "plus a short data specification", size=FS_S, color=MUTED, anchor="middle",
     italic=True)

# ------------------------------------------------------------------------------------------
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:g} {H:g}" width="{W:g}mm" height="{H:g}mm" '
       f'font-family="Inter, Helvetica, Arial, sans-serif">' + "".join(parts) + "</svg>")
(OUT / "fig_overview.svg").write_text(svg)
html = (f'<!doctype html><html><head><meta charset="utf-8"><style>@page{{size:{W:g}mm {H:g}mm;margin:0}}'
        f'html,body{{margin:0;padding:0;background:#fff}}svg{{display:block}}</style></head><body>{svg}</body></html>')
(OUT / "fig_overview.html").write_text(html)

import subprocess  # noqa: E402

subprocess.run(["node", str(HERE / "render_svg.js"), str(OUT / "fig_overview.html"), str(OUT / "fig_overview.pdf"),
                str(OUT / "fig_overview.png"), f"{W:g}", f"{H:g}"], check=True)
print("wrote", OUT / "fig_overview.pdf")
