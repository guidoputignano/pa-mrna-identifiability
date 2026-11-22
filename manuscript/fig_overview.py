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
defs = []


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
panel(BX, 0, BW, R1H, "b", "From biology to the model")
GRAD = {
    "gAmber": ("#fff3cf", "#f4c95d"), "gMito": ("#e9f7f2", "#bfe6d9"), "gCell": ("#fffdfa", "#f7f1e8"),
    "gVessel": ("#f9dcd6", "#efb3a8"), "gKidney": ("#d98f7c", "#b9604d"), "gPill": ("#fdebdc", "#f6c9a3"),
    "gFluid": ("#fde7a8", "#f2c14e"), "gLNP": ("#e3edfd", "#b7cdf6"),
}
for gid, (c0, c1) in GRAD.items():
    defs.append(f'<linearGradient id="{gid}" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="{c0}"/>'
                f'<stop offset="1" stop-color="{c1}"/></linearGradient>')
PCC_OWN, PCC_DEL = "#16866f", DEL
SMALL = 1.85

# --- sources, outside the cell
sx = BX + 3.0
src = [("amino acids", 15.0), ("odd-chain\nfatty acids", 24.0), ("gut bacteria", 34.5)]
# amino-acid chain
for j in range(4):
    circle(sx + 1.0 + j * 1.5, src[0][1] - 1.6, 0.6, "#e7b43c")
    if j:
        line(sx + 1.0 + (j - 1) * 1.5 + 0.6, src[0][1] - 1.6, sx + 1.0 + j * 1.5 - 0.6, src[0][1] - 1.6, "#c99a2e", 0.25)
# fatty acid zigzag with head
zz = " ".join(f"{sx + 1.6 + k * 1.0:.2f},{src[1][1] - 2.6 + (0.6 if k % 2 else 0):.2f}" for k in range(7))
add(f'<polyline points="{zz}" fill="none" stroke="#c99a2e" stroke-width="0.35" stroke-linejoin="round"/>')
circle(sx + 1.0, src[1][1] - 2.3, 0.7, "#e7b43c")
# bacterium
rect(sx + 0.8, src[2][1] - 2.6, 4.4, 1.8, "#e8d7a6", "#c99a2e", 0.25, rx=0.9)
path(f"M{sx + 5.2},{src[2][1] - 1.7} q0.6,-0.6 1.2,0 t1.2,0", stroke="#c99a2e", w=0.22)
for lab, yy0 in src:
    text(sx, yy0 + 0.6, lab, size=SMALL, color=INK)

# --- liver cell
cx0, cy0, cx1, cy1 = BX + 17.5, 9.8, BX + 68.0, 47.5
path(f"M{cx0 + 4},{cy0} H{cx1 - 5} Q{cx1},{cy0} {cx1},{cy0 + 5} V{cy1 - 5} Q{cx1},{cy1} {cx1 - 5},{cy1} "
     f"H{cx0 + 5} Q{cx0},{cy1} {cx0},{cy1 - 5} V{cy0 + 4} Q{cx0},{cy0} {cx0 + 4},{cy0} Z",
     fill="url(#gCell)", stroke="#d9cfc2", w=0.35)
text(cx0 + 2.2, cy0 + 3.6, "liver cell", size=SMALL, color=MUTED, italic=True)
# nucleus
circle(cx0 + 6.5, cy1 - 8.6, 3.4, "#eeeaf6", "#cfc8e2", 0.25)
circle(cx0 + 7.2, cy1 - 8.0, 1.0, "#d9d2ec")

# propionyl-CoA node
nx, ny = BX + 32.0, 27.5
circle(nx, ny, 5.0, "url(#gAmber)", MET, 0.35)
text(nx, ny - 0.3, "propionyl-", size=SMALL, anchor="middle", weight=700, color="#6b4e00")
text(nx, ny + 2.0, "CoA", size=SMALL, anchor="middle", weight=700, color="#6b4e00")
for _, yy0 in src:
    path(f"M{sx + 9.3},{yy0 - 1.2} C{sx + 12.5},{yy0 - 1.2} {nx - 8.5},{ny + (yy0 - ny) * 0.25} {nx - 5.6},{ny + (yy0 - ny) * 0.12}",
         stroke=MET, w=0.4)
    ang = math.atan2(ny + (yy0 - ny) * 0.12 - (ny + (yy0 - ny) * 0.25), (nx - 5.6) - (nx - 8.5))
    hx, hy = nx - 5.4, ny + (yy0 - ny) * 0.12
    add(f'<polygon points="{hx:.2f},{hy:.2f} {hx - 1.2 * math.cos(ang) - 0.5 * math.sin(ang):.2f},'
        f'{hy - 1.2 * math.sin(ang) + 0.5 * math.cos(ang):.2f} {hx - 1.2 * math.cos(ang) + 0.5 * math.sin(ang):.2f},'
        f'{hy - 1.2 * math.sin(ang) - 0.5 * math.cos(ang):.2f}" fill="{MET}"/>')

# mitochondrion with PCC (own = green, delivered = blue)
mx0, my0, mw, mh = BX + 44.0, 14.0, 21.0, 10.4
rect(mx0, my0, mw, mh, "url(#gMito)", ENZ, 0.4, rx=mh / 2)
rect(mx0 + 0.9, my0 + 0.9, mw - 1.8, mh - 1.8, "none", "#9ad3c2", 0.3, rx=(mh - 1.8) / 2)
for k in range(5):
    cxk = mx0 + 4.6 + k * 3.0
    if k % 2 == 0:
        path(f"M{cxk},{my0 + 0.9} V{my0 + mh * 0.58}", stroke="#9ad3c2", w=0.55)
    else:
        path(f"M{cxk},{my0 + mh - 0.9} V{my0 + mh * 0.42}", stroke="#9ad3c2", w=0.55)
blobs = [(PCC_OWN, 3.0, 6.6), (PCC_DEL, 6.1, 3.4), (PCC_DEL, 9.1, 6.8), (PCC_OWN, 12.1, 3.3), (PCC_DEL, 15.1, 6.9),
         (PCC_DEL, 17.9, 4.0)]
for col, dx, dy in blobs:
    circle(mx0 + dx, my0 + dy, 0.85, col, "#ffffff", 0.2)
text(mx0 + mw / 2, my0 - 1.2, "mitochondrion: PCC", size=SMALL, anchor="middle", weight=700, color=ENZ)
# own vs delivered legend
circle(mx0 + 1.2, my0 + mh + 2.2, 0.7, PCC_OWN)
text(mx0 + 2.3, my0 + mh + 2.85, "own θ", size=SMALL, color=PCC_OWN, weight=600)
circle(mx0 + 10.2, my0 + mh + 2.2, 0.7, PCC_DEL)
text(mx0 + 11.3, my0 + mh + 2.85, "delivered Δ", size=SMALL, color=PCC_DEL, weight=600)
arrow(nx + 4.4, ny - 2.6, mx0 - 0.3, my0 + mh * 0.62, ENZ, 0.6, 1.5)
text(mx0 + mw - 0.4, my0 + mh + 6.0, "→ energy (TCA cycle)", size=SMALL, color=ENZ, anchor="end")
# mRNA-LNP entering the cell, translated into the delivered enzyme
lx, ly = BX + 37.6, 8.4
circle(lx, ly, 2.0, "url(#gLNP)", DEL, 0.3)
for j in range(10):
    a = 2 * math.pi * j / 10
    circle(lx + 2.0 * math.cos(a), ly + 2.0 * math.sin(a), 0.32, DEL)
path(f"M{lx - 0.9},{ly} q0.45,-0.7 0.9,0 t0.9,0", stroke=DEL, w=0.3)
text(lx - 2.8, ly + 0.7, "mRNA-LNP", size=SMALL, color=DEL, anchor="end", weight=700)
path(f"M{lx + 1.6},{ly + 1.8} C{lx + 3.2},{ly + 4.0} {lx + 4.0},{ly + 4.6} {mx0 + 1.6},{my0 + 2.2}",
     stroke=DEL, w=0.35, extra='stroke-dasharray="0.8 0.6"')
# carnitine conjugation
px0, py0, pw, ph = BX + 45.0, 34.4, 20.5, 6.0
rect(px0, py0, pw, ph, "url(#gPill)", OTH, 0.4, rx=3.0)
text(px0 + pw / 2, py0 + 3.8, "propionylcarnitine", size=SMALL, anchor="middle", weight=700, color="#9c4a12")
arrow(nx + 4.0, ny + 3.0, px0 - 0.3, py0 + ph / 2, OTH, 0.6, 1.5)
text(nx + 2.0, ny + 11.0, "+ carnitine", size=SMALL, color="#9c4a12", weight=600)

# blood vessel with red cells
vx0, vy0, vx1, vh = BX + 3.0, 50.6, BX + 68.0, 5.0
rect(vx0, vy0, vx1 - vx0, vh, "url(#gVessel)", "#e19a8e", 0.25, rx=vh / 2)
for j in range(9):
    ex, ey = vx0 + 22.0 + j * 4.9, vy0 + (1.6 if j % 2 else 3.3)
    add(f'<ellipse cx="{ex:.2f}" cy="{ey:.2f}" rx="1.05" ry="0.6" fill="#d4544a" opacity="0.75"/>')
text(vx0 + 2.0, vy0 + 3.3, "blood", size=SMALL, color="#9b2b22", weight=700)
# marker into blood and the blood draw
path(f"M{nx - 1.5},{ny + 5.0} V{vy0 - 0.4}", stroke=MET, w=0.4, extra='stroke-dasharray="0.9 0.6"')
add(f'<polygon points="{nx - 1.5:.2f},{vy0 - 0.2:.2f} {nx - 2.1:.2f},{vy0 - 1.4:.2f} {nx - 0.9:.2f},{vy0 - 1.4:.2f}" fill="{MET}"/>')
text(nx + 0.6, cy1 - 3.6, "3-hydroxypropionate", size=SMALL, color="#8a6400", weight=600)
text(nx + 0.6, cy1 - 1.4, "(the marker)", size=SMALL, color=MUTED)
tube(nx - 1.5, vy0 + vh + 1.4, 0.9)
text(nx + 0.4, vy0 + vh + 4.6, "blood draw", size=SMALL, color=BLOOD, weight=600)
arrow(px0 + pw / 2 + 4.0, py0 + ph, px0 + pw / 2 + 4.0, vy0 - 0.2, OTH, 0.45, 1.2)
# kidney and urine
kx, ky = BX + 61.5, 62.6
path(f"M{kx},{ky - 3.3} C{kx + 3.1},{ky - 3.5} {kx + 3.5},{ky + 3.4} {kx},{ky + 3.3} C{kx - 1.8},{ky + 3.3} "
     f"{kx - 2.3},{ky + 1.7} {kx - 1.1},{ky + 0.8} C{kx - 0.5},{ky + 0.3} {kx - 0.5},{ky - 0.3} {kx - 1.1},{ky - 0.8} "
     f"C{kx - 2.3},{ky - 1.7} {kx - 1.8},{ky - 3.3} {kx},{ky - 3.3} Z", fill="url(#gKidney)", stroke="none")
path(f"M{kx - 0.6},{ky} C{kx - 2.0},{ky} {kx - 2.4},{ky + 1.6} {kx - 2.6},{ky + 4.6}", stroke="#b9604d", w=0.45)
arrow(kx + 1.2, vy0 + vh + 0.1, kx + 1.2, ky - 3.2, OTH, 0.4, 1.0)
text(kx - 3.6, ky - 1.0, "kidney → urine", size=SMALL, color="#9c4a12", anchor="end", weight=600)
text(kx - 3.6, ky + 1.3, "acylcarnitines compete with", size=SMALL, color=MUTED, anchor="end")
text(kx - 3.6, ky + 3.6, "carnitine for reabsorption", size=SMALL, color=MUTED, anchor="end")

# --- 'modelled as' connector
gx = BX + 71.2
line(gx, 11.0, gx, 66.0, "#d3d9e0", 0.3, dash="0.8 0.8")
circle(gx, 29.0, 2.6, "#ffffff", "#c6ced8", 0.3)
path(f"M{gx - 0.7},{27.6} l1.4,1.4 l-1.4,1.4", stroke=MUTED, w=0.45)
add(f'<text x="{gx - 1.2:.2f}" y="{44.0:.2f}" text-anchor="middle" transform="rotate(-90 {gx - 1.2:.2f} 44.0)" '
    f'style="font-size:{SMALL}px;fill:{MUTED};font-style:italic">modelled as</text>')

# --- the model (right)
tx0, tw, ty0, th_ = BX + 77.5, 15.0, 18.5, 21.0
text(BX + 75.0, 10.8, "one balance, two exits", size=SMALL, color=MUTED, italic=True)
arrow(tx0 + tw / 2, 12.4, tx0 + tw / 2, ty0 + 3.0, INK, 0.5, 1.4)
text(tx0 + tw / 2 + 1.4, 15.0, "load L(t) + L₀", size=SMALL, color=INK, weight=600)
fy = ty0 + 7.0
path(f"M{tx0 + 0.25},{fy} H{tx0 + tw - 0.25} V{ty0 + th_ - 2.2} Q{tx0 + tw - 0.25},{ty0 + th_ - 0.25} "
     f"{tx0 + tw - 2.2},{ty0 + th_ - 0.25} H{tx0 + 2.2} Q{tx0 + 0.25},{ty0 + th_ - 0.25} {tx0 + 0.25},{ty0 + th_ - 2.2} Z",
     fill="url(#gFluid)", stroke="none")
line(tx0 + 0.4, fy, tx0 + tw - 0.4, fy, "#e0a92a", 0.3)
path(f"M{tx0},{ty0} V{ty0 + th_ - 2.4} Q{tx0},{ty0 + th_} {tx0 + 2.4},{ty0 + th_} H{tx0 + tw - 2.4} "
     f"Q{tx0 + tw},{ty0 + th_} {tx0 + tw},{ty0 + th_ - 2.4} V{ty0}", stroke=INK, w=0.45)
text(tx0 + tw / 2, fy + 6.4, "m", size=4.2, anchor="middle", italic=True, weight=600, color="#6b4e00")
text(tx0 + tw / 2, fy + 10.2, "metabolite", size=SMALL, anchor="middle", color="#6b4e00")
# exits
for yy0, col, lab1, lab2 in ((fy + 2.6, ENZ, "enzyme", "(θ + Δ(t)) g(m)"), (fy + 9.6, OTH, "other route", "α(m)")):
    rect(tx0 + tw, yy0 - 0.9, 5.2, 1.8, col, rx=0.3)
    add(f'<polygon points="{tx0 + tw + 7.4:.2f},{yy0:.2f} {tx0 + tw + 5.0:.2f},{yy0 - 1.7:.2f} {tx0 + tw + 5.0:.2f},{yy0 + 1.7:.2f}" fill="{col}"/>')
    text(tx0 + tw + 8.4, yy0 - 0.5, lab1, size=FS, color=col, weight=700)
    text(tx0 + tw + 8.4, yy0 + 2.1, lab2, size=SMALL, color=col, weight=600)
# three hypotheses for alpha(m)
hy = 48.0
text(BX + 75.0, hy - 1.6, "three hypotheses for α(m)", size=SMALL, color=MUTED, italic=True)
shapes = [
    ("carnitine", "rises, then falls", lambda u: 2.4 * u * math.exp(-1.5 * u) / 0.65, ""),
    ("linear", "straight line", lambda u: 0.95 * u, 'stroke-dasharray="0.3 0.75"'),
    ("saturable", "levels off", lambda u: 1.25 * u / (0.28 + u), 'stroke-dasharray="1.1 0.7"'),
]
for i, (lab, sub, f, dash) in enumerate(shapes):
    ox, oy, ow, oh = BX + 75.0 + 14.6 * i, hy, 12.8, 9.0
    rect(ox, oy, ow, oh, "#ffffff", "#d3d9e0", 0.25, rx=1.2)
    line(ox + 1.2, oy + oh - 1.1, ox + ow - 1.0, oy + oh - 1.1, "#c6ced8", 0.2)
    line(ox + 1.2, oy + oh - 1.1, ox + 1.2, oy + 1.0, "#c6ced8", 0.2)
    pts = [f"{ox + 1.4 + (ow - 2.6) * (k / 40):.2f},{oy + oh - 1.3 - (oh - 2.6) * min(f(k / 40), 1.0):.2f}" for k in range(41)]
    path("M" + " L".join(pts), stroke=OTH, w=0.55, extra=dash)
    text(ox + ow / 2, oy + oh + 2.7, lab, size=SMALL, anchor="middle", weight=700)
    text(ox + ow / 2, oy + oh + 5.0, sub, size=SMALL, anchor="middle", color=MUTED)

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
       f'font-family="Inter, Helvetica, Arial, sans-serif"><defs>' + "".join(defs) + "</defs>" + "".join(parts) + "</svg>")
(OUT / "fig_overview.svg").write_text(svg)
html = (f'<!doctype html><html><head><meta charset="utf-8"><style>@page{{size:{W:g}mm {H:g}mm;margin:0}}'
        f'html,body{{margin:0;padding:0;background:#fff}}svg{{display:block}}</style></head><body>{svg}</body></html>')
(OUT / "fig_overview.html").write_text(html)

import subprocess  # noqa: E402

subprocess.run(["node", str(HERE / "render_svg.js"), str(OUT / "fig_overview.html"), str(OUT / "fig_overview.pdf"),
                str(OUT / "fig_overview.png"), f"{W:g}", f"{H:g}"], check=True)
print("wrote", OUT / "fig_overview.pdf")
