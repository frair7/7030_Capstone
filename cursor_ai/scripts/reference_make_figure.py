"""
Reproduces the DMD/Dp427m exon-map + patient-deletion figure.

Exon "puzzle piece" shapes
---------------------------
dp427m_phase_shapes.csv encodes, for every exon, the shape of its
5' (left) and 3' (right) edge: Flat / Round / Point / N/A -- the
classic jigsaw convention for reading-frame phase at each exon
junction. Each exon's right edge gets a convex bump (Round ->
semicircular, Point -> triangular) and the following exon's left
edge gets the matching concave notch, so consecutive exons
interlock like puzzle pieces. Flat / N/A edges are plain straight
lines.

Everything else (domain map, hinge positions, patient deletion
bars, %Dys(WB) heat-map, MW column, phenotype boxes) is reproduced
from the published figure layout. Domain exon boundaries are
approximate -- edit DOMAIN_MAP if you have the precise coordinates.
"""

"""
Reproduces the DMD/Dp427m exon-map + patient-deletion figure.

Exon "puzzle piece" shapes
---------------------------
Each exon's 5' (left) and 3' (right) edge shape (Flat / Round / Point / N/A)
follows the classic jigsaw convention for reading-frame phase at each exon
junction. Each exon's right edge gets a convex bump (Round -> semicircular,
Point -> triangular) and the following exon's left edge gets the matching
concave notch, so consecutive exons interlock like puzzle pieces. Flat / N/A
edges are plain straight lines.

Per-exon domain coloring
-------------------------
EXON_TABLE below is derived from the authoritative Ensembl REST API exon
table (GRCh38, ENST00000357033.9) plus the supplied domain/feature nucleotide
coordinates. Exon width is proportional to true CDS coding length (not mRNA
length including UTR). Some exons contain a domain junction partway through
(e.g. exon 8 is part Actin Binding Domain / part H1; exon 17 is
R3 -> H2 -> R4) -- `parts` gives each exon's 1-3 fractional [f0,f1) sub
segments with the exact domain fill/border color for that segment, verified
to match the given per-exon "# of Parts" and color columns exactly for all
79 exons.

Everything else (domain map, hinge positions, patient deletion bars,
%Dys(WB) heat-map, MW column, phenotype boxes) is reproduced from the
published figure layout.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path

# ------------------------------------------------------------------
# 1. Authoritative per-exon table (shapes, true CDS bp, domain sub-segments)
# ------------------------------------------------------------------
EXON_TABLE = [
    {'n': 1, 'bp': 31, 'five_prime': 'N/A', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 2, 'bp': 62, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 3, 'bp': 93, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 4, 'bp': 78, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 5, 'bp': 93, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 6, 'bp': 173, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 7, 'bp': 119, 'five_prime': 'Point', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#FFCCCC', '#EB6F6F')]},
    {'n': 8, 'bp': 182, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 0.3901, '#FFCCCC', '#EB6F6F'), (0.3901, 1.0000, '#AEDCCA', '#4BAF89')]},
    {'n': 9, 'bp': 129, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#AEDCCA', '#4BAF89')]},
    {'n': 10, 'bp': 189, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 0.2540, '#AEDCCA', '#4BAF89'), (0.2540, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 11, 'bp': 182, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 12, 'bp': 151, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 13, 'bp': 120, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 14, 'bp': 102, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 15, 'bp': 108, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 16, 'bp': 180, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 17, 'bp': 176, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 0.0511, '#FFFFCC', '#F0C266'), (0.0511, 0.9034, '#AEDCCA', '#4BAF89'), (0.9034, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 18, 'bp': 124, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 19, 'bp': 88, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 20, 'bp': 242, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 21, 'bp': 181, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 22, 'bp': 146, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 23, 'bp': 213, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 24, 'bp': 114, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 25, 'bp': 156, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 26, 'bp': 171, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 27, 'bp': 183, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 28, 'bp': 135, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 29, 'bp': 150, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 30, 'bp': 162, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 31, 'bp': 111, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 32, 'bp': 174, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 33, 'bp': 156, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 34, 'bp': 171, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 35, 'bp': 180, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 36, 'bp': 129, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 37, 'bp': 171, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 38, 'bp': 123, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 39, 'bp': 138, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 40, 'bp': 153, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 41, 'bp': 183, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 42, 'bp': 195, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 43, 'bp': 173, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 44, 'bp': 148, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 45, 'bp': 176, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 46, 'bp': 148, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 47, 'bp': 150, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 48, 'bp': 186, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 49, 'bp': 102, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 50, 'bp': 109, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 0.6330, '#FFFFCC', '#F0C266'), (0.6330, 1.0000, '#AEDCCA', '#4BAF89')]},
    {'n': 51, 'bp': 233, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 0.4335, '#AEDCCA', '#4BAF89'), (0.4335, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 52, 'bp': 118, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 53, 'bp': 212, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 54, 'bp': 155, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 55, 'bp': 190, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 56, 'bp': 173, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 57, 'bp': 157, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 58, 'bp': 121, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 59, 'bp': 269, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 60, 'bp': 147, 'five_prime': 'Flat', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#FFFFCC', '#F0C266')]},
    {'n': 61, 'bp': 79, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 0.4557, '#FFFFCC', '#F0C266'), (0.4557, 1.0000, '#AEDCCA', '#4BAF89')]},
    {'n': 62, 'bp': 61, 'five_prime': 'Round', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#AEDCCA', '#4BAF89')]},
    {'n': 63, 'bp': 62, 'five_prime': 'Point', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#AEDCCA', '#4BAF89')]},
    {'n': 64, 'bp': 75, 'five_prime': 'Round', 'three_prime': 'Round', 'parts': [(0.0000, 0.6667, '#AEDCCA', '#4BAF89'), (0.6667, 1.0000, '#DAE9F8', '#4D93D9')]},
    {'n': 65, 'bp': 202, 'five_prime': 'Round', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#DAE9F8', '#4D93D9')]},
    {'n': 66, 'bp': 86, 'five_prime': 'Point', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#DAE9F8', '#4D93D9')]},
    {'n': 67, 'bp': 158, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#DAE9F8', '#4D93D9')]},
    {'n': 68, 'bp': 167, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#DAE9F8', '#4D93D9')]},
    {'n': 69, 'bp': 112, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#DAE9F8', '#4D93D9')]},
    {'n': 70, 'bp': 137, 'five_prime': 'Flat', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 71, 'bp': 39, 'five_prime': 'Point', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 72, 'bp': 66, 'five_prime': 'Point', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 73, 'bp': 66, 'five_prime': 'Point', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 74, 'bp': 159, 'five_prime': 'Point', 'three_prime': 'Point', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 75, 'bp': 244, 'five_prime': 'Point', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 76, 'bp': 124, 'five_prime': 'Flat', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 77, 'bp': 93, 'five_prime': 'Round', 'three_prime': 'Round', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 78, 'bp': 32, 'five_prime': 'Round', 'three_prime': 'Flat', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
    {'n': 79, 'bp': 12, 'five_prime': 'Flat', 'three_prime': 'N/A', 'parts': [(0.0000, 1.0000, '#F0E3F4', '#D86DCD')]},
]
N_EXONS = len(EXON_TABLE)
exon_shapes = {e['n']: (e['five_prime'], e['three_prime']) for e in EXON_TABLE}

# ------------------------------------------------------------------
# 2. Domain map (for the summary Domain row -- whole-exon ranges)
# ------------------------------------------------------------------
# Exon ranges derived from the authoritative nucleotide-coordinate domain/feature
# table (5'UTR = 244 nt, CDS starts at nt 245), mapped onto exon boundaries via
# the transcript's exon base-pair lengths, then resolved sequentially so hinges
# get priority as slim distinct segments at exon junctions (no gaps/overlaps
# across all 79 exons). Colors are taken directly from the supplied table.
DOMAIN_MAP = [
    (1, 8, "Actin Binding Domain", "#FFCCCC"),
    (9, 10, "H1", "#AEDCCA"),
    (11, 12, "R1", "#FFFFCC"),
    (13, 14, "R2", "#FFFFCC"),
    (15, 17, "R3", "#FFFFCC"),
    (18, 18, "H2", "#AEDCCA"),
    (19, 20, "R4", "#FFFFCC"),
    (21, 21, "R5", "#FFFFCC"),
    (22, 23, "R6", "#FFFFCC"),
    (24, 26, "R7", "#FFFFCC"),
    (27, 28, "R8", "#FFFFCC"),
    (29, 30, "R9", "#FFFFCC"),
    (31, 32, "R10", "#FFFFCC"),
    (33, 34, "R11", "#FFFFCC"),
    (35, 36, "R12", "#FFFFCC"),
    (37, 38, "R13", "#FFFFCC"),
    (39, 40, "R14", "#FFFFCC"),
    (41, 41, "R15", "#FFFFCC"),
    (42, 44, "R16", "#FFFFCC"),
    (45, 46, "R17", "#FFFFCC"),
    (47, 48, "R18", "#FFFFCC"),
    (49, 50, "R19", "#FFFFCC"),
    (51, 51, "H3", "#AEDCCA"),
    (52, 53, "R20", "#FFFFCC"),
    (54, 55, "R21", "#FFFFCC"),
    (56, 57, "R22", "#FFFFCC"),
    (58, 59, "R23", "#FFFFCC"),
    (60, 61, "R24", "#FFFFCC"),
    (62, 64, "H4", "#AEDCCA"),
    (65, 68, "CR Domain", "#DAE9F8"),
    (69, 79, "C-terminal Domain", "#F0E3F4"),
]


# ------------------------------------------------------------------
# 3. Geometry helpers: build a "puzzle piece" exon polygon
# ------------------------------------------------------------------
BUMP_FRAC = 0.22  # bump depth as a fraction of the exon's own width


def _side_points(x, y0, y1, shape, direction, bump):
    if shape not in ("Round", "Point"):
        return [(x, y0), (x, y1)]
    ymid = (y0 + y1) / 2.0
    h = (y1 - y0) * 0.28
    dx = direction * bump
    if shape == "Point":
        return [(x, y0), (x, ymid - h), (x + dx, ymid), (x, ymid + h), (x, y1)]
    pts = [(x, y0), (x, ymid - h)]
    n = 8
    for i in range(1, n):
        theta = np.pi * i / n
        px = x + dx * np.sin(theta)
        py = (ymid - h) + (2 * h) * (i / n)
        pts.append((px, py))
    pts += [(x, ymid + h), (x, y1)]
    return pts


def exon_polygon(x, width, y, height, left_shape, right_shape):
    x0, x1 = x, x + width
    y0, y1 = y, y + height
    bump = BUMP_FRAC * width
    left_pts = _side_points(x0, y0, y1, left_shape, direction=+1, bump=bump)
    right_pts = _side_points(x1, y1, y0, right_shape, direction=+1, bump=bump)
    verts = left_pts + right_pts
    verts.append(verts[0])
    return verts


def add_exon_patch(ax, x, width, y, height, left_shape, right_shape, parts, lw=0.8):
    """Draw one exon as its 1-3 domain-colored sub-segments (from `parts`,
    fractions of the exon's own width), clipped to the jigsaw puzzle-piece
    outline, with a single unified border on top."""
    verts = exon_polygon(x, width, y, height, left_shape, right_shape)
    codes = [Path.MOVETO] + [Path.LINETO] * (len(verts) - 2) + [Path.CLOSEPOLY]
    outline_path = Path(verts, codes)
    clip_patch = mpatches.PathPatch(outline_path, facecolor="none", edgecolor="none")
    ax.add_patch(clip_patch)

    for f0, f1, fill, _border in parts:
        seg = mpatches.Rectangle((x + f0 * width, y), max(1e-6, (f1 - f0) * width), height,
                                  facecolor=fill, edgecolor="none", zorder=3)
        seg.set_clip_path(clip_patch)
        ax.add_patch(seg)

    # thin dividers between sub-parts
    for f0, f1, fill, border in parts[1:]:
        dx = x + f0 * width
        ax.plot([dx, dx], [y, y + height], color=border, lw=0.9, zorder=4, alpha=0.8)

    outline_color = parts[-1][3]
    ax.add_patch(mpatches.PathPatch(outline_path, facecolor="none", edgecolor=outline_color,
                                     lw=lw, zorder=5))


# ------------------------------------------------------------------
# 4. Patient / deletion data
# ------------------------------------------------------------------
patients = [
    ("A", "P-1", 423.3, 12.42, 3, 3, "point", ""),
    ("A", "P-2", 423.3, 12.04, 5, 5, "point", ""),
    ("A", "P-3", 423.3, 6.72, 5, 5, "point", "\u2021"),
    ("A", "P-4", 420.4, 20.28, 3, 4, "point", ""),

    ("B", "P-5", 284.6, 55.66, 3, 27, "gray", ""),
    ("B", "P-6", 273.4, 49.16, 3, 29, "black", ""),
    ("B", "P-7", 256.3, 0.26, 3, 32, "black", ""),
    ("B", "P-8", 375.6, 0.92, 5, 13, "black", ""),
    ("B", "P-9", 211.6, 12.65, 6, 41, "gray", ""),

    ("C", "P-10", 375.6, 49.27, 10, 18, "gray", ""),
    ("C", "P-11", 227.2, 78.08, 10, 42, "gray", ""),
    ("C", "P-12", 214.1, 266.60, 10, 44, "gray", "\u2021"),
    ("C", "P-13", 234.6, 300.70, 12, 43, "gray", ""),
    ("C", "P-14", 395.6, 135.60, 13, 18, "gray", ""),
    ("C", "P-15", 326.6, 334.00, 13, 29, "gray", ""),
    ("C", "P-16", 284.9, 266.00, 13, 36, "gray", ""),
    ("C", "P-17", 238.6, 122.90, 17, 44, "gray", ""),
    ("C", "P-18", 235.4, 113.90, 20, 50, "gray", ""),
    ("C", "P-19", 378.6, 155.60, 30, 37, "white", ""),

    ("D", "P-20", 408.7, 28.84, 45, 47, "gray", ""),
    ("D", "P-21", 408.7, 19.87, 45, 47, "gray", ""),
    ("D", "P-22", 408.7, 37.20, 45, 47, "gray", ""),
    ("D", "P-23", 401.4, 74.57, 45, 48, "gray", ""),
    ("D", "P-24", 384.9, 81.17, 45, 51, "gray", ""),
    ("D", "P-25", 371.9, 69.83, 45, 53, "gray", ""),
    ("D", "P-26", 371.9, 23.77, 45, 53, "gray", ""),

    ("E", "P-27", 403.0, 76.63, 48, 51, "white", "*"),
    ("E", "P-28", 410.2, 75.99, 49, 51, "white", ""),
    ("E", "P-29", 387.5, 165.70, 50, 55, "gray", ""),
]

# ------------------------------------------------------------------
# 5. Layout constants
# ------------------------------------------------------------------
ROW_H = 1.0
EXON_ROW_Y = 0
DOMAIN_ROW_Y = -1.3
TABLE_TOP_Y = -3.0
PATIENT_ROW_H = 0.85
PATIENT_GAP = 0.0
ROW_STEP = PATIENT_ROW_H + PATIENT_GAP

LEFT_MARGIN = 6.0

# Table columns, positioned to sit directly left of the transcript track
PCT_X1 = LEFT_MARGIN - 0.4
PCT_X0 = PCT_X1 - 1.8
MW_X1 = PCT_X0
MW_X0 = MW_X1 - 1.4
PART_X1 = MW_X0
PART_X0 = PART_X1 - 2.2
GROUP_X1 = PART_X0
GROUP_X0 = GROUP_X1 - 0.9
TABLE_X0 = GROUP_X0

n_patients = len(patients)
table_bottom_y = TABLE_TOP_Y - n_patients * ROW_STEP

fig_width = LEFT_MARGIN + N_EXONS * 1.0 + 4
fig_height_data = table_bottom_y - 1.5

fig, ax = plt.subplots(figsize=(28, 15))
ax.set_xlim(TABLE_X0 - 0.5, fig_width)
ax.set_ylim(fig_height_data, 3.7)
ax.axis("off")

# ------------------------------------------------------------------
# 6. Exon track -- widths proportional to true CDS coding length, with a
# minimum width floor so short exons (e.g. 71-74) stay readable, matching
# the interactive HTML version's layout algorithm.
# ------------------------------------------------------------------
TOTAL_TRACK_WIDTH = N_EXONS * 1.0
MIN_EXON_W = 0.35

total_bp = sum(e["bp"] for e in EXON_TABLE)
raw_w = [e["bp"] / total_bp * TOTAL_TRACK_WIDTH for e in EXON_TABLE]
is_fixed = [w < MIN_EXON_W for w in raw_w]
fixed_total = sum(MIN_EXON_W for f in is_fixed if f)
raw_large_total = sum(w for w, f in zip(raw_w, is_fixed) if not f)
remaining = max(10.0, TOTAL_TRACK_WIDTH - fixed_total)
exon_widths = [MIN_EXON_W if f else (w / raw_large_total * remaining)
               for w, f in zip(raw_w, is_fixed)]

exon_x = {}
x = LEFT_MARGIN
for e, w in zip(EXON_TABLE, exon_widths):
    n = e["n"]
    add_exon_patch(ax, x, w, EXON_ROW_Y, ROW_H, e["five_prime"], e["three_prime"], e["parts"])
    fs = 5 if w >= 0.7 else 4
    ax.text(x + w / 2, EXON_ROW_Y + ROW_H / 2, str(n), ha="center", va="center",
            fontsize=fs, zorder=6)
    exon_x[n] = (x, x + w)
    x += w
transcript_right = x

ax.text(LEFT_MARGIN - 0.3, EXON_ROW_Y + ROW_H / 2, "Transcript \u2192",
        ha="right", va="center", fontsize=9, fontweight="bold")

for label, exon in [("Dp427m,b,c,l", 1), ("Dp260", 30), ("Dp140", 56), ("Dp71", 62)]:
    xa = exon_x[exon][0]
    ax.annotate(label, xy=(xa, EXON_ROW_Y + ROW_H + 0.03),
                xytext=(xa, EXON_ROW_Y + ROW_H + 0.55),
                ha="center", fontsize=8,
                arrowprops=dict(arrowstyle="-|>", lw=1))

# ------------------------------------------------------------------
# 7. Domain row
# ------------------------------------------------------------------
ax.text(LEFT_MARGIN - 0.3, DOMAIN_ROW_Y + ROW_H / 2, "Domain \u2192",
        ha="right", va="center", fontsize=9, fontweight="bold")

for s, e, label, color in DOMAIN_MAP:
    x0, x1 = exon_x[s][0], exon_x[e][1]
    ax.add_patch(mpatches.FancyBboxPatch(
        (x0, DOMAIN_ROW_Y), x1 - x0, ROW_H,
        boxstyle="round,pad=0,rounding_size=0.12",
        facecolor=color, edgecolor="black", lw=0.8, zorder=3))
    fs = 8 if (x1 - x0) > 1.5 else 6
    ax.text((x0 + x1) / 2, DOMAIN_ROW_Y + ROW_H / 2, label, ha="center", va="center",
            fontsize=fs, fontweight="bold", zorder=4)

# vertical dotted guide lines at H2 / H3 boundaries, through the whole chart
for exon_bound in (17, 18, 50, 51):
    # boundaries around the (now single-exon) H2 (18) and H3 (51) hinges
    xg = exon_x[exon_bound][0]
    ax.plot([xg, xg], [DOMAIN_ROW_Y, table_bottom_y], linestyle=":", color="gray",
            lw=0.8, zorder=1)

# ------------------------------------------------------------------
# 8. Patient table (Group | Participant No. | MW | %Dys(WB) | deletion map)
# ------------------------------------------------------------------
cmap = plt.get_cmap("RdYlGn")
vmin, vmax = 0.0, 340.0


def pct_color(v):
    return cmap(min(max(v, vmin), vmax) / vmax)


STYLE_FACE = {"gray": "#9a9a9a", "black": "black", "white": "white", "point": "#c9c9c9"}

# header row
header_h = 0.9
header_y = TABLE_TOP_Y + 0.05
col_headers = [
    (GROUP_X0, GROUP_X1, "Group"),
    (PART_X0, PART_X1, "Participant No."),
    (MW_X0, MW_X1, "MW"),
    (PCT_X0, PCT_X1, "%Dys(WB)"),
]
for x0, x1, label in col_headers:
    ax.add_patch(mpatches.Rectangle((x0, header_y), x1 - x0, header_h,
                                     facecolor="white", edgecolor="black", lw=1.0, zorder=3))
    ax.text((x0 + x1) / 2, header_y + header_h / 2, label, ha="center", va="center",
            fontsize=8.5, fontweight="bold", zorder=4)

# body rows
y = TABLE_TOP_Y
group_bounds = {}
for grp, pid, mw, pct, s, e, style, note in patients:
    group_bounds.setdefault(grp, [y, y])
    group_bounds[grp][1] = y - PATIENT_ROW_H

    # Participant No. cell
    ax.add_patch(mpatches.Rectangle((PART_X0, y - PATIENT_ROW_H), PART_X1 - PART_X0,
                                     PATIENT_ROW_H, facecolor="white", edgecolor="black",
                                     lw=0.6, zorder=2))
    ax.text((PART_X0 + PART_X1) / 2, y - PATIENT_ROW_H / 2, pid, ha="center", va="center",
            fontsize=8)

    # MW cell
    ax.add_patch(mpatches.Rectangle((MW_X0, y - PATIENT_ROW_H), MW_X1 - MW_X0,
                                     PATIENT_ROW_H, facecolor="white", edgecolor="black",
                                     lw=0.6, zorder=2))
    ax.text((MW_X0 + MW_X1) / 2, y - PATIENT_ROW_H / 2, f"{mw:.1f}", ha="center",
            va="center", fontsize=8)

    # %Dys(WB) heat cell
    ax.add_patch(mpatches.Rectangle((PCT_X0, y - PATIENT_ROW_H), PCT_X1 - PCT_X0,
                                     PATIENT_ROW_H, facecolor=pct_color(pct),
                                     edgecolor="black", lw=0.6, zorder=2))
    ax.text((PCT_X0 + PCT_X1) / 2, y - PATIENT_ROW_H / 2, f"{pct:.2f}", ha="center",
            va="center", fontsize=7.5, fontweight="bold")

    # deletion bar aligned to transcript exon coordinates
    x0b, x1b = exon_x[s][0], exon_x[e][1]
    if style == "point":
        ax.add_patch(mpatches.Rectangle((x0b, y - PATIENT_ROW_H + 0.08), x1b - x0b,
                                         PATIENT_ROW_H - 0.16, facecolor=STYLE_FACE["point"],
                                         edgecolor="black", lw=0.8, zorder=3))
        lbl = f"{s}{note}" if s == e else f"{s}-{e}"
        ax.text((x0b + x1b) / 2, y - PATIENT_ROW_H / 2, lbl, ha="center", va="center",
                fontsize=7, fontweight="bold", zorder=4)
    else:
        face = STYLE_FACE[style]
        text_color = "white" if style == "black" else "black"
        ax.add_patch(mpatches.Rectangle((x0b, y - PATIENT_ROW_H + 0.05), x1b - x0b,
                                         PATIENT_ROW_H - 0.10, facecolor=face,
                                         edgecolor="black", lw=0.8, zorder=3))
        ax.text((x0b + x1b) / 2, y - PATIENT_ROW_H / 2, f"{s}-{e}{note}", ha="center",
                va="center", fontsize=7, color=text_color, fontweight="bold", zorder=4)

    y -= PATIENT_ROW_H

# Group column: one merged, shaded cell per group with a bold letter
for grp, (ytop, ybot) in group_bounds.items():
    ax.add_patch(mpatches.Rectangle((GROUP_X0, ybot), GROUP_X1 - GROUP_X0, ytop - ybot,
                                     facecolor="#d9d9d9", edgecolor="black", lw=1.0, zorder=2))
    ax.text((GROUP_X0 + GROUP_X1) / 2, (ytop + ybot) / 2, grp, ha="center", va="center",
            fontsize=15, fontweight="bold", zorder=4)
    # thicker separator line above each group block (except the very first)
    ax.plot([TABLE_X0, transcript_right], [ytop, ytop], color="black", lw=1.4, zorder=5)

# outer border for the whole table (group+participant+MW+%Dys columns)
ax.add_patch(mpatches.Rectangle((TABLE_X0, table_bottom_y), PCT_X1 - TABLE_X0,
                                 (header_y + header_h) - table_bottom_y,
                                 facecolor="none", edgecolor="black", lw=1.6, zorder=6))

# ------------------------------------------------------------------
# 9. Phenotype boxes for group D / E (right-hand side)
# ------------------------------------------------------------------
d_top, d_bot = group_bounds["D"]
pheno_x0 = transcript_right + 1.0
pheno_x1 = pheno_x0 + 13.0
pheno_labels = ["Asymptomatic/Paucisymptomatic*", "BMD/IMD\u2021", "DMD"]
n = len(pheno_labels)
pheno_h = (d_top - d_bot) / n
ax.add_patch(mpatches.Rectangle((pheno_x0, d_bot), pheno_x1 - pheno_x0, d_top - d_bot,
                                 facecolor="none", edgecolor="black", lw=1.4, zorder=6))
yy = d_top
for lbl in pheno_labels:
    ax.plot([pheno_x0, pheno_x1], [yy - pheno_h, yy - pheno_h], color="black", lw=0.8)
    ax.text((pheno_x0 + pheno_x1) / 2, yy - pheno_h / 2, lbl, ha="center", va="center",
            fontsize=9)
    yy -= pheno_h

plt.tight_layout()
plt.savefig("data/dp427m_figure.png", dpi=200, bbox_inches="tight")
print("Saved figure.")
