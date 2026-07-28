"""
Reading-frame exon shapes with per-part domain fills (transcript row).

Each EXON_TABLE part is drawn as a solid rectangle at its exact absolute x
position, clipped to the puzzle-piece outline — matching domain bubble edges.
"""

from __future__ import annotations

import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath

from src.puzzle_geometry import exon_polygon_verts


def add_exon_patch(
    ax,
    x: float,
    width: float,
    y: float,
    height: float,
    left_shape: str,
    right_shape: str,
    parts: list[tuple[float, float, str, str]],
    lw: float = 0.9,
    *,
    bump_depth: float = 0.0,
) -> None:
    """Draw one transcript exon: puzzle outline + solid per-part fills."""
    verts = exon_polygon_verts(
        x, width, y, height, left_shape, right_shape,
        bump_depth=bump_depth,
    )
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(verts) - 2) + [MplPath.CLOSEPOLY]
    outline_path = MplPath(verts, codes)
    clip_patch = mpatches.PathPatch(outline_path, facecolor="none", edgecolor="none")
    ax.add_patch(clip_patch)

    for f0, f1, fill, _border in parts:
        seg = mpatches.Rectangle(
            (x + f0 * width, y),
            max(1e-6, (f1 - f0) * width),
            height,
            facecolor=fill,
            edgecolor="none",
            zorder=3,
        )
        seg.set_clip_path(clip_patch)
        ax.add_patch(seg)

    # Thin blend at internal part boundaries (dual-domain exons)
    for i in range(len(parts) - 1):
        f1 = parts[i][1]
        bx = x + f1 * width
        blend_w = min(width * 0.06, 0.08)
        left_c = parts[i][2]
        right_c = parts[i + 1][2]
        blend = mpatches.Rectangle(
            (bx - blend_w / 2, y), blend_w, height,
            facecolor=right_c, edgecolor="none", alpha=0.45, zorder=4,
        )
        blend.set_clip_path(clip_patch)
        ax.add_patch(blend)

    outline_color = parts[-1][3]
    ax.add_patch(mpatches.PathPatch(
        outline_path, facecolor="none", edgecolor=outline_color, lw=lw, zorder=5,
    ))
