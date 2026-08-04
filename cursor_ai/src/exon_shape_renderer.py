"""
Transcript row rendering: domain-aligned colour strips with puzzle outlines on top.

Colours and borders use absolute x-coordinates shared with the protein/domain
row (coding-position alignment). Puzzle-piece geometry is independent of colour.
"""

from __future__ import annotations

from typing import Any

import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Bbox

from src.exon_gradients import build_horizontal_gradient
from src.part_segments import (
    PartSegment,
    build_part_segments,
    exon_draw_xrange,
    part_border_color,
    segments_for_exon,
    segments_in_xrange,
)
from src.puzzle_geometry import calculate_bump_depth, exon_polygon_verts

EXON_H_REF = 26.0


def draw_transcript_color_band(
    ax,
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
    widths: list[float],
    y: float,
    height: float,
    *,
    zorder: int = 2,
) -> tuple[list[PartSegment], list[mpatches.PathPatch]]:
    """Domain-aligned vertical strips with puzzle-clipped bump fills."""
    segments = build_part_segments(exon_table, exon_x)
    puzzle_clips: list[mpatches.PathPatch] = []
    n_exons = len(exon_table)
    scale = height / EXON_H_REF

    for i, (e, w) in enumerate(zip(exon_table, widths)):
        n = e["n"]
        x0, x1 = exon_x[n]
        bump_depth = calculate_bump_depth(w * scale, height) * scale
        draw_x0, draw_x1 = exon_draw_xrange(i, n_exons, x0, x1, bump_depth)

        verts = exon_polygon_verts(
            x0, w, y, height,
            e["five_prime"], e["three_prime"],
            bump_depth=bump_depth,
        )
        codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(verts) - 2) + [MplPath.CLOSEPOLY]
        puzzle_clip = mpatches.PathPatch(
            MplPath(verts, codes), facecolor="none", edgecolor="none",
        )
        ax.add_patch(puzzle_clip)
        puzzle_clips.append(puzzle_clip)

        parts = e["parts"]
        if len(parts) <= 1:
            fill = parts[0][2] if parts else "#cccccc"
            clipped = mpatches.Rectangle(
                (draw_x0, y), draw_x1 - draw_x0, height,
                facecolor=fill, edgecolor="none", zorder=zorder + 1,
            )
            clipped.set_clip_path(puzzle_clip)
            ax.add_patch(clipped)
        else:
            rgba = build_horizontal_gradient(
                parts,
                n_samples=256,
                fade_frac=0.12,
                draw_x0=draw_x0,
                draw_x1=draw_x1,
                exon_x0=x0,
                exon_x1=x1,
            )
            extent = [draw_x0, draw_x1, y, y + height]
            im = ax.imshow(
                rgba.reshape(1, -1, 4),
                extent=extent,
                origin="lower",
                aspect="auto",
                zorder=zorder + 1,
                interpolation="bilinear",
            )
            im.set_clip_path(puzzle_clip)

    return segments, puzzle_clips


def add_exon_outline(
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
    extend_right_bump: bool = True,
    segments: list[PartSegment] | None = None,
    exon_index: int = 0,
    n_exons: int = 1,
    exon_n: int | None = None,
) -> None:
    """Draw one transcript exon outline over the colour band."""
    x0, x1 = x, x + width
    draw_x0, draw_x1 = exon_draw_xrange(
        exon_index, n_exons, x0, x1, bump_depth,
    )

    verts = exon_polygon_verts(
        x, width, y, height, left_shape, right_shape, bump_depth=bump_depth,
    )
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(verts) - 2) + [MplPath.CLOSEPOLY]
    outline_path = MplPath(verts, codes)
    puzzle_clip = mpatches.PathPatch(outline_path, facecolor="none", edgecolor="none")
    ax.add_patch(puzzle_clip)

    if segments is None or exon_n is None:
        return

    border = part_border_color(parts)
    ax.add_patch(mpatches.PathPatch(
        outline_path,
        facecolor="none",
        edgecolor=border,
        lw=lw,
        zorder=5,
    ))


# Backward-compatible alias used by older call sites.
add_exon_patch = add_exon_outline


def draw_domain_row(
    ax,
    resolved_domains: list,
    segments: list[PartSegment],
    y: float,
    height: float,
    *,
    lw: float = 0.75,
    zorder: int = 3,
    font_scale: float = 1.0,
) -> None:
    """Segment-aligned fills and borders inside rounded domain bubbles."""
    from src.domain_alignment import domain_corner_radius, domain_font_size

    for domain in resolved_domains:
        w = domain.draw_width
        rx = domain_corner_radius(w, height)
        clip_patch = mpatches.FancyBboxPatch(
            (domain.draw_x0, y),
            w,
            height,
            boxstyle=f"round,pad=0,rounding_size={rx}",
            facecolor="none",
            edgecolor="none",
            zorder=zorder,
        )
        ax.add_patch(clip_patch)

        local = segments_in_xrange(segments, domain.draw_x0, domain.draw_x1)
        for seg in local:
            ix0 = max(seg.x0, domain.draw_x0)
            ix1 = min(seg.x1, domain.draw_x1)
            if ix1 <= ix0:
                continue
            fill = mpatches.Rectangle(
                (ix0, y), ix1 - ix0, height,
                facecolor=seg.fill, edgecolor="none", zorder=zorder + 1,
            )
            fill.set_clip_path(clip_patch)
            ax.add_patch(fill)

        for seg in local:
            px0 = max(seg.x0, domain.draw_x0)
            px1 = min(seg.x1, domain.draw_x1)
            if px1 <= px0:
                continue
            outline = mpatches.FancyBboxPatch(
                (domain.draw_x0, y),
                w,
                height,
                boxstyle=f"round,pad=0,rounding_size={rx}",
                facecolor="none",
                edgecolor=seg.border,
                lw=lw,
                zorder=zorder + 2,
            )
            outline.set_clip_path(clip_patch)
            outline.set_clip_box(Bbox.from_bounds(px0, y - 0.01, px1 - px0, height + 0.02))
            ax.add_patch(outline)

        fs = domain_font_size(w, height) * font_scale
        ax.text(
            (domain.draw_x0 + domain.draw_x1) / 2,
            y + height / 2,
            domain.label,
            ha="center",
            va="center",
            fontsize=fs,
            fontweight="bold",
            zorder=zorder + 3,
        )
