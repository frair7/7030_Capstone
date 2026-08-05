"""
Interactive SVG cohort mutation map for Streamlit (GeneCards-style hover).

Transcript row: puzzle-clipped fills (gradients for multi-domain exons) with
one outer outline per exon. Domain row: rounded bubbles with segment fills.
"""

from __future__ import annotations

import html
from dataclasses import dataclass, field
from typing import Any, Optional

from src.domain_alignment import (
    DOMAIN_GAP_PX,
    ResolvedDomain,
    SHOW_ALIGNMENT_GUIDES,
    domain_corner_radius,
    domain_font_size,
    hinge_guide_x_positions,
    resolve_domains,
)
from src.dp427m_exon_data import (
    ISOFORM_ANNOTATIONS,
    bar_style_for_phenotype,
    get_domain_map,
    get_exon_table,
)
from src.exon_display_config import pct_dys_color
from src.exon_gradients import svg_linear_gradient_def
from src.mutation_viz import exon_range_from_row, sort_mutations_by_exon_range
from src.part_segments import (
    PartSegment,
    build_part_segments,
    exon_draw_xrange,
    part_border_color,
    segments_in_xrange,
)
from src.puzzle_geometry import (
    build_exon_path,
    build_junction_path,
    calculate_bump_depth,
    domain_bubble_svg_d,
)

# Layout constants (pixels)
# Left metadata panel ends at TABLE_RIGHT; transcript/protein track begins at TRACK_LEFT.
_REF_GROUP, _REF_PART, _REF_MW, _REF_PCT = 0.9, 2.2, 1.4, 1.8
_REF_TABLE_GAP = 0.4
_TABLE_UNIT_PX = 30

TABLE_LEFT = 4
GROUP_W = _REF_GROUP * _TABLE_UNIT_PX
PART_W = _REF_PART * _TABLE_UNIT_PX
MW_W = _REF_MW * _TABLE_UNIT_PX
PCT_W = _REF_PCT * _TABLE_UNIT_PX
TABLE_TRACK_GAP = _REF_TABLE_GAP * _TABLE_UNIT_PX

COL_GROUP = (TABLE_LEFT, TABLE_LEFT + GROUP_W)
COL_PART = (COL_GROUP[1], COL_GROUP[1] + PART_W)
COL_MW = (COL_PART[1], COL_PART[1] + MW_W)
COL_PCT = (COL_MW[1], COL_MW[1] + PCT_W)
TABLE_RIGHT = COL_PCT[1]
TRACK_LEFT = TABLE_RIGHT + TABLE_TRACK_GAP
LEFT = TRACK_LEFT  # exon / mutation track origin

TITLE_Y = 15
TRANSCRIPT_Y = 52
EXON_H = 26
TRANSCRIPT_DOMAIN_GAP = 12  # white break between transcript and domain rows
DOMAIN_Y = TRANSCRIPT_Y + EXON_H + TRANSCRIPT_DOMAIN_GAP
DOMAIN_H = 26
ROW_GAP = 18
TABLE_TOP = DOMAIN_Y + DOMAIN_H + ROW_GAP + 28
PATIENT_ROW_H = 28
MIN_EXON_PX = 14
N_EXONS = 79

# Mutation bar insets — white margin around bars (not gray flanking boxes)
MUT_BAR_V_PAD = 6
MUT_BAR_H_PAD = 4

# Skip-candidate outline colors (combo 1 red, 2 green, 3 blue)
SKIP_COMBO_COLORS = ("#ff1a1a", "#00c853", "#2196f3")
SKIP_DELETED_EXON_FILL = "#4a4a4a"
SKIP_DELETED_EXON_BORDER = "#2b2b2b"


@dataclass
class SkipSchematicOverlay:
    """Exon-skipping overlay layers drawn on the shared transcript schematic."""

    mutation_deleted_exons: set[int] = field(default_factory=set)
    skip_candidate_exons: list[set[int]] = field(default_factory=list)
    skip_target_exons: set[int] = field(default_factory=set)
    junction_upstream: Optional[int] = None
    junction_downstream: Optional[int] = None


def _exon_puzzle_path(
    exon_table: list[dict[str, Any]],
    widths: list[float],
    exon_x: dict[int, tuple[float, float]],
    exon_num: int,
) -> str:
    """Return SVG path for one transcript exon's puzzle shape."""
    n_exons = len(exon_table)
    for i, (e, w) in enumerate(zip(exon_table, widths)):
        if e["n"] != exon_num:
            continue
        x0, x1 = exon_x[exon_num]
        bump_depth = calculate_bump_depth(w, EXON_H)
        return build_exon_path(
            x0, x1, TRANSCRIPT_Y, EXON_H,
            e["five_prime"], e["three_prime"], bump_depth,
        )
    return ""


def _skip_overlay_svgs(
    overlay: SkipSchematicOverlay,
    exon_x: dict[int, tuple[float, float]],
    exon_table: list[dict[str, Any]],
    widths: list[float],
) -> list[str]:
    """Render skip-analysis overlays: colored outlines for proposed skip exons."""
    svgs: list[str] = []
    y = TRANSCRIPT_Y
    h = EXON_H
    bar_y = y + h + 6

    skip_sets: list[tuple[set[int], str]] = []
    if overlay.skip_candidate_exons:
        for i, exon_set in enumerate(overlay.skip_candidate_exons[:3]):
            if exon_set:
                skip_sets.append((exon_set, SKIP_COMBO_COLORS[i]))
    elif overlay.skip_target_exons:
        skip_sets.append((overlay.skip_target_exons, SKIP_COMBO_COLORS[0]))

    for exon_set, color in skip_sets:
        for ex in sorted(exon_set):
            puzzle_d = _exon_puzzle_path(exon_table, widths, exon_x, ex)
            if not puzzle_d:
                continue
            svgs.append(
                f'<path d="{puzzle_d}" fill="rgba(255,255,255,0)" stroke="{color}" '
                f'stroke-width="3.2" pointer-events="none"/>'
            )

    removed = sorted(overlay.mutation_deleted_exons)
    for skip_set, _ in skip_sets:
        removed = sorted(set(removed) | skip_set)

    def _contiguous_segments(exons: list[int]) -> list[list[int]]:
        if not exons:
            return []
        segments: list[list[int]] = [[exons[0]]]
        for ex in exons[1:]:
            if ex == segments[-1][-1] + 1:
                segments[-1].append(ex)
            else:
                segments.append([ex])
        return segments

    mutation_segments = _contiguous_segments(sorted(overlay.mutation_deleted_exons))
    skip_only = sorted(
        set().union(*(s for s, _ in skip_sets)) if skip_sets else overlay.skip_target_exons
    )
    skip_segments = _contiguous_segments(skip_only)

    bar_y = y + h + 6
    for segment in mutation_segments + skip_segments:
        if not segment:
            continue
        x0 = exon_x[segment[0]][0]
        x1 = exon_x[segment[-1]][1]
        svgs.append(
            f'<line x1="{x0:.1f}" y1="{bar_y:.1f}" x2="{x1:.1f}" y2="{bar_y:.1f}" '
            f'stroke="#111" stroke-width="1.4" pointer-events="none"/>'
        )
        svgs.append(
            f'<line x1="{x0:.1f}" y1="{bar_y - 4:.1f}" x2="{x0:.1f}" y2="{bar_y + 4:.1f}" '
            f'stroke="#111" stroke-width="1.4" pointer-events="none"/>'
        )
        svgs.append(
            f'<line x1="{x1:.1f}" y1="{bar_y - 4:.1f}" x2="{x1:.1f}" y2="{bar_y + 4:.1f}" '
            f'stroke="#111" stroke-width="1.4" pointer-events="none"/>'
        )

    if overlay.junction_upstream and overlay.junction_downstream:
        if overlay.junction_upstream in exon_x:
            xj = exon_x[overlay.junction_upstream][1]
            svgs.append(
                f'<line x1="{xj:.1f}" y1="{y - 6:.1f}" x2="{xj:.1f}" y2="{bar_y + 10:.1f}" '
                f'stroke="#0f766e" stroke-width="2.5" pointer-events="none"/>'
            )
            svgs.append(
                f'<text x="{xj:.1f}" y="{y - 10:.1f}" text-anchor="middle" font-size="10" '
                f'font-weight="700" fill="#0f766e">upstream|downstream</text>'
            )

    return svgs


def _row_label_x() -> float:
    """Right-aligned labels in the gap between table and track."""
    return TRACK_LEFT - 6


def _domain_row_svgs(
    resolved: list[ResolvedDomain],
    segments: list[PartSegment],
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Segment-aligned fills and borders inside rounded domain bubbles."""
    defs: list[str] = []
    fills: list[str] = []
    outlines: list[str] = []
    labels: list[str] = []

    for domain in resolved:
        x, y = domain.draw_x0, DOMAIN_Y
        w, h = domain.draw_width, DOMAIN_H
        rx = domain_corner_radius(w, h)
        path_d = domain_bubble_svg_d(x, y, w, h, radius=rx)
        clip_id = f"clip-dom-{domain.label}"
        defs.append(f'<clipPath id="{clip_id}"><path d="{path_d}"/></clipPath>')

        local = segments_in_xrange(segments, domain.draw_x0, domain.draw_x1)
        fill_parts: list[str] = []
        for seg in local:
            ix0 = max(seg.x0, domain.draw_x0)
            ix1 = min(seg.x1, domain.draw_x1)
            if ix1 <= ix0:
                continue
            fill_parts.append(
                f'<rect x="{ix0:.2f}" y="{y:.2f}" width="{(ix1 - ix0):.2f}" '
                f'height="{h:.2f}" fill="{seg.fill}" stroke="none"/>'
            )
        if fill_parts:
            fills.append(f'<g clip-path="url(#{clip_id})">{"".join(fill_parts)}</g>')

        for k, seg in enumerate(local):
            px0 = max(seg.x0, domain.draw_x0)
            px1 = min(seg.x1, domain.draw_x1)
            if px1 <= px0:
                continue
            part_clip_id = f"clip-dom-{domain.label}-s{k}"
            defs.append(
                f'<clipPath id="{part_clip_id}">'
                f'<rect x="{px0:.2f}" y="{y - 1:.2f}" '
                f'width="{(px1 - px0):.2f}" height="{h + 2:.2f}"/>'
                f'</clipPath>'
            )
            outlines.append(
                f'<g clip-path="url(#{clip_id})">'
                f'<g clip-path="url(#{part_clip_id})">'
                f'<path d="{path_d}" fill="none" stroke="{seg.border}" stroke-width="0.75"/>'
                f'</g></g>'
            )

        labels.append(
            f'<text x="{(domain.draw_x0 + domain.draw_x1) / 2:.1f}" '
            f'y="{DOMAIN_Y + DOMAIN_H / 2 + 4:.1f}" text-anchor="middle" '
            f'font-size="{domain_font_size(w, h):.1f}" font-weight="700" pointer-events="none">'
            f'{html.escape(domain.label)}</text>'
        )

    return defs, fills, outlines, labels


def _alignment_guide_svgs(
    domains: list[ResolvedDomain],
    y_top: float,
    y_bot: float,
) -> list[str]:
    guides: list[str] = []
    for domain in domains:
        for xg in (domain.biological_x0, domain.biological_x1):
            guides.append(
                f'<line x1="{xg:.2f}" y1="{y_top:.2f}" x2="{xg:.2f}" y2="{y_bot:.2f}" '
                f'stroke="#c44" stroke-width="0.5" opacity="0.65"/>'
            )
    return guides


def _compute_layout(chart_width: float) -> tuple[list[float], dict[int, tuple[float, float]]]:
    exon_table = get_exon_table()
    track_width = max(400.0, chart_width - LEFT)
    total_bp = sum(e["bp"] for e in exon_table)
    raw = [e["bp"] / total_bp * track_width for e in exon_table]
    fixed = [w < MIN_EXON_PX for w in raw]
    fixed_total = sum(MIN_EXON_PX for f in fixed if f)
    raw_large = sum(w for w, f in zip(raw, fixed) if not f)
    remaining = max(120.0, track_width - fixed_total)
    widths = [
        MIN_EXON_PX if f else (w / raw_large * remaining)
        for w, f in zip(raw, fixed)
    ]
    exon_x: dict[int, tuple[float, float]] = {}
    x = LEFT
    for e, w in zip(exon_table, widths):
        exon_x[e["n"]] = (x, x + w)
        x += w
    return widths, exon_x


def _bar_style(row: dict[str, Any]) -> dict[str, str]:
    style = bar_style_for_phenotype(str(row.get("phenotype", "")))
    if style == "black":
        return {"fill": "#111", "stroke": "#111", "text": "#fff"}
    if style == "outline":
        return {"fill": "#fff", "stroke": "#111", "text": "#111"}
    return {"fill": "#9a9a9a", "stroke": "#111", "text": "#111"}


def _mutation_bar_box(
    x0b: float,
    x1b: float,
    y: float,
) -> tuple[float, float, float, float]:
    """Bar rect with white margin on all sides (no gray flanking boxes)."""
    span = max(x1b - x0b, 1.0)
    h_pad = min(MUT_BAR_H_PAD, span * 0.12)
    bx0 = x0b + h_pad
    bx1 = x1b - h_pad
    by = y + MUT_BAR_V_PAD
    bh = PATIENT_ROW_H - 2 * MUT_BAR_V_PAD
    return bx0, by, bx1 - bx0, bh


def _color_band_svgs(
    segments: list[PartSegment],
    *,
    y: float,
    height: float,
) -> list[str]:
    """Domain-aligned vertical colour strips at a given row y."""
    svgs: list[str] = []
    for seg in segments:
        w = seg.x1 - seg.x0
        if w <= 0:
            continue
        svgs.append(
            f'<rect x="{seg.x0:.2f}" y="{y:.2f}" width="{w:.2f}" '
            f'height="{height:.2f}" fill="{seg.fill}" stroke="none"/>'
        )
    return svgs


def _transcript_color_band_svgs(segments: list[PartSegment]) -> list[str]:
    return _color_band_svgs(segments, y=TRANSCRIPT_Y, height=EXON_H)


def _transcript_exon_svgs(
    exon_table: list[dict[str, Any]],
    widths: list[float],
    exon_x: dict[int, tuple[float, float]],
    *,
    visible_exon_range: tuple[int, int] | None = None,
    gray_exons: Optional[set[int]] = None,
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Clipped exon fills (with gradients for multi-part) and one outline per exon."""
    defs: list[str] = []
    fills: list[str] = []
    outlines: list[str] = []
    junction_overlays: list[str] = []
    n_exons = len(exon_table)
    lo, hi = visible_exon_range if visible_exon_range else (1, N_EXONS)

    for i, (e, w) in enumerate(zip(exon_table, widths)):
        n = e["n"]
        if n < lo or n > hi:
            continue
        x0, x1 = exon_x[n]
        bump_depth = calculate_bump_depth(w, EXON_H)
        draw_x0, draw_x1 = exon_draw_xrange(i, n_exons, x0, x1, bump_depth)
        parts = e["parts"]

        puzzle_d = build_exon_path(
            x0, x1, TRANSCRIPT_Y, EXON_H,
            e["five_prime"], e["three_prime"], bump_depth,
        )
        clip_id = f"clip-e{n}"
        defs.append(f'<clipPath id="{clip_id}"><path d="{puzzle_d}"/></clipPath>')

        is_deleted = gray_exons is not None and n in gray_exons
        if is_deleted:
            fill_color = SKIP_DELETED_EXON_FILL
            border = SKIP_DELETED_EXON_BORDER
        elif len(parts) <= 1:
            fill_color = parts[0][2] if parts else "#cccccc"
            border = part_border_color(parts)
        else:
            fill_color = None
            border = part_border_color(parts)

        if is_deleted or len(parts) <= 1:
            fills.append(
                f'<g clip-path="url(#{clip_id})">'
                f'<rect x="{draw_x0:.2f}" y="{TRANSCRIPT_Y:.2f}" '
                f'width="{(draw_x1 - draw_x0):.2f}" height="{EXON_H:.2f}" '
                f'fill="{fill_color}" stroke="none"/></g>'
            )
        else:
            grad_id = f"grad-e{n}"
            defs.append(svg_linear_gradient_def(
                grad_id, draw_x0, draw_x1, TRANSCRIPT_Y, parts,
                exon_x0=x0, exon_x1=x1,
            ))
            fills.append(
                f'<g clip-path="url(#{clip_id})">'
                f'<rect x="{draw_x0:.2f}" y="{TRANSCRIPT_Y:.2f}" '
                f'width="{(draw_x1 - draw_x0):.2f}" height="{EXON_H:.2f}" '
                f'fill="url(#{grad_id})" stroke="none"/></g>'
            )

        outlines.append(
            f'<path d="{puzzle_d}" fill="none" stroke="{border}" '
            f'stroke-width="0.75" pointer-events="none"/>'
        )

        tip = (
            f"Exon {n} · {e['bp']} bp CDS\\n"
            f"5′ {e['five_prime']} · 3′ {e['three_prime']}"
            + (f"\\n{len(parts)} domain segments" if len(parts) > 1 else "")
        )
        outlines.append(
            f'<path class="exon" data-tip="{html.escape(tip)}" d="{puzzle_d}" '
            f'fill="none" stroke="transparent" stroke-width="10"/>'
        )

        if i < n_exons - 1:
            next_e = exon_table[i + 1]
            if lo <= next_e["n"] <= hi:
                next_w = widths[i + 1]
                j_depth = min(bump_depth, calculate_bump_depth(next_w, EXON_H))
                junction_d = build_junction_path(
                    x1, TRANSCRIPT_Y, EXON_H,
                    e["three_prime"], next_e["five_prime"],
                    j_depth,
                )
                junction_overlays.append(
                    f'<path d="{junction_d}" fill="none" stroke="{border}" '
                    f'stroke-width="0.75"/>'
                )

        fs = 11 if w >= 16 else (9 if w >= 10 else 7)
        if w >= 6:
            label_fill = "#ffffff" if is_deleted else "#111111"
            outlines.append(
                f'<text x="{x0 + w/2:.1f}" y="{TRANSCRIPT_Y + EXON_H/2 + 4:.1f}" '
                f'text-anchor="middle" font-size="{fs}" font-weight="600" '
                f'fill="{label_fill}" pointer-events="none">{n}</text>'
            )

    return defs, fills, outlines, junction_overlays


def _vertical_table_header_svg(
    x0: float,
    x1: float,
    header_y: float,
    header_h: float,
    label: str,
) -> list[str]:
    """Rotated table header label (reads bottom-to-top)."""
    cx = (x0 + x1) / 2
    cy = header_y + header_h / 2
    return [
        f'<rect x="{x0}" y="{header_y}" width="{x1 - x0}" height="{header_h}" '
        f'fill="#ffffff" stroke="#000" stroke-width="1"/>',
        f'<text class="table-header-vertical" x="{cx:.1f}" y="{cy:.1f}" '
        f'transform="rotate(-90 {cx:.1f} {cy:.1f})" '
        f'text-anchor="middle" dominant-baseline="middle" '
        f'font-size="11" font-weight="700">{html.escape(label)}</text>',
    ]


def _segment_boundary_lines_svgs(
    segments: list[PartSegment],
    exon_table: list[dict[str, Any]],
    widths: list[float],
    exon_x: dict[int, tuple[float, float]],
    *,
    y: float,
    height: float,
    row_height: float = EXON_H,
) -> list[str]:
    """Vertical dividers at internal part boundaries within multi-part exons."""
    svgs: list[str] = []
    boundaries: set[float] = set()
    for i in range(len(segments) - 1):
        if abs(segments[i].x1 - segments[i + 1].x0) < 1e-6:
            boundaries.add(segments[i].x1)
    n_exons = len(exon_table)
    for i, (e, w) in enumerate(zip(exon_table, widths)):
        if len(e["parts"]) < 2:
            continue
        x0, x1 = exon_x[e["n"]]
        bump_depth = calculate_bump_depth(w, row_height)
        _, draw_x1 = exon_draw_xrange(i, n_exons, x0, x1, bump_depth)
        for bx in boundaries:
            if x0 < bx < draw_x1:
                border = next(
                    (s.border for s in segments if abs(s.x1 - bx) < 1e-6),
                    "#333",
                )
                svgs.append(
                    f'<line x1="{bx:.2f}" y1="{y:.2f}" x2="{bx:.2f}" '
                    f'y2="{y + height:.2f}" stroke="{border}" '
                    f'stroke-width="0.9" opacity="0.8"/>'
                )
    return svgs


def map_iframe_height(n_patient_rows: int, *, chart_width: float = 2000) -> int:
    """Streamlit iframe height — matches intrinsic SVG pixel height."""
    n_rows = max(n_patient_rows, 1)
    return int(TABLE_TOP + n_rows * PATIENT_ROW_H + 64)


def _wrap_map_html(svg_body: str, *, end_x: float, height: float, title: str) -> str:
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"/>
<style>
  body {{ margin:0; background:#fff; }}
  .map-scroll {{ overflow-x:auto; overflow-y:hidden; width:100%; }}
  svg {{ display:block; font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif; }}
  .exon, .mut-bar {{ cursor: pointer; }}
  .exon:hover, .mut-bar:hover {{ filter: brightness(0.94); }}
  #tip {{
    position: fixed; display:none; background:#1a2332; color:#fff;
    padding:8px 12px; border-radius:6px; font-size:13px; line-height:1.45;
    max-width:280px; pointer-events:none; z-index:9999;
    box-shadow:0 4px 14px rgba(0,0,0,.22); white-space:pre-line;
  }}
</style></head><body>
<div id="tip"></div>
<div class="map-scroll">
<svg viewBox="0 0 {end_x:.0f} {height:.0f}" width="{end_x:.0f}" height="{height:.0f}"
     xmlns="http://www.w3.org/2000/svg">
  <text x="{end_x/2:.0f}" y="{TITLE_Y}" text-anchor="middle" font-size="16" font-weight="700">{html.escape(title)}</text>
  {svg_body}
</svg>
</div>
<script>
const tip = document.getElementById('tip');
document.querySelectorAll('[data-tip]').forEach(el => {{
  el.addEventListener('mousemove', e => {{
    tip.style.display = 'block';
    tip.textContent = el.getAttribute('data-tip').replace(/\\\\n/g, '\\n');
    tip.style.left = (e.clientX + 14) + 'px';
    tip.style.top = (e.clientY + 14) + 'px';
  }});
  el.addEventListener('mouseleave', () => {{ tip.style.display = 'none'; }});
}});
</script></body></html>"""


def _build_map_layers(
    rows: list[dict[str, Any]],
    *,
    chart_width: float,
    show_patient_table: bool = True,
    skip_overlay: Optional[SkipSchematicOverlay] = None,
) -> tuple[str, float, float]:
    """Assemble SVG body; returns (body, end_x, height)."""
    rows = sort_mutations_by_exon_range(rows)
    widths, exon_x = _compute_layout(chart_width)
    exon_table = get_exon_table()
    end_x = exon_x[N_EXONS][1] + 24

    if show_patient_table:
        n_rows = max(len(rows), 1)
        height = TABLE_TOP + n_rows * PATIENT_ROW_H + 40
        y_bot_guides = TABLE_TOP + len(rows) * PATIENT_ROW_H
    else:
        height = DOMAIN_Y + DOMAIN_H + 48
        y_bot_guides = DOMAIN_Y + DOMAIN_H

    domain_map = get_domain_map()
    resolved = resolve_domains(domain_map, exon_table, exon_x, gap=DOMAIN_GAP_PX)
    part_segments = build_part_segments(exon_table, exon_x)

    domain_defs, domain_fill_svgs, domain_outline_svgs, domain_label_svgs = _domain_row_svgs(
        resolved, part_segments,
    )
    exon_defs, exon_fill_svgs, exon_outline_svgs, junction_overlays = _transcript_exon_svgs(
        exon_table, widths, exon_x,
        gray_exons=skip_overlay.mutation_deleted_exons if skip_overlay else None,
    )

    iso_svgs: list[str] = []
    for exon, label in ISOFORM_ANNOTATIONS:
        xa = exon_x[exon][0] + 2
        iso_svgs.append(
            f'<path d="M {xa:.1f} {TRANSCRIPT_Y:.1f} L {xa:.1f} {TRANSCRIPT_Y - 14:.1f} '
            f'L {xa + 30:.1f} {TRANSCRIPT_Y - 14:.1f}" fill="none" stroke="#222" stroke-width="1.2"/>'
        )
        iso_svgs.append(
            f'<polygon points="{xa+30:.1f},{TRANSCRIPT_Y-17:.1f} {xa+37:.1f},{TRANSCRIPT_Y-14:.1f} '
            f'{xa+30:.1f},{TRANSCRIPT_Y-11:.1f}" fill="#222"/>'
        )
        iso_svgs.append(
            f'<text x="{xa+3:.1f}" y="{TRANSCRIPT_Y - 18:.1f}" font-size="11" font-weight="700">'
            f'{html.escape(label)}</text>'
        )

    guide_svgs: list[str] = []
    y_top = DOMAIN_Y
    for xg in hinge_guide_x_positions(domain_map, exon_table, exon_x):
        guide_svgs.append(
            f'<line class="hinge-guide" x1="{xg:.2f}" y1="{y_top:.2f}" x2="{xg:.2f}" y2="{y_bot_guides:.2f}" '
            f'stroke="#666" stroke-width="1.2" stroke-dasharray="5,4"/>'
        )
    if SHOW_ALIGNMENT_GUIDES:
        guide_svgs.extend(_alignment_guide_svgs(resolved, y_top, y_bot_guides))

    table_svgs: list[str] = []
    if show_patient_table:
        header_h = 72
        header_y = TABLE_TOP - header_h - 4
        for x0, x1, label in [
            (COL_GROUP[0], COL_GROUP[1], "Group"),
            (COL_PART[0], COL_PART[1], "Participant"),
            (COL_MW[0], COL_MW[1], "MW"),
            (COL_PCT[0], COL_PCT[1], "%Dys(WB)"),
        ]:
            table_svgs.extend(_vertical_table_header_svg(x0, x1, header_y, header_h, label))

        if not rows:
            table_svgs.append(
                f'<text x="{LEFT + 200:.0f}" y="{TABLE_TOP + 20:.0f}" font-size="14" fill="#666">'
                f'Select mutations below and click Plot selected on map</text>'
            )
        else:
            group_bounds: dict[str, list[float]] = {}
            y = TABLE_TOP
            for row in rows:
                pid = html.escape(str(row.get("id", "")))
                grp = str(row.get("group", "") or "—")
                kda = row.get("expected_protein_size_kda", "")
                try:
                    mw = f"{float(kda):.1f}"
                except (TypeError, ValueError):
                    mw = str(kda) if kda else "—"
                pct_raw = row.get("pct_dys_wb", "")
                try:
                    pct = float(pct_raw)
                    pct_text = f"{pct:.2f}"
                    pct_fill = pct_dys_color(pct)
                except (TypeError, ValueError):
                    pct_text = str(pct_raw) if pct_raw else ""
                    pct_fill = "#ffffff"

                group_bounds.setdefault(grp, [y, y + PATIENT_ROW_H])
                group_bounds[grp][1] = y + PATIENT_ROW_H

                for x0, x1, text in [
                    (COL_PART[0], COL_PART[1], pid),
                    (COL_MW[0], COL_MW[1], mw),
                ]:
                    table_svgs.append(
                        f'<rect x="{x0}" y="{y}" width="{x1 - x0}" height="{PATIENT_ROW_H}" '
                        f'fill="#ffffff" stroke="#000" stroke-width="0.6"/>'
                    )
                    table_svgs.append(
                        f'<text x="{(x0 + x1) / 2:.1f}" y="{y + PATIENT_ROW_H / 2 + 4:.1f}" '
                        f'text-anchor="middle" font-size="11">{text}</text>'
                    )

                table_svgs.append(
                    f'<rect x="{COL_PCT[0]}" y="{y}" width="{COL_PCT[1] - COL_PCT[0]}" '
                    f'height="{PATIENT_ROW_H}" fill="{pct_fill}" stroke="#000" stroke-width="0.6"/>'
                )
                if pct_text:
                    table_svgs.append(
                        f'<text x="{(COL_PCT[0] + COL_PCT[1]) / 2:.1f}" y="{y + PATIENT_ROW_H / 2 + 4:.1f}" '
                        f'text-anchor="middle" font-size="10" font-weight="700">{pct_text}</text>'
                    )

                table_svgs.append(
                    f'<rect x="{TABLE_RIGHT:.1f}" y="{y}" width="{TRACK_LEFT - TABLE_RIGHT:.1f}" '
                    f'height="{PATIENT_ROW_H}" fill="#ffffff" stroke="none"/>'
                )
                table_svgs.append(
                    f'<rect x="{LEFT}" y="{y}" width="{end_x - LEFT:.1f}" height="{PATIENT_ROW_H}" '
                    f'fill="#ffffff" stroke="none"/>'
                )

                rng = exon_range_from_row(row)
                if rng:
                    first, last = rng
                    x0b, x1b = exon_x[first][0], exon_x[last][1]
                    sty = _bar_style(row)
                    lbl = str(first) if first == last else f"{first}-{last}"
                    bx, by, bw, bh = _mutation_bar_box(x0b, x1b, y)
                    table_svgs.append(
                        f'<rect class="mut-bar" data-tip="{html.escape(pid)}: exons {lbl}" '
                        f'x="{bx:.1f}" y="{by:.1f}" width="{bw:.1f}" height="{bh:.1f}" '
                        f'fill="{sty["fill"]}" stroke="{sty["stroke"]}" stroke-width="1.2" rx="1"/>'
                    )
                    table_svgs.append(
                        f'<text x="{(x0b + x1b) / 2:.1f}" y="{y + PATIENT_ROW_H / 2 + 4:.1f}" '
                        f'text-anchor="middle" font-size="11" font-weight="700" fill="{sty["text"]}">'
                        f'{lbl}</text>'
                    )
                y += PATIENT_ROW_H

            for grp, (ytop, ybot) in group_bounds.items():
                table_svgs.append(
                    f'<rect x="{COL_GROUP[0]}" y="{ytop}" width="{COL_GROUP[1] - COL_GROUP[0]}" '
                    f'height="{ybot - ytop}" fill="#ffffff" stroke="#000" stroke-width="1"/>'
                )
                table_svgs.append(
                    f'<text x="{(COL_GROUP[0] + COL_GROUP[1]) / 2:.1f}" '
                    f'y="{(ytop + ybot) / 2 + 5:.1f}" text-anchor="middle" '
                    f'font-size="18" font-weight="700">{html.escape(grp)}</text>'
                )
                table_svgs.append(
                    f'<line x1="{TABLE_LEFT}" y1="{ytop:.1f}" x2="{end_x:.1f}" y2="{ytop:.1f}" '
                    f'stroke="#000" stroke-width="1.4"/>'
                )

            table_bottom = TABLE_TOP + len(rows) * PATIENT_ROW_H
            table_svgs.append(
                f'<rect x="{TABLE_LEFT}" y="{header_y}" width="{TABLE_RIGHT - TABLE_LEFT}" '
                f'height="{table_bottom - header_y}" fill="none" stroke="#000" stroke-width="1.6"/>'
            )

    track_width = end_x - LEFT
    track_layers = [
        f'<defs>{"".join(exon_defs)}{"".join(domain_defs)}</defs>',
        *exon_fill_svgs,
        *domain_fill_svgs,
        *domain_outline_svgs,
        *domain_label_svgs,
        f'<rect x="{LEFT:.0f}" y="{TRANSCRIPT_Y + EXON_H:.2f}" '
        f'width="{track_width:.0f}" height="{TRANSCRIPT_DOMAIN_GAP:.2f}" '
        f'fill="#ffffff" stroke="none"/>',
        *exon_outline_svgs,
        *junction_overlays,
        *iso_svgs,
    ]
    label_layers = [
        f'<text x="{_row_label_x():.1f}" y="{TRANSCRIPT_Y + EXON_H/2 + 4}" '
        f'text-anchor="end" font-size="13" font-weight="700">Transcript →</text>',
        f'<text x="{_row_label_x():.1f}" y="{DOMAIN_Y + DOMAIN_H/2 + 4}" '
        f'text-anchor="end" font-size="13" font-weight="700">Domain →</text>',
    ]
    overlay_layers = (
        _skip_overlay_svgs(skip_overlay, exon_x, exon_table, widths) if skip_overlay else []
    )
    svg_body = "\n".join(track_layers + overlay_layers + label_layers + table_svgs + guide_svgs)
    return svg_body, end_x, height


def build_interactive_map_html(
    rows: list[dict[str, Any]],
    *,
    chart_width: float = 2000,
    title: str = "DMD mutation map",
    skip_overlay: Optional[SkipSchematicOverlay] = None,
) -> str:
    svg_body, end_x, height = _build_map_layers(
        rows, chart_width=chart_width, skip_overlay=skip_overlay,
    )
    return _wrap_map_html(svg_body, end_x=end_x, height=height, title=title)


def build_skip_target_map_html(
    target_exons: list[int],
    *,
    chart_width: float = 2000,
    title: str = "Proposed skip target",
) -> str:
    """Schematic of a skip target without a catalog mutation row."""
    overlay = SkipSchematicOverlay(skip_target_exons=set(target_exons))
    svg_body, end_x, height = _build_map_layers(
        [], chart_width=chart_width, show_patient_table=False, skip_overlay=overlay,
    )
    return _wrap_map_html(svg_body, end_x=end_x, height=height, title=title)
