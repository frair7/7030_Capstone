"""
Interactive SVG cohort mutation map for Streamlit (GeneCards-style hover).

Transcript row: interlocking puzzle-piece exons with domain-aligned gradients.
Domain row: independent rounded bubbles (one color per domain segment).
"""

from __future__ import annotations

import html
from typing import Any

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
from src.mutation_viz import exon_range_from_row
from src.part_segments import part_rects_for_exon
from src.puzzle_geometry import (
    build_exon_path,
    build_exon_stroke_path,
    build_junction_path,
    calculate_bump_depth,
)

# Layout constants (pixels)
LEFT = 148
TRANSCRIPT_Y = 42
EXON_H = 26
ROW_GAP_PX = 2
DOMAIN_Y = TRANSCRIPT_Y + EXON_H + ROW_GAP_PX
DOMAIN_H = 26
ROW_GAP = 18
TABLE_TOP = DOMAIN_Y + DOMAIN_H + ROW_GAP + 28
PATIENT_ROW_H = 28
MIN_EXON_PX = 14
N_EXONS = 79


def _draw_domain_box_svg(domain: ResolvedDomain) -> list[str]:
    """Rounded rect at draw coords; biological coords used for alignment only."""
    w = domain.draw_width
    h = DOMAIN_H
    rx = domain_corner_radius(w, h)
    svgs = [
        f'<rect x="{domain.draw_x0:.2f}" y="{DOMAIN_Y:.2f}" width="{w:.2f}" '
        f'height="{h:.2f}" rx="{rx:.2f}" ry="{rx:.2f}" fill="{domain.fill}" '
        f'stroke="#222" stroke-width="0.9"/>',
        f'<text x="{(domain.draw_x0 + domain.draw_x1) / 2:.1f}" '
        f'y="{DOMAIN_Y + DOMAIN_H / 2 + 4:.1f}" text-anchor="middle" '
        f'font-size="{domain_font_size(w, h):.1f}" font-weight="700" pointer-events="none">'
        f'{html.escape(domain.label)}</text>',
    ]
    return svgs


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
    total_bp = sum(e["bp"] for e in exon_table)
    raw = [e["bp"] / total_bp * chart_width for e in exon_table]
    fixed = [w < MIN_EXON_PX for w in raw]
    fixed_total = sum(MIN_EXON_PX for f in fixed if f)
    raw_large = sum(w for w, f in zip(raw, fixed) if not f)
    remaining = max(120.0, chart_width - fixed_total)
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
    return {"fill": "#9a9a9a", "stroke": "#333", "text": "#111"}


def build_interactive_map_html(
    rows: list[dict[str, Any]],
    *,
    chart_width: float = 2000,
    title: str = "DMD mutation map",
) -> str:
    widths, exon_x = _compute_layout(chart_width)
    exon_table = get_exon_table()
    end_x = exon_x[N_EXONS][1] + 24
    n_rows = max(len(rows), 1)
    height = TABLE_TOP + n_rows * PATIENT_ROW_H + 40

    defs: list[str] = []
    exon_svgs: list[str] = []
    domain_svgs: list[str] = []

    domain_map = get_domain_map()
    resolved = resolve_domains(domain_map, exon_table, exon_x, gap=DOMAIN_GAP_PX)
    for domain in resolved:
        domain_svgs.extend(_draw_domain_box_svg(domain))

    # Transcript: clip puzzle outline + solid per-part rects at absolute x
    junction_overlays: list[str] = []
    n_exons = len(exon_table)
    for i, (e, w) in enumerate(zip(exon_table, widths)):
        n = e["n"]
        x0 = exon_x[n][0]
        x1 = exon_x[n][1]
        bump_depth = calculate_bump_depth(w, EXON_H)
        clip_id = f"clip-e{n}"
        path_d = build_exon_path(
            x0, x1, TRANSCRIPT_Y, EXON_H,
            e["five_prime"], e["three_prime"],
            bump_depth,
        )
        defs.append(f'<clipPath id="{clip_id}"><path d="{path_d}"/></clipPath>')

        part_rects = part_rects_for_exon(n, exon_x, e["parts"])
        fill_svgs: list[str] = []
        for j, (rx0, rx1, fill, _border) in enumerate(part_rects):
            rw = rx1 - rx0
            fill_svgs.append(
                f'<rect x="{rx0:.2f}" y="{TRANSCRIPT_Y:.2f}" width="{rw:.2f}" '
                f'height="{EXON_H:.2f}" fill="{fill}" stroke="none"/>'
            )
            if j < len(part_rects) - 1:
                blend_w = min(w * 0.06, 6.0)
                bx = rx1
                fill_svgs.append(
                    f'<rect x="{(bx - blend_w/2):.2f}" y="{TRANSCRIPT_Y:.2f}" '
                    f'width="{blend_w:.2f}" height="{EXON_H:.2f}" fill="{part_rects[j+1][2]}" '
                    f'opacity="0.4" stroke="none"/>'
                )

        tip = (
            f"Exon {n} · {e['bp']} bp CDS\\n"
            f"5′ {e['five_prime']} · 3′ {e['three_prime']}"
            + (f"\\n{len(e['parts'])} domain segments" if len(e["parts"]) > 1 else "")
        )
        exon_svgs.append(f'<g clip-path="url(#{clip_id})">{"".join(fill_svgs)}</g>')

        stroke_d = build_exon_stroke_path(
            x0, x1, TRANSCRIPT_Y, EXON_H,
            e["five_prime"], e["three_prime"], bump_depth,
            stroke_left=(i == 0),
            stroke_right=(i == n_exons - 1),
        )
        exon_svgs.append(
            f'<path class="exon" data-tip="{html.escape(tip)}" d="{stroke_d}" '
            f'fill="none" stroke="#4a4a4a" stroke-width="0.6"/>'
        )

        if i < n_exons - 1:
            next_e = exon_table[i + 1]
            next_w = widths[i + 1]
            j_depth = min(bump_depth, calculate_bump_depth(next_w, EXON_H))
            junction_d = build_junction_path(
                x1, TRANSCRIPT_Y, EXON_H,
                e["three_prime"], next_e["five_prime"],
                j_depth,
            )
            junction_overlays.append(
                f'<path d="{junction_d}" fill="none" stroke="#4a4a4a" stroke-width="0.55"/>'
            )
        fs = 11 if w >= 16 else (9 if w >= 10 else 7)
        exon_svgs.append(
            f'<text x="{x0 + w/2:.1f}" y="{TRANSCRIPT_Y + EXON_H/2 + 4:.1f}" '
            f'text-anchor="middle" font-size="{fs}" font-weight="600" pointer-events="none">{n}</text>'
        )

    # Isoform arrows
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

    # Hinge guides — dashed lines at H1–H4 biological start boundaries
    guide_svgs: list[str] = []
    y_top = TRANSCRIPT_Y
    y_bot = TABLE_TOP + len(rows) * PATIENT_ROW_H
    for xg in hinge_guide_x_positions(domain_map, exon_table, exon_x):
        guide_svgs.append(
            f'<line x1="{xg:.2f}" y1="{y_top:.2f}" x2="{xg:.2f}" y2="{y_bot:.2f}" '
            f'stroke="#888" stroke-width="1" stroke-dasharray="5,4"/>'
        )
    if SHOW_ALIGNMENT_GUIDES:
        guide_svgs.extend(_alignment_guide_svgs(resolved, y_top, y_bot))

    # Patient table + bars
    table_svgs: list[str] = []
    col_group, col_part, col_mw, col_pct = 12, 48, 98, 128
    table_svgs.append(
        f'<text x="{col_group}" y="{TABLE_TOP - 8}" font-size="12" font-weight="700">Group</text>'
        f'<text x="{col_part}" y="{TABLE_TOP - 8}" font-size="12" font-weight="700">Participant</text>'
        f'<text x="{col_mw}" y="{TABLE_TOP - 8}" font-size="12" font-weight="700">MW</text>'
        f'<text x="{col_pct}" y="{TABLE_TOP - 8}" font-size="12" font-weight="700">%Dys(WB)</text>'
    )

    if not rows:
        table_svgs.append(
            f'<text x="{LEFT + 200:.0f}" y="{TABLE_TOP + 20:.0f}" font-size="14" fill="#666">'
            f'Select mutations below and click Plot selected on map</text>'
        )
    else:
        y = TABLE_TOP
        for row in rows:
            pid = html.escape(str(row.get("id", "")))
            grp = html.escape(str(row.get("group", "") or "—"))
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
                pct_fill = "#eee"

            table_svgs.append(f'<text x="{col_group}" y="{y + 18}" font-size="12" font-weight="700">{grp}</text>')
            table_svgs.append(f'<text x="{col_part}" y="{y + 18}" font-size="12">{pid}</text>')
            table_svgs.append(f'<text x="{col_mw}" y="{y + 18}" font-size="12">{mw}</text>')
            table_svgs.append(
                f'<rect x="{col_pct - 4}" y="{y + 2}" width="52" height="22" fill="{pct_fill}" '
                f'stroke="#333" stroke-width="0.8"/>'
            )
            if pct_text:
                table_svgs.append(
                    f'<text x="{col_pct + 22}" y="{y + 18}" text-anchor="middle" '
                    f'font-size="11" font-weight="700">{pct_text}</text>'
                )

            rng = exon_range_from_row(row)
            if rng:
                first, last = rng
                x0b, x1b = exon_x[first][0], exon_x[last][1]
                sty = _bar_style(row)
                lbl = str(first) if first == last else f"{first}-{last}"
                table_svgs.append(
                    f'<rect class="mut-bar" data-tip="{html.escape(pid)}: exons {lbl}" '
                    f'x="{x0b:.1f}" y="{y + 4:.1f}" width="{x1b - x0b:.1f}" height="{PATIENT_ROW_H - 8:.1f}" '
                    f'fill="{sty["fill"]}" stroke="{sty["stroke"]}" stroke-width="1.2" rx="1"/>'
                )
                table_svgs.append(
                    f'<text x="{(x0b+x1b)/2:.1f}" y="{y + 18}" text-anchor="middle" '
                    f'font-size="11" font-weight="700" fill="{sty["text"]}">{lbl}</text>'
                )
            y += PATIENT_ROW_H

    svg_body = "\n".join([
        f'<defs>{"".join(defs)}</defs>',
        *guide_svgs,
        *domain_svgs,
        *exon_svgs,
        *junction_overlays,
        *iso_svgs,
        f'<text x="10" y="{TRANSCRIPT_Y + EXON_H/2 + 4}" font-size="13" font-weight="700">Transcript →</text>',
        f'<text x="10" y="{DOMAIN_Y + DOMAIN_H/2 + 4}" font-size="13" font-weight="700">Domain →</text>',
        *table_svgs,
    ])

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"/>
<style>
  body {{ margin:0; background:#fff; }}
  svg {{ font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif; }}
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
<svg viewBox="0 0 {end_x:.0f} {height:.0f}" width="100%" xmlns="http://www.w3.org/2000/svg">
  <text x="{end_x/2:.0f}" y="22" text-anchor="middle" font-size="16" font-weight="700">{html.escape(title)}</text>
  {svg_body}
</svg>
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
