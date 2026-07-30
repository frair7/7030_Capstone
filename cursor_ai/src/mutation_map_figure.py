"""
Publication-style DMD cohort mutation map (matplotlib).

Ported from the user's authoritative Claude make_figure.py reference:
  - Proportional exon widths from true CDS bp (min-width floor)
  - Puzzle-piece shapes with per-exon domain sub-segments
  - Domain row from nucleotide-coordinate feature table
  - Group | Participant | MW | %Dys(WB) patient table
"""

from __future__ import annotations

from typing import Any, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.figure import Figure

from src.dp427m_exon_data import (
    ISOFORM_ANNOTATIONS,
    bar_style_for_phenotype,
    get_domain_map,
    get_exon_table,
)
from src.domain_alignment import (
    DOMAIN_GAP_DATA,
    ResolvedDomain,
    SHOW_ALIGNMENT_GUIDES,
    hinge_guide_x_positions,
    resolve_domains,
)
from src.exon_display_config import (
    COLOR_MUTATION_BAR,
    COLOR_MUTATION_BAR_BLACK,
    COLOR_MUTATION_BAR_OUTLINE,
    pct_dys_color,
)
from src.exon_shape_renderer import add_exon_outline, draw_domain_row, draw_transcript_color_band
from src.part_segments import build_part_segments
from src.models import ExonRecord
from src.mutation_viz import exon_range_from_row
from src.puzzle_geometry import calculate_bump_depth

# Layout — transcript and domain rows share identical exon x-positions
ROW_H = 1.0
ROW_GAP = 12 / 26          # white break between transcript and domain rows (matches SVG)
EXON_ROW_Y = 0.0
DOMAIN_ROW_Y = EXON_ROW_Y - ROW_H - ROW_GAP
TABLE_TOP_Y = DOMAIN_ROW_Y - ROW_H - 1.65
PATIENT_ROW_H = 0.85
ROW_STEP = PATIENT_ROW_H
LEFT_MARGIN = 6.0
N_EXONS = 79
TOTAL_TRACK_WIDTH = N_EXONS * 1.0
MIN_EXON_W = 0.35
EXON_H = 26.0  # pixel reference for bump-depth scaling
INCHES_PER_DATA_UNIT = 0.38
FONT_SCALE = 1.45

MUT_BAR_V_PAD = 0.10
MUT_BAR_H_PAD = 0.06


def _draw_domains(
    ax: plt.Axes,
    exon_x: dict[int, tuple[float, float]],
) -> list[ResolvedDomain]:
    ax.text(
        LEFT_MARGIN - 0.3, DOMAIN_ROW_Y + ROW_H / 2, "Domain \u2192",
        ha="right", va="center", fontsize=9 * FONT_SCALE / 1.1, fontweight="bold",
    )
    exon_table = get_exon_table()
    segments = build_part_segments(exon_table, exon_x)
    resolved = resolve_domains(
        get_domain_map(), exon_table, exon_x, gap=DOMAIN_GAP_DATA,
    )
    draw_domain_row(
        ax, resolved, segments, DOMAIN_ROW_Y, ROW_H,
        zorder=3, font_scale=FONT_SCALE / 1.1,
    )
    return resolved


def _draw_alignment_guides_mpl(
    ax: plt.Axes,
    domains: list[ResolvedDomain],
    y_top: float,
    y_bottom: float,
) -> None:
    for domain in domains:
        for xg in (domain.biological_x0, domain.biological_x1):
            ax.plot(
                [xg, xg], [y_top, y_bottom],
                linestyle="-", color="#c44", lw=0.5, alpha=0.65, zorder=0,
            )

STYLE_FACE = {
    "grey": COLOR_MUTATION_BAR,
    "black": COLOR_MUTATION_BAR_BLACK,
    "outline": COLOR_MUTATION_BAR_OUTLINE,
}


def _mutation_bar_box(
    x0b: float,
    x1b: float,
    y: float,
) -> tuple[float, float, float, float]:
    """Bar rect with white margin on all sides."""
    span = max(x1b - x0b, 0.01)
    h_pad = min(MUT_BAR_H_PAD, span * 0.12)
    bx0 = x0b + h_pad
    bx1 = x1b - h_pad
    by = y - PATIENT_ROW_H + MUT_BAR_V_PAD
    bh = PATIENT_ROW_H - 2 * MUT_BAR_V_PAD
    return bx0, by, bx1 - bx0, bh


def _compute_exon_layout() -> tuple[list[float], dict[int, tuple[float, float]]]:
    exon_table = get_exon_table()
    total_bp = sum(e["bp"] for e in exon_table)
    raw_w = [e["bp"] / total_bp * TOTAL_TRACK_WIDTH for e in exon_table]
    is_fixed = [w < MIN_EXON_W for w in raw_w]
    fixed_total = sum(MIN_EXON_W for f in is_fixed if f)
    raw_large_total = sum(w for w, f in zip(raw_w, is_fixed) if not f)
    remaining = max(10.0, TOTAL_TRACK_WIDTH - fixed_total)
    widths = [
        MIN_EXON_W if f else (w / raw_large_total * remaining)
        for w, f in zip(raw_w, is_fixed)
    ]
    exon_x: dict[int, tuple[float, float]] = {}
    x = LEFT_MARGIN
    for e, w in zip(exon_table, widths):
        exon_x[e["n"]] = (x, x + w)
        x += w
    return widths, exon_x


def _table_columns() -> tuple[float, float, float, float, float, float, float, float]:
    pct_x1 = LEFT_MARGIN - 0.4
    pct_x0 = pct_x1 - 1.8
    mw_x1 = pct_x0
    mw_x0 = mw_x1 - 1.4
    part_x1 = mw_x0
    part_x0 = part_x1 - 2.2
    group_x1 = part_x0
    group_x0 = group_x1 - 0.9
    return group_x0, group_x1, part_x0, part_x1, mw_x0, mw_x1, pct_x0, pct_x1


def _bar_label(row: dict[str, Any]) -> str:
    s = str(row.get("start_region", ""))
    t = str(row.get("stop_region", ""))
    if not s or not t:
        return "?"
    s_num = s.replace("e", "").replace("i", "i")
    t_num = t.replace("e", "").replace("i", "i")
    return s_num if s_num == t_num else f"{s_num}-{t_num}"


def _draw_transcript(ax: plt.Axes, widths: list[float], exon_x: dict[int, tuple[float, float]]) -> None:
    exon_table = get_exon_table()
    n_exons = len(exon_table)
    segments, _ = draw_transcript_color_band(
        ax, exon_table, exon_x, widths, EXON_ROW_Y, ROW_H,
    )
    for i, (e, w) in enumerate(zip(exon_table, widths)):
        x0 = exon_x[e["n"]][0]
        bump_depth = calculate_bump_depth(w * EXON_H, EXON_H) * (ROW_H / EXON_H)
        add_exon_outline(
            ax, x0, w, EXON_ROW_Y, ROW_H,
            e["five_prime"], e["three_prime"], e["parts"],
            bump_depth=bump_depth,
            extend_right_bump=(i < n_exons - 1),
            segments=segments,
            exon_index=i,
            n_exons=n_exons,
            exon_n=e["n"],
        )
        fs = (7 if w >= 0.7 else 6) * FONT_SCALE / 1.2
        ax.text(
            x0 + w / 2, EXON_ROW_Y + ROW_H / 2, str(e["n"]),
            ha="center", va="center", fontsize=fs, fontweight="bold", zorder=6,
        )

    ax.text(
        LEFT_MARGIN - 0.3, EXON_ROW_Y + ROW_H / 2, "Transcript \u2192",
        ha="right", va="center", fontsize=9 * FONT_SCALE / 1.1, fontweight="bold",
    )

    for exon, label in ISOFORM_ANNOTATIONS:
        xa = exon_x[exon][0]
        ax.annotate(
            label, xy=(xa, EXON_ROW_Y + ROW_H + 0.02),
            xytext=(xa, EXON_ROW_Y + ROW_H + 0.55),
            ha="center", fontsize=8 * FONT_SCALE / 1.1,
            arrowprops=dict(arrowstyle="-|>", lw=1),
        )


def _draw_hinge_guides(
    ax: plt.Axes,
    exon_x: dict[int, tuple[float, float]],
    y_bottom: float,
) -> None:
    y_top = DOMAIN_ROW_Y + ROW_H
    for xg in hinge_guide_x_positions(get_domain_map(), get_exon_table(), exon_x):
        ax.plot(
            [xg, xg], [y_top, y_bottom],
            linestyle=(0, (5, 4)), color="#666", lw=1.2, zorder=20,
        )


def _draw_table_headers(ax: plt.Axes, header_y: float) -> None:
    gx0, gx1, px0, px1, mx0, mx1, cx0, cx1 = _table_columns()
    header_h = 0.9
    for x0, x1, label in [
        (gx0, gx1, "Group"),
        (px0, px1, "Participant No."),
        (mx0, mx1, "MW"),
        (cx0, cx1, "%Dys(WB)"),
    ]:
        ax.add_patch(mpatches.Rectangle(
            (x0, header_y), x1 - x0, header_h,
            facecolor="white", edgecolor="black", lw=1.0, zorder=3,
        ))
        ax.text(
            (x0 + x1) / 2, header_y + header_h / 2, label,
            ha="center", va="center", fontsize=8.5 * FONT_SCALE / 1.1, fontweight="bold", zorder=4,
        )


def _draw_mutation_row(
    ax: plt.Axes,
    row: dict[str, Any],
    y: float,
    exon_x: dict[int, tuple[float, float]],
    *,
    bar_style: str,
) -> None:
    gx0, gx1, px0, px1, mx0, mx1, cx0, cx1 = _table_columns()
    table_right = cx1
    track_left = LEFT_MARGIN

    pid = str(row.get("id", ""))
    kda = row.get("expected_protein_size_kda", "")
    try:
        mw_text = f"{float(kda):.1f}"
    except (TypeError, ValueError):
        mw_text = str(kda) if kda else "—"

    pct_raw = row.get("pct_dys_wb", "")
    try:
        pct = float(pct_raw)
        pct_text = f"{pct:.2f}"
        pct_fill = pct_dys_color(pct)
    except (TypeError, ValueError):
        pct_text = str(pct_raw) if pct_raw else ""
        pct_fill = "#eeeeee"

    for x0, x1, text in [(px0, px1, pid), (mx0, mx1, mw_text)]:
        ax.add_patch(mpatches.Rectangle(
            (x0, y - PATIENT_ROW_H), x1 - x0, PATIENT_ROW_H,
            facecolor="white", edgecolor="black", lw=0.6, zorder=2,
        ))
        ax.text((x0 + x1) / 2, y - PATIENT_ROW_H / 2, text, ha="center", va="center", fontsize=8 * FONT_SCALE / 1.1)

    ax.add_patch(mpatches.Rectangle(
        (cx0, y - PATIENT_ROW_H), cx1 - cx0, PATIENT_ROW_H,
        facecolor=pct_fill, edgecolor="black", lw=0.6, zorder=2,
    ))
    if pct_text:
        ax.text(
            (cx0 + cx1) / 2, y - PATIENT_ROW_H / 2, pct_text,
            ha="center", va="center", fontsize=7.5 * FONT_SCALE / 1.1, fontweight="bold",
        )

    transcript_right = exon_x[N_EXONS][1]
    ax.add_patch(mpatches.Rectangle(
        (table_right, y - PATIENT_ROW_H), track_left - table_right, PATIENT_ROW_H,
        facecolor="white", edgecolor="none", zorder=1,
    ))
    ax.add_patch(mpatches.Rectangle(
        (track_left, y - PATIENT_ROW_H), transcript_right - track_left, PATIENT_ROW_H,
        facecolor="white", edgecolor="none", zorder=1,
    ))

    rng = exon_range_from_row(row)
    if not rng:
        return
    first, last = rng
    x0b, x1b = exon_x[first][0], exon_x[last][1]
    label = _bar_label(row)

    bx, by, bw, bh = _mutation_bar_box(x0b, x1b, y)
    face = STYLE_FACE[bar_style]
    tc = "white" if bar_style == "black" else "black"
    ax.add_patch(mpatches.Rectangle(
        (bx, by), bw, bh,
        facecolor=face, edgecolor="black", lw=0.8, zorder=3,
    ))
    ax.text(
        (x0b + x1b) / 2, y - PATIENT_ROW_H / 2, label,
        ha="center", va="center", fontsize=7 * FONT_SCALE / 1.1, color=tc, fontweight="bold", zorder=4,
    )


def _draw_group_column(
    ax: plt.Axes,
    rows: list[dict[str, Any]],
    y_start: float,
    transcript_right: float,
) -> None:
    gx0, gx1, *_ = _table_columns()
    table_x0 = gx0
    group_bounds: dict[str, list[float]] = {}
    y = y_start
    for row in rows:
        grp = str(row.get("group", "") or "—")
        group_bounds.setdefault(grp, [y, y - PATIENT_ROW_H])
        group_bounds[grp][1] = y - PATIENT_ROW_H
        y -= PATIENT_ROW_H

    for grp, (ytop, ybot) in group_bounds.items():
        ax.add_patch(mpatches.Rectangle(
            (gx0, ybot), gx1 - gx0, ytop - ybot,
            facecolor="white", edgecolor="black", lw=1.0, zorder=2,
        ))
        ax.text(
            (gx0 + gx1) / 2, (ytop + ybot) / 2, grp,
            ha="center", va="center", fontsize=15 * FONT_SCALE / 1.1, fontweight="bold", zorder=4,
        )
        ax.plot([table_x0, transcript_right], [ytop, ytop], color="black", lw=1.4, zorder=5)


def _draw_phenotype_boxes(
    ax: plt.Axes,
    group_bounds: dict[str, list[float]],
    transcript_right: float,
) -> None:
    """Right-hand phenotype legend spanning groups D+E (reference figure)."""
    if "D" not in group_bounds:
        return
    d_top, d_bot = group_bounds["D"]
    if "E" in group_bounds:
        d_bot = group_bounds["E"][1]
    pheno_x0 = transcript_right + 1.0
    pheno_x1 = pheno_x0 + 13.0
    labels = ["Asymptomatic/Paucisymptomatic*", "BMD/IMD\u2021", "DMD"]
    pheno_h = (d_top - d_bot) / len(labels)
    ax.add_patch(mpatches.Rectangle(
        (pheno_x0, d_bot), pheno_x1 - pheno_x0, d_top - d_bot,
        facecolor="none", edgecolor="black", lw=1.4, zorder=6,
    ))
    yy = d_top
    for lbl in labels:
        ax.plot([pheno_x0, pheno_x1], [yy - pheno_h, yy - pheno_h], color="black", lw=0.8)
        ax.text(
            (pheno_x0 + pheno_x1) / 2, yy - pheno_h / 2, lbl,
            ha="center", va="center", fontsize=9 * FONT_SCALE / 1.1,
        )
        yy -= pheno_h


def create_cohort_mutation_map(
    selected_rows: list[dict[str, Any]],
    exons: list[ExonRecord],
    *,
    preview_row: Optional[dict[str, Any]] = None,
    title: str = "DMD mutation map",
) -> Figure:
    rows = list(selected_rows)
    if preview_row and preview_row.get("start_region"):
        rows = [preview_row] + [r for r in rows if r.get("id") != preview_row.get("id")]

    widths, exon_x = _compute_exon_layout()
    transcript_right = exon_x[N_EXONS][1]
    n_patients = len(rows)
    table_bottom_y = TABLE_TOP_Y - n_patients * ROW_STEP
    table_x0 = _table_columns()[0]

    data_w = transcript_right - table_x0 + 14
    data_h = (EXON_ROW_Y + ROW_H + 1.4) - table_bottom_y + 1.0
    fig_w = max(28, data_w * INCHES_PER_DATA_UNIT)
    fig_h = max(10, min(15, data_h * INCHES_PER_DATA_UNIT + 2))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
    ax.set_xlim(table_x0 - 0.5, transcript_right + 14)
    ax.set_ylim(table_bottom_y - 1.0, EXON_ROW_Y + ROW_H + 1.4)
    ax.axis("off")

    _draw_transcript(ax, widths, exon_x)
    resolved_domains = _draw_domains(ax, exon_x)
    ax.plot(
        [LEFT_MARGIN, LEFT_MARGIN],
        [DOMAIN_ROW_Y, EXON_ROW_Y + ROW_H],
        color="black", lw=1.2, zorder=6,
    )
    if SHOW_ALIGNMENT_GUIDES:
        _draw_alignment_guides_mpl(ax, resolved_domains, EXON_ROW_Y, table_bottom_y)

    group_bounds: dict[str, list[float]] = {}
    if rows:
        _draw_table_headers(ax, TABLE_TOP_Y + 0.05)
        y = TABLE_TOP_Y
        for row in rows:
            style = bar_style_for_phenotype(str(row.get("phenotype", "")))
            _draw_mutation_row(ax, row, y, exon_x, bar_style=style)
            grp = str(row.get("group", "") or "—")
            group_bounds.setdefault(grp, [y, y - PATIENT_ROW_H])
            group_bounds[grp][1] = y - PATIENT_ROW_H
            y -= PATIENT_ROW_H
        _draw_group_column(ax, rows, TABLE_TOP_Y, transcript_right)

        header_y = TABLE_TOP_Y + 0.05
        header_h = 0.9
        table_x0 = _table_columns()[0]
        _, _, _, _, _, _, cx0, cx1 = _table_columns()
        table_bottom_y_rows = TABLE_TOP_Y - len(rows) * ROW_STEP
        ax.add_patch(mpatches.Rectangle(
            (table_x0, table_bottom_y_rows), cx1 - table_x0, (header_y + header_h) - table_bottom_y_rows,
            facecolor="none", edgecolor="black", lw=1.6, zorder=6,
        ))
        ax.plot(
            [LEFT_MARGIN, LEFT_MARGIN],
            [DOMAIN_ROW_Y, table_bottom_y_rows],
            color="black", lw=1.2, zorder=6,
        )
        _draw_phenotype_boxes(ax, group_bounds, transcript_right)
    else:
        ax.text(
            LEFT_MARGIN + 30, TABLE_TOP_Y - 0.5,
            "Select mutations below and click Plot selected on map",
            ha="center", fontsize=10, color="#666",
        )

    _draw_hinge_guides(ax, exon_x, table_bottom_y)

    ax.set_title(title, fontsize=12 * FONT_SCALE / 1.1, fontweight="bold", pad=10)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.94, bottom=0.04)
    return fig
