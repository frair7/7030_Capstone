"""
Plotly-based DMD transcript and protein-domain visualizations.

Provides:
  - Exon-order schematic (default): every exon equally visible
  - Genomic/transcript-scale view: proportional intron representation
  - Mutation, deletion, duplication, and skip-candidate highlighting
  - Aligned protein-domain track
  - Accessible data table for the same information
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.coordinate_mapper import build_exon_records
from src.models import ExonRecord, SkipCandidate
from src.reference_data import load_domain_table

# Colour palette — original, restrained scientific style
COLORS = {
    "retained": "#4C78A8",
    "default": "#9DB4C8",
    "deleted": "#D62728",
    "duplicated": "#F28E2B",
    "skip_candidate": "#9467BD",
    "mutation": "#EECA3B",
    "intron": "#E8E8E8",
    "domain_bg": "#F7F7F7",
}

DOMAIN_COLORS = {
    "Actin-binding domain": "#1f4e79",
    "Central rod domain": "#6baed6",
    "Hinge regions": "#fd8d3c",
    "Cysteine-rich domain": "#DAE9F8",
    "Carboxy-terminal domain": "#9e9ac8",
}


class ViewMode(str, Enum):
    """Transcript map display mode."""

    SCHEMATIC = "schematic"
    GENOMIC = "genomic"


class ExonState(str, Enum):
    """Visual state for an exon rectangle."""

    DEFAULT = "default"
    RETAINED = "retained"
    DELETED = "deleted"
    DUPLICATED = "duplicated"
    SKIP_CANDIDATE = "skip_candidate"
    MUTATION = "mutation"


@dataclass
class ExonLayout:
    """Computed layout for one exon rectangle."""

    exon_number: int
    x_start: float
    x_end: float
    y_center: float = 0.5
    height: float = 0.7


@dataclass
class VisualizationState:
    """Highlighting state passed to the visualization layer."""

    deleted_exons: set[int] = field(default_factory=set)
    duplicated_exons: set[int] = field(default_factory=set)
    skip_candidate_exons: set[int] = field(default_factory=set)
    mutation_exons: set[int] = field(default_factory=set)
    retained_exons: set[int] = field(default_factory=set)

    def state_for_exon(self, exon_number: int) -> ExonState:
        if exon_number in self.mutation_exons:
            return ExonState.MUTATION
        if exon_number in self.deleted_exons:
            return ExonState.DELETED
        if exon_number in self.duplicated_exons:
            return ExonState.DUPLICATED
        if exon_number in self.skip_candidate_exons:
            return ExonState.SKIP_CANDIDATE
        if exon_number in self.retained_exons:
            return ExonState.RETAINED
        return ExonState.DEFAULT


def _exon_state_color(state: ExonState) -> str:
    return COLORS.get(state.value, COLORS["default"])


def compute_exon_layouts(
    exons: list[ExonRecord],
    view_mode: ViewMode = ViewMode.SCHEMATIC,
    *,
    schematic_width: float = 1.0,
    schematic_gap: float = 0.15,
) -> list[ExonLayout]:
    """
    Compute x-axis positions for exons in transcript order.

    Schematic mode assigns equal width to every exon.
    Genomic mode uses transcript coordinates (includes compressed intron gaps).
    """
    layouts: list[ExonLayout] = []

    if view_mode == ViewMode.SCHEMATIC:
        cursor = 0.0
        for exon in exons:
            x_start = cursor
            x_end = cursor + schematic_width
            layouts.append(
                ExonLayout(
                    exon_number=exon.exon_number,
                    x_start=x_start,
                    x_end=x_end,
                )
            )
            cursor = x_end + schematic_gap
        return layouts

    # Genomic / transcript-scale: use transcript coordinates from reference table
    for exon in exons:
        layouts.append(
            ExonLayout(
                exon_number=exon.exon_number,
                x_start=float(exon.transcript_exon_start),
                x_end=float(exon.transcript_exon_end),
            )
        )
    return layouts


def _hover_text(exon: ExonRecord) -> str:
    return (
        f"<b>Exon {exon.exon_number}</b><br>"
        f"Genomic: chr{exon.chromosome}:{exon.genomic_start:,}-"
        f"{exon.genomic_end:,}<br>"
        f"Transcript: {exon.transcript_exon_start}-{exon.transcript_exon_end}<br>"
        f"CDS: {exon.cds_start or '—'}-{exon.cds_end or '—'}<br>"
        f"Coding length: {exon.coding_length_bp} bp<br>"
        f"Phases (5′/3′): {exon.splice_phase_5prime}/"
        f"{exon.splice_phase_3prime}"
    )


def build_exon_table(
    exons: list[ExonRecord],
    viz_state: Optional[VisualizationState] = None,
) -> pd.DataFrame:
    """Accessible table mirroring the graphic."""
    state = viz_state or VisualizationState()
    rows = []
    for exon in exons:
        rows.append(
            {
                "exon_number": exon.exon_number,
                "genomic_start": exon.genomic_start,
                "genomic_end": exon.genomic_end,
                "transcript_start": exon.transcript_exon_start,
                "transcript_end": exon.transcript_exon_end,
                "cds_start": exon.cds_start,
                "cds_end": exon.cds_end,
                "coding_length_bp": exon.coding_length_bp,
                "phase_5prime": exon.splice_phase_5prime,
                "phase_3prime": exon.splice_phase_3prime,
                "visual_state": state.state_for_exon(exon.exon_number).value,
            }
        )
    return pd.DataFrame(rows)


def create_transcript_figure(
    exons: Optional[list[ExonRecord]] = None,
    viz_state: Optional[VisualizationState] = None,
    view_mode: ViewMode = ViewMode.SCHEMATIC,
    *,
    title: str = "DMD Dp427m transcript (exons 1–79, 5′ → 3′)",
    show_exon_labels: bool = True,
    height: int = 220,
) -> go.Figure:
    """Build the exon map as a standalone Plotly figure."""
    records = exons if exons is not None else build_exon_records()
    state = viz_state or VisualizationState()
    layouts = compute_exon_layouts(records, view_mode)
    index = {e.exon_number: e for e in records}

    fig = go.Figure()
    hover_x: list[float] = []
    hover_y: list[float] = []
    hover_text: list[str] = []
    hover_custom: list[str] = []

    for layout in layouts:
        exon = index[layout.exon_number]
        exon_state = state.state_for_exon(exon.exon_number)
        color = _exon_state_color(exon_state)
        line_color = "#333333"
        line_width = 1.5
        if exon_state == ExonState.MUTATION:
            line_color = COLORS["mutation"]
            line_width = 3

        y0 = layout.y_center - layout.height / 2
        y1 = layout.y_center + layout.height / 2
        fig.add_shape(
            type="rect",
            x0=layout.x_start,
            x1=layout.x_end,
            y0=y0,
            y1=y1,
            fillcolor=color,
            line=dict(color=line_color, width=line_width),
            layer="below",
        )

        if show_exon_labels and (view_mode == ViewMode.SCHEMATIC or exon.exon_number % 5 == 0 or exon.exon_number in {1, 79}):
            fig.add_annotation(
                x=(layout.x_start + layout.x_end) / 2,
                y=layout.y_center,
                text=str(exon.exon_number),
                showarrow=False,
                font=dict(size=7, color="#1a1a1a"),
            )

        cx = (layout.x_start + layout.x_end) / 2
        hover_x.append(cx)
        hover_y.append(layout.y_center)
        hover_text.append(_hover_text(exon))
        hover_custom.append(exon_state.value)

        # Draw intron connector in schematic mode
        if view_mode == ViewMode.SCHEMATIC and layout.exon_number < len(layouts):
            next_layout = layouts[layout.exon_number]  # exon_number is 1-based
            fig.add_shape(
                type="line",
                x0=layout.x_end,
                x1=next_layout.x_start,
                y0=layout.y_center,
                y1=layout.y_center,
                line=dict(color=COLORS["intron"], width=2),
                layer="below",
            )

    fig.add_trace(
        go.Scatter(
            x=hover_x,
            y=hover_y,
            mode="markers",
            marker=dict(size=12, opacity=0),
            text=hover_text,
            customdata=hover_custom,
            hovertemplate="%{text}<extra>State: %{customdata}</extra>",
            name="Exons",
        )
    )

    x_title = (
        "Exon order (5′ → 3′, schematic)"
        if view_mode == ViewMode.SCHEMATIC
        else "Transcript position (bp, incl. introns)"
    )
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        height=height,
        margin=dict(l=40, r=20, t=50, b=40),
        xaxis=dict(title=x_title, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, range=[0, 1]),
        plot_bgcolor="white",
        showlegend=False,
        dragmode="pan",
    )
    return fig


def create_domain_figure(
    domains: Optional[list[dict[str, Any]]] = None,
    *,
    protein_length: int = 3685,
    title: str = "Dystrophin protein domains (amino-acid coordinates)",
    height: int = 160,
) -> go.Figure:
    """Build aligned protein-domain track."""
    domain_rows = domains if domains is not None else load_domain_table()
    fig = go.Figure()

    y_positions = {
        "Actin-binding domain": 0.75,
        "Central rod domain": 0.55,
        "Hinge regions": 0.35,
        "Cysteine-rich domain": 0.55,
        "Carboxy-terminal domain": 0.75,
    }

    for row in domain_rows:
        name = row["domain_name"]
        aa_start = int(row["protein_start_aa"])
        aa_end = int(row["protein_end_aa"])
        color = DOMAIN_COLORS.get(name, "#888888")
        y = y_positions.get(name, 0.5)
        fig.add_shape(
            type="rect",
            x0=aa_start,
            x1=aa_end,
            y0=y - 0.12,
            y1=y + 0.12,
            fillcolor=color,
            opacity=0.85,
            line=dict(color="#333", width=1),
            layer="below",
        )
        fig.add_annotation(
            x=(aa_start + aa_end) / 2,
            y=y,
            text=name.replace(" domain", ""),
            showarrow=False,
            font=dict(size=8, color="white"),
        )

    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        height=height,
        margin=dict(l=40, r=20, t=50, b=40),
        xaxis=dict(
            title="Amino-acid position",
            range=[0, protein_length + 50],
            showgrid=True,
            gridcolor="#eee",
        ),
        yaxis=dict(visible=False, range=[0, 1]),
        plot_bgcolor=COLORS["domain_bg"],
        showlegend=False,
        dragmode="pan",
    )
    return fig


def create_combined_figure(
    exons: Optional[list[ExonRecord]] = None,
    viz_state: Optional[VisualizationState] = None,
    view_mode: ViewMode = ViewMode.SCHEMATIC,
    *,
    domains: Optional[list[dict[str, Any]]] = None,
    title: str = "DMD Mutation and Exon Map",
) -> go.Figure:
    """
    Combined transcript + domain figure with two aligned rows.

    In genomic view mode the top panel uses transcript coordinates; the bottom
    panel uses amino-acid coordinates (separate scales, vertically aligned).
    """
    records = exons if exons is not None else build_exon_records()
    state = viz_state or VisualizationState()
    domain_rows = domains if domains is not None else load_domain_table()
    layouts = compute_exon_layouts(records, view_mode)
    index = {e.exon_number: e for e in records}

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=False,
        row_heights=[0.55, 0.45],
        vertical_spacing=0.12,
        subplot_titles=(
            "Dp427m transcript (exons 1–79, 5′ → 3′)",
            "Protein domains (amino acids)",
        ),
    )

    # --- Row 1: exons ---
    for layout in layouts:
        exon = index[layout.exon_number]
        exon_state = state.state_for_exon(exon.exon_number)
        color = _exon_state_color(exon_state)
        y0, y1 = 0.15, 0.85
        fig.add_shape(
            type="rect",
            x0=layout.x_start,
            x1=layout.x_end,
            y0=y0,
            y1=y1,
            fillcolor=color,
            line=dict(color="#333", width=1),
            layer="below",
            row=1,
            col=1,
        )
        if view_mode == ViewMode.SCHEMATIC or exon.exon_number % 10 == 0 or exon.exon_number in {1, 79}:
            fig.add_annotation(
                x=(layout.x_start + layout.x_end) / 2,
                y=0.5,
                text=str(exon.exon_number),
                showarrow=False,
                font=dict(size=7),
                row=1,
                col=1,
            )
        fig.add_trace(
            go.Scatter(
                x=[(layout.x_start + layout.x_end) / 2],
                y=[0.5],
                mode="markers",
                marker=dict(size=10, opacity=0),
                text=[_hover_text(exon)],
                customdata=[exon_state.value],
                hovertemplate="%{text}<extra>State: %{customdata}</extra>",
                showlegend=False,
            ),
            row=1,
            col=1,
        )

    # --- Row 2: domains ---
    for row in domain_rows:
        name = row["domain_name"]
        aa_start = int(row["protein_start_aa"])
        aa_end = int(row["protein_end_aa"])
        color = DOMAIN_COLORS.get(name, "#888")
        fig.add_shape(
            type="rect",
            x0=aa_start,
            x1=aa_end,
            y0=0.35,
            y1=0.65,
            fillcolor=color,
            opacity=0.85,
            line=dict(color="#333", width=1),
            layer="below",
            row=2,
            col=1,
        )

    x1_title = (
        "Exon order (schematic)"
        if view_mode == ViewMode.SCHEMATIC
        else "Transcript position (bp)"
    )
    fig.update_layout(
        title=dict(text=title, font=dict(size=15)),
        height=480,
        margin=dict(l=50, r=20, t=80, b=40),
        plot_bgcolor="white",
        showlegend=False,
        dragmode="pan",
    )
    fig.update_xaxes(title_text=x1_title, row=1, col=1, showgrid=False)
    fig.update_xaxes(title_text="Amino-acid position", row=2, col=1, showgrid=True)
    fig.update_yaxes(visible=False, row=1, col=1)
    fig.update_yaxes(visible=False, row=2, col=1)

    return fig


def visualization_state_from_variant(
    first_exon: Optional[int],
    last_exon: Optional[int],
    *,
    variant_type: str = "deletion",
    skip_candidate: Optional[SkipCandidate] = None,
) -> VisualizationState:
    """Build a VisualizationState from a parsed variant and optional skip candidate."""
    state = VisualizationState()
    if first_exon is None or last_exon is None:
        return state

    affected = set(range(first_exon, last_exon + 1))
    if variant_type == "deletion":
        state.deleted_exons = affected
    elif variant_type == "duplication":
        state.duplicated_exons = affected
    else:
        state.mutation_exons = affected

    if skip_candidate:
        state.skip_candidate_exons = set(skip_candidate.additional_skipped_exons)
        # Retained = not deleted by mutation or additional skip
        all_removed = set(skip_candidate.deleted_exons) | set(
            skip_candidate.additional_skipped_exons
        )
        state.retained_exons = {
            e.exon_number
            for e in build_exon_records()
            if e.exon_number not in all_removed
        }

    return state


def legend_items() -> list[dict[str, str]]:
    """Return legend swatches for the UI."""
    return [
        {"label": "Retained exon", "color": COLORS["retained"]},
        {"label": "Deleted exon", "color": COLORS["deleted"]},
        {"label": "Duplicated exon", "color": COLORS["duplicated"]},
        {"label": "Proposed skip (computational)", "color": COLORS["skip_candidate"]},
        {"label": "Mutation highlight", "color": COLORS["mutation"]},
    ]
