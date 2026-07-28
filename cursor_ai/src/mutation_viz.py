"""Map catalog mutation rows to exon visualization state."""

from __future__ import annotations

from typing import Any, Optional

from src.region_parser import parse_region_text
from src.visualization import VisualizationState


def exon_range_from_row(row: dict[str, Any]) -> Optional[tuple[int, int]]:
    """Extract inclusive exon range when both bounds are exonic."""
    start = parse_region_text(str(row.get("start_region", "")))
    stop = parse_region_text(str(row.get("stop_region", "")))
    first = start.start_exon or stop.start_exon
    last = stop.stop_exon or start.stop_exon or start.stop_exon
    if first is None or last is None:
        # Try single exon from either field
        if start.start_exon:
            return start.start_exon, start.start_exon
        if stop.stop_exon:
            return stop.stop_exon, stop.stop_exon
        return None
    if first > last:
        first, last = last, first
    return first, last


def viz_state_from_row(row: dict[str, Any]) -> VisualizationState:
    """Build highlight state for one catalog row."""
    state = VisualizationState()
    rng = exon_range_from_row(row)
    if not rng:
        return state

    first, last = rng
    affected = set(range(first, last + 1))
    mclass = str(row.get("mutation_class", "")).lower()
    msub = str(row.get("mutation_subclass", "")).lower()

    if msub == "duplication":
        state.duplicated_exons = affected
    elif msub == "deletion":
        state.deleted_exons = affected
    else:
        state.mutation_exons = affected

    return state


def viz_state_from_rows(rows: list[dict[str, Any]]) -> VisualizationState:
    """Merge highlight state from multiple selected rows."""
    merged = VisualizationState()
    for row in rows:
        partial = viz_state_from_row(row)
        merged.deleted_exons |= partial.deleted_exons
        merged.duplicated_exons |= partial.duplicated_exons
        merged.mutation_exons |= partial.mutation_exons
        merged.skip_candidate_exons |= partial.skip_candidate_exons
        merged.retained_exons |= partial.retained_exons
    return merged
