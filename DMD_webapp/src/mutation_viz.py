"""Map catalog mutation rows to exon visualization state."""

from __future__ import annotations

from typing import Any, Optional

import pandas as pd

from src.config import REFERENCE
from src.region_parser import parse_region_text
from src.visualization import VisualizationState

N_EXONS = REFERENCE.coding_exon_count
_UNPARSEABLE_SORT_KEY = (N_EXONS + 1, N_EXONS + 1)


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


def exon_sort_key_from_row(row: dict[str, Any]) -> tuple[int, int]:
    """Sort key (start_exon, end_exon) for catalog rows and map plotting."""
    rng = exon_range_from_row(row)
    if rng is None:
        return _UNPARSEABLE_SORT_KEY
    return rng


def sort_mutations_by_exon_range(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return mutations sorted by (start_exon, end_exon) ascending."""
    return sorted(rows, key=exon_sort_key_from_row)


def sort_catalog_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Sort catalog rows by (start_exon, end_exon) ascending."""
    if df.empty:
        return df
    out = df.copy()
    keys = out.apply(lambda row: exon_sort_key_from_row(row.to_dict()), axis=1)
    out["_sort_start"] = [k[0] for k in keys]
    out["_sort_end"] = [k[1] for k in keys]
    out = out.sort_values(["_sort_start", "_sort_end"], kind="stable")
    return out.drop(columns=["_sort_start", "_sort_end"]).reset_index(drop=True)


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
