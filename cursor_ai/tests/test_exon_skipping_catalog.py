"""Tests for catalog-aware exon-skipping analysis."""

from __future__ import annotations

import pandas as pd
import pytest

from src.coordinate_mapper import build_exon_records
from src.exon_skipping import find_skip_candidates
from src.exon_skipping_catalog import (
    analyze_mutations_for_skip_target,
    format_exon_list,
    format_skip_candidate,
    rank_skip_candidates_for_row,
    validate_skip_target,
)
from src.frame_analysis import assess_deletion_frame
from src.models import FrameStatus
from src.mutation_catalog import CATALOG_COLUMNS, ensure_record_ids
from src.mutation_map_svg import _compute_layout, build_interactive_map_html, build_skip_target_map_html
from src.mutation_viz import exon_sort_key_from_row, sort_catalog_dataframe


def _catalog_row(
    *,
    participant: str = "S1",
    msub: str = "Deletion",
    start: str = "e45",
    stop: str = "e50",
) -> dict[str, str]:
    row = {col: "" for col in CATALOG_COLUMNS}
    row.update({
        "id": participant,
        "mutation_class": "Exonic",
        "mutation_subclass": msub,
        "start_region": start,
        "stop_region": stop,
    })
    return row


def _catalog_df(*rows: dict[str, str]) -> pd.DataFrame:
    df, _ = ensure_record_ids(pd.DataFrame(list(rows), columns=CATALOG_COLUMNS))
    return df


def test_single_exon_target_validation() -> None:
    exons, err = validate_skip_target(45, 45)
    assert err is None
    assert exons == [45]


def test_multi_exon_target_validation() -> None:
    exons, err = validate_skip_target(45, 47)
    assert err is None
    assert exons == [45, 46, 47]


def test_invalid_target_rejected() -> None:
    _, err = validate_skip_target(80, 80)
    assert err is not None
    _, err2 = validate_skip_target(50, 45)
    assert err2 is not None


def test_out_of_frame_deletion_restored_by_single_skip() -> None:
    exons = build_exon_records()
    frame = assess_deletion_frame(45, 50, exons)
    assert frame.status == FrameStatus.OUT_OF_FRAME
    candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
    assert candidates
    assert candidates[0].restores_frame


def test_in_frame_deletion_not_in_target_results_by_default() -> None:
    exons = build_exon_records()
    df = _catalog_df(_catalog_row(participant="S1", start="e3", stop="e5"))
    restored = analyze_mutations_for_skip_target(df, [2], exons)
    assert all(r.original_frame != "In Frame" for r in restored)


def test_catalog_table_sorts_by_exon_range() -> None:
    df = _catalog_df(
        _catalog_row(participant="B", start="e10", stop="e10"),
        _catalog_row(participant="A", start="e3", stop="e3"),
    )
    sorted_df = sort_catalog_dataframe(df)
    keys = [exon_sort_key_from_row(r.to_dict()) for _, r in sorted_df.iterrows()]
    assert keys == [(3, 3), (10, 10)]


def test_record_id_trailing_in_display_columns() -> None:
    from src.mutation_catalog import DISPLAY_COLUMNS

    assert DISPLAY_COLUMNS[-1] == "mutation_record_id"


def test_single_exon_duplication_forces_candidate_one() -> None:
    exons = build_exon_records()
    row = _catalog_row(msub="Duplication", start="e2", stop="e2")
    frame, candidates = rank_skip_candidates_for_row(row, exons, max_additional_skips=3)
    assert frame.status == FrameStatus.OUT_OF_FRAME
    assert candidates
    assert format_skip_candidate(candidates[0]) == "Exon 2"


def test_multi_exon_duplication_no_forced_first_candidate() -> None:
    exons = build_exon_records()
    row = _catalog_row(msub="Duplication", start="e8", stop="e11")
    frame, candidates = rank_skip_candidates_for_row(row, exons, max_additional_skips=3)
    if candidates and frame.status == FrameStatus.OUT_OF_FRAME:
        first = format_skip_candidate(candidates[0])
        assert first != "Exon 8" or len(candidates) == 1


def test_deletion_and_duplication_3_7_use_different_candidates() -> None:
    exons = build_exon_records()
    duplication = _catalog_row(
        participant="S42",
        msub="Duplication",
        start="e3",
        stop="e7",
    )
    deletion = _catalog_row(
        participant="S79",
        msub="Deletion",
        start="e3",
        stop="e7",
    )

    _, duplication_candidates = rank_skip_candidates_for_row(
        duplication,
        exons,
        max_additional_skips=3,
    )
    _, deletion_candidates = rank_skip_candidates_for_row(
        deletion,
        exons,
        max_additional_skips=3,
    )

    assert duplication_candidates
    assert deletion_candidates
    duplicated_interval = set(range(3, 8))
    assert all(
        duplicated_interval.issubset(candidate.additional_skipped_exons)
        for candidate in duplication_candidates
    )
    assert duplication_candidates[0].additional_skipped_exons == [3, 4, 5, 6, 7]
    assert any(
        candidate.additional_skipped_exons == [8, 9]
        for candidate in deletion_candidates
    )
    assert duplication_candidates[0].additional_skipped_exons != (
        deletion_candidates[0].additional_skipped_exons
    )


def test_empty_target_results_message_data() -> None:
    exons = build_exon_records()
    df = _catalog_df(_catalog_row(participant="S99", start="e1", stop="e1"))
    restored = analyze_mutations_for_skip_target(df, [78], exons)
    assert restored == []


def test_target_results_separate_id_and_preserve_candidate_rank() -> None:
    exons = build_exon_records()
    df = _catalog_df(
        _catalog_row(
            participant="S7",
            start="e7",
            stop="e7",
        )
    )

    restored = analyze_mutations_for_skip_target(
        df,
        [6],
        exons,
        target_mode="multi",
    )

    assert restored
    assert all(item.participant_id == "S7" for item in restored)
    assert all("S7" not in item.mutation_label for item in restored)
    assert [item.skip_combination_rank for item in restored[:3]] == [1, 2, 3]
    assert [item.skipped_exons for item in restored[:3]] == [
        [6, 8],
        [5, 6, 8],
        [6, 8, 9],
    ]


def test_schematic_reuses_shared_layout() -> None:
    exons = build_exon_records()
    _ = exons
    widths_a, layout_a = _compute_layout(2000)
    html = build_interactive_map_html([], chart_width=2000)
    html_target = build_skip_target_map_html([45], chart_width=2000)
    assert "viewBox" in html
    assert "viewBox" in html_target
    widths_b, layout_b = _compute_layout(2000)
    assert layout_a[45] == layout_b[45]
    assert widths_a == widths_b


def test_format_exon_list() -> None:
    assert format_exon_list([45]) == "Exon 45"
    assert format_exon_list([45, 46, 47]) == "Exons 45-47"
