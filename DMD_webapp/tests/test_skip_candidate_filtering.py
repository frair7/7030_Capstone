"""Tests for frame-restoring-only skip candidate filtering and ranking."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import build_exon_records
from src.exon_skipping import (
    MAX_DISPLAY_CANDIDATES,
    evaluate_candidate,
    find_skip_candidates,
    ranking_key,
)
from src.exon_skipping_catalog import (
    overlay_from_row_and_candidate,
    rank_skip_candidates_for_row,
)
from src.models import FrameStatus
from src.mutation_catalog import CATALOG_COLUMNS
from src.reporting import skip_candidates_dataframe
from src.transcript_reconstruction import (
    evaluate_all_strategies,
    reconstruct_deletion_skip,
)


@pytest.fixture(scope="module")
def exons():
    return build_exon_records()


def _catalog_row() -> dict[str, str]:
    row = {col: "" for col in CATALOG_COLUMNS}
    row.update({
        "id": "S1",
        "mutation_class": "Exonic",
        "mutation_subclass": "Deletion",
        "start_region": "e45",
        "stop_region": "e50",
    })
    return row


class TestFrameRestoringFilter:
    def test_non_restoring_excluded_from_final_list(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        assert all(c.additional_skipped_exons != [44] for c in candidates)
        transcript = reconstruct_deletion_skip(45, 50, (44,), exons, require_frame=False)
        assert transcript is not None
        assert not transcript.restores_frame

    def test_non_restoring_receives_no_combo_number(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        df = skip_candidates_dataframe(candidates)
        labels = set(df["candidate"].tolist())
        assert "Skip combo 4" not in labels
        assert evaluate_candidate(45, 50, (44,), exons, "del45-50") is None

    def test_non_restoring_not_in_outline(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        overlay = overlay_from_row_and_candidate(
            _catalog_row(), candidates[0], all_candidates=candidates,
        )
        assert 44 not in overlay.skip_candidate_exons[0]
        for combo_set in overlay.skip_candidate_exons:
            assert combo_set != {44}

    def test_candidate_count_only_frame_restoring(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        strategies = evaluate_all_strategies(45, 50, exons, max_additional_skips=3)
        restoring_count = sum(1 for s in strategies if s.restores_frame)
        assert len(candidates) == min(restoring_count, MAX_DISPLAY_CANDIDATES)
        assert all(c.restores_frame for c in candidates)

    def test_valid_candidates_renumbered_consecutively(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        df = skip_candidates_dataframe(candidates)
        expected = [f"Skip combo {i + 1}" for i in range(len(candidates))]
        assert df["candidate"].tolist() == expected
        assert candidates[0].additional_skipped_exons == [51]

    def test_ranking_after_frame_filter(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        for i in range(len(candidates) - 1):
            assert ranking_key(candidates[i]) <= ranking_key(candidates[i + 1])
        non_restoring = evaluate_candidate(45, 50, (44,), exons, "del45-50")
        assert non_restoring is None

    def test_fewer_than_three_candidates_allowed(self, exons) -> None:
        candidates = find_skip_candidates(
            45, 50, exons, max_additional_skips=1, max_results=3,
        )
        assert 0 < len(candidates) <= 1

    def test_zero_candidates_message_condition(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=0)
        assert candidates == []

    def test_max_additional_skips_is_search_depth(self, exons) -> None:
        shallow = find_skip_candidates(45, 50, exons, max_additional_skips=1)
        deep = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        assert len(deep) >= len(shallow)
        assert len(deep) <= MAX_DISPLAY_CANDIDATES

    def test_every_displayed_candidate_restores_frame(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        assert all(c.restores_frame is True for c in candidates)

    def test_every_displayed_candidate_structurally_valid(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        assert all(c.is_structurally_valid for c in candidates)

    def test_retained_exons_absent_from_display_table(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        df = skip_candidates_dataframe(candidates)
        assert "retained_exons" not in df.columns
        assert candidates[0].retained_exons

    def test_labels_follow_post_filter_ranking(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        overlay = overlay_from_row_and_candidate(
            _catalog_row(), candidates[0], all_candidates=candidates,
        )
        assert overlay.skip_candidate_exons[0] == {51}
        if len(candidates) > 1:
            assert overlay.skip_candidate_exons[1] == set(candidates[1].additional_skipped_exons)

    def test_advanced_noncontiguous_must_restore_frame(self, exons) -> None:
        candidates = find_skip_candidates(
            45, 50, exons, max_additional_skips=3, advanced=True,
        )
        assert all(c.restores_frame for c in candidates)
        assert all(not c.is_advanced_noncontiguous or c.restores_frame for c in candidates)

    def test_duplication_forced_candidate_restores_frame(self, exons) -> None:
        row = {col: "" for col in CATALOG_COLUMNS}
        row.update({
            "id": "S2",
            "mutation_class": "Exonic",
            "mutation_subclass": "Duplication",
            "start_region": "e2",
            "stop_region": "e2",
        })
        frame, candidates = rank_skip_candidates_for_row(row, exons, max_additional_skips=3)
        assert frame.status == FrameStatus.OUT_OF_FRAME
        assert candidates
        assert candidates[0].restores_frame
        assert candidates[0].additional_skipped_exons == [2]


class TestDeletion4550Regression:
    def test_strategy_report_and_final_list(self, exons) -> None:
        strategies = evaluate_all_strategies(45, 50, exons, max_additional_skips=3)
        assert len(strategies) == 9

        retained = []
        for s in strategies:
            pj = s.principal_junction
            if not pj and s.new_junctions:
                pj = s.new_junctions[0]
            status = "retained" if s.restores_frame else "excluded"
            retained.append((list(s.additional_skipped_exons), pj, s.restores_frame, status))

        assert any(s[2] for s in retained)
        assert any(not s[2] for s in retained)

        final = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        final_keys = [tuple(c.additional_skipped_exons) for c in final]
        restoring_keys = [
            tuple(s.additional_skipped_exons)
            for s in strategies
            if s.restores_frame
        ]
        for key in final_keys:
            assert key in restoring_keys
        assert all(c.restores_frame for c in final)
        assert final[0].additional_skipped_exons == [51]
