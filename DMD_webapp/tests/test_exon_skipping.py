"""Tests for exon_skipping module."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import build_exon_records
from src.exon_skipping import find_skip_candidates
from src.frame_analysis import assess_deletion_frame
from src.models import FrameStatus


@pytest.fixture(scope="module")
def exons():
    return build_exon_records()


class TestSkipSearch:
    def test_del45_50_finds_candidates(self, exons) -> None:
        """del45-50 is OOF; skipping exon 51 (233 bp, 233 mod 3 = 2) restores frame."""
        base = assess_deletion_frame(45, 50, exons)
        assert base.status == FrameStatus.OUT_OF_FRAME

        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        assert len(candidates) >= 1
        best = candidates[0]
        assert best.restores_frame
        assert 51 in best.additional_skipped_exons or len(best.additional_skipped_exons) >= 1

    def test_in_frame_deletion_no_candidates(self, exons) -> None:
        """In-frame deletions should not need additional skipping."""
        # Find an in-frame single exon or small deletion
        for first in range(2, 78):
            for last in range(first, first + 2):
                base = assess_deletion_frame(first, last, exons)
                if base.status == FrameStatus.IN_FRAME:
                    candidates = find_skip_candidates(first, last, exons)
                    assert candidates == []
                    return
        pytest.fail("No in-frame deletion found for negative test")

    def test_ranking_prefers_fewer_skips(self, exons) -> None:
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3)
        if len(candidates) >= 2:
            assert candidates[0].rank_score <= candidates[1].rank_score

    def test_no_candidates_within_tight_limit(self, exons) -> None:
        """With max_additional_skips=0, no candidates can be found."""
        candidates = find_skip_candidates(45, 50, exons, max_additional_skips=0)
        assert candidates == []
