"""Tests for frame_analysis module."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import build_exon_records
from src.frame_analysis import (
    assess_deletion_frame,
    assess_duplication_frame,
    coding_bp_for_exon_range,
)
from src.models import FrameStatus, InputMode
from src.variant_parser import parse_variant
from src.frame_analysis import assess_variant_frame


@pytest.fixture(scope="module")
def exons():
    return build_exon_records()


class TestDeletionFrame:
    def test_del45_50_out_of_frame(self, exons) -> None:
        """
        del45-50 removes 871 coding bp (871 mod 3 = 1).

        Verified against Ensembl ENST00000357033.9 exon table
        (DMD_webapp/data/dmd_exons_grch38.csv).
        """
        result = assess_deletion_frame(45, 50, exons)
        assert result.status == FrameStatus.OUT_OF_FRAME
        assert result.coding_bases_removed == 871
        assert result.upstream_exon == 44
        assert result.downstream_exon == 51

    def test_single_exon_del_phase_math(self, exons) -> None:
        """Exon 5 deletion: 93 bp removed (divisible by 3) → in frame."""
        removed = coding_bp_for_exon_range(exons, 5, 5)
        result = assess_deletion_frame(5, 5, exons)
        assert removed == 93
        assert removed % 3 == 0
        assert result.status == FrameStatus.IN_FRAME

    def test_in_frame_triplet_deletion(self, exons) -> None:
        """Find a 3-exon range whose total coding length is divisible by 3."""
        for first in range(1, 77):
            for last in range(first, min(first + 3, 79)):
                bp = coding_bp_for_exon_range(exons, first, last)
                if bp % 3 == 0 and first > 1 and last < 79:
                    result = assess_deletion_frame(first, last, exons)
                    assert result.status == FrameStatus.IN_FRAME
                    return
        pytest.fail("No in-frame deletion test case found in reference data")


class TestDuplicationFrame:
    def test_dup2(self, exons) -> None:
        """Duplication of exon 2 adds 62 bp (62 mod 3 = 2) → out of frame."""
        result = assess_duplication_frame(2, 2, exons)
        assert result.status == FrameStatus.OUT_OF_FRAME
        assert result.coding_bases_added == 62

    def test_dup3_4_in_frame(self, exons) -> None:
        """
        Duplication of exons 3-4 adds 171 bp (171 mod 3 = 0).

        Parent-repo DMD_Mutations_clean lists dup e3-e4 as IF.
        """
        added = coding_bp_for_exon_range(exons, 3, 4)
        assert added % 3 == 0
        result = assess_duplication_frame(3, 4, exons)
        assert result.status == FrameStatus.IN_FRAME


class TestVariantDispatch:
    def test_snv_cannot_determine(self, exons) -> None:
        variant = parse_variant("c.5287C>T", InputMode.HGVS_CODING)
        result = assess_variant_frame(variant, exons)
        assert result.status == FrameStatus.CANNOT_DETERMINE
