"""Tests for variant_parser module."""

from __future__ import annotations

import pytest

from src.models import InputMode, VariantType
from src.variant_parser import parse_variant


class TestExonDeletions:
    def test_single_exon_deletion(self) -> None:
        result = parse_variant("del45", InputMode.EXON_DELETION)
        assert result.is_valid
        assert result.first_exon == 45
        assert result.last_exon == 45
        assert result.variant_type == VariantType.DELETION

    def test_multi_exon_deletion(self) -> None:
        result = parse_variant("del45-50", InputMode.EXON_DELETION)
        assert result.is_valid
        assert result.first_exon == 45
        assert result.last_exon == 50

    def test_deletion_verbose(self) -> None:
        result = parse_variant("deletion of exons 45 through 50", InputMode.EXON_DELETION)
        assert result.is_valid
        assert result.first_exon == 45
        assert result.last_exon == 50


class TestDuplications:
    def test_single_duplication(self) -> None:
        result = parse_variant("dup2", InputMode.EXON_DUPLICATION)
        assert result.is_valid
        assert result.first_exon == 2
        assert result.last_exon == 2
        assert result.variant_type == VariantType.DUPLICATION

    def test_multi_duplication(self) -> None:
        result = parse_variant("dup3-7", InputMode.EXON_DUPLICATION)
        assert result.is_valid
        assert result.first_exon == 3
        assert result.last_exon == 7


class TestInvalidRanges:
    def test_exon_zero(self) -> None:
        result = parse_variant("del0", InputMode.EXON_DELETION)
        assert not result.is_valid

    def test_exon_over_79(self) -> None:
        result = parse_variant("del80", InputMode.EXON_DELETION)
        assert not result.is_valid

    def test_end_before_start(self) -> None:
        result = parse_variant("del50-45", InputMode.EXON_DELETION)
        assert not result.is_valid


class TestHgvs:
    def test_snv(self) -> None:
        result = parse_variant("c.5287C>T", InputMode.HGVS_CODING)
        assert result.is_valid
        assert result.cds_position == 5287
        assert result.variant_type == VariantType.SUBSTITUTION

    def test_splice(self) -> None:
        result = parse_variant("c.1812+1G>A", InputMode.HGVS_CODING)
        assert result.is_valid
        assert result.variant_type == VariantType.SPLICE

    def test_invalid_hgvs(self) -> None:
        result = parse_variant("not-valid-hgvs", InputMode.HGVS_CODING)
        assert not result.is_valid


class TestGenomic:
    def test_single_position(self) -> None:
        result = parse_variant("X:31140000", InputMode.GENOMIC)
        assert result.is_valid
        assert result.genomic_start == 31140000

    def test_interval(self) -> None:
        result = parse_variant("chrX:31140000-31150000", InputMode.GENOMIC)
        assert result.is_valid
        assert result.genomic_start == 31140000
        assert result.genomic_end == 31150000

    def test_grch37_rejected(self) -> None:
        result = parse_variant("X:31140000", InputMode.GENOMIC, assembly="GRCh37")
        assert not result.is_valid
