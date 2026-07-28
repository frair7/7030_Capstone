"""Tests for coordinate_mapper module."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import (
    build_exon_records,
    map_genomic_position,
    map_parsed_variant,
)
from src.models import InputMode
from src.variant_parser import parse_variant


@pytest.fixture(scope="module")
def exons():
    return build_exon_records()


class TestGenomicMapping:
    def test_position_inside_exon(self, exons) -> None:
        # Exon 79 start from reference table
        exon79 = next(e for e in exons if e.exon_number == 79)
        pos = exon79.genomic_start + 10
        mapping = map_genomic_position(pos, exons)
        assert mapping.location_type == "exon"
        assert mapping.exon_number == 79

    def test_position_in_intron(self, exons) -> None:
        # Intron between exon 1 and 2 (minus strand)
        e1 = next(e for e in exons if e.exon_number == 1)
        e2 = next(e for e in exons if e.exon_number == 2)
        intron_pos = e2.genomic_end + 100
        assert intron_pos < e1.genomic_start
        mapping = map_genomic_position(intron_pos, exons)
        assert mapping.location_type == "intron"
        assert mapping.upstream_exon == 1
        assert mapping.downstream_exon == 2

    def test_outside_gene_downstream(self, exons) -> None:
        mapping = map_genomic_position(30_000_000, exons)
        assert mapping.location_type == "downstream"


class TestTranscriptOrder:
    def test_exon1_higher_genomic_than_exon79(self, exons) -> None:
        e1 = next(e for e in exons if e.exon_number == 1)
        e79 = next(e for e in exons if e.exon_number == 79)
        assert e1.genomic_start > e79.genomic_end


class TestParsedVariantMapping:
    def test_deletion_maps_to_exons(self, exons) -> None:
        variant = parse_variant("del45-50", InputMode.EXON_DELETION)
        mapping = map_parsed_variant(variant, exons)
        assert mapping.location_type == "exon"
        assert mapping.upstream_exon == 45
