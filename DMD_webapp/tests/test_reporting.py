"""Tests for reporting module."""

from __future__ import annotations

from src.frame_analysis import assess_deletion_frame
from src.models import InputMode
from src.reporting import build_text_report, skip_candidates_dataframe
from src.variant_parser import parse_variant
from src.exon_skipping import find_skip_candidates
from src.coordinate_mapper import build_exon_records, map_parsed_variant


def test_text_report_contains_disclaimer() -> None:
    exons = build_exon_records()
    variant = parse_variant("del45-50", InputMode.EXON_DELETION)
    frame = assess_deletion_frame(45, 50, exons)
    mapping = map_parsed_variant(variant, exons)
    skips = find_skip_candidates(45, 50, exons)
    report = build_text_report(variant, frame, mapping.message, skips)
    assert "NOT FOR CLINICAL" in report
    assert "del45-50" in report or "45" in report


def test_skip_dataframe_columns() -> None:
    exons = build_exon_records()
    skips = find_skip_candidates(45, 50, exons)
    df = skip_candidates_dataframe(skips)
    if not skips:
        return
    assert "additional_skipped_exons" in df.columns
    assert len(df) >= 1
