"""Tests for mutation catalog and autofill."""

from __future__ import annotations

from src.mutation_autofill import autofill_mutation_fields
from src.mutation_catalog import filter_catalog, import_legacy_catalog
from src.region_parser import parse_quick_entry, parse_region_text


def test_parse_region_range() -> None:
    r = parse_region_text("e45-e55")
    assert r.start_exon == 45
    assert r.stop_exon == 55


def test_quick_entry() -> None:
    parsed = parse_quick_entry("Exonic Deletion e45-e55")
    assert parsed["mutation_class"] == "Exonic"
    assert parsed["mutation_subclass"] == "Deletion"
    assert parsed["start_region"] == "e45"


def test_autofill_deletion_frame() -> None:
    filled = autofill_mutation_fields("Exonic", "Deletion", "e45", "e50")
    assert filled["frame"] == "Out-of-frame"
    assert filled["expected_protein_size_kda"]


def test_legacy_import_has_ids() -> None:
    df = import_legacy_catalog()
    assert len(df) == 101
    assert df.iloc[0]["id"] == "S1"


def test_filter_catalog() -> None:
    df = import_legacy_catalog()
    out = filter_catalog(df, mutation_classes=["Exonic"])
    assert all(out["mutation_class"] == "Exonic")
