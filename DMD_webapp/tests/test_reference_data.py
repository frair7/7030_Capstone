"""Tests for reference data loading edge cases."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from src.reference_data import load_exon_table
from scripts.validate_reference_data import validate_exons


def test_missing_reference_file(tmp_path: Path) -> None:
    missing = tmp_path / "no_exons.csv"
    with pytest.raises(FileNotFoundError):
        load_exon_table(missing)


def test_inconsistent_reference_data(tmp_path: Path, exon_csv: Path) -> None:
    """Validation should fail when cumulative CDS is corrupted."""
    rows = load_exon_table(exon_csv)
    rows[10]["cumulative_cds_end"] = "99999"
    bad_path = tmp_path / "bad_exons.csv"
    with bad_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    errors = validate_exons(load_exon_table(bad_path))
    assert len(errors) > 0


@pytest.fixture
def exon_csv() -> Path:
    from src.reference_data import EXONS_CSV

    return EXONS_CSV
