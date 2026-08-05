"""
Load and validate DMD reference tables from local CSV files.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from src.config import REFERENCE

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = REPO_ROOT / "reference_tables"
DATA_DIR = Path(__file__).resolve().parents[1] / "data"

EXONS_CSV = REFERENCE_DIR / "dmd_exons_grch38.csv"
MAP_STYLES_CSV = REFERENCE_DIR / "dmd_exon_map_styles.csv"
PROTEIN_DOMAINS_CSV = REFERENCE_DIR / "dmd_protein_domains_dp427m.csv"
DOMAINS_CSV = DATA_DIR / "dmd_domains.csv"


def load_exon_table(path: Path | None = None) -> list[dict[str, Any]]:
    """Load the DMD exon reference CSV as a list of row dicts."""
    csv_path = path or EXONS_CSV
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Exon reference file not found: {csv_path}. "
            "Run: python scripts/fetch_reference_data.py"
        )
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_domain_table(path: Path | None = None) -> list[dict[str, Any]]:
    """Load the broad dystrophin domain summary CSV (used for map visualization)."""
    csv_path = path or DOMAINS_CSV
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Domain reference file not found: {csv_path}. "
            "Run: python scripts/fetch_reference_data.py"
        )
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_map_styles_table(path: Path | None = None) -> list[dict[str, Any]]:
    """Load per-exon map styling (colors, segment counts, edge directions)."""
    csv_path = path or MAP_STYLES_CSV
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Exon map styles file not found: {csv_path}. "
            "See reference_tables/README.md."
        )
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_protein_domains_table(path: Path | None = None) -> list[dict[str, Any]]:
    """Load detailed Dp427m protein domain and feature annotations."""
    csv_path = path or PROTEIN_DOMAINS_CSV
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Protein domains file not found: {csv_path}. "
            "See reference_tables/README.md."
        )
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def exon_by_number(rows: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    """Index exon rows by exon_number."""
    return {int(row["exon_number"]): row for row in rows}


def genomic_span(rows: list[dict[str, Any]]) -> tuple[int, int]:
    """Return (min_start, max_end) genomic span across all exons."""
    starts = [int(r["genomic_start_grch38"]) for r in rows]
    ends = [int(r["genomic_end_grch38"]) for r in rows]
    return min(starts), max(ends)
