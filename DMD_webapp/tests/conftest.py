"""Shared pytest fixtures."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import build_exon_records


@pytest.fixture(scope="module")
def exon_coding_lengths() -> dict[int, int]:
    """Canonical Dp427m coding-length map from the project reference table."""
    exons = build_exon_records()
    return {exon.exon_number: exon.coding_length_bp for exon in exons}
