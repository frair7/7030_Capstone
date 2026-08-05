"""Tests for mutation visualization helpers."""

from __future__ import annotations

from src.mutation_viz import exon_sort_key_from_row, sort_mutations_by_exon_range


def _row(label: str, start: str, stop: str = "") -> dict[str, str]:
    return {
        "id": label,
        "start_region": start,
        "stop_region": stop or start,
        "mutation_subclass": "Deletion",
    }


def test_sort_mutations_by_exon_range_example() -> None:
    mutations = [
        _row("Del3", "e3"),
        _row("Del3-4", "e3", "e4"),
        _row("Del5", "e5"),
        _row("Del20-50", "e20", "e50"),
        _row("Del16", "e16"),
    ]
    sorted_ids = [m["id"] for m in sort_mutations_by_exon_range(mutations)]
    assert sorted_ids == ["Del3", "Del3-4", "Del5", "Del16", "Del20-50"]


def test_exon_sort_key_single_exon() -> None:
    assert exon_sort_key_from_row(_row("Del3", "e3")) == (3, 3)


def test_exon_sort_key_range() -> None:
    assert exon_sort_key_from_row(_row("Del3-4", "e3", "e4")) == (3, 4)
