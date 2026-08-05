"""Tests for visualization module."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import build_exon_records
from src.visualization import (
    ViewMode,
    VisualizationState,
    build_exon_table,
    compute_exon_layouts,
    create_combined_figure,
    create_transcript_figure,
    visualization_state_from_variant,
)


@pytest.fixture(scope="module")
def exons():
    return build_exon_records()


class TestLayouts:
    def test_schematic_has_79_layouts(self, exons) -> None:
        layouts = compute_exon_layouts(exons, ViewMode.SCHEMATIC)
        assert len(layouts) == 79

    def test_schematic_equal_widths(self, exons) -> None:
        layouts = compute_exon_layouts(exons, ViewMode.SCHEMATIC)
        widths = [l.x_end - l.x_start for l in layouts]
        assert len(set(round(w, 4) for w in widths)) == 1

    def test_genomic_uses_transcript_coords(self, exons) -> None:
        layouts = compute_exon_layouts(exons, ViewMode.GENOMIC)
        e1 = layouts[0]
        assert e1.x_start == exons[0].transcript_exon_start


class TestFigures:
    def test_transcript_figure_builds(self, exons) -> None:
        fig = create_transcript_figure(exons)
        assert len(fig.data) >= 1
        assert fig.layout.height == 220

    def test_combined_figure_builds(self, exons) -> None:
        state = visualization_state_from_variant(45, 50, variant_type="deletion")
        fig = create_combined_figure(exons, state, ViewMode.SCHEMATIC)
        assert fig.layout.height == 480


class TestHighlightState:
    def test_deletion_highlight(self) -> None:
        state = visualization_state_from_variant(45, 50, variant_type="deletion")
        assert 45 in state.deleted_exons
        assert 50 in state.deleted_exons
        assert 44 not in state.deleted_exons


class TestAccessibleTable:
    def test_table_has_79_rows(self, exons) -> None:
        df = build_exon_table(exons)
        assert len(df) == 79
        assert "visual_state" in df.columns
