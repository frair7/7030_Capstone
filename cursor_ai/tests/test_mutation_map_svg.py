"""Tests for gradient fills and interactive SVG map."""

from __future__ import annotations

from src.exon_gradients import build_horizontal_gradient, gradient_stops_for_svg
from src.mutation_map_svg import build_interactive_map_html


def test_gradient_single_part() -> None:
    parts = [(0.0, 1.0, "#FFCCCC", "#EB6F6F")]
    arr = build_horizontal_gradient(parts, n_samples=32)
    assert arr.shape == (32, 4)


def test_gradient_dual_domain_blend() -> None:
    parts = [
        (0.0, 0.39, "#FFCCCC", "#EB6F6F"),
        (0.39, 1.0, "#AEDCCA", "#4BAF89"),
    ]
    stops = gradient_stops_for_svg(parts)
    assert len(stops) >= 3
    arr = build_horizontal_gradient(parts, n_samples=64)
    mid = arr[32, :3]
    assert 0.4 < mid[0] < 1.0  # blended, not pure green or pure pink


def test_interactive_map_html_builds() -> None:
    rows = [
        {
            "id": "P-6",
            "group": "B",
            "start_region": "e3",
            "stop_region": "e29",
            "phenotype": "DMD",
            "expected_protein_size_kda": "273.4",
            "pct_dys_wb": "49.16",
        }
    ]
    html = build_interactive_map_html(rows)
    assert "data-tip" in html
    assert 'fill="#FFCCCC"' in html
    assert 'stroke="#EB6F6F"' in html
    assert "clip-e8" in html
    assert html.count('class="hinge-guide"') == 8
    assert html.rindex('class="hinge-guide"') > html.rindex("</text>")
