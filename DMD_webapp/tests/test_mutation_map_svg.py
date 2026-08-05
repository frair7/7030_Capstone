"""Tests for gradient fills and interactive SVG map."""

from __future__ import annotations

from src.dp427m_exon_data import get_exon_table
from src.exon_gradients import build_horizontal_gradient, gradient_stops_for_svg, svg_linear_gradient_def
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
    stops = gradient_stops_for_svg(
        parts, draw_x0=0.0, draw_x1=10.0, exon_x0=0.0, exon_x1=8.0,
    )
    assert len(stops) >= 3
    arr = build_horizontal_gradient(
        parts, n_samples=64, draw_x0=0.0, draw_x1=10.0, exon_x0=0.0, exon_x1=8.0,
    )
    mid = arr[32, :3]
    assert 0.4 < mid[0] < 1.0


def test_svg_gradient_uses_user_space() -> None:
    parts = [
        (0.0, 0.633, "#FFFFCC", "#F0C266"),
        (0.633, 1.0, "#AEDCCA", "#4BAF89"),
    ]
    grad = svg_linear_gradient_def(
        "grad-test", 100.0, 122.0, 52.0, parts, exon_x0=100.0, exon_x1=117.0,
    )
    assert 'gradientUnits="userSpaceOnUse"' in grad
    assert 'x1="100.00"' in grad
    assert 'x2="122.00"' in grad
    # boundary at 63.3% of 17px exon ≈ 48.9% along 22px draw span
    assert 'offset="48.' in grad or 'offset="49.' in grad


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
    assert 'fill="#FFCCCC"' in html or 'stop-color="#FFCCCC"' in html
    assert "clip-e8" in html
    assert html.count('class="hinge-guide"') == 8
    assert html.rindex('class="hinge-guide"') > html.rindex("</text>")
    assert 'class="table-header-vertical"' in html
    assert 'transform="rotate(-90' in html
    assert 'gradientUnits="userSpaceOnUse"' in html
    assert 'id="grad-e17"' in html
    assert 'id="grad-e8"' in html
    import re
    transcript_internal = re.findall(
        r'<line x1="[^"]+" y1="52\.00"[^>]+opacity="0\.8"',
        html,
    )
    assert transcript_internal == []
    # No vertical track separator at TRACK_LEFT
    assert 'y1="42.0" x2="205.0" y2="106.0"' not in html


def test_multicolor_exons_use_single_outline() -> None:
    html = build_interactive_map_html([])
    exon_table = get_exon_table()
    for n in (17, 50, 61, 64):
        e = next(x for x in exon_table if x["n"] == n)
        assert len(e["parts"]) > 1
        assert f'id="grad-e{n}"' in html
        assert f"clip-e{n}-s0" not in html


def test_single_color_exons_unchanged() -> None:
    html = build_interactive_map_html([])
    assert 'id="grad-e6"' not in html
    assert 'id="grad-e43"' not in html
    assert "clip-e6" in html
    assert "clip-e43" in html
