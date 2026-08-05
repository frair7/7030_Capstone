"""Tests for interlocking puzzle geometry."""

from __future__ import annotations

import re

import pytest

from src.puzzle_geometry import (
    build_exon_path,
    build_exon_stroke_path,
    calculate_bump_depth,
    left_edge_path,
    normalize_edge_shape,
    right_edge_path,
)

TEST_EXONS: list[tuple[EdgeShape, EdgeShape]] = [
    ("Flat", "Flat"),
    ("Flat", "Point"),
    ("Point", "Round"),
    ("Round", "Flat"),
    ("Flat", "Round"),
    ("Round", "Point"),
    ("Point", "Flat"),
]

X0, X1, Y_TOP, HEIGHT, DEPTH = 10.0, 40.0, 5.0, 26.0, 4.0


def _parse_points(path: str) -> list[tuple[float, float]]:
    pts: list[tuple[float, float]] = []
    for cmd in re.findall(r"[ML]\s*([-\d.]+)\s+([-\d.]+)", path):
        pts.append((float(cmd[0]), float(cmd[1])))
    return pts


def test_normalize_edge_shape() -> None:
    assert normalize_edge_shape("point") == "Point"
    assert normalize_edge_shape("ROUND") == "Round"
    assert normalize_edge_shape("n/a") == "N/A"
    assert normalize_edge_shape(None) == "N/A"
    with pytest.raises(ValueError, match="Unknown transcript edge shape"):
        normalize_edge_shape("unknown")


def test_calculate_bump_depth_bounds() -> None:
    assert calculate_bump_depth(100.0, 26.0) == pytest.approx(5.0)
    assert calculate_bump_depth(10.0, 26.0) == pytest.approx(3.8)
    assert calculate_bump_depth(3.0, 26.0) == pytest.approx(1.25)
    assert calculate_bump_depth(3.0, 4.0) == pytest.approx(1.25)


@pytest.mark.parametrize("shape", ["Flat", "Point", "Round", "N/A"])
def test_right_edge_path_commands(shape: str) -> None:
    path = " ".join(right_edge_path(X1, Y_TOP, Y_TOP + HEIGHT, shape, DEPTH))
    assert path.startswith("L ") or path.startswith("C ")
    assert "NaN" not in path
    if shape == "Flat":
        assert path == f"L {X1:.3f} {Y_TOP + HEIGHT:.3f}"
    elif shape == "Point":
        y_mid = Y_TOP + HEIGHT / 2.0
        assert f"L {X1 + DEPTH:.3f} {y_mid:.3f}" in path
        assert path.endswith(f"L {X1:.3f} {Y_TOP + HEIGHT:.3f}")
    elif shape == "Round":
        assert "C " in path
        assert f"{X1 + DEPTH:.3f}" in path


@pytest.mark.parametrize("shape", ["Flat", "Point", "Round", "N/A"])
def test_left_edge_path_commands(shape: str) -> None:
    path = " ".join(left_edge_path(X0, Y_TOP, Y_TOP + HEIGHT, shape, DEPTH))
    assert path.startswith("L ") or path.startswith("C ")
    assert "NaN" not in path
    if shape == "Flat":
        assert path == f"L {X0:.3f} {Y_TOP:.3f}"
    elif shape == "Point":
        y_mid = Y_TOP + HEIGHT / 2.0
        assert f"L {X0 + DEPTH:.3f} {y_mid:.3f}" in path
        assert path.endswith(f"L {X0:.3f} {Y_TOP:.3f}")
    elif shape == "Round":
        assert "C " in path
        assert f"{X0 + DEPTH:.3f}" in path


@pytest.mark.parametrize("five_prime,three_prime", TEST_EXONS)
def test_build_exon_path_closed_without_diagonal(
    five_prime: EdgeShape,
    three_prime: EdgeShape,
) -> None:
    path = build_exon_path(
        X0, X1, Y_TOP, HEIGHT,
        five_prime, three_prime, DEPTH,
    )
    assert path.startswith(f"M {X0:.3f} {Y_TOP:.3f}")
    assert path.endswith(f"L {X0:.3f} {Y_TOP:.3f} Z")
    assert "NaN" not in path

    pts = _parse_points(path)
    assert pts[0] == pytest.approx((X0, Y_TOP))
    assert pts[-1] == pytest.approx((X0, Y_TOP))

    y_bottom = Y_TOP + HEIGHT
    for (xa, ya), (xb, yb) in zip(pts, pts[1:]):
        if abs(xb - xa) > 0.01 and abs(yb - ya) > 0.01:
            near_right = abs(xa - X1) < 0.01 or abs(xb - X1) < 0.01
            near_left = abs(xa - X0) < 0.01 or abs(xb - X0) < 0.01
            near_protrusion = (
                abs(xa - (X0 + DEPTH)) < 0.01 or abs(xb - (X0 + DEPTH)) < 0.01
                or abs(xa - (X1 + DEPTH)) < 0.01 or abs(xb - (X1 + DEPTH)) < 0.01
            )
            assert near_right or near_left or near_protrusion, (
                f"diagonal across body: {(xa, ya)} -> {(xb, yb)}"
            )
            assert X0 - 0.01 <= xa <= X1 + DEPTH + 0.01
        assert Y_TOP - 0.01 <= ya <= y_bottom + 0.01


def test_build_exon_stroke_path_skips_junction_edges() -> None:
    """Interior exons must not draw flat vertical lines on shared junction edges."""
    path = build_exon_stroke_path(
        X0, X1, Y_TOP, HEIGHT,
        "Point", "Point", DEPTH,
        stroke_left=False,
        stroke_right=False,
    )
    assert f"L {X1:.3f} {Y_TOP + HEIGHT:.3f}" not in path
    assert f"L {X0:.3f} {Y_TOP:.3f}" not in path
    assert f"M {X1:.3f} {Y_TOP + HEIGHT:.3f}" in path
    assert f"M {X0:.3f} {Y_TOP:.3f}" in path
