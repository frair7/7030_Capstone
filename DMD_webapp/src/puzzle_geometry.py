"""
Diagram.net-style interlocking exon puzzle geometry.

Adjacent exons share a nominal junction x; the left exon's 3′ edge protrudes
right and the right exon's 5′ edge follows the same boundary. Paths always
return explicitly to the upper-left corner before closing.
"""

from __future__ import annotations

from typing import Literal

EdgeShape = Literal["Flat", "Point", "Round", "N/A"]


def normalize_edge_shape(value: str | None) -> EdgeShape:
    if value is None:
        return "N/A"

    cleaned = str(value).strip().lower()

    mapping: dict[str, EdgeShape] = {
        "flat": "Flat",
        "point": "Point",
        "round": "Round",
        "n/a": "N/A",
        "na": "N/A",
        "none": "N/A",
    }

    if cleaned not in mapping:
        raise ValueError(f"Unknown transcript edge shape: {value!r}")

    return mapping[cleaned]


def calculate_bump_depth(exon_width: float, exon_height: float) -> float:
    """Keep the edge shape visible without consuming most of a narrow exon."""
    return max(
        1.25,
        min(
            exon_height * 0.30,
            exon_width * 0.38,
            5.0,
        ),
    )


def right_edge_path(
    x: float,
    y_top: float,
    y_bottom: float,
    shape: EdgeShape | str,
    depth: float,
) -> str:
    """Draw the 3′ edge from top to bottom. Point and Round protrude toward +x."""
    edge = normalize_edge_shape(shape if isinstance(shape, str) else None)
    y_mid = (y_top + y_bottom) / 2.0

    if edge in {"Flat", "N/A"}:
        return [f"L {x:.3f} {y_bottom:.3f}"]

    if edge == "Point":
        return [
            f"L {x + depth:.3f} {y_mid:.3f}",
            f"L {x:.3f} {y_bottom:.3f}",
        ]

    if edge == "Round":
        return [
            (
                f"C {x + depth:.3f} {y_top:.3f}, "
                f"{x + depth:.3f} {y_bottom:.3f}, "
                f"{x:.3f} {y_bottom:.3f}"
            ),
        ]

    raise ValueError(f"Unsupported right edge shape: {edge}")


def left_edge_path(
    x: float,
    y_top: float,
    y_bottom: float,
    shape: EdgeShape | str,
    depth: float,
) -> str:
    """Draw the 5′ edge from bottom to top along the shared junction boundary."""
    edge = normalize_edge_shape(shape if isinstance(shape, str) else None)
    y_mid = (y_top + y_bottom) / 2.0

    if edge in {"Flat", "N/A"}:
        return [f"L {x:.3f} {y_top:.3f}"]

    if edge == "Point":
        return [
            f"L {x + depth:.3f} {y_mid:.3f}",
            f"L {x:.3f} {y_top:.3f}",
        ]

    if edge == "Round":
        return [
            (
                f"C {x + depth:.3f} {y_bottom:.3f}, "
                f"{x + depth:.3f} {y_top:.3f}, "
                f"{x:.3f} {y_top:.3f}"
            ),
        ]

    raise ValueError(f"Unsupported left edge shape: {edge}")


def build_exon_path(
    x0: float,
    x1: float,
    y_top: float,
    exon_height: float,
    five_prime_shape: EdgeShape | str,
    three_prime_shape: EdgeShape | str,
    bump_depth: float,
) -> str:
    """Build one closed exon puzzle-piece outline."""
    y_bottom = y_top + exon_height

    five_shape = normalize_edge_shape(
        five_prime_shape if isinstance(five_prime_shape, str) else None
    )
    three_shape = normalize_edge_shape(
        three_prime_shape if isinstance(three_prime_shape, str) else None
    )

    commands: list[str] = [
        f"M {x0:.3f} {y_top:.3f}",
        f"L {x1:.3f} {y_top:.3f}",
    ]

    commands.extend(
        right_edge_path(
            x=x1,
            y_top=y_top,
            y_bottom=y_bottom,
            shape=three_shape,
            depth=bump_depth,
        )
    )

    commands.append(f"L {x0:.3f} {y_bottom:.3f}")

    commands.extend(
        left_edge_path(
            x=x0,
            y_top=y_top,
            y_bottom=y_bottom,
            shape=five_shape,
            depth=bump_depth,
        )
    )

    commands.append(f"L {x0:.3f} {y_top:.3f}")
    commands.append("Z")

    return " ".join(commands)


def _right_edge_points(
    x: float,
    y_top: float,
    y_bottom: float,
    shape: EdgeShape | str,
    depth: float,
    *,
    steps: int = 12,
) -> list[tuple[float, float]]:
    """Vertex samples for 3′ edge (top → bottom), excluding start point."""
    edge = normalize_edge_shape(shape if isinstance(shape, str) else None)
    if edge in {"Flat", "N/A"} or depth <= 0:
        return [(x, y_bottom)]

    y_mid = (y_top + y_bottom) / 2.0
    if edge == "Point":
        return [(x + depth, y_mid), (x, y_bottom)]

    if edge == "Round":
        return _sample_cubic(
            (x, y_top),
            (x + depth, y_top),
            (x + depth, y_bottom),
            (x, y_bottom),
            steps,
        )

    raise ValueError(f"Unsupported right edge shape: {edge}")


def _left_edge_points(
    x: float,
    y_top: float,
    y_bottom: float,
    shape: EdgeShape | str,
    depth: float,
    *,
    steps: int = 12,
) -> list[tuple[float, float]]:
    """Vertex samples for 5′ edge (bottom → top), excluding start point."""
    edge = normalize_edge_shape(shape if isinstance(shape, str) else None)
    if edge in {"Flat", "N/A"} or depth <= 0:
        return [(x, y_top)]

    y_mid = (y_top + y_bottom) / 2.0
    if edge == "Point":
        return [(x + depth, y_mid), (x, y_top)]

    if edge == "Round":
        return _sample_cubic(
            (x, y_bottom),
            (x + depth, y_bottom),
            (x + depth, y_top),
            (x, y_top),
            steps,
        )

    raise ValueError(f"Unsupported left edge shape: {edge}")


def _sample_cubic(
    p0: tuple[float, float],
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
    steps: int,
) -> list[tuple[float, float]]:
    pts: list[tuple[float, float]] = []
    for i in range(1, steps + 1):
        t = i / steps
        u = 1.0 - t
        x = (
            u ** 3 * p0[0]
            + 3 * u ** 2 * t * p1[0]
            + 3 * u * t ** 2 * p2[0]
            + t ** 3 * p3[0]
        )
        y = (
            u ** 3 * p0[1]
            + 3 * u ** 2 * t * p1[1]
            + 3 * u * t ** 2 * p2[1]
            + t ** 3 * p3[1]
        )
        pts.append((x, y))
    return pts


def exon_polygon_verts(
    x: float,
    width: float,
    y: float,
    height: float,
    left_shape: str,
    right_shape: str,
    *,
    bump_depth: float = 0.0,
) -> list[tuple[float, float]]:
    """Closed polygon vertices for matplotlib rendering."""
    x0, x1 = x, x + width
    y_top, y_bottom = y, y + height
    verts: list[tuple[float, float]] = [(x0, y_top), (x1, y_top)]
    verts.extend(_right_edge_points(x1, y_top, y_bottom, right_shape, bump_depth))
    verts.append((x0, y_bottom))
    verts.extend(_left_edge_points(x0, y_top, y_bottom, left_shape, bump_depth))
    verts.append((x0, y_top))
    return verts


def build_junction_path(
    x: float,
    y_top: float,
    exon_height: float,
    three_prime_shape: EdgeShape | str,
    five_prime_shape: EdgeShape | str,
    depth: float,
) -> str:
    """Single-stroke path for one internal exon–exon junction."""
    y_bottom = y_top + exon_height
    return " ".join([
        f"M {x:.3f} {y_top:.3f}",
        *right_edge_path(x, y_top, y_bottom, three_prime_shape, depth),
        *left_edge_path(x, y_top, y_bottom, five_prime_shape, depth),
    ])


def build_exon_stroke_path(
    x0: float,
    x1: float,
    y_top: float,
    exon_height: float,
    five_prime_shape: EdgeShape | str,
    three_prime_shape: EdgeShape | str,
    bump_depth: float,
    *,
    stroke_left: bool,
    stroke_right: bool,
) -> str:
    """Open outline path for one exon, omitting edges drawn by shared junctions."""
    y_bottom = y_top + exon_height
    parts: list[str] = [f"M {x0:.3f} {y_top:.3f}", f"L {x1:.3f} {y_top:.3f}"]

    if stroke_right:
        parts.extend(
            right_edge_path(x1, y_top, y_bottom, three_prime_shape, bump_depth)
        )
    else:
        # Junction path owns this edge — move without drawing a flat vertical line.
        parts.append(f"M {x1:.3f} {y_bottom:.3f}")

    parts.append(f"L {x0:.3f} {y_bottom:.3f}")

    if stroke_left:
        parts.extend(
            left_edge_path(x0, y_top, y_bottom, five_prime_shape, bump_depth)
        )
    else:
        parts.append(f"M {x0:.3f} {y_top:.3f}")

    return " ".join(parts)


def domain_bubble_svg_d(
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    radius: float = 5.0,
) -> str:
    """Rounded rectangle path for one domain bubble."""
    r = min(radius, width / 4, height / 2)
    x2, y2 = x + width, y + height
    return (
        f"M {x + r:.2f} {y:.2f} "
        f"L {x2 - r:.2f} {y:.2f} Q {x2:.2f} {y:.2f} {x2:.2f} {y + r:.2f} "
        f"L {x2:.2f} {y2 - r:.2f} Q {x2:.2f} {y2:.2f} {x2 - r:.2f} {y2:.2f} "
        f"L {x + r:.2f} {y2:.2f} Q {x:.2f} {y2:.2f} {x:.2f} {y2 - r:.2f} "
        f"L {x:.2f} {y + r:.2f} Q {x:.2f} {y:.2f} {x + r:.2f} {y:.2f} Z"
    )
