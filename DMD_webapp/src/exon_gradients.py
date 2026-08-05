"""Color gradient helpers for dual-domain exons."""

from __future__ import annotations

import numpy as np


def _hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    h = hex_color.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def _rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    return "#{:02x}{:02x}{:02x}".format(
        int(max(0, min(1, rgb[0])) * 255),
        int(max(0, min(1, rgb[1])) * 255),
        int(max(0, min(1, rgb[2])) * 255),
    )


def _lerp(c0: tuple[float, float, float], c1: tuple[float, float, float], t: float) -> tuple[float, float, float]:
    return tuple(c0[i] + (c1[i] - c0[i]) * t for i in range(3))


def blend_fraction(*, exon_width: float, fade_frac: float = 0.12) -> float:
    """Fade width as a fraction of exon width (min 2 px, max 8 px)."""
    if exon_width <= 0:
        return fade_frac
    blend_px = min(max(exon_width * fade_frac, 2.0), 8.0)
    return blend_px / exon_width


def _exon_fraction_to_offset(
    fraction: float,
    *,
    draw_x0: float,
    draw_x1: float,
    exon_x0: float,
    exon_x1: float,
) -> float:
    """Map an exon-local fraction [0,1] to a gradient stop offset % along draw span."""
    span = draw_x1 - draw_x0
    if span <= 0:
        return 0.0
    ax = exon_x0 + fraction * (exon_x1 - exon_x0)
    return (ax - draw_x0) / span * 100.0


def gradient_stops_for_svg(
    parts: list[tuple[float, float, str, str]],
    *,
    draw_x0: float,
    draw_x1: float,
    exon_x0: float,
    exon_x1: float,
    fade_frac: float = 0.12,
) -> list[tuple[float, str]]:
    """SVG gradient stops at absolute positions mapped onto the draw span."""
    if len(parts) <= 1:
        fill = parts[0][2] if parts else "#cccccc"
        return [(0.0, fill), (100.0, fill)]

    exon_width = exon_x1 - exon_x0
    fade_frac = blend_fraction(exon_width=exon_width, fade_frac=fade_frac)
    half = fade_frac / 2.0

    def off(f: float) -> float:
        return _exon_fraction_to_offset(
            f, draw_x0=draw_x0, draw_x1=draw_x1, exon_x0=exon_x0, exon_x1=exon_x1,
        )

    stops: list[tuple[float, str]] = [(0.0, parts[0][2])]
    for idx in range(len(parts) - 1):
        _f0a, f1a, fill_a, _ = parts[idx]
        _f0b, _f1b, fill_b, _ = parts[idx + 1]
        stops.append((off(max(0.0, f1a - half)), fill_a))
        mid = _rgb_to_hex(_lerp(_hex_to_rgb(fill_a), _hex_to_rgb(fill_b), 0.5))
        stops.append((off(f1a), mid))
        stops.append((off(min(1.0, f1a + half)), fill_b))
    stops.append((100.0, parts[-1][2]))

    deduped: list[tuple[float, str]] = []
    for pct, color in sorted(stops, key=lambda s: s[0]):
        if deduped and abs(deduped[-1][0] - pct) < 1e-6:
            deduped[-1] = (pct, color)
        else:
            deduped.append((pct, color))
    return deduped


def svg_linear_gradient_def(
    gradient_id: str,
    draw_x0: float,
    draw_x1: float,
    y: float,
    parts: list[tuple[float, float, str, str]],
    *,
    exon_x0: float,
    exon_x1: float,
    fade_frac: float = 0.12,
) -> str:
    """userSpaceOnUse horizontal gradient for one exon's internal colour segments."""
    stops = gradient_stops_for_svg(
        parts,
        draw_x0=draw_x0,
        draw_x1=draw_x1,
        exon_x0=exon_x0,
        exon_x1=exon_x1,
        fade_frac=fade_frac,
    )
    stop_elems = "".join(
        f'<stop offset="{pct:.3f}%" stop-color="{color}"/>'
        for pct, color in stops
    )
    return (
        f'<linearGradient id="{gradient_id}" gradientUnits="userSpaceOnUse" '
        f'x1="{draw_x0:.2f}" y1="{y:.2f}" x2="{draw_x1:.2f}" y2="{y:.2f}">'
        f"{stop_elems}</linearGradient>"
    )


def build_horizontal_gradient(
    parts: list[tuple[float, float, str, str]],
    *,
    n_samples: int = 256,
    fade_frac: float = 0.12,
    draw_x0: float = 0.0,
    draw_x1: float = 1.0,
    exon_x0: float = 0.0,
    exon_x1: float = 1.0,
) -> np.ndarray:
    """
    Build an (n_samples, 4) RGBA array with smooth blends at part boundaries.

    Samples are spaced along [draw_x0, draw_x1] but colours follow exon-local
    fractions defined on [exon_x0, exon_x1].
    """
    rgba = np.ones((n_samples, 4), dtype=float)
    if not parts:
        return rgba

    exon_w = exon_x1 - exon_x0
    fade_frac = blend_fraction(exon_width=exon_w, fade_frac=fade_frac) if exon_w > 0 else fade_frac
    abs_xs = np.linspace(draw_x0, draw_x1, n_samples)

    for i, ax in enumerate(abs_xs):
        x = (ax - exon_x0) / exon_w if exon_w > 0 else 0.0
        c = _hex_to_rgb(parts[-1][2])
        for f0, f1, fill, _ in parts:
            if f0 <= x <= f1 + 1e-9:
                c = _hex_to_rgb(fill)
                break
        if len(parts) > 1 and fade_frac > 0:
            for j in range(len(parts) - 1):
                b = parts[j][1]
                half = fade_frac / 2
                if abs(x - b) <= half:
                    c = _lerp(
                        _hex_to_rgb(parts[j][2]),
                        _hex_to_rgb(parts[j + 1][2]),
                        min(1.0, max(0.0, (x - (b - half)) / max(2 * half, 1e-6))),
                    )
                    break
        rgba[i, :3] = c
    return rgba
