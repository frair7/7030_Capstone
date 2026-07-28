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


def build_horizontal_gradient(
    parts: list[tuple[float, float, str, str]],
    *,
    n_samples: int = 256,
    fade_frac: float = 0.12,
) -> np.ndarray:
    """
    Build an (n_samples, 4) RGBA array with smooth blends at part boundaries.

    parts: list of (f0, f1, fill_hex, border_hex) fractions along exon width.
    """
    rgba = np.ones((n_samples, 4), dtype=float)
    if not parts:
        return rgba

    xs = np.linspace(0.0, 1.0, n_samples)
    for i, x in enumerate(xs):
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


def gradient_stops_for_svg(
    parts: list[tuple[float, float, str, str]],
    *,
    fade_frac: float = 0.08,
) -> list[tuple[float, str]]:
    """SVG gradient stops centred on domain breakpoints within the exon."""
    if len(parts) <= 1:
        fill = parts[0][2] if parts else "#cccccc"
        return [(0, fill), (100, fill)]

    stops: list[tuple[float, str]] = [(0, parts[0][2])]
    for idx in range(len(parts) - 1):
        _f0a, f1a, fill_a, _ = parts[idx]
        f0b, _f1b, fill_b, _ = parts[idx + 1]
        boundary = f1a * 100  # == f0b * 100 when data is consistent
        half = fade_frac * 50
        stops.append((max(0, boundary - half), fill_a))
        mid = _rgb_to_hex(_lerp(_hex_to_rgb(fill_a), _hex_to_rgb(fill_b), 0.5))
        stops.append((boundary, mid))
        stops.append((min(100, boundary + half), fill_b))
    stops.append((100, parts[-1][2]))
    return sorted(stops, key=lambda s: s[0])
