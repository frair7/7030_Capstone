"""
Horizontal part segments from EXON_TABLE — shared x-coordinates for
transcript fills and domain bubble alignment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class PartSegment:
    x0: float
    x1: float
    fill: str
    border: str
    exon: int
    label: str | None = None


def _norm_hex(color: str) -> str:
    return color.strip().upper()


def build_part_segments(
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
) -> list[PartSegment]:
    """Absolute x-ranges for each fractional part in every exon."""
    segments: list[PartSegment] = []
    for e in exon_table:
        n = e["n"]
        xa, xb = exon_x[n]
        w = xb - xa
        for f0, f1, fill, border in e["parts"]:
            segments.append(PartSegment(
                x0=xa + f0 * w,
                x1=xa + f1 * w,
                fill=fill,
                border=border,
                exon=n,
            ))
    return segments


def find_part_segments(
    label: str,
    domain_map: list[tuple[int, int, str, str]],
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
) -> list[PartSegment]:
    """
    Resolve absolute part segments for one domain label from EXON_TABLE.parts.

    Each domain is matched by colour within its DOMAIN_MAP exon span. When the
    previous domain has a different colour, the trailing partial segment in
    exon (start - 1) is included (e.g. H1 in exon 8, R1 in exon 10).
    """
    entry = next((d for d in domain_map if d[2] == label), None)
    if entry is None:
        return []

    start, end, _, color = entry
    target = _norm_hex(color)
    idx = next(i for i, d in enumerate(domain_map) if d[2] == label)
    prev_color = domain_map[idx - 1][3] if idx > 0 else None
    same_color_as_prev = (
        prev_color is not None and _norm_hex(prev_color) == target
    )
    extra_exon = None if same_color_as_prev or start <= 1 else start - 1

    segments: list[PartSegment] = []
    for e in exon_table:
        n = e["n"]
        in_main = start <= n <= end
        in_extra = extra_exon is not None and n == extra_exon
        if not in_main and not in_extra:
            continue
        xa, xb = exon_x[n]
        w = xb - xa
        for f0, f1, fill, border in e["parts"]:
            if _norm_hex(fill) != target:
                continue
            segments.append(PartSegment(
                x0=xa + f0 * w,
                x1=xa + f1 * w,
                fill=fill,
                border=border,
                exon=n,
                label=label,
            ))
    return segments


def part_rects_for_exon(
    exon_n: int,
    exon_x: dict[int, tuple[float, float]],
    parts: list[tuple[float, float, str, str]],
) -> list[tuple[float, float, str, str]]:
    """Return (x0, x1, fill, border) in absolute coordinates for one exon."""
    xa, xb = exon_x[exon_n]
    w = xb - xa
    return [(xa + f0 * w, xa + f1 * w, fill, border) for f0, f1, fill, border in parts]
