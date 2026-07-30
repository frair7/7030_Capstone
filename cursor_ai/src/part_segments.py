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


def _contiguous_parts_from_5prime(
    parts: list[tuple[float, float, str, str]],
    target: str,
) -> list[tuple[float, float, str, str]]:
    """Colour-matched parts contiguous from the 5′ edge."""
    kept: list[tuple[float, float, str, str]] = []
    for part in parts:
        if _norm_hex(part[2]) != target:
            break
        kept.append(part)
    return kept


def _contiguous_parts_from_3prime(
    parts: list[tuple[float, float, str, str]],
    target: str,
) -> list[tuple[float, float, str, str]]:
    """Colour-matched parts contiguous from the 3′ edge."""
    kept: list[tuple[float, float, str, str]] = []
    for part in reversed(parts):
        if _norm_hex(part[2]) != target:
            break
        kept.insert(0, part)
    return kept


def _interior_color_parts(
    parts: list[tuple[float, float, str, str]],
    target: str,
) -> list[tuple[float, float, str, str]]:
    """Matching parts sandwiched between other colours (e.g. H2 in exon 17)."""
    leading = _contiguous_parts_from_5prime(parts, target)
    trailing = _contiguous_parts_from_3prime(parts, target)
    lead_n = len(leading)
    trail_n = len(trailing)
    if lead_n + trail_n >= len(parts):
        return []
    return [
        p for i, p in enumerate(parts)
        if lead_n <= i < len(parts) - trail_n and _norm_hex(p[2]) == target
    ]


def _parts_for_extra_exon(
    parts: list[tuple[float, float, str, str]],
    target: str,
) -> list[tuple[float, float, str, str]]:
    """Parts for the exon immediately 5′ of a domain's main span."""
    if len(parts) <= 1:
        return [p for p in parts if _norm_hex(p[2]) == target]

    interior = _interior_color_parts(parts, target)
    if interior:
        return interior

    trailing = _contiguous_parts_from_3prime(parts, target)
    if trailing and len(trailing) < len(parts):
        return trailing

    leading = _contiguous_parts_from_5prime(parts, target)
    if leading and len(leading) < len(parts):
        return leading

    return [p for p in parts if _norm_hex(p[2]) == target]


def _parts_after_leading_other(
    parts: list[tuple[float, float, str, str]],
    target: str,
) -> list[tuple[float, float, str, str]]:
    """Parts matching target that follow a leading block of another colour (e.g. H1 in exon 8)."""
    seen_other = False
    kept: list[tuple[float, float, str, str]] = []
    for part in parts:
        if _norm_hex(part[2]) == target:
            if seen_other:
                kept.append(part)
        else:
            seen_other = True
    return kept


def _parts_for_domain_in_exon(
    parts: list[tuple[float, float, str, str]],
    target: str,
    *,
    is_start_exon: bool,
    is_end_exon: bool,
) -> list[tuple[float, float, str, str]]:
    """Select fractional parts belonging to one domain within a possibly split exon."""
    if len(parts) <= 1:
        return [p for p in parts if _norm_hex(p[2]) == target]

    if is_start_exon and is_end_exon:
        return [p for p in parts if _norm_hex(p[2]) == target]

    if is_end_exon:
        return _contiguous_parts_from_5prime(parts, target)

    if is_start_exon:
        if _norm_hex(parts[0][2]) != target:
            return _parts_after_leading_other(parts, target)
        return _contiguous_parts_from_5prime(parts, target)

    return [p for p in parts if _norm_hex(p[2]) == target]


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
    bridge_exon = start - 2 if start >= 2 else None

    segments: list[PartSegment] = []
    for e in exon_table:
        n = e["n"]
        in_main = start <= n <= end
        in_extra = extra_exon is not None and n == extra_exon
        in_bridge = (
            bridge_exon is not None
            and n == bridge_exon
            and not same_color_as_prev
        )
        if not in_main and not in_extra and not in_bridge:
            continue
        xa, xb = exon_x[n]
        w = xb - xa

        if in_bridge:
            part_list = _contiguous_parts_from_3prime(e["parts"], target)
            if len(part_list) >= len(e["parts"]):
                part_list = []
        elif in_extra:
            part_list = _parts_for_extra_exon(e["parts"], target)
        else:
            part_list = _parts_for_domain_in_exon(
                e["parts"], target,
                is_start_exon=(n == start),
                is_end_exon=(n == end),
            )

        for f0, f1, fill, border in part_list:
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


def segments_in_xrange(
    segments: list[PartSegment],
    x0: float,
    x1: float,
) -> list[PartSegment]:
    """Part segments whose absolute x-range overlaps [x0, x1]."""
    return [s for s in segments if s.x1 > x0 and s.x0 < x1]


def segments_for_exon(segments: list[PartSegment], exon_n: int) -> list[PartSegment]:
    """Part segments belonging to one exon (avoids neighbour bleed into bump padding)."""
    return [s for s in segments if s.exon == exon_n]


def exon_draw_xrange(
    exon_index: int,
    n_exons: int,
    x0: float,
    x1: float,
    bump_depth: float,
) -> tuple[float, float]:
    """Exon x-range including 3′ bump protrusion for colour and stroke clipping."""
    pad_r = bump_depth if exon_index < n_exons - 1 else 0.0
    return x0, x1 + pad_r


def part_border_color(parts: list[tuple[float, float, str, str]], *, index: int = -1) -> str:
    """Border stroke colour for a part segment (falls back to fill if needed)."""
    _f0, _f1, fill, border = parts[index]
    return border or fill
