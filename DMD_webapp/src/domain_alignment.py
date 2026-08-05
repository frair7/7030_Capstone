"""
Align domain-row bubbles to transcript exon part boundaries.

Biological x-coordinates come from EXON_TABLE.parts via find_part_segments.
Draw coordinates apply a visual inset only (DOMAIN_GAP) — never shift biology.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any

from src.part_segments import find_part_segments

DOMAIN_GAP_PX = 0.75
DOMAIN_GAP_DATA = 0.04

SHOW_ALIGNMENT_GUIDES = os.environ.get("SHOW_ALIGNMENT_GUIDES", "").lower() in (
    "1", "true", "yes",
)


@dataclass(frozen=True)
class DomainBubble:
    """Legacy bubble — x0/x1 are biological (alignment) coordinates."""

    x0: float
    x1: float
    label: str
    fill: str
    border: str
    is_repeat: bool = False


@dataclass(frozen=True)
class ResolvedDomain:
    label: str
    fill: str
    border: str
    biological_x0: float
    biological_x1: float
    draw_x0: float
    draw_x1: float
    is_repeat: bool = False

    @property
    def biological_width(self) -> float:
        return self.biological_x1 - self.biological_x0

    @property
    def draw_width(self) -> float:
        return self.draw_x1 - self.draw_x0


def _norm_hex(color: str) -> str:
    return color.strip().upper()


def _is_repeat_label(label: str) -> bool:
    return len(label) >= 2 and label[0] == "R" and label[1:].isdigit()


def draw_coords(
    biological_x0: float,
    biological_x1: float,
    gap: float,
) -> tuple[float, float]:
    """Visual inset for domain boxes; skip inset when the span is too narrow."""
    width = biological_x1 - biological_x0
    if width <= gap:
        return biological_x0, biological_x1
    return biological_x0 + gap / 2, biological_x1 - gap / 2


def domain_corner_radius(width: float, height: float) -> float:
    return min(4.0, height * 0.18, width * 0.18)


def domain_font_size(width: float, height: float) -> float:
    return max(4.0, min(12.0, min(width, height) * 0.35))


def resolve_domain(
    label: str,
    domain_map: list[tuple[int, int, str, str]],
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
    *,
    gap: float = DOMAIN_GAP_PX,
) -> ResolvedDomain | None:
    segments = find_part_segments(label, domain_map, exon_table, exon_x)
    if not segments:
        return None

    entry = next(d for d in domain_map if d[2] == label)
    _, _, _, color = entry
    biological_x0 = min(s.x0 for s in segments)
    biological_x1 = max(s.x1 for s in segments)
    if biological_x1 <= biological_x0:
        return None

    draw_x0, draw_x1 = draw_coords(biological_x0, biological_x1, gap)
    border = segments[0].border
    return ResolvedDomain(
        label=label,
        fill=color,
        border=border,
        biological_x0=biological_x0,
        biological_x1=biological_x1,
        draw_x0=draw_x0,
        draw_x1=draw_x1,
        is_repeat=_is_repeat_label(label),
    )


def resolve_domains(
    domain_map: list[tuple[int, int, str, str]],
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
    *,
    gap: float = DOMAIN_GAP_PX,
) -> list[ResolvedDomain]:
    domains: list[ResolvedDomain] = []
    for _, _, label, _ in domain_map:
        resolved = resolve_domain(label, domain_map, exon_table, exon_x, gap=gap)
        if resolved is not None:
            domains.append(resolved)
    return domains


def build_aligned_domain_bubbles(
    domain_map: list[tuple[int, int, str, str]],
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
    *,
    repeat_gap: float = 0.0,  # noqa: ARG001 — kept for API compat; no longer used
) -> list[DomainBubble]:
    """One bubble per DOMAIN_MAP entry with biological x extents from part segments."""
    del repeat_gap
    return [
        DomainBubble(
            x0=d.biological_x0,
            x1=d.biological_x1,
            label=d.label,
            fill=d.fill,
            border=d.border,
            is_repeat=d.is_repeat,
        )
        for d in resolve_domains(domain_map, exon_table, exon_x, gap=0.0)
    ]


def hinge_guide_x_positions(
    domain_map: list[tuple[int, int, str, str]],
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
) -> list[float]:
    """Dashed guide x-positions at the 5′ start and 3′ end of each hinge (H1–H4)."""
    domains = resolve_domains(domain_map, exon_table, exon_x, gap=0.0)
    by_label = {d.label: d for d in domains}
    positions: list[float] = []
    for label in ("H1", "H2", "H3", "H4"):
        if label in by_label:
            d = by_label[label]
            positions.append(d.biological_x0)
            positions.append(d.biological_x1)
    return positions


def validate_domain_alignment(
    domain_map: list[tuple[int, int, str, str]],
    exon_table: list[dict[str, Any]],
    exon_x: dict[int, tuple[float, float]],
) -> None:
    """Assert H1 and R1–R24 boundaries match part segments without overlap."""
    domains = resolve_domains(domain_map, exon_table, exon_x, gap=0.0)
    by_label = {d.label: d for d in domains}

    h1 = by_label.get("H1")
    if h1 is None:
        raise AssertionError("H1 domain not resolved")

    xa8, xb8 = exon_x[8]
    w8 = xb8 - xa8
    h1_e8_x0 = xa8 + 0.3901 * w8
    if abs(h1.biological_x0 - h1_e8_x0) > 1e-6:
        raise AssertionError(
            f"H1 x0 {h1.biological_x0} != exon 8 part start {h1_e8_x0}"
        )

    xa10, xb10 = exon_x[10]
    w10 = xb10 - xa10
    h1_e10_x1 = xa10 + 0.2540 * w10
    if abs(h1.biological_x1 - h1_e10_x1) > 1e-6:
        raise AssertionError(
            f"H1 x1 {h1.biological_x1} != exon 10 part end {h1_e10_x1}"
        )

    repeats = [d for d in domains if d.is_repeat]
    if len(repeats) != 24:
        raise AssertionError(f"expected 24 repeat domains, got {len(repeats)}")

    for prev, curr in zip(repeats, repeats[1:]):
        if prev.biological_x1 > curr.biological_x0 + 1e-6:
            raise AssertionError(
                f"{prev.label} x1 {prev.biological_x1} overlaps "
                f"{curr.label} x0 {curr.biological_x0}"
            )

    r1 = by_label.get("R1")
    if r1 is not None and abs(h1.biological_x1 - r1.biological_x0) > 1e-6:
        raise AssertionError("H1 end does not meet R1 start")
