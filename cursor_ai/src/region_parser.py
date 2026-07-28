"""
Parse exon/intron region notation (e#, i#, ranges, junctional).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

from src.config import REFERENCE

EXON_COUNT = REFERENCE.coding_exon_count

# e45, e45-e55, e45-ei45, i7-e8, i11
_REGION_TOKEN = re.compile(r"^(e|i)(\d+)$", re.I)
_RANGE_PATTERN = re.compile(
    r"^\s*(e|i)?(\d+)\s*[-–]\s*(e|i)?(\d+)\s*$", re.I
)
_SINGLE_PATTERN = re.compile(r"^\s*(e|i)(\d+)\s*$", re.I)
_JUNCTION_PATTERN = re.compile(
    r"^\s*(e|i)(\d+)\s*[-–]\s*(e|i)(\d+)\s*$", re.I
)


@dataclass
class ParsedRegion:
    """Parsed start/stop region tokens."""

    start_region: str
    stop_region: str
    start_exon: Optional[int] = None
    stop_exon: Optional[int] = None
    start_is_intron: bool = False
    stop_is_intron: bool = False
    is_junctional: bool = False
    warnings: list[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.warnings is None:
            self.warnings = []


def _format_token(prefix: str, number: int) -> str:
    return f"{prefix.lower()}{number}"


def parse_region_text(text: str) -> ParsedRegion:
    """
    Parse region strings like e45, e45-e55, i11, e45-ei45, i7-e8.
    """
    raw = (text or "").strip()
    if not raw:
        return ParsedRegion(start_region="", stop_region="", warnings=["Empty region"])

    # Junctional: e45-ei45 or i7-e8 (mixed e/i markers)
    junction = re.match(r"^\s*(e|i)(\d+)\s*[-–]\s*(e|i)(\d+)\s*$", raw, re.I)
    if junction:
        p1, n1, p2, n2 = junction.groups()
        start_tok = _format_token(p1, int(n1))
        stop_tok = _format_token(p2, int(n2))
        is_junctional = p1.lower() != p2.lower() or True
        return ParsedRegion(
            start_region=start_tok,
            stop_region=stop_tok,
            start_exon=int(n1) if p1.lower() == "e" else None,
            stop_exon=int(n2) if p2.lower() == "e" else None,
            start_is_intron=p1.lower() == "i",
            stop_is_intron=p2.lower() == "i",
            is_junctional=is_junctional,
            warnings=["Junctional region — frame may require manual review"],
        )

    # Range: e45-e55 or 45-55
    rng = re.match(r"^\s*(e)?(\d+)\s*[-–]\s*(e)?(\d+)\s*$", raw, re.I)
    if rng:
        _, n1, _, n2 = rng.groups()
        first, last = int(n1), int(n2)
        if last < first:
            first, last = last, first
        return ParsedRegion(
            start_region=_format_token("e", first),
            stop_region=_format_token("e", last),
            start_exon=first,
            stop_exon=last,
        )

    # Single: e45 or i11
    single = _SINGLE_PATTERN.match(raw)
    if single:
        prefix, num = single.groups()
        n = int(num)
        tok = _format_token(prefix, n)
        is_intron = prefix.lower() == "i"
        return ParsedRegion(
            start_region=tok,
            stop_region=tok,
            start_exon=n if not is_intron else None,
            stop_exon=n if not is_intron else None,
            start_is_intron=is_intron,
            stop_is_intron=is_intron,
        )

    return ParsedRegion(
        start_region=raw,
        stop_region=raw,
        warnings=[f"Could not fully parse region: {raw}"],
    )


def parse_quick_entry(text: str) -> dict[str, str]:
    """
    Parse quick-entry strings like 'exonic deletion e45-e55'.

    Returns partial field dict for the intake form.
    """
    raw = (text or "").strip()
    result: dict[str, str] = {}
    if not raw:
        return result

    lower = raw.lower()
    # Mutation class
    if "exonic" in lower:
        result["mutation_class"] = "Exonic"
    elif "subexonic" in lower or "sub-exonic" in lower:
        result["mutation_class"] = "Subexonic"
    elif "intronic" in lower:
        result["mutation_class"] = "Intronic"

    # Subclass
    subclasses = [
        "Duplication",
        "Deletion",
        "Nonsense",
        "Frameshift",
        "Splice",
        "Pseudoexon",
        "Exon/Intron Junctional Deletion",
    ]
    for sub in subclasses:
        if sub.lower() in lower:
            result["mutation_subclass"] = sub
            break

    # Region token(s)
    region_match = re.search(
        r"((?:e|i)\d+(?:\s*[-–]\s*(?:e|i)?\d+)?)", raw, re.I
    )
    if region_match:
        region = parse_region_text(region_match.group(1))
        result["start_region"] = region.start_region
        result["stop_region"] = region.stop_region

    return result
