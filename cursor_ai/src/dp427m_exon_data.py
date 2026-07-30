"""
Authoritative Dp427m exon table and domain map.

Sourced from the user's Claude make_figure.py reference (scripts/reference_make_figure.py),
which encodes per-exon CDS bp, puzzle-piece shapes, and fractional domain sub-segments
derived from the supplied domain/feature nucleotide coordinate table.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

REFERENCE_SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "reference_make_figure.py"
)

ISOFORM_ANNOTATIONS: list[tuple[int, str]] = [
    (1, "Dp427m,b,c,l"),
    (30, "Dp260"),
    (56, "Dp140"),
    (62, "Dp71"),
]

HINGE_GUIDE_AT_EXONS = [17, 18, 50, 51]

PHENOTYPE_OPTIONS = [
    "",
    "Healthy control",
    "Asymptomatic",
    "Paucisymptomatic",
    "IMD",
    "BMD",
    "DMD",
]


@lru_cache(maxsize=2)
def load_dp427m_data(mtime: float) -> tuple[list[dict[str, Any]], list[tuple[int, int, str, str]]]:
    del mtime  # cache key only — reload when reference script file changes
    src = REFERENCE_SCRIPT.read_text(encoding="utf-8")
    start = src.index("EXON_TABLE = [")
    end = src.index("# 3. Geometry helpers")
    ns: dict[str, Any] = {}
    exec(src[start:end], ns)  # noqa: S102 — trusted local reference constants
    return ns["EXON_TABLE"], ns["DOMAIN_MAP"]


def _reference_mtime() -> float:
    return REFERENCE_SCRIPT.stat().st_mtime


def get_exon_table() -> list[dict[str, Any]]:
    return load_dp427m_data(_reference_mtime())[0]


def get_domain_map() -> list[tuple[int, int, str, str]]:
    return load_dp427m_data(_reference_mtime())[1]


def bar_style_for_phenotype(phenotype: str) -> str:
    p = (phenotype or "").strip().lower()
    if p == "dmd":
        return "black"
    if p in ("asymptomatic", "paucisymptomatic", "healthy control"):
        return "outline"
    return "grey"
