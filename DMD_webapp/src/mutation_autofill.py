"""
Auto-fill mutation catalog fields from class, subclass, and region notation.
"""

from __future__ import annotations

from typing import Any, Optional

from src.coordinate_mapper import build_exon_records
from src.frame_analysis import (
    assess_deletion_frame,
    assess_duplication_frame,
    coding_bp_for_exon_range,
)
from src.models import FrameStatus
from src.reference_data import load_domain_table
from src.region_parser import ParsedRegion, parse_region_text

# Average protein MW estimate: ~110 Da per amino acid
DA_PER_AA = 110.0
TOTAL_CODING_BP = 11058  # Dp427m CDS incl. stop


def _frame_label(status: FrameStatus) -> str:
    mapping = {
        FrameStatus.IN_FRAME: "In-frame",
        FrameStatus.OUT_OF_FRAME: "Out-of-frame",
        FrameStatus.CANNOT_DETERMINE: "Cannot determine",
    }
    return mapping.get(status, "Cannot determine")


def _domains_for_exon_range(first_exon: int, last_exon: int) -> str:
    """Return comma-separated domain names overlapping deleted exon CDS."""
    try:
        exons = build_exon_records()
        index = {e.exon_number: e for e in exons}
        cds_start = index[first_exon].cds_start
        cds_end = index[last_exon].cds_end
        if cds_start is None or cds_end is None:
            return ""
        domains = load_domain_table()
        names: list[str] = []
        for row in domains:
            d_start = int(row["cds_start_bp"])
            d_end = int(row["cds_end_bp"])
            if not (cds_end < d_start or cds_start > d_end):
                names.append(row["domain_name"])
        return ", ".join(names)
    except Exception:
        return ""


def _estimate_protein_kda(remaining_coding_bp: int) -> str:
    aa = remaining_coding_bp // 3
    kda = (aa * DA_PER_AA) / 1000.0
    return f"{kda:.1f}"


def _molecular_consequence(
    mutation_class: str,
    mutation_subclass: str,
    frame: str,
) -> str:
    parts = [mutation_class, mutation_subclass]
    if frame and frame not in ("", "Cannot determine", "N/A"):
        parts.append(f"({frame.lower()})")
    return " ".join(p for p in parts if p)


def autofill_mutation_fields(
    mutation_class: str = "",
    mutation_subclass: str = "",
    start_region: str = "",
    stop_region: str = "",
    *,
    existing: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """
    Compute auto-fillable catalog fields from structured intake values.

    Only overwrites fields that can be derived; preserves manual entries
    in ``existing`` for fields that cannot be computed.
    """
    base = dict(existing or {})
    result: dict[str, Any] = {
        "mutation_class": mutation_class or base.get("mutation_class", ""),
        "mutation_subclass": mutation_subclass or base.get("mutation_subclass", ""),
        "start_region": start_region or base.get("start_region", ""),
        "stop_region": stop_region or base.get("stop_region", ""),
    }

    # Parse regions
    start_parsed = parse_region_text(result["start_region"])
    stop_parsed = parse_region_text(result["stop_region"])

    # If only start provided as range, split it
    if start_parsed.start_exon and not stop_parsed.stop_exon:
        if start_parsed.stop_exon:
            result["stop_region"] = start_parsed.stop_region
            stop_parsed = start_parsed

    region = start_parsed
    if stop_parsed.stop_exon:
        region = ParsedRegion(
            start_region=start_parsed.start_region or stop_parsed.start_region,
            stop_region=stop_parsed.stop_region,
            start_exon=start_parsed.start_exon or stop_parsed.start_exon,
            stop_exon=stop_parsed.stop_exon or start_parsed.stop_exon,
            is_junctional=start_parsed.is_junctional or stop_parsed.is_junctional,
            warnings=(start_parsed.warnings or []) + (stop_parsed.warnings or []),
        )

    mclass = (result["mutation_class"] or "").strip()
    msub = (result["mutation_subclass"] or "").strip()

    # Frame + protein size for whole-exon events
    if (
        mclass.lower() == "exonic"
        and region.start_exon
        and region.stop_exon
        and not region.is_junctional
        and not region.start_is_intron
    ):
        first, last = region.start_exon, region.stop_exon
        if first > last:
            first, last = last, first
        exons = build_exon_records()

        if msub.lower() == "deletion":
            frame_result = assess_deletion_frame(first, last, exons)
            result["frame"] = _frame_label(frame_result.status)
            if frame_result.remaining_coding_bp:
                result["expected_protein_size_kda"] = _estimate_protein_kda(
                    frame_result.remaining_coding_bp
                )
            result["domains_affected"] = _domains_for_exon_range(first, last)
            result["molecular_consequence"] = _molecular_consequence(
                mclass, msub, result["frame"]
            )
        elif msub.lower() == "duplication":
            frame_result = assess_duplication_frame(first, last, exons)
            result["frame"] = _frame_label(frame_result.status)
            if frame_result.remaining_coding_bp:
                result["expected_protein_size_kda"] = _estimate_protein_kda(
                    frame_result.remaining_coding_bp
                )
            added = coding_bp_for_exon_range(exons, first, last)
            result["domains_affected"] = _domains_for_exon_range(first, last)
            result["molecular_consequence"] = _molecular_consequence(
                mclass, msub, result["frame"]
            )
            _ = added
    elif msub.lower() in {"nonsense", "frameshift", "splice"}:
        result["frame"] = "Out-of-frame"
        result["molecular_consequence"] = _molecular_consequence(
            mclass, msub, result.get("frame", "")
        )
    elif region.is_junctional or region.start_is_intron:
        result["frame"] = result.get("frame") or "Cannot determine"
        result["molecular_consequence"] = _molecular_consequence(mclass, msub, "")

    # Preserve manual fields not computed
    for key in (
        "cdna",
        "rna",
        "protein",
        "experimental_research_category",
        "tissue_category_comments",
        "general_comments",
    ):
        if key in base and base[key]:
            result[key] = base[key]

    return result
