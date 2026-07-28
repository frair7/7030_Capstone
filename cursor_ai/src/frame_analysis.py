"""
Reading-frame analysis for whole-exon deletions and tandem duplications.

Calculations use CDS lengths and splice phases from the reference table,
not hard-coded in-frame/out-of-frame lists.
"""

from __future__ import annotations

from typing import Optional

from src.config import REFERENCE
from src.models import ExonRecord, FrameResult, FrameStatus, ParsedVariant, VariantType
from src.coordinate_mapper import build_exon_records


def _exon_index(exons: list[ExonRecord]) -> dict[int, ExonRecord]:
    return {e.exon_number: e for e in exons}


def total_coding_bp(exons: list[ExonRecord]) -> int:
    """Sum coding bases across all exons."""
    return sum(e.coding_length_bp for e in exons)


def coding_bp_for_exon_range(
    exons: list[ExonRecord], first: int, last: int
) -> int:
    """Sum coding bases for inclusive exon range."""
    index = _exon_index(exons)
    return sum(index[i].coding_length_bp for i in range(first, last + 1))


def _junction_compatible(
    upstream_phase_3: int,
    removed_or_added_bp: int,
    downstream_phase_5: int,
) -> bool:
    """Check splice-phase compatibility at a novel exon junction."""
    return (upstream_phase_3 + removed_or_added_bp) % 3 == downstream_phase_5


def assess_deletion_frame(
    first_exon: int,
    last_exon: int,
    exons: Optional[list[ExonRecord]] = None,
) -> FrameResult:
    """
    Assess reading frame after deleting exons ``first_exon`` through ``last_exon``.

    The junction joins the upstream flanking exon (first-1) to the downstream
    flanking exon (last+1).
    """
    records = exons if exons is not None else build_exon_records()
    index = _exon_index(records)
    assumptions = [
        "Whole-exon deletion with canonical splice junctions assumed.",
        "Partial exon deletions are not evaluated by this function.",
    ]

    if first_exon < 1 or last_exon > REFERENCE.coding_exon_count or last_exon < first_exon:
        return FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation="Invalid exon range for deletion analysis.",
            assumptions=assumptions,
        )

    removed_bp = coding_bp_for_exon_range(records, first_exon, last_exon)
    upstream = first_exon - 1
    downstream = last_exon + 1

    if upstream < 1 or downstream > REFERENCE.coding_exon_count:
        return FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation=(
                f"Deletion of exons {first_exon}-{last_exon} removes the "
                f"{'5′' if upstream < 1 else '3′'} end of the transcript; "
                "frame cannot be assessed as a simple internal junction."
            ),
            coding_bases_removed=removed_bp,
            assumptions=assumptions,
        )

    up_exon = index[upstream]
    down_exon = index[downstream]
    if up_exon.splice_phase_3prime is None or down_exon.splice_phase_5prime is None:
        return FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation="Missing splice phase data for flanking exons.",
            assumptions=assumptions,
        )

    phase_up = up_exon.splice_phase_3prime
    phase_down = down_exon.splice_phase_5prime
    junction_phase = (phase_up + removed_bp) % 3
    remaining = total_coding_bp(records) - removed_bp
    in_frame = removed_bp % 3 == 0 and junction_phase == phase_down

    if in_frame:
        status = FrameStatus.IN_FRAME
        explanation = (
            f"Deleting exons {first_exon}-{last_exon} removes {removed_bp} coding "
            f"bases ({removed_bp % 3} mod 3 remainder). The new junction between "
            f"exon {upstream} and exon {downstream} preserves the reading frame "
            f"(phase {phase_up} → {junction_phase}, matching exon {downstream} "
            f"5′ phase {phase_down})."
        )
    else:
        status = FrameStatus.OUT_OF_FRAME
        reasons = []
        if removed_bp % 3 != 0:
            reasons.append(
                f"{removed_bp} deleted coding bases is not divisible by 3 "
                f"(remainder {removed_bp % 3})"
            )
        if junction_phase != phase_down:
            reasons.append(
                f"junction phase {junction_phase} ≠ downstream 5′ phase {phase_down}"
            )
        explanation = (
            f"Deleting exons {first_exon}-{last_exon} is predicted OUT OF FRAME. "
            f"Junction exon {upstream}|{downstream}: " + "; ".join(reasons) + "."
        )

    return FrameResult(
        status=status,
        explanation=explanation,
        upstream_exon=upstream,
        downstream_exon=downstream,
        upstream_phase_3prime=phase_up,
        downstream_phase_5prime=phase_down,
        junction_phase=junction_phase,
        coding_bases_removed=removed_bp,
        remaining_coding_bp=remaining,
        estimated_protein_aa=remaining // 3,
        assumptions=assumptions,
    )


def assess_duplication_frame(
    first_exon: int,
    last_exon: int,
    exons: Optional[list[ExonRecord]] = None,
) -> FrameResult:
    """
    Assess reading frame after tandem duplication of exons first..last.

    Models an extra copy inserted immediately downstream of the original block
    in transcript order. Exact duplication structure may vary in real genomes.
    """
    records = exons if exons is not None else build_exon_records()
    index = _exon_index(records)
    assumptions = [
        "Tandem duplication immediately downstream of the original exon block.",
        "Duplicated exons retain internal splice phases from the reference model.",
        "Inverted or complex rearrangements are not modeled.",
    ]

    if first_exon < 1 or last_exon > REFERENCE.coding_exon_count or last_exon < first_exon:
        return FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation="Invalid exon range for duplication analysis.",
            assumptions=assumptions,
        )

    added_bp = coding_bp_for_exon_range(records, first_exon, last_exon)
    upstream = first_exon - 1 if first_exon > 1 else None
    downstream = last_exon + 1 if last_exon < REFERENCE.coding_exon_count else None

    total = total_coding_bp(records) + added_bp
    in_frame = added_bp % 3 == 0

    if upstream and downstream:
        up_exon = index[upstream]
        block_end = index[last_exon]
        down_exon = index[downstream]
        if (
            up_exon.splice_phase_3prime is not None
            and block_end.splice_phase_3prime is not None
            and down_exon.splice_phase_5prime is not None
        ):
            # Junction at end of duplicated block joining to downstream exon
            dup_block_phase = (block_end.splice_phase_3prime + added_bp) % 3
            in_frame = in_frame and dup_block_phase == down_exon.splice_phase_5prime

    if in_frame:
        status = FrameStatus.IN_FRAME
        explanation = (
            f"Tandem duplication of exons {first_exon}-{last_exon} adds "
            f"{added_bp} coding bases. The net CDS length change is divisible "
            f"by 3, preserving the reading frame under the tandem model."
        )
    else:
        status = FrameStatus.OUT_OF_FRAME
        explanation = (
            f"Tandem duplication of exons {first_exon}-{last_exon} adds "
            f"{added_bp} coding bases ({added_bp % 3} mod 3 remainder), "
            "predicted OUT OF FRAME under the tandem duplication model."
        )

    return FrameResult(
        status=status,
        explanation=explanation,
        upstream_exon=upstream,
        downstream_exon=downstream,
        coding_bases_added=added_bp,
        remaining_coding_bp=total,
        estimated_protein_aa=total // 3,
        assumptions=assumptions,
    )


def assess_variant_frame(
    variant: ParsedVariant,
    exons: Optional[list[ExonRecord]] = None,
) -> FrameResult:
    """Dispatch frame analysis based on parsed variant type."""
    if not variant.is_valid:
        return FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation="Cannot analyze frame: variant input is invalid.",
        )

    if variant.variant_type == VariantType.DELETION and variant.exon_range:
        first, last = variant.exon_range
        return assess_deletion_frame(first, last, exons)

    if variant.variant_type == VariantType.DUPLICATION and variant.exon_range:
        first, last = variant.exon_range
        return assess_duplication_frame(first, last, exons)

    if variant.variant_type in {VariantType.SUBSTITUTION, VariantType.SPLICE}:
        return FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation=(
                "Single-nucleotide or splice-site variants require defined "
                "sequence consequences; whole-exon frame rules do not apply."
            ),
            assumptions=["RNA-level consequence may be required."],
        )

    return FrameResult(
        status=FrameStatus.CANNOT_DETERMINE,
        explanation="Frame effect cannot be calculated for this variant type.",
    )
