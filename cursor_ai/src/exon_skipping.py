"""
Search for computational exon-skipping candidates that restore reading frame.

Candidates are ranked by transparent, non-clinical criteria only.
"""

from __future__ import annotations

import itertools
from typing import Optional

from src.config import REFERENCE
from src.frame_analysis import assess_deletion_frame, coding_bp_for_exon_range, total_coding_bp
from src.models import ExonRecord, FrameResult, FrameStatus, SkipCandidate, SkipEvidence
from src.coordinate_mapper import build_exon_records


def _domains_for_cds_range(cds_start: int, cds_end: int) -> list[str]:
    """Return domain names overlapping a CDS range (requires domain table)."""
    from src.reference_data import load_domain_table

    affected: list[str] = []
    try:
        domains = load_domain_table()
    except FileNotFoundError:
        return affected
    for row in domains:
        d_start = int(row["cds_start_bp"])
        d_end = int(row["cds_end_bp"])
        if not (cds_end < d_start or cds_start > d_end):
            affected.append(row["domain_name"])
    return affected


def _evaluate_skip_set(
    mutation_first: int,
    mutation_last: int,
    additional_skips: tuple[int, ...],
    exons: list[ExonRecord],
    mutation_label: str,
) -> Optional[SkipCandidate]:
    """Evaluate one additional-skip combination."""
    index = {e.exon_number: e for e in exons}
    deleted = set(range(mutation_first, mutation_last + 1))
    skip_set = set(additional_skips)

    if skip_set & deleted:
        return None  # cannot skip already-deleted exons

    all_removed = sorted(deleted | skip_set)
    if not all_removed:
        return None

    first_removed = all_removed[0]
    last_removed = all_removed[-1]
    upstream = first_removed - 1
    downstream = last_removed + 1

    if upstream < 1 or downstream > REFERENCE.coding_exon_count:
        return None

    removed_bp = sum(index[e].coding_length_bp for e in all_removed)
    additional_bp = sum(index[e].coding_length_bp for e in skip_set)

    up_exon = index[upstream]
    down_exon = index[downstream]
    if up_exon.splice_phase_3prime is None or down_exon.splice_phase_5prime is None:
        return None

    junction_phase = (up_exon.splice_phase_3prime + removed_bp) % 3
    restores = removed_bp % 3 == 0 and junction_phase == down_exon.splice_phase_5prime

    if not restores:
        return None

    total_bp = total_coding_bp(exons)
    remaining = total_bp - removed_bp
    contiguous = int(
        skip_set == set(range(min(skip_set), max(skip_set) + 1)) if skip_set else 1
    )
    # Lower rank_score is better: (additional exons, additional bp, non-contiguous penalty)
    rank = (len(skip_set), additional_bp, 0 if contiguous else 1)

    cds_removed_start = up_exon.cumulative_cds_end + 1 if up_exon.cumulative_cds_end else 1
    cds_removed_end = down_exon.cumulative_cds_start - 1 if down_exon.cumulative_cds_start else remaining

    return SkipCandidate(
        original_mutation=mutation_label,
        deleted_exons=list(range(mutation_first, mutation_last + 1)),
        additional_skipped_exons=sorted(skip_set),
        final_upstream_exon=upstream,
        final_downstream_exon=downstream,
        total_coding_bases_removed=removed_bp,
        additional_coding_bases_removed=additional_bp,
        restores_frame=True,
        estimated_remaining_coding_bp=remaining,
        estimated_protein_aa=remaining // 3,
        evidence_class=SkipEvidence.COMPUTATIONAL,
        rank_score=rank,
        assumptions=[
            "Computational frame restoration only — not therapeutic evidence.",
            "Assumes additional exons can be skipped without disrupting splice regulation.",
        ],
        affected_domains=_domains_for_cds_range(cds_removed_start, cds_removed_end),
    )


def find_skip_candidates(
    mutation_first: int,
    mutation_last: int,
    exons: Optional[list[ExonRecord]] = None,
    *,
    max_additional_skips: int = 3,
    mutation_label: str = "",
    advanced: bool = False,
) -> list[SkipCandidate]:
    """
    Search for additional exon skips that restore frame after a deletion.

    By default searches adjacent/contiguous strategies near deletion boundaries.
    Set ``advanced=True`` to include non-contiguous multi-exon combinations.
    """
    records = exons if exons is not None else build_exon_records()
    label = mutation_label or f"del{mutation_first}-{mutation_last}"

    base_frame: FrameResult = assess_deletion_frame(mutation_first, mutation_last, records)
    if base_frame.status == FrameStatus.IN_FRAME:
        return []

    candidates: list[SkipCandidate] = []

    # Candidate pool: retained exons adjacent to deletion boundaries, expanding outward
    adjacent_pool: list[int] = []
    for offset in range(1, 6):
        up = mutation_first - offset
        down = mutation_last + offset
        if up >= 1:
            adjacent_pool.append(up)
        if down <= REFERENCE.coding_exon_count:
            adjacent_pool.append(down)

    search_sizes = range(1, max_additional_skips + 1)
    for size in search_sizes:
        for combo in itertools.combinations(adjacent_pool, size):
            if not advanced and size > 1:
                # Require contiguous additional skips in basic mode
                sorted_combo = sorted(combo)
                if sorted_combo != list(range(sorted_combo[0], sorted_combo[-1] + 1)):
                    continue
            candidate = _evaluate_skip_set(
                mutation_first, mutation_last, combo, records, label
            )
            if candidate:
                candidates.append(candidate)

    if advanced:
        # Broader search across all retained exons (still bounded by max_additional_skips)
        retained = [
            e.exon_number
            for e in records
            if e.exon_number < mutation_first or e.exon_number > mutation_last
        ]
        for size in search_sizes:
            for combo in itertools.combinations(retained, size):
                if combo in [tuple(c.additional_skipped_exons) for c in candidates]:
                    continue
                candidate = _evaluate_skip_set(
                    mutation_first, mutation_last, combo, records, label
                )
                if candidate:
                    candidates.append(candidate)

    # Deduplicate by additional skip set
    seen: set[tuple[int, ...]] = set()
    unique: list[SkipCandidate] = []
    for cand in sorted(candidates, key=lambda c: c.rank_score):
        key = tuple(cand.additional_skipped_exons)
        if key not in seen:
            seen.add(key)
            unique.append(cand)

    return unique
