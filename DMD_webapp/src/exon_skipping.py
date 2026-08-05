"""
Search for computational exon-skipping candidates that restore reading frame.

Only frame-restoring strategies are returned as displayed skip candidates.
"""

from __future__ import annotations

import itertools
from typing import Optional

from src.config import REFERENCE
from src.frame_analysis import assess_deletion_frame
from src.models import ExonRecord, FrameResult, FrameStatus, SkipCandidate, SkipEvidence
from src.coordinate_mapper import build_exon_records
from src.transcript_reconstruction import (
    ReconstructedTranscript,
    _contiguous_blocks,
    _extension_counts,
    generate_boundary_extension_candidates,
    reconstruct_deletion_skip,
)

MAX_DISPLAY_CANDIDATES = 3


def _domains_for_cds_range(cds_start: int, cds_end: int) -> list[str]:
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


def ranking_key(candidate: SkipCandidate) -> tuple:
    """Rank only frame-restoring candidates (restores_frame already filtered)."""
    return candidate.rank_score


def _build_rank_score(
    transcript: ReconstructedTranscript,
    mutation_first: int,
    mutation_last: int,
    *,
    is_advanced: bool,
) -> tuple:
    additional = transcript.additional_skipped_exons
    up_count, down_count = _extension_counts(
        mutation_first, mutation_last, set(additional),
    )
    block_count = len(_contiguous_blocks(additional))
    return (
        len(additional),
        block_count,
        transcript.additional_coding_bases_removed,
        -transcript.estimated_protein_aa,
        1 if is_advanced else 0,
        up_count,
        down_count,
        tuple(additional),
    )


def _candidate_from_transcript(
    transcript: ReconstructedTranscript,
    mutation_label: str,
    exons: list[ExonRecord],
    mutation_first: int,
    mutation_last: int,
    *,
    is_advanced: bool = False,
) -> SkipCandidate:
    index = {e.exon_number: e for e in exons}
    principal = transcript.principal_junction
    if principal:
        final_up, final_down = principal
    elif transcript.new_junctions:
        final_up, final_down = transcript.new_junctions[0]
    else:
        final_up, final_down = 0, 0

    cds_start = 1
    if transcript.new_junctions:
        up = transcript.new_junctions[0][0]
        if up in index and index[up].cumulative_cds_end:
            cds_start = index[up].cumulative_cds_end + 1

    additional = transcript.additional_skipped_exons
    up_count, down_count = _extension_counts(
        mutation_first, mutation_last, set(additional),
    )

    return SkipCandidate(
        original_mutation=mutation_label,
        deleted_exons=transcript.original_mutation_exons,
        additional_skipped_exons=additional,
        all_removed_exons=transcript.all_removed_exons,
        retained_exons=transcript.retained_exons,
        new_junctions=transcript.new_junctions,
        mutation_boundary_junction=transcript.mutation_boundary_junction,
        principal_junction=principal,
        final_upstream_exon=final_up,
        final_downstream_exon=final_down,
        total_coding_bases_removed=transcript.total_coding_bases_removed,
        additional_coding_bases_removed=transcript.additional_coding_bases_removed,
        restores_frame=transcript.restores_frame,
        estimated_remaining_coding_bp=transcript.estimated_remaining_coding_bp,
        estimated_protein_aa=transcript.estimated_protein_aa,
        evidence_class=SkipEvidence.COMPUTATIONAL,
        rank_score=_build_rank_score(
            transcript, mutation_first, mutation_last, is_advanced=is_advanced,
        ),
        is_boundary_adjacent=transcript.is_boundary_adjacent,
        is_structurally_valid=transcript.is_structurally_valid,
        is_advanced_noncontiguous=is_advanced,
        skip_block_count=len(_contiguous_blocks(additional)),
        upstream_extension_count=up_count,
        downstream_extension_count=down_count,
        assumptions=[
            "Computational frame restoration only — not therapeutic evidence.",
            "Assumes additional exons can be skipped without disrupting splice regulation.",
            "Frame assessed at every novel junction in the reconstructed transcript.",
        ],
        affected_domains=_domains_for_cds_range(
            cds_start, transcript.estimated_remaining_coding_bp,
        ),
    )


def evaluate_candidate(
    mutation_first: int,
    mutation_last: int,
    additional_skips: tuple[int, ...],
    exons: list[ExonRecord],
    mutation_label: str,
    *,
    is_advanced: bool = False,
) -> Optional[SkipCandidate]:
    """Reconstruct and return a candidate only when frame is restored."""
    transcript = reconstruct_deletion_skip(
        mutation_first,
        mutation_last,
        additional_skips,
        exons,
        require_frame=False,
    )
    if not transcript or not transcript.is_structurally_valid:
        return None
    if not transcript.restores_frame:
        return None
    return _candidate_from_transcript(
        transcript, mutation_label, exons, mutation_first, mutation_last,
        is_advanced=is_advanced,
    )


def _filter_valid_candidates(candidates: list[SkipCandidate]) -> list[SkipCandidate]:
    """Keep only structurally valid, frame-restoring candidates."""
    return [
        c for c in candidates
        if c is not None and c.is_structurally_valid and c.restores_frame is True
    ]


def find_skip_candidates(
    mutation_first: int,
    mutation_last: int,
    exons: Optional[list[ExonRecord]] = None,
    *,
    max_additional_skips: int = 3,
    mutation_label: str = "",
    advanced: bool = False,
    max_results: int = MAX_DISPLAY_CANDIDATES,
) -> list[SkipCandidate]:
    """
    Search boundary-extension strategies; return only frame-restoring candidates.

  ``max_additional_skips`` controls search depth, not how many results are shown.
    """
    records = exons if exons is not None else build_exon_records()
    label = mutation_label or f"del{mutation_first}-{mutation_last}"

    base_frame: FrameResult = assess_deletion_frame(mutation_first, mutation_last, records)
    if base_frame.status == FrameStatus.IN_FRAME:
        return []

    generated = generate_boundary_extension_candidates(
        mutation_first, mutation_last, max_additional_skips=max_additional_skips,
    )

    evaluated: list[SkipCandidate] = []
    for combo in generated:
        candidate = evaluate_candidate(
            mutation_first, mutation_last, combo, records, label, is_advanced=False,
        )
        if candidate:
            evaluated.append(candidate)

    if advanced:
        retained = [
            e.exon_number for e in records
            if e.exon_number < mutation_first or e.exon_number > mutation_last
        ]
        seen = {tuple(c.additional_skipped_exons) for c in evaluated}
        base_combos = set(generated)
        for size in range(1, max_additional_skips + 1):
            for combo in itertools.combinations(retained, size):
                if combo in seen or combo in base_combos:
                    continue
                candidate = evaluate_candidate(
                    mutation_first, mutation_last, combo, records, label,
                    is_advanced=True,
                )
                if candidate:
                    candidate.assumptions.append(
                        "Advanced disconnected skip strategy — not a single boundary extension."
                    )
                    evaluated.append(candidate)
                    seen.add(combo)

    valid = _filter_valid_candidates(evaluated)
    ranked = sorted(valid, key=ranking_key)

    unique: list[SkipCandidate] = []
    seen_keys: set[tuple[int, ...]] = set()
    for cand in ranked:
        key = tuple(cand.additional_skipped_exons)
        if key not in seen_keys:
            seen_keys.add(key)
            unique.append(cand)

    return unique[:max_results]


# Backward-compatible aliases
_evaluate_skip_set = evaluate_candidate
_candidate_from_transcript_export = _candidate_from_transcript
