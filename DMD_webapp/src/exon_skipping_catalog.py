"""Catalog-aware exon-skipping analysis for the Streamlit page."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from src.config import REFERENCE
from src.coordinate_mapper import build_exon_records
from src.exon_skipping import MAX_DISPLAY_CANDIDATES
from src.exon_skipping_analysis import (
    ExonMutation as AnalysisMutation,
    MutationClass as AnalysisMutationClass,
    filter_candidates_by_target_exons,
    find_frame_restoring_candidates,
)
from src.frame_analysis import assess_deletion_frame, assess_duplication_frame
from src.models import ExonRecord, FrameResult, FrameStatus, SkipCandidate, SkipEvidence
from src.mutation_catalog import active_catalog
from src.mutation_map_svg import SkipSchematicOverlay
from src.mutation_viz import exon_range_from_row, sort_catalog_dataframe

MUTATION_CLASS_FILTERS = {
    "Deletion": "deletion",
    "Duplication": "duplication",
    "Nonsense": "nonsense",
    "Frameshift": "frameshift",
    "Splice-altering": "splice",
}


def format_exon_list(exons: list[int]) -> str:
    """Format exon number(s) for display."""
    if not exons:
        return "—"
    exons = sorted(exons)
    if len(exons) == 1:
        return f"Exon {exons[0]}"
    if exons == list(range(exons[0], exons[-1] + 1)):
        return f"Exons {exons[0]}-{exons[-1]}"
    return "Exons " + ", ".join(str(e) for e in exons)


def format_skip_candidate(candidate: Optional[SkipCandidate]) -> str:
    """Format a skip candidate's additional skipped exons for table display."""
    if candidate is None:
        return "—"
    if not candidate.additional_skipped_exons:
        return "—"
    return format_exon_list(candidate.additional_skipped_exons)


def format_skip_target_only(candidate: Optional[SkipCandidate]) -> str:
    """Format only the additional/target exons skipped (not mutation-deleted exons)."""
    if candidate is None:
        return "—"
    if not candidate.additional_skipped_exons:
        return "—"
    return format_exon_list(candidate.additional_skipped_exons)


def mutation_description(row: dict[str, Any]) -> str:
    """Human-readable mutation label from a catalog row."""
    pid = str(row.get("id", "")).strip()
    mclass = str(row.get("mutation_class", "")).strip()
    msub = str(row.get("mutation_subclass", "")).strip()
    start = str(row.get("start_region", "")).strip()
    stop = str(row.get("stop_region", "")).strip()
    region = start if start == stop else f"{start}-{stop}"
    parts = [p for p in (pid, mclass, msub, region) if p]
    return " / ".join(parts) if parts else "Unknown mutation"


def validate_skip_target(start_exon: int, end_exon: int) -> tuple[list[int], Optional[str]]:
    """Validate a contiguous skip-target exon range."""
    n = REFERENCE.coding_exon_count
    if start_exon < 1 or end_exon > n:
        return [], f"Exon numbers must be between 1 and {n}."
    if end_exon < start_exon:
        return [], "End exon must be greater than or equal to the start exon."
    target = list(range(start_exon, end_exon + 1))
    if target != list(range(target[0], target[-1] + 1)):
        return [], "Skip target must be a contiguous exon range."
    return target, None


def frame_result_for_catalog_row(
    row: dict[str, Any],
    exons: Optional[list[ExonRecord]] = None,
) -> tuple[Optional[tuple[int, int]], FrameResult]:
    """Return exon range and frame assessment for a catalog row."""
    records = exons if exons is not None else build_exon_records()
    rng = exon_range_from_row(row)
    if not rng:
        return None, FrameResult(
            status=FrameStatus.CANNOT_DETERMINE,
            explanation="Could not parse exon range from catalog row.",
        )

    first, last = rng
    msub = str(row.get("mutation_subclass", "")).strip().lower()
    if msub == "deletion":
        return rng, assess_deletion_frame(first, last, records)
    if msub == "duplication":
        return rng, assess_duplication_frame(first, last, records)
    return rng, FrameResult(
        status=FrameStatus.CANNOT_DETERMINE,
        explanation=f"Frame rules do not apply to subclass '{msub}'.",
    )


def _row_matches_class_filter(row: dict[str, Any], class_filter: Optional[set[str]]) -> bool:
    if not class_filter:
        return True
    msub = str(row.get("mutation_subclass", "")).strip().lower()
    for label in class_filter:
        key = MUTATION_CLASS_FILTERS.get(label, label.lower())
        if key == "deletion" and msub == "deletion":
            return True
        if key == "duplication" and msub == "duplication":
            return True
        if key == "nonsense" and msub == "nonsense":
            return True
        if key == "frameshift" and "frame" in msub:
            return True
        if key == "splice" and msub in {"splice", "splice-altering"}:
            return True
    return False


def rank_skip_candidates_for_row(
    row: dict[str, Any],
    exons: Optional[list[ExonRecord]] = None,
    *,
    max_additional_skips: int = 3,
    advanced: bool = False,
) -> tuple[FrameResult, list[SkipCandidate]]:
    """Rank mutation-specific, frame-restoring candidates for one row."""
    records = exons if exons is not None else build_exon_records()
    rng, frame = frame_result_for_catalog_row(row, records)
    label = mutation_description(row)

    if not rng or frame.status != FrameStatus.OUT_OF_FRAME:
        return frame, []

    first, last = rng
    msub = str(row.get("mutation_subclass", "")).strip().lower()
    analysis_mutation = _analysis_mutation_from_row(first, last, msub)
    if analysis_mutation is None:
        return frame, []

    affected_count = len(analysis_mutation.affected_exons)
    maximum_total_targets = max_additional_skips
    if analysis_mutation.mutation_class == AnalysisMutationClass.DUPLICATION:
        maximum_total_targets += affected_count

    analysis_candidates = find_frame_restoring_candidates(
        mutation=analysis_mutation,
        exon_coding_lengths=_coding_length_map(records),
        maximum_upstream_targets=max_additional_skips,
        maximum_downstream_targets=max_additional_skips,
        maximum_total_targets=maximum_total_targets,
    )

    # ``advanced`` is retained for API compatibility. Candidate generation is
    # always restricted to structured, mutation-adjacent blocks.
    _ = advanced
    candidates = [
        candidate
        for analysis_candidate in analysis_candidates
        if (
            candidate := _legacy_candidate_from_analysis(
                analysis_candidate,
                label,
                records,
            )
        )
    ]
    return frame, candidates[:MAX_DISPLAY_CANDIDATES]


@dataclass
class RestoredMutationRow:
    """One catalog mutation restored by a proposed skip target."""

    row: dict[str, Any]
    participant_id: str
    mutation_record_id: str
    description: str
    mutation_label: str
    mutation_class: str
    original_frame: str
    resulting_frame: str
    skip_combination_rank: int
    skipped_exons: list[int]
    candidate: SkipCandidate


def _coding_length_map(exons: list[ExonRecord]) -> dict[int, int]:
    return {exon.exon_number: exon.coding_length_bp for exon in exons}


def mutation_label(row: dict[str, Any]) -> str:
    """Describe a mutation without embedding the participant ID."""
    mclass = str(row.get("mutation_class", "")).strip()
    msub = str(row.get("mutation_subclass", "")).strip()
    start = str(row.get("start_region", "")).strip()
    stop = str(row.get("stop_region", "")).strip()
    region = start if start == stop else f"{start}-{stop}"
    parts = [part for part in (mclass, msub, region) if part]
    return " / ".join(parts) if parts else "Unknown mutation"


def _analysis_mutation_from_row(
    first: int,
    last: int,
    msub: str,
) -> Optional[AnalysisMutation]:
    if msub == "deletion":
        mutation_class = AnalysisMutationClass.DELETION
    elif msub == "duplication":
        mutation_class = AnalysisMutationClass.DUPLICATION
    else:
        return None

    return AnalysisMutation(
        mutation_class=mutation_class,
        first_exon=first,
        last_exon=last,
    )


def _legacy_candidate_from_analysis(
    analysis_candidate,
    label: str,
    records: list[ExonRecord],
) -> Optional[SkipCandidate]:
    """Adapt authoritative analysis output to the existing UI model."""
    mutation = analysis_candidate.mutation
    skipped = tuple(analysis_candidate.skipped_exons)
    affected = set(mutation.affected_exons)
    skipped_set = set(skipped)
    transcript_exons = set(range(1, REFERENCE.coding_exon_count + 1))

    if mutation.mutation_class == AnalysisMutationClass.DELETION:
        removed_from_reference = affected | skipped_set
        mutation_boundary = (mutation.first_exon - 1, mutation.last_exon + 1)
    else:
        # One duplicated copy is removed therapeutically; only extension
        # targets remove sequence from the canonical transcript.
        removed_from_reference = skipped_set - affected
        mutation_boundary = None

    def contiguous_blocks(exons: set[int]) -> list[tuple[int, int]]:
        if not exons:
            return []
        ordered = sorted(exons)
        blocks: list[tuple[int, int]] = []
        start = previous = ordered[0]
        for exon in ordered[1:]:
            if exon != previous + 1:
                blocks.append((start, previous))
                start = exon
            previous = exon
        blocks.append((start, previous))
        return blocks

    removed_blocks = contiguous_blocks(removed_from_reference)
    new_junctions = [
        (start - 1, end + 1)
        for start, end in removed_blocks
        if start > 1 and end < REFERENCE.coding_exon_count
    ]
    principal_junction = new_junctions[0] if new_junctions else None
    final_upstream, final_downstream = principal_junction or (0, 0)
    target_blocks = contiguous_blocks(skipped_set)
    total_coding_bp = sum(record.coding_length_bp for record in records)
    estimated_remaining = total_coding_bp + analysis_candidate.final_coding_delta

    return SkipCandidate(
        original_mutation=label,
        deleted_exons=list(mutation.affected_exons),
        additional_skipped_exons=list(skipped),
        all_removed_exons=sorted(removed_from_reference),
        retained_exons=sorted(transcript_exons - removed_from_reference),
        new_junctions=new_junctions,
        mutation_boundary_junction=mutation_boundary,
        final_upstream_exon=final_upstream,
        final_downstream_exon=final_downstream,
        total_coding_bases_removed=analysis_candidate.skipped_coding_bases,
        additional_coding_bases_removed=analysis_candidate.skipped_coding_bases,
        restores_frame=analysis_candidate.restores_frame,
        estimated_remaining_coding_bp=estimated_remaining,
        estimated_protein_aa=estimated_remaining // 3,
        evidence_class=SkipEvidence.COMPUTATIONAL,
        rank_score=(
            analysis_candidate.target_count,
            analysis_candidate.additional_exons_skipped,
            analysis_candidate.skipped_coding_bases,
            skipped,
        ),
        is_boundary_adjacent=True,
        is_structurally_valid=True,
        is_advanced_noncontiguous=False,
        skip_block_count=len(target_blocks),
        upstream_extension_count=sum(
            exon < mutation.first_exon for exon in skipped_set
        ),
        downstream_extension_count=sum(
            exon > mutation.last_exon for exon in skipped_set
        ),
        principal_junction=principal_junction,
        assumptions=[
            "Computational frame restoration only — not therapeutic evidence.",
            f"Mutation-specific strategy: {analysis_candidate.strategy}.",
            (
                f"Final coding delta: {analysis_candidate.final_coding_delta} bp "
                f"(remainder {analysis_candidate.frame_remainder})."
            ),
        ],
    )


def analyze_mutations_for_skip_target(
    catalog_df: pd.DataFrame,
    target_exons: list[int],
    exons: Optional[list[ExonRecord]] = None,
    *,
    include_in_frame: bool = False,
    class_filter: Optional[set[str]] = None,
    target_mode: str = "all",
) -> list[RestoredMutationRow]:
    """
    List catalog mutations for which skipping ``target_exons`` is part of a
    frame-restoring candidate strategy.

    Each mutation is evaluated independently. Target exons are a display
    filter applied only after valid candidates are calculated.
    """
    records = exons if exons is not None else build_exon_records()
    coding_lengths = _coding_length_map(records)
    active = sort_catalog_dataframe(active_catalog(catalog_df))
    results: list[RestoredMutationRow] = []

    for _, raw in active.iterrows():
        row = raw.to_dict()
        if not _row_matches_class_filter(row, class_filter):
            continue

        rng, frame = frame_result_for_catalog_row(row, records)
        if not rng:
            continue
        if frame.status == FrameStatus.IN_FRAME and not include_in_frame:
            continue
        if frame.status == FrameStatus.CANNOT_DETERMINE:
            continue

        first, last = rng
        msub = str(row.get("mutation_subclass", "")).strip().lower()
        label = mutation_description(row)
        analysis_mutation = _analysis_mutation_from_row(first, last, msub)
        if analysis_mutation is None:
            continue

        valid_candidates = find_frame_restoring_candidates(
            mutation=analysis_mutation,
            exon_coding_lengths=coding_lengths,
            maximum_upstream_targets=8,
            maximum_downstream_targets=8,
            maximum_total_targets=10,
        )

        displayed_candidates = filter_candidates_by_target_exons(
            candidates=valid_candidates,
            target_exons=target_exons,
            target_mode=target_mode,
        )
        candidate_ranks = {
            candidate: rank
            for rank, candidate in enumerate(valid_candidates, start=1)
        }

        for analysis_candidate in displayed_candidates:
            legacy_candidate = _legacy_candidate_from_analysis(
                analysis_candidate,
                label,
                records,
            )
            if not legacy_candidate or not legacy_candidate.restores_frame:
                continue

            results.append(
                RestoredMutationRow(
                    row=row,
                    participant_id=str(row.get("id", "")).strip(),
                    mutation_record_id=str(row.get("mutation_record_id", "")),
                    description=label,
                    mutation_label=mutation_label(row),
                    mutation_class=str(row.get("mutation_subclass", "")),
                    original_frame=frame.status.value.replace("_", " ").title(),
                    resulting_frame="In Frame",
                    skip_combination_rank=candidate_ranks[analysis_candidate],
                    skipped_exons=list(analysis_candidate.skipped_exons),
                    candidate=legacy_candidate,
                )
            )

    return results


def build_catalog_options(catalog_df: pd.DataFrame) -> list[tuple[str, str]]:
    """Build selectbox options (label, mutation_record_id) sorted by exon range."""
    active = sort_catalog_dataframe(active_catalog(catalog_df))
    options: list[tuple[str, str]] = []
    for _, raw in active.iterrows():
        row = raw.to_dict()
        rec_id = str(row.get("mutation_record_id", ""))
        options.append((mutation_description(row), rec_id))
    return options


def catalog_row_to_map_row(row: dict[str, Any]) -> dict[str, Any]:
    """Ensure catalog row has fields needed by the SVG map renderer."""
    return {k: str(v) if v is not None else "" for k, v in row.items()}


def skip_target_preview_row(target_exons: list[int]) -> dict[str, Any]:
    """Synthetic catalog row for target-only schematic (no patient mutation)."""
    if not target_exons:
        return {"id": "Skip target", "mutation_class": "", "mutation_subclass": "", "start_region": "", "stop_region": ""}
    first, last = target_exons[0], target_exons[-1]
    return {
        "id": "Skip target",
        "mutation_class": "Target",
        "mutation_subclass": "Skip",
        "start_region": f"e{first}",
        "stop_region": f"e{last}",
        "phenotype": "",
        "group": "",
    }


def overlay_from_row_and_candidate(
    row: dict[str, Any],
    candidate: Optional[SkipCandidate],
    *,
    all_candidates: Optional[list[SkipCandidate]] = None,
) -> SkipSchematicOverlay:
    """Build schematic overlay state for mutation + up to three skip candidates."""
    overlay = SkipSchematicOverlay()
    rng = exon_range_from_row(row)
    msub = str(row.get("mutation_subclass", "")).strip().lower()
    if rng and msub == "deletion":
        first, last = rng
        overlay.mutation_deleted_exons = set(range(first, last + 1))

    source = all_candidates if all_candidates is not None else ([candidate] if candidate else [])
    overlay.skip_candidate_exons = [
        set(c.additional_skipped_exons) for c in source[:3] if c
    ]

    active = candidate or (source[0] if source else None)
    if active and active.new_junctions:
        overlay.junction_upstream = active.new_junctions[0][0]
        overlay.junction_downstream = active.new_junctions[0][1]
    elif active and active.mutation_boundary_junction:
        overlay.junction_upstream = active.mutation_boundary_junction[0]
        overlay.junction_downstream = active.mutation_boundary_junction[1]
    return overlay


def catalog_row_variant_text(row: dict[str, Any]) -> str:
    """Convert a catalog row to sidebar-style variant text for ensure_analysis."""
    rng = exon_range_from_row(row)
    if not rng:
        return ""
    first, last = rng
    msub = str(row.get("mutation_subclass", "")).strip().lower()
    if msub == "deletion":
        return f"del{first}" if first == last else f"del{first}-{last}"
    if msub == "duplication":
        return f"dup{first}" if first == last else f"dup{first}-{last}"
    return ""


def frame_status_label(status: FrameStatus) -> str:
    return status.value.replace("_", " ").title()

