"""
Explicit transcript reconstruction for exon-skipping candidate analysis.

Never collapses noncontiguous removed exons into a continuous interval.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from src.config import REFERENCE
from src.coordinate_mapper import build_exon_records
from src.frame_analysis import coding_bp_for_exon_range, total_coding_bp
from src.models import ExonRecord


@dataclass
class ReconstructedTranscript:
    """Ordered retained exons and junctions after mutation + additional skips."""

    original_mutation_exons: list[int]
    additional_skipped_exons: list[int]
    all_removed_exons: list[int]
    retained_exons: list[int]
    new_junctions: list[tuple[int, int]]
    mutation_boundary_junction: Optional[tuple[int, int]]
    principal_junction: Optional[tuple[int, int]]
    total_coding_bases_removed: int
    additional_coding_bases_removed: int
    estimated_remaining_coding_bp: int
    estimated_protein_aa: int
    restores_frame: bool
    is_boundary_adjacent: bool
    is_structurally_valid: bool = True


def _contiguous_blocks(exons: list[int]) -> list[list[int]]:
    if not exons:
        return []
    blocks: list[list[int]] = [[exons[0]]]
    for exon in exons[1:]:
        if exon == blocks[-1][-1] + 1:
            blocks[-1].append(exon)
        else:
            blocks.append([exon])
    return blocks


def _extension_counts(
    mutation_first: int,
    mutation_last: int,
    additional: set[int],
) -> tuple[int, int]:
    upstream = sum(1 for e in additional if e < mutation_first)
    downstream = sum(1 for e in additional if e > mutation_last)
    return upstream, downstream


def generate_boundary_extension_candidates(
    mutation_first: int,
    mutation_last: int,
    *,
    max_additional_skips: int,
) -> list[tuple[int, ...]]:
    """
    Sequential upstream, downstream, and combined boundary-extension strategies.

    Generates combinations in increasing order of additional exon count.
    """
    n = REFERENCE.coding_exon_count
    seen: set[tuple[int, ...]] = set()
    combos: list[tuple[int, ...]] = []

    def add(combo: tuple[int, ...]) -> None:
        if not combo or len(combo) > max_additional_skips:
            return
        if any(e < 1 or e > n for e in combo):
            return
        if combo not in seen:
            seen.add(combo)
            combos.append(combo)

    for total in range(1, max_additional_skips + 1):
        for down_len in range(0, total + 1):
            up_len = total - down_len
            if down_len > 0:
                downstream = tuple(range(mutation_last + 1, mutation_last + 1 + down_len))
                if downstream[-1] <= n:
                    add(downstream)
            if up_len > 0:
                upstream = tuple(range(mutation_first - up_len, mutation_first))
                if upstream[0] >= 1:
                    add(upstream)
            if down_len > 0 and up_len > 0:
                downstream = tuple(range(mutation_last + 1, mutation_last + 1 + down_len))
                upstream = tuple(range(mutation_first - up_len, mutation_first))
                if downstream[-1] <= n and upstream[0] >= 1:
                    combined = tuple(sorted(set(upstream) | set(downstream)))
                    if len(combined) == total:
                        add(combined)

    return combos


def boundary_adjacent_skip_combos(
    mutation_first: int,
    mutation_last: int,
    *,
    max_additional_skips: int,
) -> list[tuple[int, ...]]:
    """Backward-compatible alias for boundary extension candidate generation."""
    return generate_boundary_extension_candidates(
        mutation_first, mutation_last, max_additional_skips=max_additional_skips,
    )


def _novel_junctions(retained: list[int]) -> list[tuple[int, int]]:
    """Junctions where consecutive retained exons are not adjacent in the reference."""
    junctions: list[tuple[int, int]] = []
    for i in range(len(retained) - 1):
        upstream, downstream = retained[i], retained[i + 1]
        if downstream > upstream + 1:
            junctions.append((upstream, downstream))
    return junctions


def _junction_in_frame(
    upstream: int,
    downstream: int,
    removed: set[int],
    index: dict[int, ExonRecord],
) -> bool:
    """True when splice phases match across a novel junction spanning removed exons."""
    if upstream < 1 or downstream > REFERENCE.coding_exon_count:
        return False
    up_exon = index.get(upstream)
    down_exon = index.get(downstream)
    if not up_exon or not down_exon:
        return False
    if up_exon.splice_phase_3prime is None or down_exon.splice_phase_5prime is None:
        return False
    between = set(range(upstream + 1, downstream))
    if between - removed:
        return False
    removed_bp = sum(index[e].coding_length_bp for e in between if e in index)
    junction_phase = (up_exon.splice_phase_3prime + removed_bp) % 3
    return junction_phase == down_exon.splice_phase_5prime


def _is_boundary_adjacent_deletion_skip(
    mutation_first: int,
    mutation_last: int,
    additional: set[int],
) -> bool:
    """Additional skips form one contiguous block touching a deletion boundary."""
    if not additional:
        return True
    ordered = sorted(additional)
    if ordered != list(range(ordered[0], ordered[-1] + 1)):
        return False
    touches_downstream = ordered[0] == mutation_last + 1
    touches_upstream = ordered[-1] == mutation_first - 1
    return touches_downstream or touches_upstream


def _principal_junction(
    mutation_first: int,
    mutation_last: int,
    new_junctions: list[tuple[int, int]],
) -> Optional[tuple[int, int]]:
    """Single extension junction for continuous boundary-adjacent strategies."""
    boundary_up = mutation_first - 1
    for upstream, downstream in new_junctions:
        if upstream == boundary_up and downstream > mutation_last + 1:
            return (upstream, downstream)
    if len(new_junctions) == 1:
        return new_junctions[0]
    return None


def reconstruct_deletion_skip(
    mutation_first: int,
    mutation_last: int,
    additional_skips: tuple[int, ...],
    exons: Optional[list[ExonRecord]] = None,
    *,
    require_frame: bool = True,
) -> Optional[ReconstructedTranscript]:
    """
    Reconstruct transcript after a whole-exon deletion plus additional skipped exons.
    """
    records = exons if exons is not None else build_exon_records()
    index = {e.exon_number: e for e in records}
    n = REFERENCE.coding_exon_count

    if mutation_first < 1 or mutation_last > n or mutation_last < mutation_first:
        return None

    original_deleted = set(range(mutation_first, mutation_last + 1))
    additional = set(additional_skips)
    if additional & original_deleted:
        return None

    all_removed_set = original_deleted | additional
    retained = [e for e in range(1, n + 1) if e not in all_removed_set]
    if not retained:
        return None

    new_junctions = _novel_junctions(retained)
    mutation_boundary: Optional[tuple[int, int]] = None
    if mutation_first > 1 and mutation_last < n:
        mutation_boundary = (mutation_first - 1, mutation_last + 1)

    restores = True
    if all_removed_set:
        total_removed_bp = sum(index[e].coding_length_bp for e in all_removed_set)
        if total_removed_bp % 3 != 0:
            restores = False
    else:
        total_removed_bp = 0

    for junction in new_junctions:
        if not _junction_in_frame(junction[0], junction[1], all_removed_set, index):
            restores = False
            break

    additional_bp = sum(index[e].coding_length_bp for e in additional if e in index)
    total_bp = total_coding_bp(records)
    remaining = total_bp - total_removed_bp

    is_adjacent = _is_boundary_adjacent_deletion_skip(
        mutation_first, mutation_last, additional,
    )
    principal = _principal_junction(mutation_first, mutation_last, new_junctions)
    structurally_valid = bool(retained)

    if require_frame and not restores:
        return None

    return ReconstructedTranscript(
        original_mutation_exons=sorted(original_deleted),
        additional_skipped_exons=sorted(additional),
        all_removed_exons=sorted(all_removed_set),
        retained_exons=retained,
        new_junctions=new_junctions,
        mutation_boundary_junction=mutation_boundary,
        principal_junction=principal,
        total_coding_bases_removed=total_removed_bp,
        additional_coding_bases_removed=additional_bp,
        estimated_remaining_coding_bp=remaining,
        estimated_protein_aa=remaining // 3,
        restores_frame=restores,
        is_boundary_adjacent=is_adjacent,
        is_structurally_valid=structurally_valid,
    )


@dataclass
class StrategyEvaluation:
    """Result of evaluating one internally tested skip strategy."""

    additional_skipped_exons: tuple[int, ...]
    restores_frame: bool
    is_structurally_valid: bool
    retained: bool
    new_junctions: list[tuple[int, int]]
    principal_junction: Optional[tuple[int, int]]


def evaluate_all_strategies(
    mutation_first: int,
    mutation_last: int,
    exons: Optional[list[ExonRecord]] = None,
    *,
    max_additional_skips: int = 3,
    advanced: bool = False,
) -> list[StrategyEvaluation]:
    """Evaluate every generated strategy (for regression reporting)."""
    records = exons if exons is not None else build_exon_records()
    combos = list(generate_boundary_extension_candidates(
        mutation_first, mutation_last, max_additional_skips=max_additional_skips,
    ))
    if advanced:
        import itertools
        retained = [
            e.exon_number for e in records
            if e.exon_number < mutation_first or e.exon_number > mutation_last
        ]
        base = set(combos)
        for size in range(1, max_additional_skips + 1):
            for combo in itertools.combinations(retained, size):
                if combo not in base:
                    combos.append(combo)

    results: list[StrategyEvaluation] = []
    for combo in combos:
        transcript = reconstruct_deletion_skip(
            mutation_first, mutation_last, combo, records, require_frame=False,
        )
        if transcript is None:
            results.append(StrategyEvaluation(
                additional_skipped_exons=combo,
                restores_frame=False,
                is_structurally_valid=False,
                retained=False,
                new_junctions=[],
                principal_junction=None,
            ))
        else:
            results.append(StrategyEvaluation(
                additional_skipped_exons=combo,
                restores_frame=transcript.restores_frame,
                is_structurally_valid=transcript.is_structurally_valid,
                retained=True,
                new_junctions=transcript.new_junctions,
                principal_junction=transcript.principal_junction,
            ))
    return results
