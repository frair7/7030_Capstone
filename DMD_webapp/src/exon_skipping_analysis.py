"""Mutation-specific exon-skipping candidate analysis.

This module determines whether an exon or exon combination restores the
DMD reading frame for a particular whole-exon deletion or duplication.

Important:
    Candidate generation and frame validation are mutation-specific.
    A selected target exon is only a display filter applied after valid
    candidates have been calculated.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from itertools import product
from typing import Iterable, Mapping, Sequence


class MutationClass(str, Enum):
    DELETION = "deletion"
    DUPLICATION = "duplication"


@dataclass(frozen=True)
class ExonRecord:
    """Coding annotation for one transcript exon."""

    exon_number: int
    coding_length: int

    def __post_init__(self) -> None:
        if self.exon_number < 1:
            raise ValueError("exon_number must be at least 1")

        if self.coding_length < 0:
            raise ValueError("coding_length cannot be negative")


@dataclass(frozen=True)
class ExonMutation:
    """A whole-exon deletion or tandem duplication."""

    mutation_class: MutationClass
    first_exon: int
    last_exon: int

    def __post_init__(self) -> None:
        if self.first_exon < 1:
            raise ValueError("first_exon must be at least 1")

        if self.last_exon < self.first_exon:
            raise ValueError("last_exon must be >= first_exon")

    @property
    def affected_exons(self) -> tuple[int, ...]:
        return tuple(range(self.first_exon, self.last_exon + 1))

    @property
    def label(self) -> str:
        if self.first_exon == self.last_exon:
            exon_text = str(self.first_exon)
        else:
            exon_text = f"{self.first_exon}-{self.last_exon}"

        return f"{self.mutation_class.value} {exon_text}"


@dataclass(frozen=True)
class SkipCandidate:
    """A validated exon-skipping candidate."""

    mutation: ExonMutation
    skipped_exons: tuple[int, ...]
    mutation_coding_delta: int
    skipped_coding_bases: int
    final_coding_delta: int
    restores_frame: bool
    target_count: int
    additional_exons_skipped: int
    strategy: str

    @property
    def skipped_exons_label(self) -> str:
        return ", ".join(str(exon) for exon in self.skipped_exons)

    @property
    def frame_remainder(self) -> int:
        return self.final_coding_delta % 3


def build_exon_length_map(
    exon_records: Iterable[ExonRecord],
) -> dict[int, int]:
    """Convert exon records into an exon-number-to-coding-length dictionary."""

    exon_lengths: dict[int, int] = {}

    for record in exon_records:
        if record.exon_number in exon_lengths:
            raise ValueError(
                f"Duplicate exon annotation for exon {record.exon_number}"
            )

        exon_lengths[record.exon_number] = record.coding_length

    if not exon_lengths:
        raise ValueError("No exon records were supplied")

    return exon_lengths


def coding_length_for_exons(
    exons: Iterable[int],
    exon_coding_lengths: Mapping[int, int],
) -> int:
    """Return the total coding length for the requested exons."""

    total = 0

    for exon_number in exons:
        if exon_number not in exon_coding_lengths:
            raise KeyError(
                f"Missing coding-length annotation for exon {exon_number}"
            )

        total += int(exon_coding_lengths[exon_number])

    return total


def mutation_coding_delta(
    mutation: ExonMutation,
    exon_coding_lengths: Mapping[int, int],
) -> int:
    """Return coding bases added or removed by the mutation.

    A deletion has a negative delta.
    A duplication has a positive delta.
    """

    affected_length = coding_length_for_exons(
        mutation.affected_exons,
        exon_coding_lengths,
    )

    if mutation.mutation_class == MutationClass.DELETION:
        return -affected_length

    if mutation.mutation_class == MutationClass.DUPLICATION:
        return affected_length

    raise ValueError(
        f"Unsupported mutation class: {mutation.mutation_class}"
    )


def candidate_is_allowed(
    mutation: ExonMutation,
    skipped_exons: Sequence[int],
) -> bool:
    """Return whether a skip set is structurally allowed for the mutation."""
    affected = set(mutation.affected_exons)
    skipped = set(int(exon) for exon in skipped_exons)

    if mutation.mutation_class == MutationClass.DELETION:
        return not affected.intersection(skipped)

    if mutation.mutation_class == MutationClass.DUPLICATION:
        return affected.issubset(skipped)

    return False


def evaluate_skip_set(
    mutation: ExonMutation,
    skipped_exons: Sequence[int],
    exon_coding_lengths: Mapping[int, int],
    strategy: str,
) -> SkipCandidate:
    """Evaluate whether a skip set restores the reading frame."""

    normalized_skip_set = tuple(sorted(set(int(exon) for exon in skipped_exons)))

    if not normalized_skip_set:
        raise ValueError("At least one skipped exon is required")

    transcript_exons = set(exon_coding_lengths)

    unknown_exons = set(normalized_skip_set) - transcript_exons
    if unknown_exons:
        raise ValueError(
            f"Unknown skipped exon(s): {sorted(unknown_exons)}"
        )

    affected_set = set(mutation.affected_exons)

    if not candidate_is_allowed(mutation, normalized_skip_set):
        if mutation.mutation_class == MutationClass.DELETION:
            invalid_targets = affected_set.intersection(normalized_skip_set)
            raise ValueError(
                "Deleted exons cannot also be exon-skipping targets: "
                f"{sorted(invalid_targets)}"
            )

        missing_targets = affected_set.difference(normalized_skip_set)
        raise ValueError(
            "Duplication candidates must include the duplicated interval; "
            f"missing exon(s): {sorted(missing_targets)}"
        )

    mutation_delta = mutation_coding_delta(
        mutation,
        exon_coding_lengths,
    )

    skipped_bases = coding_length_for_exons(
        normalized_skip_set,
        exon_coding_lengths,
    )

    # Skipping removes coding sequence, so it always contributes a
    # negative coding-length delta.
    final_delta = mutation_delta - skipped_bases
    restores_frame = final_delta % 3 == 0

    if mutation.mutation_class == MutationClass.DUPLICATION:
        additional_exons = len(
            set(normalized_skip_set) - affected_set
        )
    else:
        additional_exons = len(normalized_skip_set)

    return SkipCandidate(
        mutation=mutation,
        skipped_exons=normalized_skip_set,
        mutation_coding_delta=mutation_delta,
        skipped_coding_bases=skipped_bases,
        final_coding_delta=final_delta,
        restores_frame=restores_frame,
        target_count=len(normalized_skip_set),
        additional_exons_skipped=additional_exons,
        strategy=strategy,
    )


def _contiguous_upstream_blocks(
    boundary_exon: int,
    minimum_exon: int,
    maximum_block_size: int,
) -> list[tuple[int, ...]]:
    """Return contiguous upstream blocks ending at boundary_exon - 1."""

    blocks: list[tuple[int, ...]] = [tuple()]

    available = boundary_exon - minimum_exon
    max_size = min(maximum_block_size, max(available, 0))

    for size in range(1, max_size + 1):
        start = boundary_exon - size
        stop = boundary_exon
        blocks.append(tuple(range(start, stop)))

    return blocks


def _contiguous_downstream_blocks(
    boundary_exon: int,
    maximum_exon: int,
    maximum_block_size: int,
) -> list[tuple[int, ...]]:
    """Return contiguous downstream blocks starting at boundary_exon + 1."""

    blocks: list[tuple[int, ...]] = [tuple()]

    available = maximum_exon - boundary_exon
    max_size = min(maximum_block_size, max(available, 0))

    for size in range(1, max_size + 1):
        start = boundary_exon + 1
        stop = boundary_exon + size + 1
        blocks.append(tuple(range(start, stop)))

    return blocks


def generate_candidate_skip_sets(
    mutation: ExonMutation,
    exon_coding_lengths: Mapping[int, int],
    maximum_upstream_targets: int = 8,
    maximum_downstream_targets: int = 8,
    maximum_total_targets: int = 10,
) -> list[tuple[tuple[int, ...], str]]:
    """Generate biologically structured exon-skipping candidates.

    Deletions:
        The deleted exons are absent and cannot be targeted. Candidates are
        generated from contiguous upstream and downstream blocks bordering
        the deletion.

        Example for deletion exon 7:
            upstream block: 6
            downstream block: 8
            combined candidate: 6, 8

    Duplications:
        The duplicated interval is the base correction target. Contiguous
        upstream and downstream extensions are then considered.

        Example for duplication exons 3-4:
            base candidate: 3, 4
            extended candidate: 3, 4, 5
            extended candidate: 3, 4, 5, 6, 7, 8

    This function generates possibilities. Frame restoration is evaluated
    separately and invalid candidates are discarded.
    """

    all_exons = sorted(exon_coding_lengths)

    if not all_exons:
        return []

    minimum_exon = min(all_exons)
    maximum_exon = max(all_exons)

    upstream_blocks = _contiguous_upstream_blocks(
        boundary_exon=mutation.first_exon,
        minimum_exon=minimum_exon,
        maximum_block_size=maximum_upstream_targets,
    )

    downstream_blocks = _contiguous_downstream_blocks(
        boundary_exon=mutation.last_exon,
        maximum_exon=maximum_exon,
        maximum_block_size=maximum_downstream_targets,
    )

    if mutation.mutation_class == MutationClass.DUPLICATION:
        required_targets = mutation.affected_exons
    else:
        required_targets = tuple()

    generated: dict[tuple[int, ...], str] = {}

    for upstream, downstream in product(
        upstream_blocks,
        downstream_blocks,
    ):
        combined = tuple(
            sorted(
                set(required_targets)
                | set(upstream)
                | set(downstream)
            )
        )

        if not combined:
            continue

        if len(combined) > maximum_total_targets:
            continue

        if mutation.mutation_class == MutationClass.DELETION:
            if set(combined).intersection(mutation.affected_exons):
                continue

        if required_targets and combined == required_targets:
            strategy = "duplicated interval only"
        elif upstream and downstream:
            strategy = "upstream and downstream"
        elif upstream:
            strategy = "upstream"
        elif downstream:
            strategy = "downstream"
        else:
            strategy = "mutation interval"

        generated[combined] = strategy

    return sorted(
        generated.items(),
        key=lambda item: (
            len(item[0]),
            item[0],
        ),
    )


def find_frame_restoring_candidates(
    mutation: ExonMutation,
    exon_coding_lengths: Mapping[int, int],
    maximum_upstream_targets: int = 8,
    maximum_downstream_targets: int = 8,
    maximum_total_targets: int = 10,
) -> list[SkipCandidate]:
    """Return only candidates that restore the coding frame."""

    candidate_sets = generate_candidate_skip_sets(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
        maximum_upstream_targets=maximum_upstream_targets,
        maximum_downstream_targets=maximum_downstream_targets,
        maximum_total_targets=maximum_total_targets,
    )

    valid_candidates: list[SkipCandidate] = []

    for skipped_exons, strategy in candidate_sets:
        candidate = evaluate_skip_set(
            mutation=mutation,
            skipped_exons=skipped_exons,
            exon_coding_lengths=exon_coding_lengths,
            strategy=strategy,
        )

        if candidate.restores_frame:
            valid_candidates.append(candidate)

    valid_candidates.sort(
        key=lambda candidate: (
            candidate.target_count,
            candidate.additional_exons_skipped,
            candidate.skipped_coding_bases,
            candidate.skipped_exons,
        )
    )

    return valid_candidates


def filter_candidates_by_target_exon(
    candidates: Iterable[SkipCandidate],
    target_exon: int,
    target_mode: str,
) -> list[SkipCandidate]:
    """Filter valid candidates using the selected therapeutic target exon.

    This is a display filter only. It must never be used to determine
    whether a candidate restores the frame.
    """

    normalized_mode = target_mode.strip().lower()

    if normalized_mode not in {"single", "multi", "all"}:
        raise ValueError(
            "target_mode must be 'single', 'multi', or 'all'"
        )

    filtered: list[SkipCandidate] = []

    for candidate in candidates:
        if target_exon not in candidate.skipped_exons:
            continue

        if normalized_mode == "single" and candidate.target_count != 1:
            continue

        if normalized_mode == "multi" and candidate.target_count <= 1:
            continue

        filtered.append(candidate)

    return filtered


def filter_candidates_by_target_exons(
    candidates: Iterable[SkipCandidate],
    target_exons: list[int],
    target_mode: str,
) -> list[SkipCandidate]:
    """Filter candidates by one or more target exons and skipping strategy."""

    if not target_exons:
        return []

    normalized_mode = target_mode.strip().lower()
    target_set = set(target_exons)

    if len(target_set) == 1:
        return filter_candidates_by_target_exon(
            candidates,
            next(iter(target_set)),
            normalized_mode,
        )

    filtered: list[SkipCandidate] = []

    for candidate in candidates:
        if not target_set.issubset(set(candidate.skipped_exons)):
            continue

        if normalized_mode == "single" and candidate.target_count != 1:
            continue

        if normalized_mode == "multi" and candidate.target_count <= 1:
            continue

        filtered.append(candidate)

    return filtered


def candidate_to_dict(candidate: SkipCandidate) -> dict[str, object]:
    """Convert a candidate into a table-ready dictionary."""

    return {
        "mutation_class": candidate.mutation.mutation_class.value,
        "mutation": candidate.mutation.label,
        "first_exon": candidate.mutation.first_exon,
        "last_exon": candidate.mutation.last_exon,
        "skipped_exons": candidate.skipped_exons_label,
        "target_count": candidate.target_count,
        "strategy": candidate.strategy,
        "mutation_coding_delta": candidate.mutation_coding_delta,
        "skipped_coding_bases": candidate.skipped_coding_bases,
        "final_coding_delta": candidate.final_coding_delta,
        "frame_remainder": candidate.frame_remainder,
        "restores_frame": candidate.restores_frame,
    }
