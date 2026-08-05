"""Tests for mutation-specific exon-skipping candidate analysis."""

from __future__ import annotations

import pytest

from src.exon_skipping_analysis import (
    ExonMutation,
    MutationClass,
    candidate_is_allowed,
    evaluate_skip_set,
    filter_candidates_by_target_exon,
    filter_candidates_by_target_exons,
    find_frame_restoring_candidates,
)


def candidate_sets(candidates):
    return {
        candidate.skipped_exons
        for candidate in candidates
    }


def test_deletion_45_single_exon_44_restores_frame(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(
        mutation_class=MutationClass.DELETION,
        first_exon=45,
        last_exon=45,
    )

    candidate = evaluate_skip_set(
        mutation=mutation,
        skipped_exons=(44,),
        exon_coding_lengths=exon_coding_lengths,
        strategy="upstream",
    )

    assert candidate.restores_frame is True
    assert candidate.frame_remainder == 0


def test_deletion_7_expected_multi_exon_candidates(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(
        mutation_class=MutationClass.DELETION,
        first_exon=7,
        last_exon=7,
    )

    candidates = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
        maximum_upstream_targets=6,
        maximum_downstream_targets=3,
        maximum_total_targets=6,
    )

    sets = candidate_sets(candidates)

    assert (6, 8) in sets
    assert (6, 8, 9) in sets
    assert (5, 6, 8) in sets
    assert (2, 3, 4, 5, 6) in sets

    assert all(
        7 not in candidate.skipped_exons
        for candidate in candidates
    )

    assert not any(
        candidate.target_count == 1
        for candidate in candidates
    )


def test_duplication_3_4_expected_candidates(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(
        mutation_class=MutationClass.DUPLICATION,
        first_exon=3,
        last_exon=4,
    )

    candidates = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
        maximum_upstream_targets=2,
        maximum_downstream_targets=6,
        maximum_total_targets=8,
    )

    sets = candidate_sets(candidates)

    assert (3, 4) in sets
    assert (3, 4, 5) in sets
    assert (3, 4, 5, 6, 7, 8) in sets

    assert not any(
        candidate.target_count == 1
        for candidate in candidates
    )


def test_duplication_2_is_corrected_by_targeting_exon_2(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(
        mutation_class=MutationClass.DUPLICATION,
        first_exon=2,
        last_exon=2,
    )

    candidate = evaluate_skip_set(
        mutation=mutation,
        skipped_exons=(2,),
        exon_coding_lengths=exon_coding_lengths,
        strategy="duplicated interval only",
    )

    assert candidate.restores_frame is True


def test_unrelated_deletion_does_not_appear_for_exon_45(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(
        mutation_class=MutationClass.DELETION,
        first_exon=3,
        last_exon=7,
    )

    candidates = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
    )

    exon_45_candidates = [
        candidate
        for candidate in candidates
        if 45 in candidate.skipped_exons
    ]

    assert exon_45_candidates == []


def test_every_candidate_restores_frame(
    exon_coding_lengths,
) -> None:
    mutations = [
        ExonMutation(MutationClass.DELETION, 7, 7),
        ExonMutation(MutationClass.DELETION, 45, 45),
        ExonMutation(MutationClass.DUPLICATION, 2, 2),
        ExonMutation(MutationClass.DUPLICATION, 3, 4),
    ]

    for mutation in mutations:
        candidates = find_frame_restoring_candidates(
            mutation=mutation,
            exon_coding_lengths=exon_coding_lengths,
        )
        for candidate in candidates:
            assert candidate.final_coding_delta % 3 == 0
            assert candidate.restores_frame is True


def test_no_deleted_exon_in_skip_targets(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(MutationClass.DELETION, 7, 7)
    candidates = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
    )

    for candidate in candidates:
        assert 7 not in candidate.skipped_exons


def test_duplication_candidates_include_duplicated_interval(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(MutationClass.DUPLICATION, 3, 4)
    candidates = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
    )

    for candidate in candidates:
        assert set(mutation.affected_exons).issubset(set(candidate.skipped_exons))


def test_duplication_rejects_downstream_only_skip_set(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(MutationClass.DUPLICATION, 3, 7)

    assert candidate_is_allowed(mutation, (8, 9)) is False
    with pytest.raises(
        ValueError,
        match="must include the duplicated interval",
    ):
        evaluate_skip_set(
            mutation=mutation,
            skipped_exons=(8, 9),
            exon_coding_lengths=exon_coding_lengths,
            strategy="downstream",
        )


def test_deletion_3_7_skip_8_9_uses_coding_delta(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(MutationClass.DELETION, 3, 7)
    candidate = evaluate_skip_set(
        mutation=mutation,
        skipped_exons=(8, 9),
        exon_coding_lengths=exon_coding_lengths,
        strategy="downstream",
    )

    affected_bases = sum(exon_coding_lengths[exon] for exon in range(3, 8))
    skipped_bases = exon_coding_lengths[8] + exon_coding_lengths[9]
    assert candidate.final_coding_delta == -affected_bases - skipped_bases
    assert candidate.frame_remainder == 0
    assert candidate.restores_frame is True


def test_target_filter_never_creates_candidates(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(MutationClass.DELETION, 45, 45)
    candidates = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
    )

    filtered = filter_candidates_by_target_exon(
        candidates,
        target_exon=44,
        target_mode="single",
    )

    assert len(filtered) <= len(candidates)
    assert all(c.target_count == 1 for c in filtered)
    assert all(44 in c.skipped_exons for c in filtered)


def test_target_exon_45_does_not_surface_unrelated_mutations(
    exon_coding_lengths,
) -> None:
    unrelated = ExonMutation(MutationClass.DELETION, 3, 7)
    candidates = find_frame_restoring_candidates(
        mutation=unrelated,
        exon_coding_lengths=exon_coding_lengths,
    )

    filtered = filter_candidates_by_target_exons(
        candidates,
        target_exons=[45],
        target_mode="all",
    )

    assert filtered == []


def test_candidate_order_is_deterministic(
    exon_coding_lengths,
) -> None:
    mutation = ExonMutation(MutationClass.DELETION, 7, 7)

    first = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
    )
    second = find_frame_restoring_candidates(
        mutation=mutation,
        exon_coding_lengths=exon_coding_lengths,
    )

    assert first == second
