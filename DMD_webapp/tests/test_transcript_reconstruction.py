"""Regression tests for explicit transcript reconstruction in exon skipping."""

from __future__ import annotations

import pytest

from src.coordinate_mapper import build_exon_records
from src.exon_skipping import find_skip_candidates
from src.transcript_reconstruction import reconstruct_deletion_skip


@pytest.fixture(scope="module")
def exons():
    return build_exon_records()


class TestDeletion4550Reconstruction:
  def test_skip_53_structure(self, exons) -> None:
      t = reconstruct_deletion_skip(45, 50, (53,), exons, require_frame=False)
      assert t is not None
      assert t.original_mutation_exons == [45, 46, 47, 48, 49, 50]
      assert t.additional_skipped_exons == [53]
      assert t.all_removed_exons == [45, 46, 47, 48, 49, 50, 53]
      assert 51 in t.retained_exons
      assert 52 in t.retained_exons
      assert (44, 51) in t.new_junctions
      assert (52, 54) in t.new_junctions
      assert (44, 54) not in t.new_junctions

  def test_skip_51_53_structure(self, exons) -> None:
      t = reconstruct_deletion_skip(45, 50, (51, 52, 53), exons, require_frame=False)
      assert t is not None
      assert t.all_removed_exons == list(range(45, 54))
      assert (44, 54) in t.new_junctions
      assert t.additional_skipped_exons == [51, 52, 53]

  def test_skip_51_structure(self, exons) -> None:
      t = reconstruct_deletion_skip(45, 50, (51,), exons, require_frame=False)
      assert t is not None
      assert (44, 52) in t.new_junctions
      assert t.additional_skipped_exons == [51]

  def test_no_interval_filling(self, exons) -> None:
      t = reconstruct_deletion_skip(45, 50, (53,), exons, require_frame=False)
      assert t is not None
      assert 51 in t.retained_exons
      assert 52 in t.retained_exons

  def test_coding_bases_only_removed_exons(self, exons) -> None:
      t = reconstruct_deletion_skip(45, 50, (53,), exons, require_frame=False)
      assert t is not None
      index = {e.exon_number: e for e in exons}
      expected = sum(index[e].coding_length_bp for e in t.all_removed_exons)
      assert t.total_coding_bases_removed == expected
      assert 51 not in t.all_removed_exons
      assert 52 not in t.all_removed_exons

  def test_isolated_skip_53_not_in_default_candidates(self, exons) -> None:
      candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3, advanced=False)
      assert all(c.additional_skipped_exons != [53] for c in candidates)

  def test_boundary_skip_51_in_candidates(self, exons) -> None:
      candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3, advanced=False)
      assert any(c.additional_skipped_exons == [51] for c in candidates)

  def test_candidate_fields_match_reconstruction(self, exons) -> None:
      candidates = find_skip_candidates(45, 50, exons, max_additional_skips=3, advanced=False)
      assert candidates
      best = candidates[0]
      t = reconstruct_deletion_skip(
          45, 50, tuple(best.additional_skipped_exons), exons, require_frame=False,
      )
      assert t is not None
      assert best.all_removed_exons == t.all_removed_exons
      assert best.retained_exons == t.retained_exons
      assert best.new_junctions == t.new_junctions
      assert best.total_coding_bases_removed == t.total_coding_bases_removed
