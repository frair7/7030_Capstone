"""
Data models for variants, mapping results, frame analysis, and skip candidates.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class InputMode(str, Enum):
    """Supported variant input modes."""

    EXON_DELETION = "exon_deletion"
    EXON_DUPLICATION = "exon_duplication"
    HGVS_CODING = "hgvs_coding"
    GENOMIC = "genomic"
    EXON_SELECTION = "exon_selection"


class VariantType(str, Enum):
    """High-level variant consequence category."""

    DELETION = "deletion"
    DUPLICATION = "duplication"
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    SPLICE = "splice"
    COMPLEX = "complex"
    UNKNOWN = "unknown"


class FrameStatus(str, Enum):
    """Reading-frame assessment outcome."""

    IN_FRAME = "in_frame"
    OUT_OF_FRAME = "out_of_frame"
    CANNOT_DETERMINE = "cannot_determine"


class SkipEvidence(str, Enum):
    """How a skip candidate is classified (non-clinical)."""

    COMPUTATIONAL = "computationally_restores_frame"
    UNSUPPORTED = "unsupported_theoretical_candidate"


@dataclass
class ExonRecord:
    """One row from the DMD exon reference table."""

    exon_number: int
    chromosome: str
    genomic_start: int
    genomic_end: int
    strand: str
    transcript_id: str
    ensembl_exon_id: str
    transcript_exon_start: int
    transcript_exon_end: int
    cds_start: Optional[int]
    cds_end: Optional[int]
    coding_length_bp: int
    splice_phase_5prime: Optional[int]
    splice_phase_3prime: Optional[int]
    cumulative_cds_start: Optional[int]
    cumulative_cds_end: Optional[int]

    @classmethod
    def from_row(cls, row: dict) -> ExonRecord:
        """Build from a CSV DictReader row."""
        def _opt_int(key: str) -> Optional[int]:
            val = row.get(key, "")
            return int(val) if val not in ("", None) else None

        return cls(
            exon_number=int(row["exon_number"]),
            chromosome=row["chromosome"],
            genomic_start=int(row["genomic_start_grch38"]),
            genomic_end=int(row["genomic_end_grch38"]),
            strand=row["strand"],
            transcript_id=row["transcript_id"],
            ensembl_exon_id=row["ensembl_exon_id"],
            transcript_exon_start=int(row["transcript_exon_start"]),
            transcript_exon_end=int(row["transcript_exon_end"]),
            cds_start=_opt_int("cds_start"),
            cds_end=_opt_int("cds_end"),
            coding_length_bp=int(row["coding_length_bp"]),
            splice_phase_5prime=_opt_int("splice_phase_5prime"),
            splice_phase_3prime=_opt_int("splice_phase_3prime"),
            cumulative_cds_start=_opt_int("cumulative_cds_start"),
            cumulative_cds_end=_opt_int("cumulative_cds_end"),
        )


@dataclass
class ParsedVariant:
    """Normalized variant representation after parsing user input."""

    raw_input: str
    input_mode: InputMode
    variant_type: VariantType
    first_exon: Optional[int] = None
    last_exon: Optional[int] = None
    transcript_id: Optional[str] = None
    hgvs_coding: Optional[str] = None
    genomic_chromosome: Optional[str] = None
    genomic_start: Optional[int] = None
    genomic_end: Optional[int] = None
    cds_position: Optional[int] = None
    selected_exons: list[int] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def exon_range(self) -> Optional[tuple[int, int]]:
        if self.first_exon is not None and self.last_exon is not None:
            return self.first_exon, self.last_exon
        return None


@dataclass
class LocationMapping:
    """Result of mapping a coordinate to the transcript model."""

    location_type: str  # exon, intron, upstream, downstream, boundary
    exon_number: Optional[int] = None
    upstream_exon: Optional[int] = None
    downstream_exon: Optional[int] = None
    genomic_position: Optional[int] = None
    genomic_start: Optional[int] = None
    genomic_end: Optional[int] = None
    message: str = ""
    warnings: list[str] = field(default_factory=list)


@dataclass
class FrameResult:
    """Reading-frame analysis for a whole-exon variant."""

    status: FrameStatus
    explanation: str
    upstream_exon: Optional[int] = None
    downstream_exon: Optional[int] = None
    upstream_phase_3prime: Optional[int] = None
    downstream_phase_5prime: Optional[int] = None
    junction_phase: Optional[int] = None
    coding_bases_removed: Optional[int] = None
    coding_bases_added: Optional[int] = None
    remaining_coding_bp: Optional[int] = None
    estimated_protein_aa: Optional[int] = None
    assumptions: list[str] = field(default_factory=list)


@dataclass
class SkipCandidate:
    """A computational exon-skipping strategy."""

    original_mutation: str
    deleted_exons: list[int]
    additional_skipped_exons: list[int]
    final_upstream_exon: int
    final_downstream_exon: int
    total_coding_bases_removed: int
    additional_coding_bases_removed: int
    restores_frame: bool
    estimated_remaining_coding_bp: int
    estimated_protein_aa: int
    evidence_class: SkipEvidence
    rank_score: tuple[int, int, int]  # fewer skips, fewer bp, contiguous bonus
    assumptions: list[str] = field(default_factory=list)
    affected_domains: list[str] = field(default_factory=list)
