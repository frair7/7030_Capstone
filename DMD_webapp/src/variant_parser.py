"""
Parse user variant inputs into structured ParsedVariant objects.

Supported modes:
  - Exon-level deletion / duplication (del45, del45-50, dup2, etc.)
  - Basic HGVS coding notation (partial support)
  - GRCh38 genomic coordinates on chromosome X
  - Direct exon selection lists
"""

from __future__ import annotations

import re
from typing import Optional

from src.config import REFERENCE
from src.models import InputMode, ParsedVariant, VariantType

EXON_COUNT = REFERENCE.coding_exon_count
REFSEQ_ID = REFERENCE.refseq_transcript

# del45, del45-50, deletion of exons 45 through 50
_EXON_DEL_PATTERNS = [
    re.compile(
        r"^del(?:etion)?\s*(?:of\s*)?(?:exons?\s*)?(\d+)(?:\s*[-–]\s*(\d+))?$",
        re.I,
    ),
    re.compile(
        r"^deletion\s+of\s+exons?\s+(\d+)\s+through\s+(\d+)$",
        re.I,
    ),
    re.compile(r"^(?:exon\s*)?(\d+)(?:\s*[-–]\s*(\d+))?\s*del(?:etion)?$", re.I),
    re.compile(r"^e(\d+)(?:\s*[-–]\s*e?(\d+))?\s*del$", re.I),
]

_EXON_DUP_PATTERNS = [
    re.compile(r"^dup(?:lication)?\s*(?:of\s*)?(?:exons?\s*)?(\d+)(?:\s*[-–]\s*(\d+))?$", re.I),
    re.compile(r"^(?:exon\s*)?(\d+)(?:\s*[-–]\s*(\d+))?\s*dup(?:lication)?$", re.I),
    re.compile(r"^e(\d+)(?:\s*[-–]\s*e?(\d+))?\s*dup$", re.I),
]

_GENOMIC_PATTERNS = [
    re.compile(
        r"^(?:chr)?([xX])(?::|-)?(\d+)(?:\s*[-–]\s*(\d+))?$"
    ),
]

# c.5287C>T, NM_004006.3:c.5287C>T
_HGVS_SUB = re.compile(
    r"^(?:(NM_\d+\.\d+):)?c\.(\d+)([ACGTUMRYSWKVBHDN])(?:>|=)([ACGTUMRYSWKVBHDN])$",
    re.I,
)
_HGVS_DEL = re.compile(
    r"^(?:(NM_\d+\.\d+):)?c\.(\d+)(?:_(\d+))?del$",
    re.I,
)
_HGVS_SPLICE = re.compile(
    r"^(?:(NM_\d+\.\d+):)?c\.(\d+)([+-])(\d+)([ACGTUMRYSWKVBHDN])(?:>|=)([ACGTUMRYSWKVBHDN])$",
    re.I,
)


def _validate_exon_range(first: int, last: int) -> list[str]:
    """Return validation errors for an exon range."""
    errors: list[str] = []
    if first < 1 or last > EXON_COUNT:
        errors.append(f"Exon numbers must be between 1 and {EXON_COUNT}.")
    if last < first:
        errors.append("End exon cannot be less than start exon.")
    return errors


def _parse_exon_bounds(match: re.Match) -> tuple[int, int]:
    first = int(match.group(1))
    last = int(match.group(2)) if match.group(2) else first
    return first, last


def parse_exon_variant(
    text: str,
    variant_type: VariantType,
    input_mode: InputMode,
) -> ParsedVariant:
    """Parse deletion or duplication shorthand."""
    raw = text.strip()
    patterns = _EXON_DEL_PATTERNS if variant_type == VariantType.DELETION else _EXON_DUP_PATTERNS
    for pattern in patterns:
        match = pattern.match(raw)
        if match:
            first, last = _parse_exon_bounds(match)
            errors = _validate_exon_range(first, last)
            return ParsedVariant(
                raw_input=raw,
                input_mode=input_mode,
                variant_type=variant_type,
                first_exon=first,
                last_exon=last,
                errors=errors,
            )
    return ParsedVariant(
        raw_input=raw,
        input_mode=input_mode,
        variant_type=VariantType.UNKNOWN,
        errors=[f"Could not parse exon {variant_type.value}: '{raw}'"],
    )


def parse_genomic_coordinate(text: str, assembly: str = "GRCh38") -> ParsedVariant:
    """Parse GRCh38 chromosome-X coordinate strings."""
    raw = text.strip()
    result = ParsedVariant(
        raw_input=raw,
        input_mode=InputMode.GENOMIC,
        variant_type=VariantType.UNKNOWN,
    )
    if assembly.upper() not in {"GRCH38", "HG38", "GRCH37", "HG19"}:
        result.errors.append(f"Unsupported genome assembly: {assembly}")
        return result
    if assembly.upper() in {"GRCH37", "HG19"}:
        result.errors.append(
            "GRCh37/HG19 coordinates are not supported. Please use GRCh38."
        )
        return result

    for pattern in _GENOMIC_PATTERNS:
        match = pattern.match(raw)
        if match:
            chrom = match.group(1).upper()
            start = int(match.group(2))
            end = int(match.group(3)) if match.group(3) else start
            if end < start:
                result.errors.append("Genomic end coordinate is less than start.")
                return result
            result.genomic_chromosome = chrom
            result.genomic_start = start
            result.genomic_end = end
            result.variant_type = (
                VariantType.DELETION if start != end else VariantType.SUBSTITUTION
            )
            return result

    result.errors.append(
        f"Could not parse genomic coordinate: '{raw}'. "
        "Expected formats: X:position, chrX:position, X:start-end"
    )
    return result


def parse_hgvs_coding(text: str) -> ParsedVariant:
    """Parse basic HGVS coding strings (limited first-version support)."""
    raw = text.strip()
    result = ParsedVariant(
        raw_input=raw,
        input_mode=InputMode.HGVS_CODING,
        variant_type=VariantType.UNKNOWN,
        hgvs_coding=raw,
    )

    for pattern in (_HGVS_SUB, _HGVS_SPLICE, _HGVS_DEL):
        match = pattern.match(raw)
        if not match:
            continue

        if pattern is _HGVS_SUB:
            tid, pos, _ref, _alt = match.groups()
            result.transcript_id = tid or REFSEQ_ID
            result.cds_position = int(pos)
            result.variant_type = VariantType.SUBSTITUTION
            result.warnings.append(
                "SNV frame effect requires sequence-level consequence; "
                "may be reported as cannot determine unless exon mapping succeeds."
            )
            return result

        if pattern is _HGVS_SPLICE:
            tid, exon_pos, intron_dir, intron_offset, _ref, _alt = match.groups()
            result.transcript_id = tid or REFSEQ_ID
            result.cds_position = int(exon_pos)
            result.variant_type = VariantType.SPLICE
            result.warnings.append(
                f"Splice variant at c.{exon_pos}{intron_dir}{intron_offset}; "
                "exon assignment requires reference mapping."
            )
            return result

        if pattern is _HGVS_DEL:
            tid, start, end = match.groups()
            result.transcript_id = tid or REFSEQ_ID
            result.cds_position = int(start)
            result.variant_type = VariantType.DELETION
            result.warnings.append(
                "HGVS deletion not yet fully normalized to exon boundaries; "
                "coordinate mapping may be uncertain."
            )
            return result

    result.errors.append(
        f"Unsupported or ambiguous HGVS coding string: '{raw}'. "
        "Supported examples: c.5287C>T, c.1812+1G>A, c.100_200del"
    )
    return result


def parse_exon_selection(exon_list: list[int]) -> ParsedVariant:
    """Parse a direct multiselect of exon numbers."""
    if not exon_list:
        return ParsedVariant(
            raw_input=str(exon_list),
            input_mode=InputMode.EXON_SELECTION,
            variant_type=VariantType.UNKNOWN,
            errors=["No exons selected."],
        )
    sorted_exons = sorted(set(exon_list))
    first, last = sorted_exons[0], sorted_exons[-1]
    errors = _validate_exon_range(first, last)
    warnings: list[str] = []
    if sorted_exons != list(range(first, last + 1)):
        warnings.append(
            "Non-contiguous exon selection; treating as deletion of listed exons."
        )
    return ParsedVariant(
        raw_input=",".join(str(e) for e in sorted_exons),
        input_mode=InputMode.EXON_SELECTION,
        variant_type=VariantType.DELETION,
        first_exon=first,
        last_exon=last,
        selected_exons=sorted_exons,
        warnings=warnings,
        errors=errors,
    )


def parse_variant(
    text: str,
    input_mode: InputMode,
    *,
    selected_exons: Optional[list[int]] = None,
    assembly: str = "GRCh38",
) -> ParsedVariant:
    """
    Dispatch parser by input mode.

    Parameters
    ----------
    text:
        Raw user input string (ignored for EXON_SELECTION if selected_exons given).
    input_mode:
        One of the InputMode enum values.
    selected_exons:
        Exon numbers for direct selection mode.
    assembly:
        Genome assembly label for genomic coordinate mode.
    """
    if input_mode == InputMode.EXON_DELETION:
        return parse_exon_variant(text, VariantType.DELETION, InputMode.EXON_DELETION)
    if input_mode == InputMode.EXON_DUPLICATION:
        return parse_exon_variant(
            text, VariantType.DUPLICATION, InputMode.EXON_DUPLICATION
        )
    if input_mode == InputMode.HGVS_CODING:
        return parse_hgvs_coding(text)
    if input_mode == InputMode.GENOMIC:
        return parse_genomic_coordinate(text, assembly=assembly)
    if input_mode == InputMode.EXON_SELECTION:
        return parse_exon_selection(selected_exons or [])
    return ParsedVariant(
        raw_input=text,
        input_mode=input_mode,
        variant_type=VariantType.UNKNOWN,
        errors=[f"Unsupported input mode: {input_mode}"],
    )
