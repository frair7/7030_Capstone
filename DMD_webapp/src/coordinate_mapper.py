"""
Map genomic coordinates and CDS positions to DMD exon/intron locations.
"""

from __future__ import annotations

from typing import Optional

from src.config import REFERENCE
from src.models import ExonRecord, LocationMapping, ParsedVariant
from src.reference_data import exon_by_number, genomic_span, load_exon_table


def build_exon_records(rows: Optional[list[dict]] = None) -> list[ExonRecord]:
    """Load and convert exon table rows to ExonRecord objects."""
    data = rows if rows is not None else load_exon_table()
    return [ExonRecord.from_row(row) for row in data]


def map_genomic_position(
    position: int,
    exons: list[ExonRecord],
    chromosome: str = REFERENCE.chromosome,
) -> LocationMapping:
    """
    Map a single GRCh38 genomic position on chromosome X.

    For minus-strand DMD, transcript order runs opposite to ascending
    genomic coordinates.
    """
    gene_low, gene_high = genomic_span(
        [{"genomic_start_grch38": e.genomic_start, "genomic_end_grch38": e.genomic_end}
         for e in exons]
    )

    if position < gene_low:
        return LocationMapping(
            location_type="downstream",
            genomic_position=position,
            message=(
                f"Position {position} is downstream of the DMD gene body "
                f"(below {gene_low} on chr{chromosome})."
            ),
        )
    if position > gene_high:
        return LocationMapping(
            location_type="upstream",
            genomic_position=position,
            message=(
                f"Position {position} is upstream of the DMD gene body "
                f"(above {gene_high} on chr{chromosome})."
            ),
        )

    indexed = exon_by_number(
        [{"exon_number": e.exon_number,
          "genomic_start_grch38": e.genomic_start,
          "genomic_end_grch38": e.genomic_end}
         for e in exons]
    )

    for exon_num in range(1, len(exons) + 1):
        row = indexed[exon_num]
        g_start = int(row["genomic_start_grch38"])
        g_end = int(row["genomic_end_grch38"])
        if g_start <= position <= g_end:
            boundary = ""
            if position == g_start or position == g_end:
                boundary = " (exon boundary)"
            return LocationMapping(
                location_type="exon",
                exon_number=exon_num,
                genomic_position=position,
                message=f"Position {position} lies within exon {exon_num}{boundary}.",
            )

        if exon_num < len(exons):
            next_row = indexed[exon_num + 1]
            next_start = int(next_row["genomic_start_grch38"])
            next_end = int(next_row["genomic_end_grch38"])
            # Minus strand: intron between exon N and N+1 lies between
            # the lower boundary of exon N and the higher boundary of exon N+1.
            intron_low = min(g_end, next_end) + 1
            intron_high = max(g_start, next_start) - 1
            if intron_low <= position <= intron_high:
                return LocationMapping(
                    location_type="intron",
                    upstream_exon=exon_num,
                    downstream_exon=exon_num + 1,
                    genomic_position=position,
                    message=(
                        f"Position {position} lies in the intron between "
                        f"exons {exon_num} and {exon_num + 1}."
                    ),
                )

    return LocationMapping(
        location_type="unknown",
        genomic_position=position,
        message=f"Could not map position {position} within the DMD reference model.",
        warnings=["Ambiguous mapping — position within gene span but not resolved."],
    )


def map_genomic_interval(
    start: int,
    end: int,
    exons: list[ExonRecord],
) -> LocationMapping:
    """Map a genomic interval; summarize spanning exons/introns."""
    if end < start:
        return LocationMapping(
            location_type="error",
            genomic_start=start,
            genomic_end=end,
            message="Invalid interval: end < start.",
        )

    start_map = map_genomic_position(start, exons)
    end_map = map_genomic_position(end, exons)

    if start_map.location_type == end_map.location_type == "exon":
        if start_map.exon_number == end_map.exon_number:
            return LocationMapping(
                location_type="exon",
                exon_number=start_map.exon_number,
                genomic_start=start,
                genomic_end=end,
                message=(
                    f"Interval {start}-{end} lies within exon {start_map.exon_number}."
                ),
            )
        return LocationMapping(
            location_type="multi_exon",
            genomic_start=start,
            genomic_end=end,
            upstream_exon=start_map.exon_number,
            downstream_exon=end_map.exon_number,
            message=(
                f"Interval {start}-{end} spans exons "
                f"{start_map.exon_number} through {end_map.exon_number}."
            ),
            warnings=["Interval may include partial exons and introns."],
        )

    return LocationMapping(
        location_type="complex",
        genomic_start=start,
        genomic_end=end,
        message=(
            f"Interval {start}-{end}: start={start_map.location_type}, "
            f"end={end_map.location_type}."
        ),
        warnings=["Partial exon or intronic span — interpret with caution."],
    )


def map_cds_position(
    cds_position: int,
    exons: list[ExonRecord],
) -> LocationMapping:
    """Map a 1-based CDS nucleotide position to an exon."""
    for exon in exons:
        if exon.cds_start is None or exon.cds_end is None:
            continue
        if exon.cds_start <= cds_position <= exon.cds_end:
            return LocationMapping(
                location_type="exon",
                exon_number=exon.exon_number,
                message=(
                    f"CDS position c.{cds_position} lies within exon "
                    f"{exon.exon_number}."
                ),
            )
    return LocationMapping(
        location_type="unknown",
        message=f"CDS position {cds_position} is outside the coding region.",
        warnings=["Cannot assign exon for this CDS coordinate."],
    )


def map_parsed_variant(
    variant: ParsedVariant,
    exons: Optional[list[ExonRecord]] = None,
) -> LocationMapping:
    """Map a parsed variant to a LocationMapping when possible."""
    records = exons if exons is not None else build_exon_records()

    if variant.genomic_start is not None:
        if variant.genomic_end is not None and variant.genomic_end != variant.genomic_start:
            return map_genomic_interval(variant.genomic_start, variant.genomic_end, records)
        return map_genomic_position(variant.genomic_start, records)

    if variant.cds_position is not None:
        return map_cds_position(variant.cds_position, records)

    if variant.first_exon is not None:
        return LocationMapping(
            location_type="exon",
            exon_number=variant.first_exon,
            upstream_exon=variant.first_exon,
            downstream_exon=variant.last_exon,
            message=(
                f"Variant affects exons {variant.first_exon}"
                f"{f'-{variant.last_exon}' if variant.last_exon != variant.first_exon else ''}."
            ),
        )

    return LocationMapping(
        location_type="unknown",
        message="No coordinate mapping available for this variant.",
        warnings=variant.warnings,
    )
