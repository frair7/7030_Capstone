#!/usr/bin/env python3
"""
Validate cached DMD reference data.

Usage (from cursor_ai/):
    python scripts/validate_reference_data.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config import REFERENCE  # noqa: E402
from src.reference_data import EXONS_CSV, load_exon_table  # noqa: E402


def validate_exons(rows: list[dict]) -> list[str]:
    """Return a list of validation error messages (empty if valid)."""
    errors: list[str] = []

    if len(rows) != REFERENCE.coding_exon_count:
        errors.append(
            f"Expected {REFERENCE.coding_exon_count} exons, found {len(rows)}"
        )

    exon_numbers = [int(r["exon_number"]) for r in rows]
    if len(set(exon_numbers)) != len(exon_numbers):
        errors.append("Duplicate exon_number values found")
    if sorted(exon_numbers) != list(range(1, REFERENCE.coding_exon_count + 1)):
        errors.append("Exon numbers are not exactly 1..79 in order")

    cumulative = 0
    prev_transcript_end = 0
    for row in rows:
        exon_num = int(row["exon_number"])
        g_start = int(row["genomic_start_grch38"])
        g_end = int(row["genomic_end_grch38"])
        if g_start > g_end:
            errors.append(f"Exon {exon_num}: genomic_start > genomic_end")
        if g_start < 1:
            errors.append(f"Exon {exon_num}: invalid genomic_start {g_start}")

        t_start = int(row["transcript_exon_start"])
        t_end = int(row["transcript_exon_end"])
        if t_start > t_end:
            errors.append(f"Exon {exon_num}: transcript start > end")
        if t_start != prev_transcript_end + 1 and exon_num > 1:
            errors.append(
                f"Exon {exon_num}: transcript coordinates not contiguous "
                f"(expected start {prev_transcript_end + 1}, got {t_start})"
            )
        prev_transcript_end = t_end
        expected_t_len = t_end - t_start + 1
        genomic_len = g_end - g_start + 1
        if expected_t_len != genomic_len:
            errors.append(f"Exon {exon_num}: transcript length != genomic length")

        coding_len = int(row["coding_length_bp"])
        if coding_len < 0:
            errors.append(f"Exon {exon_num}: negative coding_length_bp")

        if row["strand"] != REFERENCE.strand:
            errors.append(
                f"Exon {exon_num}: strand {row['strand']} != {REFERENCE.strand}"
            )
        if row["chromosome"] != REFERENCE.chromosome:
            errors.append(
                f"Exon {exon_num}: chromosome {row['chromosome']} != "
                f"{REFERENCE.chromosome}"
            )
        if not row.get("source") or not row.get("source_version"):
            errors.append(f"Exon {exon_num}: missing source or source_version")
        if not row.get("date_accessed"):
            errors.append(f"Exon {exon_num}: missing date_accessed")

        if coding_len > 0:
            cum_start = int(row["cumulative_cds_start"])
            cum_end = int(row["cumulative_cds_end"])
            if cum_start != cumulative + 1:
                errors.append(
                    f"Exon {exon_num}: cumulative_cds_start {cum_start} != "
                    f"expected {cumulative + 1}"
                )
            if cum_end - cum_start + 1 != coding_len:
                errors.append(f"Exon {exon_num}: cumulative CDS length mismatch")
            cumulative = cum_end

            phase_5 = row["splice_phase_5prime"]
            phase_3 = row["splice_phase_3prime"]
            if str(phase_5) not in {"0", "1", "2"}:
                errors.append(f"Exon {exon_num}: invalid splice_phase_5prime {phase_5}")
            if str(phase_3) not in {"0", "1", "2"}:
                errors.append(f"Exon {exon_num}: invalid splice_phase_3prime {phase_3}")
            if int(phase_5) != (cum_start - 1) % 3:
                errors.append(f"Exon {exon_num}: splice_phase_5prime inconsistent")
            if int(phase_3) != cum_end % 3:
                errors.append(f"Exon {exon_num}: splice_phase_3prime inconsistent")

    expected_cds = REFERENCE.protein_length_aa * 3 + 3
    if cumulative != expected_cds:
        errors.append(
            f"Total CDS {cumulative} bp != expected {expected_cds} bp "
            f"({REFERENCE.protein_length_aa} aa + stop)"
        )

    # Minus-strand: exon 1 should have the highest genomic coordinates
    if rows and rows[0]["strand"] == "reverse":
        e1 = rows[0]
        e79 = rows[-1]
        if int(e1["genomic_start_grch38"]) <= int(e79["genomic_end_grch38"]):
            errors.append(
                "Transcript order check failed: exon 1 should lie 5' of exon 79 "
                "on the minus strand (higher genomic coordinates for exon 1)"
            )

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate DMD reference data.")
    parser.add_argument(
        "--exons",
        type=Path,
        default=EXONS_CSV,
        help="Path to exon CSV",
    )
    args = parser.parse_args()

    print(f"Validating {args.exons} ...")
    rows = load_exon_table(args.exons)
    errors = validate_exons(rows)

    if errors:
        print(f"FAILED — {len(errors)} error(s):")
        for err in errors:
            print(f"  - {err}")
        return 1

    total_coding = sum(int(r["coding_length_bp"]) for r in rows)
    print("PASSED")
    print(f"  Exons: {len(rows)}")
    print(f"  Total coding bp: {total_coding}")
    print(f"  Genomic span: {rows[0]['genomic_start_grch38']} – "
          f"{rows[-1]['genomic_end_grch38']} (see note: minus strand)")
    print(f"  Source: {rows[0]['source']} ({rows[0]['source_version']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
