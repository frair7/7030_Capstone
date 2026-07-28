#!/usr/bin/env python3
"""
Fetch authoritative DMD Dp427m exon reference data and cache locally.

Primary source: Ensembl REST API (GRCh38), transcript ENST00000357033.
Cross-check metadata: NCBI RefSeq NM_004006.3.

Usage (from cursor_ai/):
    python scripts/fetch_reference_data.py
    python scripts/fetch_reference_data.py --output data/dmd_exons_grch38.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date
from pathlib import Path
from typing import Any

import requests

# Allow running as a script from cursor_ai/
ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
REFERENCE_DIR = REPO_ROOT / "reference_tables"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config import REFERENCE  # noqa: E402

ENSEMBL_BASE = "https://rest.ensembl.org"
CACHE_DIR = ROOT / "data" / "cache"
DEFAULT_OUTPUT = REFERENCE_DIR / "dmd_exons_grch38.csv"
DOMAINS_OUTPUT = ROOT / "data" / "dmd_domains.csv"
EDGE_SHAPES_CSV = ROOT / "data" / "dp427m_exon_edge_shapes.csv"

CSV_COLUMNS = [
    "exon_number",
    "five_prime_shape",
    "three_prime_shape",
    "chromosome",
    "genomic_start_grch38",
    "genomic_end_grch38",
    "strand",
    "transcript_id",
    "ensembl_exon_id",
    "transcript_exon_start",
    "transcript_exon_end",
    "cds_start",
    "cds_end",
    "coding_length_bp",
    "splice_phase_5prime",
    "splice_phase_3prime",
    "cumulative_cds_start",
    "cumulative_cds_end",
    "source",
    "source_version",
    "date_accessed",
]


def _ensembl_get(path: str, params: dict[str, Any] | None = None) -> Any:
    """GET request to Ensembl REST API with basic error handling."""
    url = f"{ENSEMBL_BASE}{path}"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers, params=params or {}, timeout=120)
    response.raise_for_status()
    return response.json()


def fetch_ensembl_transcript(ensembl_id: str) -> dict[str, Any]:
    """Retrieve expanded transcript record from Ensembl."""
    ensembl_base_id = ensembl_id.split(".")[0]
    return _ensembl_get(f"/lookup/id/{ensembl_base_id}", params={"expand": 1})


def genomic_positions_in_transcript_order(
    genomic_start: int, genomic_end: int, strand: int
) -> list[int]:
    """Return genomic coordinates visited 5'→3' along the transcript."""
    if strand == 1:
        return list(range(genomic_start, genomic_end + 1))
    return list(range(genomic_end, genomic_start - 1, -1))


def build_exon_table(transcript: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Build one row per exon in transcript order (exon 1 = 5' end).

    Splice phases are computed from cumulative coding length mod 3.
    """
    translation = transcript["Translation"]
    cds_low = min(translation["start"], translation["end"])
    cds_high = max(translation["start"], translation["end"])
    strand = transcript["strand"]
    strand_label = "reverse" if strand == -1 else "forward"
    chromosome = transcript["seq_region_name"]
    transcript_id = f"{transcript['id']}.{transcript['version']}"
    ensembl_version = str(transcript["version"])
    today = date.today().isoformat()

    exons = transcript["Exon"]
    if len(exons) != REFERENCE.coding_exon_count:
        raise ValueError(
            f"Expected {REFERENCE.coding_exon_count} exons, got {len(exons)} "
            f"from Ensembl {transcript_id}"
        )

    rows: list[dict[str, Any]] = []
    transcript_pos = 0
    cumulative_cds = 0

    for exon_number, exon in enumerate(exons, start=1):
        g_start = int(exon["start"])
        g_end = int(exon["end"])
        exon_length = g_end - g_start + 1
        transcript_exon_start = transcript_pos + 1
        transcript_exon_end = transcript_pos + exon_length
        transcript_pos += exon_length

        coding_positions = [
            pos
            for pos in genomic_positions_in_transcript_order(g_start, g_end, strand)
            if cds_low <= pos <= cds_high
        ]
        coding_length = len(coding_positions)

        if coding_length > 0:
            splice_phase_5prime = cumulative_cds % 3
            cumulative_cds_start = cumulative_cds + 1
            cumulative_cds += coding_length
            cumulative_cds_end = cumulative_cds
            splice_phase_3prime = cumulative_cds % 3
            cds_start = cumulative_cds_start
            cds_end = cumulative_cds_end
        else:
            splice_phase_5prime = ""
            splice_phase_3prime = ""
            cumulative_cds_start = ""
            cumulative_cds_end = ""
            cds_start = ""
            cds_end = ""

        rows.append(
            {
                "exon_number": exon_number,
                "chromosome": chromosome,
                "genomic_start_grch38": g_start,
                "genomic_end_grch38": g_end,
                "strand": strand_label,
                "transcript_id": transcript_id,
                "ensembl_exon_id": exon["id"],
                "transcript_exon_start": transcript_exon_start,
                "transcript_exon_end": transcript_exon_end,
                "cds_start": cds_start,
                "cds_end": cds_end,
                "coding_length_bp": coding_length,
                "splice_phase_5prime": splice_phase_5prime,
                "splice_phase_3prime": splice_phase_3prime,
                "cumulative_cds_start": cumulative_cds_start,
                "cumulative_cds_end": cumulative_cds_end,
                "source": "Ensembl REST API",
                "source_version": f"GRCh38 {transcript_id}",
                "date_accessed": today,
            }
        )

    expected_cds = REFERENCE.protein_length_aa * 3 + 3  # includes stop codon
    if cumulative_cds != expected_cds:
        raise ValueError(
            f"Total CDS length {cumulative_cds} bp != expected {expected_cds} bp"
        )

    return rows


def load_edge_shapes() -> dict[int, tuple[str, str]]:
    """Load curated 5'/3' puzzle-piece shapes keyed by exon number."""
    if not EDGE_SHAPES_CSV.exists():
        return {}
    with EDGE_SHAPES_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        int(row["exon"]): (row["five_prime_shape"], row["three_prime_shape"])
        for row in rows
    }


def apply_edge_shapes(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Attach five_prime_shape and three_prime_shape to fetched exon rows."""
    shapes = load_edge_shapes()
    enriched: list[dict[str, Any]] = []
    for row in rows:
        exon_number = int(row["exon_number"])
        fp, tp = shapes.get(exon_number, ("", ""))
        enriched.append(
            {
                "exon_number": row["exon_number"],
                "five_prime_shape": fp,
                "three_prime_shape": tp,
                **{k: v for k, v in row.items() if k != "exon_number"},
            }
        )
    return enriched


def write_exon_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write exon reference table to CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def cache_raw_transcript(transcript: dict[str, Any]) -> Path:
    """Save raw Ensembl JSON for reproducibility."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / "ensembl_ENST00000357033.json"
    cache_path.write_text(json.dumps(transcript, indent=2), encoding="utf-8")
    return cache_path


def build_domains_table() -> list[dict[str, Any]]:
    """
    Broad dystrophin protein domains in amino-acid coordinates.

    Boundaries are taken from UniProt P11532 feature annotations for
    full-length dystrophin (same protein as NP_003997.2). Genomic coordinates
    are intentionally omitted; map via CDS separately when needed.
    """
    today = date.today().isoformat()
    return [
        {
            "domain_id": "ABD",
            "domain_name": "Actin-binding domain",
            "protein_start_aa": 1,
            "protein_end_aa": 246,
            "cds_start_bp": 1,
            "cds_end_bp": 738,
            "display_fill_hex": "",
            "display_border_hex": "",
            "notes": "N-terminal actin-binding / CH domain region",
            "source": "UniProt P11532",
            "source_version": "NP_003997.2",
            "date_accessed": today,
        },
        {
            "domain_id": "ROD",
            "domain_name": "Central rod domain",
            "protein_start_aa": 247,
            "protein_end_aa": 3040,
            "cds_start_bp": 739,
            "cds_end_bp": 9120,
            "display_fill_hex": "",
            "display_border_hex": "",
            "notes": "Spectrin-repeat rod region (broad; repeats not split)",
            "source": "UniProt P11532",
            "source_version": "NP_003997.2",
            "date_accessed": today,
        },
        {
            "domain_id": "HINGE",
            "domain_name": "Hinge regions",
            "protein_start_aa": 541,
            "protein_end_aa": 2503,
            "cds_start_bp": 1621,
            "cds_end_bp": 7509,
            "display_fill_hex": "",
            "display_border_hex": "",
            "notes": (
                "Broad hinge span covering hinge 1 (541-611), hinge 2 "
                "(2134-2193), hinge 3 (2444-2503); overlaps rod domain"
            ),
            "source": "UniProt P11532",
            "source_version": "NP_003997.2",
            "date_accessed": today,
        },
        {
            "domain_id": "CR",
            "domain_name": "Cysteine-rich domain",
            "protein_start_aa": 3081,
            "protein_end_aa": 3360,
            "cds_start_bp": 9241,
            "cds_end_bp": 10080,
            "display_fill_hex": "#DAE9F8",
            "display_border_hex": "#4D93D9",
            "notes": "Cysteine-rich / WW-binding region",
            "source": "UniProt P11532",
            "source_version": "NP_003997.2",
            "date_accessed": today,
        },
        {
            "domain_id": "CT",
            "domain_name": "Carboxy-terminal domain",
            "protein_start_aa": 3361,
            "protein_end_aa": 3685,
            "cds_start_bp": 10081,
            "cds_end_bp": 11055,
            "display_fill_hex": "",
            "display_border_hex": "",
            "notes": "C-terminal region preceding stop codon",
            "source": "UniProt P11532",
            "source_version": "NP_003997.2",
            "date_accessed": today,
        },
    ]


DOMAIN_COLUMNS = [
    "domain_id",
    "domain_name",
    "protein_start_aa",
    "protein_end_aa",
    "cds_start_bp",
    "cds_end_bp",
    "display_fill_hex",
    "display_border_hex",
    "notes",
    "source",
    "source_version",
    "date_accessed",
]


def write_domains_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write domain reference table to CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DOMAIN_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch DMD reference exon data.")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output CSV path for exon table",
    )
    parser.add_argument(
        "--domains-output",
        type=Path,
        default=DOMAINS_OUTPUT,
        help="Output CSV path for domain table",
    )
    parser.add_argument(
        "--ensembl-id",
        default=REFERENCE.ensembl_transcript.split(".")[0],
        help="Ensembl transcript ID (version optional)",
    )
    args = parser.parse_args()

    print(f"Fetching {args.ensembl_id} from Ensembl ({REFERENCE.genome_assembly})...")
    transcript = fetch_ensembl_transcript(args.ensembl_id)
    cache_path = cache_raw_transcript(transcript)
    print(f"Cached raw response: {cache_path}")

    rows = build_exon_table(transcript)
    rows = apply_edge_shapes(rows)
    write_exon_csv(rows, args.output)
    print(f"Wrote {len(rows)} exons to {args.output}")

    domains = build_domains_table()
    write_domains_csv(domains, args.domains_output)
    print(f"Wrote {len(domains)} domain rows to {args.domains_output}")

    total_coding = sum(int(r["coding_length_bp"]) for r in rows)
    print(
        f"Summary: {len(rows)} exons, {total_coding} coding bp, "
        f"strand={rows[0]['strand']}, chr{rows[0]['chromosome']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
