"""
Centralized reference configuration for the DMD Explorer.

Transcript and assembly versions can be updated here as new releases become
available. The Streamlit UI always displays these values to the user.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class TranscriptReference:
    """Canonical DMD transcript metadata for Dp427m."""

    gene_symbol: str = "DMD"
    protein_isoform: str = "Dp427m"
    refseq_transcript: str = "NM_004006.3"
    refseq_protein: str = "NP_003997.2"
    ensembl_transcript: str = "ENST00000357033.9"
    genome_assembly: str = "GRCh38"
    chromosome: str = "X"
    strand: str = "reverse"
    coding_exon_count: int = 79
    protein_length_aa: int = 3685


# Singleton used throughout the application
REFERENCE = TranscriptReference()


def reference_summary() -> dict[str, str | int]:
    """Return reference metadata as a plain dict for display in the UI."""
    ref = REFERENCE
    return {
        "Gene": ref.gene_symbol,
        "Protein isoform": ref.protein_isoform,
        "RefSeq transcript": ref.refseq_transcript,
        "RefSeq protein": ref.refseq_protein,
        "Ensembl transcript": ref.ensembl_transcript,
        "Genome assembly": ref.genome_assembly,
        "Chromosome": ref.chromosome,
        "Strand": ref.strand,
        "Coding exons": ref.coding_exon_count,
        "Protein length (aa)": ref.protein_length_aa,
    }
