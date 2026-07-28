"""
Generate downloadable CSV, plain-text, and HTML analysis reports.
"""

from __future__ import annotations

import io
from datetime import datetime, timezone
from typing import Any, Optional

import pandas as pd

from src.config import REFERENCE
from src.models import FrameResult, ParsedVariant, SkipCandidate
from src.visualization import VisualizationState, build_exon_table


def skip_candidates_dataframe(
    candidates: list[SkipCandidate],
) -> pd.DataFrame:
    """Tabular skip-candidate results."""
    if not candidates:
        return pd.DataFrame(
            columns=[
                "original_mutation",
                "deleted_exons",
                "additional_skipped_exons",
                "final_junction",
                "additional_bp_removed",
                "total_bp_removed",
                "restores_frame",
                "estimated_protein_aa",
                "evidence_class",
            ]
        )
    return pd.DataFrame(
        [
            {
                "original_mutation": c.original_mutation,
                "deleted_exons": str(c.deleted_exons),
                "additional_skipped_exons": str(c.additional_skipped_exons),
                "final_junction": (
                    f"exon {c.final_upstream_exon}|{c.final_downstream_exon}"
                ),
                "additional_bp_removed": c.additional_coding_bases_removed,
                "total_bp_removed": c.total_coding_bases_removed,
                "restores_frame": c.restores_frame,
                "estimated_protein_aa": c.estimated_protein_aa,
                "evidence_class": c.evidence_class.value,
            }
            for c in candidates
        ]
    )


def variant_summary_dataframe(
    variant: ParsedVariant,
    frame_result: FrameResult,
    mapping_message: str,
) -> pd.DataFrame:
    """Single-row variant summary for CSV export."""
    return pd.DataFrame(
        [
            {
                "raw_input": variant.raw_input,
                "input_mode": variant.input_mode.value,
                "variant_type": variant.variant_type.value,
                "first_exon": variant.first_exon,
                "last_exon": variant.last_exon,
                "frame_status": frame_result.status.value,
                "frame_explanation": frame_result.explanation,
                "junction_upstream": frame_result.upstream_exon,
                "junction_downstream": frame_result.downstream_exon,
                "mapping": mapping_message,
                "transcript": REFERENCE.refseq_transcript,
                "assembly": REFERENCE.genome_assembly,
            }
        ]
    )


def build_results_csv(
    variant: ParsedVariant,
    frame_result: FrameResult,
    mapping_message: str,
    skip_candidates: list[SkipCandidate],
    exons: list[Any],
    viz_state: Optional[VisualizationState] = None,
) -> bytes:
    """Combined CSV export: summary + skip candidates + exon table."""
    buffer = io.StringIO()
    variant_summary_dataframe(variant, frame_result, mapping_message).to_csv(
        buffer, index=False
    )
    buffer.write("\n")
    skip_candidates_dataframe(skip_candidates).to_csv(buffer, index=False)
    buffer.write("\n")
    build_exon_table(exons, viz_state).to_csv(buffer, index=False)
    return buffer.getvalue().encode("utf-8")


def build_text_report(
    variant: ParsedVariant,
    frame_result: FrameResult,
    mapping_message: str,
    skip_candidates: list[SkipCandidate],
) -> str:
    """Plain-text analysis summary."""
    lines = [
        "DMD Mutation and Exon-Skipping Explorer — Analysis Report",
        "=" * 60,
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "RESEARCH AND EDUCATIONAL USE ONLY — NOT FOR CLINICAL DIAGNOSIS",
        "",
        "Reference",
        f"  Transcript: {REFERENCE.refseq_transcript}",
        f"  Protein:    {REFERENCE.refseq_protein}",
        f"  Assembly:   {REFERENCE.genome_assembly}",
        f"  Strand:     {REFERENCE.strand}",
        "",
        "Variant input",
        f"  Raw:        {variant.raw_input}",
        f"  Mode:       {variant.input_mode.value}",
        f"  Type:       {variant.variant_type.value}",
    ]
    if variant.exon_range:
        lines.append(f"  Exons:      {variant.exon_range[0]}-{variant.exon_range[1]}")
    if mapping_message:
        lines.append(f"  Mapping:    {mapping_message}")
    if variant.warnings:
        lines.append(f"  Warnings:   {'; '.join(variant.warnings)}")
    if variant.errors:
        lines.append(f"  Errors:     {'; '.join(variant.errors)}")

    lines.extend(
        [
            "",
            "Reading-frame result",
            f"  Status:     {frame_result.status.value}",
            f"  Detail:     {frame_result.explanation}",
        ]
    )
    if frame_result.upstream_exon and frame_result.downstream_exon:
        lines.append(
            f"  Junction:   exon {frame_result.upstream_exon}|"
            f"{frame_result.downstream_exon}"
        )

    lines.append("")
    lines.append("Candidate exon-skipping strategies (computational only)")
    if skip_candidates:
        for i, c in enumerate(skip_candidates, 1):
            lines.append(
                f"  {i}. Additional skips {c.additional_skipped_exons} → "
                f"junction {c.final_upstream_exon}|{c.final_downstream_exon} "
                f"(+{c.additional_coding_bases_removed} bp)"
            )
    else:
        lines.append("  None found within search limits.")

    lines.extend(
        [
            "",
            "Limitations",
            "  - Computational frame restoration ≠ therapeutic feasibility",
            "  - HGVS normalization is partial in this version",
            "  - Splice-phase model assumes whole-exon events",
        ]
    )
    return "\n".join(lines)


def build_html_report(
    variant: ParsedVariant,
    frame_result: FrameResult,
    mapping_message: str,
    skip_candidates: list[SkipCandidate],
) -> str:
    """Simple HTML report for browser download."""
    text = build_text_report(
        variant, frame_result, mapping_message, skip_candidates
    )
    body = "<br>".join(text.replace("&", "&amp;").replace("<", "&lt;").splitlines())
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>DMD Explorer Report</title>
  <style>
    body {{ font-family: sans-serif; max-width: 800px; margin: 2rem auto; }}
    .warning {{ color: #b00020; font-weight: bold; border: 2px solid #b00020;
                padding: 1rem; margin-bottom: 1rem; }}
  </style>
</head>
<body>
  <div class="warning">Research and educational use only — not for clinical diagnosis.</div>
  <p>{body}</p>
</body>
</html>"""
