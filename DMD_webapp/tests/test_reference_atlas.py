"""Validation tests for the continuous-scrolling DMD Reference Map."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.config import REFERENCE
from src.dp427m_exon_data import get_exon_table
from src.reference_atlas import (
    build_detailed_protein_features,
    build_exon_reference_frame,
    create_genomic_viewer,
    create_protein_feature_viewer,
    create_protein_overview,
    create_transcript_viewer,
    genomic_reference_table,
    resource_table,
    transcript_reference_table,
)
from src.reference_data import (
    load_domain_table,
    load_exon_table,
    load_protein_domains_table,
)


@pytest.fixture(scope="module")
def atlas_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    exons = build_exon_reference_frame(load_exon_table(), load_domain_table())
    features = build_detailed_protein_features(
        load_protein_domains_table(),
        exons,
    )
    return exons, features


def test_exon_reference_uses_complete_spliced_transcript(
    atlas_data,
) -> None:
    exons, _ = atlas_data

    assert len(exons) == 79
    assert int(exons["transcript_exon_start"].min()) == 1
    assert int(exons["transcript_exon_end"].max()) == 13992
    assert int(exons["coding_length_bp"].sum()) == 11058
    assert int(exons["utr_5prime_length_bp"].sum()) == 237
    assert int(exons["utr_3prime_length_bp"].sum()) == 2697


def test_genomic_coding_bounds_are_numeric_and_length_checked(
    atlas_data,
) -> None:
    exons, _ = atlas_data
    table = genomic_reference_table(exons)

    assert (
        exons["coding_genomic_length_bp"] == exons["coding_length_bp"]
    ).all()
    assert (
        exons["coding_genomic_start"] >= exons["genomic_start_grch38"]
    ).all()
    assert (
        exons["coding_genomic_end"] <= exons["genomic_end_grch38"]
    ).all()
    assert pd.api.types.is_integer_dtype(table["Coding Genomic Start"].dtype)
    assert set(table["Genome Assembly"]) == {REFERENCE.genome_assembly}


def test_transcript_table_keeps_coordinate_systems_separate(
    atlas_data,
) -> None:
    exons, _ = atlas_data
    table = transcript_reference_table(exons)

    assert "Transcript Start" in table
    assert "CDS Start" in table
    assert "Encoded Amino-Acid Start" in table
    assert int(table.loc[0, "Transcript Start"]) == 1
    assert int(table.loc[0, "CDS Start"]) == 1
    assert int(table.loc[0, "5′ UTR Length"]) == 237


def test_detailed_features_convert_transcript_nt_to_amino_acids(
    atlas_data,
) -> None:
    _, features = atlas_data

    assert not features.empty
    assert features["Amino-Acid Start"].between(
        1,
        REFERENCE.protein_length_aa,
    ).all()
    assert features["Amino-Acid End"].between(
        1,
        REFERENCE.protein_length_aa,
    ).all()
    assert not features["Feature"].str.contains("UTR", case=False).any()
    abd = features.loc[features["Feature"] == "Actin binding domain"].iloc[0]
    assert int(abd["Source Transcript-nt Start"]) == 278
    assert int(abd["Amino-Acid Start"]) == 12
    assert "derived" in str(abd["Notes"]).lower()


def test_viewers_use_level_specific_axis_titles(atlas_data) -> None:
    exons, features = atlas_data
    genomic = create_genomic_viewer(exons)
    schematic = create_transcript_viewer(
        exons,
        get_exon_table(),
        schematic=True,
    )
    transcript = create_transcript_viewer(
        exons,
        get_exon_table(),
        schematic=False,
    )
    protein = create_protein_overview(load_domain_table())
    detailed = create_protein_feature_viewer(features)

    assert genomic.layout.xaxis.title.text == "Chromosome X genomic position (bp)"
    assert schematic.layout.xaxis.title.text.startswith("Exon order")
    assert (
        transcript.layout.xaxis.title.text
        == "Dp427m mature transcript position (bp)"
    )
    assert "amino acids" in protein.layout.xaxis.title.text
    assert detailed.layout.xaxis.title.text == protein.layout.xaxis.title.text


def test_external_resources_use_validated_local_locus(atlas_data) -> None:
    exons, _ = atlas_data
    start = int(exons["genomic_start_grch38"].min())
    end = int(exons["genomic_end_grch38"].max())
    resources = resource_table(start, end)

    assert resources["Link"].str.startswith("https://").all()
    ucsc = resources.loc[
        resources["Resource"] == "Open DMD in UCSC (GRCh38)",
        "Link",
    ].iloc[0]
    assert "db=hg38" in ucsc
    assert f"{start}-{end}" in ucsc
    blast = resources.loc[
        resources["Resource"] == "Open NCBI BLAST",
        "Purpose",
    ].iloc[0]
    assert "not a coordinate-table source" in blast


def test_reference_page_is_one_continuous_ordered_page() -> None:
    page = (
        Path(__file__).resolve().parents[1] / "pages" / "4_Reference_Map.py"
    ).read_text(encoding="utf-8")

    assert "st.tabs(" not in page
    headings = [
        '"1. Genomic Reference"',
        '"2. Transcript Reference"',
        '"3. Protein Domains and Binding Sites"',
        '"4. Data Provenance and External Resources"',
    ]
    positions = [page.index(heading) for heading in headings]
    assert positions == sorted(positions)
    assert "create_genomic_viewer" in page
    assert "create_transcript_viewer" in page
    assert "create_protein_overview" in page
