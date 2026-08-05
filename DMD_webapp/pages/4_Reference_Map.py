"""Single-page DMD genomic, transcript, protein, and provenance atlas."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from src.config import REFERENCE
from src.dp427m_exon_data import get_exon_table
from src.reference_atlas import (
    build_detailed_protein_features,
    build_exon_reference_frame,
    create_genomic_viewer,
    create_protein_feature_viewer,
    create_protein_overview,
    create_style_summary_figure,
    create_transcript_viewer,
    genomic_reference_table,
    resource_table,
    style_reference_table,
    transcript_reference_table,
)
from src.reference_data import (
    load_domain_table,
    load_exon_table,
    load_protein_domains_table,
)
from src.streamlit_theme import apply_widescreen_theme, render_nav_bar
from src.ui_helpers import check_reference_data, init_session_state, plotly_config


@st.cache_data(show_spinner=False)
def _load_reference_atlas() -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list[dict],
]:
    exon_rows = load_exon_table()
    broad_domains = load_domain_table()
    detailed_features = load_protein_domains_table()
    exons = build_exon_reference_frame(exon_rows, broad_domains)
    protein_features = build_detailed_protein_features(detailed_features, exons)
    return (
        exons,
        genomic_reference_table(exons),
        transcript_reference_table(exons),
        protein_features,
        broad_domains,
    )


def _anchor(identifier: str) -> None:
    st.markdown(f'<div id="{identifier}"></div>', unsafe_allow_html=True)


def _major_section(identifier: str, title: str, description: str) -> None:
    st.divider()
    _anchor(identifier)
    st.header(title)
    st.markdown(description)


def _download_table(label: str, frame: pd.DataFrame, filename: str) -> None:
    st.download_button(
        label,
        data=frame.to_csv(index=False).encode("utf-8"),
        file_name=filename,
        mime="text/csv",
        use_container_width=False,
    )


def _filter_exons(frame: pd.DataFrame, query: str) -> pd.DataFrame:
    query = query.strip()
    if not query:
        return frame
    return frame.loc[
        frame.astype(str).apply(
            lambda row: row.str.contains(query, case=False, regex=False).any(),
            axis=1,
        )
    ]


init_session_state()
apply_widescreen_theme()
render_nav_bar(page_title="Reference Map")
_anchor("reference-map")

st.markdown(
    "Explore the DMD gene across genomic, transcript, and protein coordinate "
    "systems. Each section uses the coordinate system appropriate for that "
    "biological level."
)

coordinate_overview = pd.DataFrame(
    [
        {
            "Biological level": "Genome",
            "Coordinate system": "Chromosome X genomic position, base pairs",
            "Main purpose": (
                "Genomic location, exon boundaries, splice sites, and sequence lookup"
            ),
        },
        {
            "Biological level": "Transcript",
            "Coordinate system": "Spliced Dp427m transcript position, base pairs",
            "Main purpose": (
                "Exon structure, transcript length, coding sequence, UTRs, "
                "and reading-frame reference"
            ),
        },
        {
            "Biological level": "Protein",
            "Coordinate system": "Dp427m amino acid position",
            "Main purpose": (
                "Protein domains, motifs, binding sites, hinges, and spectrin repeats"
            ),
        },
    ]
)

st.subheader("Coordinate Systems Overview")
st.dataframe(coordinate_overview, use_container_width=True, hide_index=True)
st.warning(
    "Genome, transcript, CDS, and protein coordinates are related but are not "
    "interchangeable. Every genomic coordinate below is explicitly labeled "
    f"{REFERENCE.genome_assembly}."
)
st.markdown(
    """
    **Jump to:** [Genomic Reference](#genomic-reference) ·
    [Transcript Reference](#transcript-reference) ·
    [Protein Domains and Binding Sites](#protein-reference) ·
    [Data Provenance](#data-provenance)
    """
)

if not check_reference_data():
    st.stop()

(
    exons,
    genomic_table,
    transcript_table,
    protein_features,
    broad_domains,
) = _load_reference_atlas()
exon_domain_rows = get_exon_table()
gene_start = int(exons["genomic_start_grch38"].min())
gene_end = int(exons["genomic_end_grch38"].max())
transcript_length = int(exons["transcript_exon_end"].max())
cds_length = int(exons["coding_length_bp"].sum())


_major_section(
    "genomic-reference",
    "1. Genomic Reference",
    "This section displays the DMD locus on chromosome X using genomic "
    "coordinates. Introns and the physical distances between exons are "
    "represented in the genomic coordinate system.",
)
assembly_col, locus_col, strand_col, transcript_col = st.columns(4)
assembly_col.metric("Genome assembly", REFERENCE.genome_assembly)
locus_col.metric("DMD locus", f"chrX:{gene_start:,}–{gene_end:,}")
strand_col.metric("Strand", "Reverse (−)")
transcript_col.metric("Transcript", REFERENCE.ensembl_transcript)

st.subheader("1.1 Genomic Viewer")
st.caption(
    "Chromosome coordinates retain standard ascending orientation. The arrow "
    "shows that DMD transcription proceeds toward decreasing coordinates."
)
st.plotly_chart(
    create_genomic_viewer(exons),
    use_container_width=True,
    config=plotly_config(),
)

st.subheader("1.2 Genomic Coding-Exon Reference Table")
st.markdown(
    f"**Genome assembly: {REFERENCE.genome_assembly}.** Primary coordinate "
    "columns report the coding portion of each exon. Full-exon genomic "
    "boundaries are retained in separately labeled columns."
)
genomic_query = st.text_input(
    "Search genomic exon reference",
    key="genomic_reference_search",
    placeholder="Search exon number, accession, chromosome, or coordinate",
)
filtered_genomic = _filter_exons(genomic_table, genomic_query)
st.dataframe(
    filtered_genomic,
    use_container_width=True,
    hide_index=True,
    height=390,
)
_download_table(
    "Download genomic coding-exon table (CSV)",
    genomic_table,
    "dmd_grch38_coding_exons.csv",
)
st.info(
    "Coordinates from different genome assemblies must not be combined without "
    "an explicit assembly conversion."
)
st.markdown("[Back to top](#reference-map)")


_major_section(
    "transcript-reference",
    "2. Transcript Reference",
    "This section displays the spliced Dp427m transcript. Introns are removed, "
    "and exon widths in transcript-scale mode are proportional to the nucleotide "
    "lengths represented in the mature transcript.",
)
st.subheader("2.1 Transcript Viewer")
transcript_view = st.radio(
    "Transcript display mode",
    ["Exon-order schematic", "Mature transcript scale"],
    horizontal=True,
    key="reference_transcript_view",
)
st.markdown(
    """
    **Exon-order schematic:** All exons are shown at a similar visual width for
    easy identification.

    **Mature transcript scale:** Exon widths reflect their nucleotide lengths
    in the spliced transcript. Introns are omitted.
    """
)
st.caption(
    f"The local Ensembl fields `transcript_exon_start/end` span the complete "
    f"{transcript_length:,}-bp spliced transcript, including {int(exons['utr_5prime_length_bp'].sum()):,} "
    f"bp of 5′ UTR and {int(exons['utr_3prime_length_bp'].sum()):,} bp of 3′ UTR. "
    "UTRs are shown in neutral gray."
)
st.plotly_chart(
    create_transcript_viewer(
        exons,
        exon_domain_rows,
        schematic=transcript_view == "Exon-order schematic",
    ),
    use_container_width=True,
    config=plotly_config(),
)

st.subheader("2.2 Exon Reference Table")
st.caption(
    "Transcript Start/End are full spliced-transcript coordinates. CDS Start/End "
    "are cumulative CDS offsets, not genomic positions or absolute transcript positions. "
    "Lengths are base pairs; encoded positions are amino acids."
)
transcript_query = st.text_input(
    "Search exon reference",
    key="transcript_reference_search",
    placeholder="Search exon, phase, or protein region",
)
filtered_transcript = _filter_exons(transcript_table, transcript_query)
st.dataframe(
    filtered_transcript,
    use_container_width=True,
    hide_index=True,
    height=420,
)
_download_table(
    "Download transcript exon table (CSV)",
    transcript_table,
    "dp427m_transcript_exons.csv",
)

st.subheader("2.3 Exon Styling Table")
st.caption(
    "One row is shown per meaningful styling class. These are visual conventions, "
    "not additional coordinate annotations."
)
style_table = style_reference_table()
st.dataframe(style_table, use_container_width=True, hide_index=True)
_download_table("Download styling classes (CSV)", style_table, "dmd_style_classes.csv")
st.markdown("[Back to top](#reference-map)")


_major_section(
    "protein-reference",
    "3. Protein Domains and Binding Sites",
    "This section displays the Dp427m protein using amino-acid coordinates. "
    "Protein coordinates are separate from genomic and transcript nucleotide coordinates.",
)
protein_col, refseq_col, uniprot_col = st.columns(3)
protein_col.metric("Dp427m protein length", f"{REFERENCE.protein_length_aa:,} aa")
refseq_col.metric("RefSeq protein", REFERENCE.refseq_protein)
uniprot_col.metric("UniProt", "P11532")

st.subheader("3.1 Protein Viewer")
st.plotly_chart(
    create_protein_overview(broad_domains),
    use_container_width=True,
    config=plotly_config(),
)

st.subheader("3.2 Protein Domain and Feature Viewer")
st.caption(
    "Overlapping structural regions, repeats, hinges, binding sites, and motifs "
    "are separated into tracks on the same Dp427m amino-acid axis."
)
st.plotly_chart(
    create_protein_feature_viewer(protein_features),
    use_container_width=True,
    config=plotly_config(),
)
st.warning(
    "Coordinate-source conflict retained for transparency: the detailed Leiden "
    "feature CSV uses transcript-nucleotide positions (5′ UTR 1–244; coding "
    "protein 245–11299), although its legacy headers say amino acids. The viewer "
    "derives amino-acid positions only for protein-coding features and preserves "
    "the original source positions in the table."
)

st.subheader("3.3 Protein Domains and Features Table")
protein_categories = sorted(protein_features["Feature Category"].unique())
selected_categories = st.multiselect(
    "Filter protein feature categories",
    protein_categories,
    default=protein_categories,
    key="protein_feature_categories",
)
protein_query = st.text_input(
    "Search protein features",
    key="protein_feature_search",
    placeholder="Search feature, source, category, or exon",
)
filtered_features = protein_features.loc[
    protein_features["Feature Category"].isin(selected_categories)
]
filtered_features = _filter_exons(filtered_features, protein_query)
st.dataframe(
    filtered_features,
    use_container_width=True,
    hide_index=True,
    height=460,
)
_download_table(
    "Download protein features (CSV)",
    protein_features,
    "dp427m_protein_features_aa.csv",
)

st.subheader("3.4 Protein Styling Summary Visualization")
st.caption(
    "Visual styling atlas only; this legend does not define scientific coordinate boundaries."
)
st.plotly_chart(
    create_style_summary_figure(),
    use_container_width=True,
    config={"displayModeBar": False},
)
st.markdown("[Back to top](#reference-map)")


_major_section(
    "data-provenance",
    "4. Data Provenance and External Resources",
    "This section documents the local coordinate sources, transformations, "
    "accessions, and external resources needed to reproduce or cross-check the atlas.",
)
date_accessed = str(exons.iloc[0]["date_accessed"])
provenance = pd.DataFrame(
    [
        ["Gene symbol", REFERENCE.gene_symbol, "NCBI Gene 1756", "", "NCBI / project config", "https://www.ncbi.nlm.nih.gov/gene/1756", date_accessed, "src/config.py", "Gene identifier only"],
        ["Chromosome", f"chr{REFERENCE.chromosome}", "", REFERENCE.genome_assembly, "Ensembl REST API", "https://rest.ensembl.org/", date_accessed, "reference_tables/dmd_exons_grch38.csv", "Standard ascending chromosome coordinates"],
        ["Genomic reference accession", "Not specified in local configuration", "Not stored", REFERENCE.genome_assembly, "—", "", date_accessed, "src/config.py", "Do not infer an accession silently"],
        ["Dp427m transcript", REFERENCE.refseq_transcript, REFERENCE.refseq_transcript, "", "NCBI RefSeq cross-reference", f"https://www.ncbi.nlm.nih.gov/nuccore/{REFERENCE.refseq_transcript}", date_accessed, "src/config.py", f"Coordinate table itself uses {REFERENCE.ensembl_transcript}"],
        ["Ensembl transcript coordinate model", REFERENCE.ensembl_transcript, REFERENCE.ensembl_transcript, REFERENCE.genome_assembly, "Ensembl REST API", "https://rest.ensembl.org/", date_accessed, "reference_tables/dmd_exons_grch38.csv", "Imported genomic/full-transcript/CDS fields"],
        ["Dp427m protein", REFERENCE.refseq_protein, f"{REFERENCE.refseq_protein}; UniProt P11532", "", "NCBI RefSeq / UniProt", "https://www.uniprot.org/uniprotkb/P11532/entry", date_accessed, "src/config.py", "Project protein references"],
        ["Transcript length", f"{transcript_length:,} bp", REFERENCE.ensembl_transcript, REFERENCE.genome_assembly, "Derived from local Ensembl exon table", "", date_accessed, "reference_tables/dmd_exons_grch38.csv", "Complete spliced transcript including UTR"],
        ["CDS length", f"{cds_length:,} bp", REFERENCE.ensembl_transcript, REFERENCE.genome_assembly, "Derived from local Ensembl exon table", "", date_accessed, "reference_tables/dmd_exons_grch38.csv", "Includes stop codon"],
        ["Protein length", f"{REFERENCE.protein_length_aa:,} aa", REFERENCE.refseq_protein, "", "Project config / RefSeq", f"https://www.ncbi.nlm.nih.gov/protein/{REFERENCE.refseq_protein}", date_accessed, "src/config.py", "Dp427m"],
        ["Exon-coordinate source", "79 exon records", REFERENCE.ensembl_transcript, REFERENCE.genome_assembly, "Ensembl REST API", "https://rest.ensembl.org/", date_accessed, "reference_tables/dmd_exons_grch38.csv", "Coding genomic bounds are derived locally from full-exon bounds, strand, and coding length"],
        ["Protein-domain source", "Broad domains", "UniProt P11532 / NP_003997.2", "", "UniProt-derived local summary", "https://www.uniprot.org/uniprotkb/P11532/entry", date_accessed, "data/dmd_domains.csv", "Direct amino-acid coordinates"],
        ["Detailed feature source", "48 transcript-position features", "Dp427m protein feature map", "", "Leiden DMD resources", "https://www.dmd.nl/DMD_home.html", date_accessed, "reference_tables/dmd_protein_domains_dp427m.csv", "Legacy headers conflict with transcript-nt values; AA coordinates derived explicitly"],
        ["Map styling", "Curated DMD palette", REFERENCE.ensembl_transcript, "", "DMD map modeling reference", "", date_accessed, "reference_tables/dmd_exon_map_styles.csv", "Manual visual curation; not a coordinate source"],
    ],
    columns=[
        "Data Element",
        "Value",
        "Accession or Version",
        "Genome Assembly",
        "Source Organization",
        "Source URL",
        "Date Accessed",
        "Local File or Module",
        "Notes",
    ],
)
st.dataframe(
    provenance,
    use_container_width=True,
    hide_index=True,
    column_config={
        "Source URL": st.column_config.LinkColumn("Source URL", display_text="Open source"),
    },
    height=475,
)
_download_table("Download provenance table (CSV)", provenance, "dmd_atlas_provenance.csv")

st.subheader("Coordinate safeguards")
st.markdown(
    f"""
    - All genomic coordinates on this page use **{REFERENCE.genome_assembly}**.
    - DMD is on the **reverse strand**; chromosome coordinates remain in standard
      ascending orientation while transcription direction is shown separately.
    - Genomic exon, full transcript-exon, CDS-offset, and protein amino-acid
      coordinates are labeled separately.
    - Protein-domain boundaries are not used to infer genomic exon boundaries.
    - BLAST is listed only for sequence validation, never as an exon-coordinate source.
    - GeneCards is an aggregated portal, not the primary source of project coordinates.
    - The project transcript accession and assembly are displayed without silent updates.
    """
)

st.subheader("External DMD and Sequence Resources")
resources = resource_table(gene_start, gene_end)
for group in resources["Group"].drop_duplicates():
    with st.expander(group, expanded=group == "DMD-Specific Resources"):
        group_rows = resources.loc[resources["Group"] == group]
        for _, resource in group_rows.iterrows():
            st.markdown(
                f"**[{resource['Resource']}]({resource['Link']})**  \n"
                f"{resource['Organization']} · {resource['Biological Level']} — "
                f"{resource['Purpose']}"
            )

st.caption(
    "NCBI BLAST is a sequence-comparison resource. It is not used as a genomic "
    "coordinate map or exon-table source."
)
st.markdown("[Back to top](#reference-map)")
