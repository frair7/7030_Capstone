"""Reference map and exon/domain tables."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from src.config import REFERENCE
from src.coordinate_mapper import build_exon_records
from src.reference_data import (
    DOMAINS_CSV,
    EXONS_CSV,
    MAP_STYLES_CSV,
    PROTEIN_DOMAINS_CSV,
    load_domain_table,
    load_exon_table,
    load_map_styles_table,
    load_protein_domains_table,
)
from src.ui_helpers import check_reference_data, init_session_state, plotly_config
from src.visualization import ViewMode, create_combined_figure, create_domain_figure

init_session_state()
st.header("Reference Map")
st.markdown(
    f"""
    Canonical **Dp427m** reference map for **{REFERENCE.refseq_transcript}**
    on **{REFERENCE.genome_assembly}**. DMD lies on the **{REFERENCE.strand}**
    strand of chromosome {REFERENCE.chromosome}; exon numbers follow transcript
    order (exon 1 = 5′ mRNA end).
    """
)

if not check_reference_data():
    st.stop()

exon_rows = load_exon_table()
domain_rows = load_domain_table()
map_style_rows = load_map_styles_table()
protein_domain_rows = load_protein_domains_table()

view = st.radio(
    "Display mode",
    ["Exon-order schematic", "Transcript-scale"],
    horizontal=True,
)
view_mode = ViewMode.SCHEMATIC if view.startswith("Exon-order") else ViewMode.GENOMIC

exons = build_exon_records()
fig = create_combined_figure(exons, view_mode=view_mode, title="DMD reference map")
st.plotly_chart(fig, use_container_width=True, config=plotly_config())

st.subheader("Exon reference table")
st.dataframe(exon_rows, use_container_width=True, height=400)

st.subheader("Exon map styling")
st.dataframe(map_style_rows, use_container_width=True, height=300)

st.subheader("Protein domain summary (visualization)")
st.dataframe(domain_rows, use_container_width=True)

st.subheader("Protein domains and features (Dp427m)")
st.dataframe(protein_domain_rows, use_container_width=True, height=300)

domain_fig = create_domain_figure(domain_rows)
st.plotly_chart(domain_fig, use_container_width=True, config=plotly_config())

st.subheader("Data provenance")
st.markdown(
    f"""
    - **Exon coordinates:** {exon_rows[0]['source']} — {exon_rows[0]['source_version']}
    - **Date accessed:** {exon_rows[0]['date_accessed']}
    - **Map styling:** {map_style_rows[0]['source']} ({map_style_rows[0]['source_version']})
    - **Domain annotations:** {protein_domain_rows[0]['source']} ({protein_domain_rows[0]['source_version']})
    - Tables are stored in `reference_tables/` at the repository root.
    """
)

col1, col2, col3, col4 = st.columns(4)
col1.download_button(
    "Exons (CSV)",
    data=EXONS_CSV.read_bytes(),
    file_name="dmd_exons_grch38.csv",
    mime="text/csv",
    use_container_width=True,
)
col2.download_button(
    "Map styles (CSV)",
    data=MAP_STYLES_CSV.read_bytes(),
    file_name="dmd_exon_map_styles.csv",
    mime="text/csv",
    use_container_width=True,
)
col3.download_button(
    "Domain summary (CSV)",
    data=DOMAINS_CSV.read_bytes(),
    file_name="dmd_domains.csv",
    mime="text/csv",
    use_container_width=True,
)
col4.download_button(
    "Protein domains (CSV)",
    data=PROTEIN_DOMAINS_CSV.read_bytes(),
    file_name="dmd_protein_domains_dp427m.csv",
    mime="text/csv",
    use_container_width=True,
)
