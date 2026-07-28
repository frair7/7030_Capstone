"""
DMD Mutation and Exon-Skipping Explorer — home page.
"""

from __future__ import annotations

import streamlit as st

from src.config import REFERENCE
from src.reference_data import EXONS_CSV
from src.ui_helpers import init_session_state, render_research_warning

st.set_page_config(
    page_title="DMD Mutation Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

init_session_state()

st.title("DMD Mutation and Exon-Skipping Explorer")
st.markdown(
    """
    Locally hosted **research and educational** web application for exploring
    mutations in the human *DMD* gene relative to the **Dp427m** transcript
    (NM_004006.3, 79 coding exons, GRCh38).
    """
)
render_research_warning()

st.info(
    f"**Active reference:** {REFERENCE.refseq_transcript} ({REFERENCE.protein_isoform}) "
    f"on **{REFERENCE.genome_assembly}**, chromosome {REFERENCE.chromosome} "
    f"({REFERENCE.strand} strand)."
)

col1, col2 = st.columns(2)
with col1:
    st.subheader("Get started")
    st.markdown(
        """
        1. Open **Mutation Explorer** in the sidebar
        2. Choose an input mode and enter a variant (e.g. `del45-50`)
        3. Click **Run analysis**
        4. Explore the transcript map, frame result, and skip candidates
        """
    )
with col2:
    st.subheader("Pages")
    st.markdown(
        """
        | Page | Description |
        |------|-------------|
        | **Mutation Explorer** | Main analysis workflow |
        | **Exon Skipping Analysis** | Candidate strategies in detail |
        | **Reference Map** | Full 79-exon reference table |
        | **Methods & Limitations** | Methods, provenance, disclaimers |
        """
    )

if EXONS_CSV.exists():
    import pandas as pd

    exons = pd.read_csv(EXONS_CSV)
    st.success(f"Reference data ready: **{len(exons)} exons** loaded.")
else:
    st.warning("Run `python scripts/fetch_reference_data.py` to fetch reference data.")

st.subheader("Supported variant formats")
st.markdown(
    """
    - **Exon deletion / duplication:** `del45-50`, `dup2`, `dup3-7`
    - **HGVS coding:** `c.5287C>T`, `c.1812+1G>A` (partial support)
    - **Genomic (GRCh38):** `X:31140000`, `chrX:start-end`
    - **Direct exon selection:** multiselect on Mutation Explorer page
    """
)

st.caption("DMD Explorer v0.3.0 — BSGP 7030 Capstone")
