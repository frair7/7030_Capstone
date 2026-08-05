"""Methods, limitations, and provenance."""

from __future__ import annotations

import streamlit as st

from src.config import REFERENCE
from src.streamlit_theme import apply_widescreen_theme, render_nav_bar
from src.ui_helpers import init_session_state

init_session_state()
apply_widescreen_theme()
render_nav_bar(page_title="Methods and Limitations")

st.error(
    "**Research and educational use only.** This application is not a clinical "
    "diagnostic tool and must not be used for patient care decisions."
)

st.subheader("Reference system")
st.markdown(
    f"""
    | Field | Value |
    |-------|-------|
    | Gene | {REFERENCE.gene_symbol} |
    | Isoform | {REFERENCE.protein_isoform} |
    | RefSeq transcript | {REFERENCE.refseq_transcript} |
    | RefSeq protein | {REFERENCE.refseq_protein} |
    | Ensembl transcript | {REFERENCE.ensembl_transcript} |
    | Assembly | {REFERENCE.genome_assembly} |
    | Chromosome / strand | {REFERENCE.chromosome} / {REFERENCE.strand} |
    | Coding exons | {REFERENCE.coding_exon_count} |
    | Protein length | {REFERENCE.protein_length_aa:,} aa |
    """
)

st.subheader("Supported variant formats")
st.markdown(
    """
    | Mode | Examples | Notes |
    |------|----------|-------|
    | Exon deletion | `del45`, `del45-50` | Whole-exon frame analysis supported |
    | Exon duplication | `dup2`, `dup3-7` | Tandem duplication model |
    | HGVS coding | `c.5287C>T`, `c.1812+1G>A` | Partial parser; may not resolve all exons |
    | Genomic (GRCh38) | `X:pos`, `chrX:start-end` | Assembly must be GRCh38 |
    | Exon selection | multiselect UI | Treated as deletion of selected exons |
    """
)

st.subheader("Reading-frame calculation")
st.markdown(
    """
    **Whole-exon deletions**

    1. Remove the selected exon(s) from the ordered transcript model.
    2. Join the upstream flanking exon (*n*−1) to the downstream flanking exon (*n*+1).
    3. Sum coding bases removed; if divisible by 3 **and** splice phases at the new
       junction are compatible, the deletion is **in frame**.
    4. Splice phases are derived from cumulative CDS length mod 3 in the reference table.

    **Tandem duplications**

    An extra copy of the exon block is modeled immediately downstream in transcript
    order. Added coding bases divisible by 3 → in frame (under stated assumptions).

    **SNVs, indels, splice variants**

    Mapped to exon/intron when possible; frame effect reported only when the
    sequence consequence is sufficiently defined.
    """
)

st.subheader("Exon-skipping search")
st.markdown(
    """
    For out-of-frame whole-exon deletions, the tool searches for **additional**
    exon skips (adjacent first, then broader in advanced mode) that restore the
    computational reading frame. Candidates are ranked by:

    1. Fewest additional exons skipped
    2. Fewest additional coding bases removed
    3. Contiguous skip preferred over non-contiguous

    **Never** ranked by claimed therapeutic efficacy.
    """
)

st.subheader("Known limitations")
st.markdown(
    """
    - Computational frame restoration ≠ therapeutic feasibility or clinical approval
    - HGVS normalization is incomplete in v1
    - Duplication structure (orientation, exact breakpoints) may be unknown
    - Partial exon deletions and deep-intronic variants may not yield frame calls
    - Domain boundaries are broad UniProt regions, not precise structural domains
    - Minus-strand display: exon 1 has the highest genomic coordinates
    """
)

st.subheader("Data provenance")
st.markdown(
    """
    - **Exon table:** Ensembl REST API, ENST00000357033.9, GRCh38
    - **Domains:** UniProt P11532 (amino-acid coordinates)
    - **Regenerate:** `python scripts/fetch_reference_data.py`
    - **Validate:** `python scripts/validate_reference_data.py`
    """
)

st.subheader("AI use statement")
st.markdown(
    "See [AI_USE.md](https://github.com) in the repository or `DMD_webapp/AI_USE.md` "
    "for documentation of Cursor/AI-assisted development."
)

st.caption("DMD Explorer — BSGP 7030 Capstone")
