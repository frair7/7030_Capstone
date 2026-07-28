"""Exon-skipping candidate analysis page."""

from __future__ import annotations

import streamlit as st

from src.models import FrameStatus, VariantType
from src.reporting import skip_candidates_dataframe
from src.ui_helpers import (
    check_reference_data,
    ensure_analysis,
    init_session_state,
    plotly_config,
    render_input_sidebar,
    render_research_warning,
)
from src.visualization import (
    create_transcript_figure,
    visualization_state_from_variant,
)

init_session_state()
st.header("Exon Skipping Analysis")
render_research_warning()

if not check_reference_data():
    st.stop()

with st.sidebar:
    settings = render_input_sidebar()

bundle = ensure_analysis(settings)
variant = bundle.variant
frame_result = bundle.frame_result
skip_candidates = bundle.skip_candidates
exons = bundle.exons

st.markdown(
    """
    This page lists **computational** exon-skipping strategies that may restore
    the reading frame after an out-of-frame whole-exon deletion. Candidates are
    ranked by fewest additional exons skipped, then fewest additional coding
    bases removed.
    """
)

if variant.variant_type != VariantType.DELETION or not variant.exon_range:
    st.info("Enter a whole-exon **deletion** on the Mutation Explorer page and run analysis.")
    st.stop()

if frame_result.status == FrameStatus.IN_FRAME:
    st.success(
        f"Deletion {variant.raw_input} is already **in frame**. "
        "No additional skipping is computationally required."
    )
    st.stop()

if frame_result.status == FrameStatus.CANNOT_DETERMINE:
    st.warning("Frame status could not be determined for this variant.")
    st.stop()

st.subheader(f"Original mutation: `{variant.raw_input}`")
st.markdown(frame_result.explanation)

if not skip_candidates:
    st.warning(
        f"No frame-restoring candidates found within "
        f"{settings['max_additional_skips']} additional exon skips."
    )
    st.stop()

# Candidate selector
idx = st.selectbox(
    "Select candidate to visualize",
    range(len(skip_candidates)),
    format_func=lambda i: (
        f"#{i+1}: skip {skip_candidates[i].additional_skipped_exons} → "
        f"junction {skip_candidates[i].final_upstream_exon}|"
        f"{skip_candidates[i].final_downstream_exon}"
    ),
)
selected = skip_candidates[idx]

c1, c2, c3 = st.columns(3)
c1.metric("Additional exons skipped", len(selected.additional_skipped_exons))
c2.metric("Additional bp removed", selected.additional_coding_bases_removed)
c3.metric("Est. protein length (aa)", selected.estimated_protein_aa)

st.dataframe(skip_candidates_dataframe(skip_candidates), use_container_width=True)

st.subheader("Selected strategy visualization")
viz_state = visualization_state_from_variant(
    variant.first_exon,
    variant.last_exon,
    variant_type="deletion",
    skip_candidate=selected,
)
fig = create_transcript_figure(
    exons,
    viz_state,
    settings["view_mode"],
    title=f"Skip candidate — {selected.additional_skipped_exons}",
)
st.plotly_chart(fig, use_container_width=True, config=plotly_config())

st.subheader("Assumptions")
for item in selected.assumptions:
    st.markdown(f"- {item}")
if selected.affected_domains:
    st.markdown("**Affected protein domains (broad):** " + ", ".join(selected.affected_domains))

st.caption(
    "Evidence class: computational frame restoration only. "
    "Not an approved or experimentally validated therapy."
)
