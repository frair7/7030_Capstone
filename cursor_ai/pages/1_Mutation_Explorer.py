"""Mutation Explorer — interactive cohort mutation map and visualization."""

from __future__ import annotations

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from src.coordinate_mapper import build_exon_records
from src.mutation_catalog import (
    CATALOG_COLUMNS,
    EXPLORER_DISPLAY_COLUMNS,
    active_catalog,
    catalog_for_plot_selection,
    filter_catalog,
    load_catalog,
)
from src.mutation_map_figure import create_cohort_mutation_map
from src.mutation_map_svg import build_interactive_map_html, map_iframe_height
from src.mutation_viz import sort_mutations_by_exon_range
from src.streamlit_theme import apply_widescreen_theme, render_nav_bar
from src.ui_helpers import check_reference_data, init_session_state, render_research_warning

init_session_state()
apply_widescreen_theme()
render_nav_bar(page_title="Mutation Explorer")

if "map_selection_rows" not in st.session_state:
    st.session_state.map_selection_rows = []
if "show_intake_preview" not in st.session_state:
    st.session_state.show_intake_preview = False

render_research_warning()

if not check_reference_data():
    st.stop()

if "mutation_catalog" not in st.session_state:
    st.session_state.mutation_catalog = load_catalog()

catalog_df = st.session_state.mutation_catalog
exons = build_exon_records()

if "phenotype" not in catalog_df.columns:
    catalog_df["phenotype"] = ""
    for col in CATALOG_COLUMNS:
        if col not in catalog_df.columns:
            catalog_df[col] = ""

active_df = active_catalog(catalog_df)


def _plot_selection_ids(selection_df: pd.DataFrame) -> set[str]:
    """Return mutation_record_id values for rows checked in the Plot column."""
    if selection_df is None or selection_df.empty:
        return set()
    if "selected" not in selection_df.columns or "mutation_record_id" not in selection_df.columns:
        return set()
    plot_mask = selection_df["selected"].fillna(False).astype(bool)
    return set(selection_df.loc[plot_mask, "mutation_record_id"].dropna().astype(str).tolist())


st.caption(
    "Visualize catalog mutations on the transcript map. "
    "Add or edit records on the **Mutation Catalog** page."
)

# =============================================================================
# COHORT MUTATION MAP
# =============================================================================
_draft = st.session_state.get("mutation_draft")
preview_row = (
    _draft
    if st.session_state.get("show_intake_preview") and _draft and _draft.get("start_region")
    else None
)
map_rows = sort_mutations_by_exon_range(st.session_state.map_selection_rows)
if preview_row:
    map_rows = [preview_row] + [r for r in map_rows if r.get("id") != preview_row.get("id")]

m1, m2, m3 = st.columns([2, 2, 6])
with m1:
    st.metric("Plotted mutations", len(map_rows))
with m2:
    if preview_row:
        st.metric("Preview", "Draft on map", help="Intake draft from Mutation Catalog")
with m3:
    st.caption(
        "Hover exons for CDS length and frame edges. Hover deletion bars for participant details. "
        "Transcript and domain rows share identical exon positions; dual-domain exons fade gradually."
    )

map_html = build_interactive_map_html(
    map_rows,
    chart_width=2000,
    title="DMD transcript & mutation alignment",
)
iframe_h = map_iframe_height(len(map_rows))
components.html(map_html, height=iframe_h, scrolling=True)

with st.expander("Download static map (PNG)", expanded=False):
    map_fig = create_cohort_mutation_map(
        st.session_state.map_selection_rows,
        exons,
        preview_row=preview_row,
        title="DMD transcript & mutation alignment",
    )
    st.pyplot(map_fig, use_container_width=True)

# =============================================================================
# PLOT ON MAP
# =============================================================================
st.subheader("Plot on map")

f1, f2, f3, f4, f5 = st.columns(5)
classes = [c for c in sorted(active_df["mutation_class"].unique()) if c]
subclasses = [s for s in sorted(active_df["mutation_subclass"].unique()) if s]
frames = [f for f in sorted(active_df["frame"].unique()) if f]
phenotypes = [p for p in sorted(active_df.get("phenotype", pd.Series(dtype=str)).unique()) if p]

with f1:
    filter_class = st.multiselect("Class", classes, key="plot_filter_class")
with f2:
    filter_subclass = st.multiselect("Subclass", subclasses, key="plot_filter_subclass")
with f3:
    filter_frame = st.multiselect("Frame", frames, key="plot_filter_frame")
with f4:
    filter_phenotype = st.multiselect("Phenotype", phenotypes, key="plot_filter_phenotype")
with f5:
    filter_search = st.text_input("Search", key="plot_filter_search")

filtered = filter_catalog(
    catalog_df,
    mutation_classes=filter_class or None,
    mutation_subclasses=filter_subclass or None,
    frames=filter_frame or None,
    search_text=filter_search,
)
if filter_phenotype:
    filtered = filtered[filtered["phenotype"].isin(filter_phenotype)]

selection_df = catalog_for_plot_selection(filtered)
if "plot_selection_df" in st.session_state:
    prev = st.session_state.plot_selection_df
    if (
        not prev.empty
        and "mutation_record_id" in prev.columns
        and "selected" in prev.columns
    ):
        prev_sel = dict(zip(
            prev["mutation_record_id"].astype(str),
            prev["selected"].fillna(False).astype(bool),
        ))
        selection_df["selected"] = selection_df["mutation_record_id"].astype(str).map(
            lambda rid: prev_sel.get(rid, False),
        )

plot_column_config = {
    "selected": st.column_config.CheckboxColumn(
        "Plot",
        help="Include this mutation on the transcript map.",
        default=False,
    ),
    "mutation_record_id": st.column_config.TextColumn(
        "Record ID", disabled=True, width="small",
    ),
    "id": st.column_config.TextColumn("Participant ID", disabled=True, width="small"),
    "mutation_class": st.column_config.TextColumn("Class", disabled=True),
    "mutation_subclass": st.column_config.TextColumn("Subclass", disabled=True),
    "start_region": st.column_config.TextColumn("Start", disabled=True),
    "stop_region": st.column_config.TextColumn("Stop", disabled=True),
    "frame": st.column_config.TextColumn("Frame", disabled=True),
    "phenotype": st.column_config.TextColumn("Phenotype", disabled=True),
    "group": st.column_config.TextColumn("Group", disabled=True),
}

edited_selection_df = st.data_editor(
    selection_df,
    key="mutation_plot_editor",
    hide_index=True,
    use_container_width=True,
    column_config=plot_column_config,
    column_order=EXPLORER_DISPLAY_COLUMNS,
    disabled=[c for c in EXPLORER_DISPLAY_COLUMNS if c != "selected"],
)
st.session_state.plot_selection_df = edited_selection_df

plot_ids = _plot_selection_ids(edited_selection_df)
st.caption(
    f"**{len(filtered)}** of **{len(active_df)}** active mutations · "
    f"**{len(plot_ids)}** selected for plotting"
)

c_plot, c_clear_map = st.columns([1, 1])
if c_plot.button("Plot selected on map", type="primary"):
    rows = filtered[filtered["mutation_record_id"].isin(plot_ids)]
    st.session_state.map_selection_rows = sort_mutations_by_exon_range(
        rows.to_dict("records"),
    )
    st.session_state.show_intake_preview = False
    st.rerun()

if c_clear_map.button("Clear map selection"):
    st.session_state.map_selection_rows = []
    st.session_state.show_intake_preview = False
    st.rerun()
