"""Mutation Explorer — catalog intake, cohort mutation map, filterable table."""

from __future__ import annotations

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from src.coordinate_mapper import build_exon_records
from src.exon_display_config import PHENOTYPE_OPTIONS
from src.mutation_autofill import autofill_mutation_fields
from src.mutation_catalog import (
    CATALOG_COLUMNS,
    DISPLAY_COLUMNS,
    catalog_for_display,
    filter_catalog,
    load_catalog,
    next_id,
    save_catalog,
)
from src.mutation_map_figure import create_cohort_mutation_map
from src.mutation_map_svg import build_interactive_map_html, map_iframe_height
from src.region_parser import parse_quick_entry
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
if "mutation_draft" not in st.session_state:
    st.session_state.mutation_draft = {col: "" for col in CATALOG_COLUMNS}

catalog_df = st.session_state.mutation_catalog
exons = build_exon_records()

if "phenotype" not in catalog_df.columns:
    catalog_df["phenotype"] = ""
    for col in CATALOG_COLUMNS:
        if col not in catalog_df.columns:
            catalog_df[col] = ""

# =============================================================================
# COHORT MUTATION MAP (primary view — full width, interactive)
# =============================================================================
preview_row = (
    st.session_state.mutation_draft
    if st.session_state.show_intake_preview
    and st.session_state.mutation_draft.get("start_region")
    else None
)
map_rows = list(st.session_state.map_selection_rows)
if preview_row:
    map_rows = [preview_row] + [r for r in map_rows if r.get("id") != preview_row.get("id")]

m1, m2, m3 = st.columns([2, 2, 6])
with m1:
    st.metric("Plotted mutations", len(map_rows))
with m2:
    if preview_row:
        st.metric("Preview", "Draft on map", help="Intake form preview is shown")
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
# INTAKE FORM
# =============================================================================
with st.expander("Enter mutation data", expanded=False):
    quick_entry = st.text_input(
        "Quick entry (optional)",
        placeholder="e.g. Exonic Deletion e45-e55",
    )

    c_apply, c_clear = st.columns(2)
    if c_apply.button("Apply autofill from quick entry", type="secondary"):
        parsed = parse_quick_entry(quick_entry)
        draft = {**st.session_state.mutation_draft, **parsed}
        filled = autofill_mutation_fields(
            draft.get("mutation_class", ""),
            draft.get("mutation_subclass", ""),
            draft.get("start_region", ""),
            draft.get("stop_region", ""),
            existing=draft,
        )
        st.session_state.mutation_draft = {**draft, **filled}
        st.rerun()

    if c_clear.button("Clear form"):
        st.session_state.mutation_draft = {col: "" for col in CATALOG_COLUMNS}
        st.session_state.show_intake_preview = False
        st.rerun()

    draft = st.session_state.mutation_draft
    _phenotype_opts = PHENOTYPE_OPTIONS
    _phenotype_index = (
        _phenotype_opts.index(draft.get("phenotype", "") or "")
        if (draft.get("phenotype", "") or "") in _phenotype_opts
        else 0
    )

    with st.form("mutation_intake_form"):
        st.caption("Fill fields manually or use quick entry + autofill. Frame and kDa compute from exon regions.")
        r1c1, r1c2, r1c3 = st.columns(3)
        with r1c1:
            draft_id = st.text_input("ID", value=draft.get("id", ""), placeholder="Auto if blank")
            mutation_class = st.selectbox(
                "Mutation class",
                ["", "Exonic", "Subexonic", "Intronic"],
                index=["", "Exonic", "Subexonic", "Intronic"].index(
                    draft.get("mutation_class", "") or ""
                )
                if (draft.get("mutation_class", "") or "") in ["", "Exonic", "Subexonic", "Intronic"]
                else 0,
            )
            mutation_subclass = st.text_input(
                "Mutation subclass", value=draft.get("mutation_subclass", ""),
            )
        with r1c2:
            start_region = st.text_input("Start region", value=draft.get("start_region", ""), placeholder="e45")
            stop_region = st.text_input("Stop region", value=draft.get("stop_region", ""), placeholder="e55")
            frame = st.selectbox(
                "Frame",
                ["", "In-frame", "Out-of-frame", "Cannot determine", "N/A"],
                index=["", "In-frame", "Out-of-frame", "Cannot determine", "N/A"].index(
                    draft.get("frame", "") or ""
                )
                if (draft.get("frame", "") or "") in ["", "In-frame", "Out-of-frame", "Cannot determine", "N/A"]
                else 0,
            )
        with r1c3:
            phenotype = st.selectbox("Phenotype", _phenotype_opts, index=_phenotype_index)
            group = st.text_input("Group", value=draft.get("group", ""), placeholder="A–E")
            pct_dys = st.text_input("%Dys(WB)", value=draft.get("pct_dys_wb", ""), placeholder="e.g. 55.66")
            cdna = st.text_input("c.DNA", value=draft.get("cdna", ""))
            rna = st.text_input("r.RNA", value=draft.get("rna", ""))
            protein = st.text_input("p.Protein", value=draft.get("protein", ""))

        r2c1, r2c2 = st.columns(2)
        with r2c1:
            domains_affected = st.text_input("Domains affected", value=draft.get("domains_affected", ""))
            expected_kda = st.text_input("Expected protein size (kDa)", value=draft.get("expected_protein_size_kda", ""))
            molecular_consequence = st.text_input("Molecular consequence", value=draft.get("molecular_consequence", ""))
        with r2c2:
            exp_category = st.text_input("Experimental research category", value=draft.get("experimental_research_category", ""))
            tissue_comments = st.text_input("Tissue category / comments", value=draft.get("tissue_category_comments", ""))
            general_comments = st.text_area("General comments", value=draft.get("general_comments", ""))

        b_autofill, b_preview, b_add = st.columns(3)
        do_autofill = b_autofill.form_submit_button("Autofill computed fields")
        do_preview = b_preview.form_submit_button("Preview on map")
        do_add = b_add.form_submit_button("Add to catalog", type="primary")

        form_values = {
            "id": draft_id, "mutation_class": mutation_class,
            "mutation_subclass": mutation_subclass,
            "start_region": start_region, "stop_region": stop_region,
            "frame": frame, "phenotype": phenotype,
            "group": group, "pct_dys_wb": pct_dys,
            "cdna": cdna, "rna": rna, "protein": protein,
            "domains_affected": domains_affected,
            "expected_protein_size_kda": expected_kda,
            "molecular_consequence": molecular_consequence,
            "experimental_research_category": exp_category,
            "tissue_category_comments": tissue_comments,
            "general_comments": general_comments,
        }

    if do_autofill or do_preview or do_add:
        filled = autofill_mutation_fields(
            form_values["mutation_class"], form_values["mutation_subclass"],
            form_values["start_region"], form_values["stop_region"],
            existing={**form_values, "phenotype": phenotype},
        )
        filled["phenotype"] = phenotype
        filled["group"] = group
        filled["pct_dys_wb"] = pct_dys
        st.session_state.mutation_draft = filled

    if do_preview:
        st.session_state.show_intake_preview = True
        st.rerun()

    if do_add:
        new_row = {**st.session_state.mutation_draft}
        if not new_row.get("id"):
            new_row["id"] = next_id(catalog_df)
        catalog_df = pd.concat([catalog_df, pd.DataFrame([new_row])], ignore_index=True)
        save_catalog(catalog_df)
        st.session_state.mutation_catalog = catalog_df
        st.session_state.mutation_draft = {col: "" for col in CATALOG_COLUMNS}
        st.session_state.show_intake_preview = False
        st.success(f"Added **{new_row['id']}** to catalog.")
        st.rerun()

# =============================================================================
# FILTERS + TABLE
# =============================================================================
st.subheader("Mutation catalog")

f1, f2, f3, f4, f5 = st.columns(5)
classes = [c for c in sorted(catalog_df["mutation_class"].unique()) if c]
subclasses = [s for s in sorted(catalog_df["mutation_subclass"].unique()) if s]
frames = [f for f in sorted(catalog_df["frame"].unique()) if f]
phenotypes = [p for p in sorted(catalog_df.get("phenotype", pd.Series(dtype=str)).unique()) if p]

with f1:
    filter_class = st.multiselect("Class", classes)
with f2:
    filter_subclass = st.multiselect("Subclass", subclasses)
with f3:
    filter_frame = st.multiselect("Frame", frames)
with f4:
    filter_phenotype = st.multiselect("Phenotype", phenotypes)
with f5:
    filter_search = st.text_input("Search")

filtered = filter_catalog(
    catalog_df,
    mutation_classes=filter_class or None,
    mutation_subclasses=filter_subclass or None,
    frames=filter_frame or None,
    search_text=filter_search,
)
if filter_phenotype:
    filtered = filtered[filtered["phenotype"].isin(filter_phenotype)]

st.caption(f"**{len(filtered)}** of **{len(catalog_df)}** mutations — select rows, then plot on map above")

display_df = catalog_for_display(filtered)

edited_df = st.data_editor(
    display_df,
    use_container_width=True,
    hide_index=True,
    num_rows="fixed",
    column_config={
        "selected": st.column_config.CheckboxColumn("Select", default=False),
        "id": st.column_config.TextColumn("ID", width="small"),
        "mutation_class": st.column_config.TextColumn("Mutation class"),
        "mutation_subclass": st.column_config.TextColumn("Subclass"),
        "start_region": st.column_config.TextColumn("Start"),
        "stop_region": st.column_config.TextColumn("Stop"),
        "frame": st.column_config.TextColumn("Frame", width="small"),
        "phenotype": st.column_config.TextColumn("Phenotype"),
        "cdna": st.column_config.TextColumn("c.DNA"),
        "rna": st.column_config.TextColumn("r.RNA"),
        "protein": st.column_config.TextColumn("p.Protein"),
        "domains_affected": st.column_config.TextColumn("Domains"),
        "expected_protein_size_kda": st.column_config.TextColumn("kDa"),
        "molecular_consequence": st.column_config.TextColumn("Consequence"),
        "experimental_research_category": st.column_config.TextColumn("Research cat."),
        "tissue_category_comments": st.column_config.TextColumn("Tissue"),
        "general_comments": st.column_config.TextColumn("Comments"),
    },
    disabled=[c for c in DISPLAY_COLUMNS if c != "selected"],
    key="catalog_editor",
)

c_plot, c_clear_map = st.columns([1, 1])
if c_plot.button("Plot selected on map", type="primary"):
    if "selected" in edited_df.columns:
        selected = edited_df[edited_df["selected"] == True]  # noqa: E712
        st.session_state.map_selection_rows = selected.drop(columns=["selected"]).to_dict("records")
        st.session_state.show_intake_preview = False
        st.rerun()
    else:
        st.warning("No rows selected.")

if c_clear_map.button("Clear map selection"):
    st.session_state.map_selection_rows = []
    st.session_state.show_intake_preview = False
    st.rerun()

with st.expander("Catalog management"):
    if st.button("Re-import legacy mutations"):
        from src.mutation_catalog import import_legacy_catalog
        st.session_state.mutation_catalog = import_legacy_catalog()
        save_catalog(st.session_state.mutation_catalog)
        st.rerun()
    st.download_button(
        "Download catalog CSV",
        data=catalog_df.to_csv(index=False).encode("utf-8"),
        file_name="mutation_catalog.csv",
        mime="text/csv",
    )
