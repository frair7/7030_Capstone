"""Mutation Catalog — data entry, editing, import/export, and audit."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from src.coordinate_mapper import build_exon_records
from src.exon_display_config import PHENOTYPE_OPTIONS
from src.mutation_autofill import autofill_mutation_fields
from src.mutation_catalog import (
    CATALOG_COLUMNS,
    CATALOG_EDITOR_COLUMNS,
    active_catalog,
    catalog_for_editor,
    ensure_record_ids,
    load_catalog,
    new_mutation_record_id,
    next_id,
    save_catalog,
)
from src.mutation_catalog_ops import (
    apply_import_preview,
    apply_table_preview,
    export_catalog_csv,
    load_audit_for_record,
    preview_csv_import,
    preview_table_edit,
    soft_delete_rows,
    validation_errors_csv,
)
from src.region_parser import parse_quick_entry
from src.streamlit_theme import apply_widescreen_theme, render_nav_bar
from src.ui_helpers import check_reference_data, init_session_state, render_research_warning

init_session_state()
apply_widescreen_theme()
render_nav_bar(page_title="Mutation Catalog")

if "show_intake_preview" not in st.session_state:
    st.session_state.show_intake_preview = False
if "catalog_table_draft" not in st.session_state:
    st.session_state.catalog_table_draft = None
if "pending_import_preview" not in st.session_state:
    st.session_state.pending_import_preview = None
if "pending_delete_ids" not in st.session_state:
    st.session_state.pending_delete_ids = []
if "pending_table_preview" not in st.session_state:
    st.session_state.pending_table_preview = None

render_research_warning()

if not check_reference_data():
    st.stop()

if "mutation_catalog" not in st.session_state:
    st.session_state.mutation_catalog = load_catalog()
if "mutation_draft" not in st.session_state:
    st.session_state.mutation_draft = {col: "" for col in CATALOG_COLUMNS}

catalog_df = st.session_state.mutation_catalog
_ = build_exon_records()

if "phenotype" not in catalog_df.columns:
    catalog_df["phenotype"] = ""
    for col in CATALOG_COLUMNS:
        if col not in catalog_df.columns:
            catalog_df[col] = ""

active_df = active_catalog(catalog_df)

st.markdown(
    "Add, edit, import, and export mutation records. "
    "Use **Mutation Explorer** to visualize selected mutations on the transcript map."
)

# =============================================================================
# INTAKE FORM
# =============================================================================
st.subheader("Enter mutation data")

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
    st.info("Draft saved. Open **Mutation Explorer** to view the preview on the map.")
    st.rerun()

if do_add:
    new_row = {**st.session_state.mutation_draft}
    if not new_row.get("id"):
        new_row["id"] = next_id(catalog_df)
    new_row["mutation_record_id"] = new_mutation_record_id()
    catalog_df = pd.concat([catalog_df, pd.DataFrame([new_row])], ignore_index=True)
    catalog_df, _ = ensure_record_ids(catalog_df)
    save_catalog(catalog_df)
    st.session_state.mutation_catalog = catalog_df
    st.session_state.catalog_table_draft = None
    st.session_state.mutation_draft = {col: "" for col in CATALOG_COLUMNS}
    st.session_state.show_intake_preview = False
    st.success(f"Added **{new_row['id']}** to catalog.")
    st.rerun()

# =============================================================================
# IMPORT / EXPORT
# =============================================================================
st.subheader("Import and export")

dm1, dm2 = st.columns(2)
with dm1:
    st.download_button(
        "Download Mutation Table CSV",
        data=export_catalog_csv(catalog_df),
        file_name="mutation_catalog.csv",
        mime="text/csv",
    )
with dm2:
    st.caption("Export includes all active records with permanent `mutation_record_id` values.")

st.markdown("#### Upload edited mutation table")
uploaded = st.file_uploader("Upload edited CSV", type=["csv"], key="mutation_csv_upload")
if uploaded is not None:
    try:
        upload_df = pd.read_csv(uploaded, dtype=str).fillna("")
        st.session_state.pending_import_preview = preview_csv_import(catalog_df, upload_df)
    except Exception as exc:
        st.error(f"Could not read CSV: {exc}")

if st.session_state.pending_import_preview is not None:
    preview = st.session_state.pending_import_preview
    st.markdown("#### Import preview")
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Total", preview.total_rows)
    c2.metric("New", preview.new_rows)
    c3.metric("Updated", preview.updated_rows)
    c4.metric("Unchanged", preview.unchanged_rows)
    c5.metric("Invalid", preview.invalid_rows)
    c6.metric("Unknown IDs", preview.unknown_ids)
    if preview.blocking_errors:
        st.error("; ".join(preview.blocking_errors))
    updated_rows = [r for r in preview.rows if r.status == "updated" and r.field_changes]
    if updated_rows:
        change_rows = []
        for item in updated_rows:
            for ch in item.field_changes:
                change_rows.append({
                    "mutation_record_id": item.mutation_record_id,
                    "field": ch["field"],
                    "existing": ch["existing"],
                    "uploaded": ch["uploaded"],
                })
        st.dataframe(pd.DataFrame(change_rows), use_container_width=True, hide_index=True)
    ic1, ic2, ic3 = st.columns(3)
    if ic1.button("Apply CSV import", type="primary", disabled=not preview.can_apply):
        try:
            updated = apply_import_preview(catalog_df, preview, source="csv_import")
            save_catalog(updated)
            st.session_state.mutation_catalog = updated
            st.session_state.catalog_table_draft = None
            st.session_state.pending_import_preview = None
            st.success("CSV import applied.")
            st.rerun()
        except ValueError as exc:
            st.error(str(exc))
    if ic2.button("Cancel CSV import"):
        st.session_state.pending_import_preview = None
        st.rerun()
    if preview.invalid_rows > 0:
        ic3.download_button(
            "Download validation report",
            data=validation_errors_csv(preview),
            file_name="mutation_import_errors.csv",
            mime="text/csv",
            key="import_errors_download",
        )

# =============================================================================
# EDITABLE TABLE
# =============================================================================
st.subheader("Edit mutation table")

if st.session_state.catalog_table_draft is None:
    st.session_state.catalog_table_draft = active_df.copy()

editor_df = catalog_for_editor(st.session_state.catalog_table_draft)

_readonly = [
    "mutation_record_id", "is_deleted", "deleted_at", "deleted_by", "deletion_reason",
]

column_config = {
    "delete_row": st.column_config.CheckboxColumn("Delete", default=False),
    "mutation_record_id": st.column_config.TextColumn(
        "Record ID", disabled=True, width="small",
    ),
    "id": st.column_config.TextColumn("Participant ID", width="small"),
    "mutation_class": st.column_config.SelectboxColumn(
        "Mutation class", options=["", "Exonic", "Subexonic", "Intronic"],
    ),
    "mutation_subclass": st.column_config.TextColumn("Subclass"),
    "start_region": st.column_config.TextColumn("Start"),
    "stop_region": st.column_config.TextColumn("Stop"),
    "frame": st.column_config.SelectboxColumn(
        "Frame", options=["", "In-frame", "Out-of-frame", "Cannot determine", "N/A"],
    ),
    "phenotype": st.column_config.SelectboxColumn("Phenotype", options=PHENOTYPE_OPTIONS),
    "is_deleted": st.column_config.TextColumn("Deleted", disabled=True),
    "deleted_at": st.column_config.TextColumn("Deleted at", disabled=True),
    "deleted_by": st.column_config.TextColumn("Deleted by", disabled=True),
    "deletion_reason": st.column_config.TextColumn("Deletion reason", disabled=True),
}

edited_df = st.data_editor(
    editor_df,
    use_container_width=True,
    hide_index=True,
    num_rows="dynamic",
    column_config=column_config,
    column_order=CATALOG_EDITOR_COLUMNS,
    disabled=_readonly + ["mutation_record_id"],
    key="mutation_catalog_editor",
)


def _editor_to_catalog(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in df.columns if c in CATALOG_COLUMNS]
    return df[cols].fillna("").astype(str)


edited_catalog = _editor_to_catalog(edited_df)
st.session_state.catalog_table_draft = edited_catalog.copy()

_stored_active = active_catalog(catalog_df)[CATALOG_COLUMNS].reset_index(drop=True)
_has_unsaved = not edited_catalog.reset_index(drop=True).equals(_stored_active)
if _has_unsaved and st.session_state.pending_table_preview is None:
    st.warning("You have unsaved table changes. Save or reset before leaving this page.")

sc1, sc2, sc3 = st.columns(3)
if sc1.button("Save mutation table changes", type="primary"):
    st.session_state.pending_table_preview = preview_table_edit(catalog_df, edited_catalog)
    st.rerun()
if sc2.button("Reset unsaved changes"):
    st.session_state.catalog_table_draft = active_catalog(catalog_df).copy()
    st.session_state.pending_table_preview = None
    st.rerun()

delete_rows = edited_df[edited_df.get("delete_row", False) == True]  # noqa: E712
if sc3.button("Delete selected mutation rows", disabled=delete_rows.empty):
    st.session_state.pending_delete_ids = delete_rows["mutation_record_id"].tolist()
    st.rerun()

if st.session_state.pending_delete_ids:
    ids = st.session_state.pending_delete_ids
    to_delete = active_df[active_df["mutation_record_id"].isin(ids)]
    st.warning(f"Delete **{len(ids)}** mutation record(s)?")
    st.dataframe(
        to_delete[["mutation_record_id", "id", "start_region", "stop_region", "mutation_subclass"]],
        use_container_width=True,
        hide_index=True,
    )
    dc1, dc2 = st.columns(2)
    if dc1.button("Confirm deletion", type="primary"):
        updated = soft_delete_rows(catalog_df, ids, reason="user confirmed table delete")
        save_catalog(updated)
        st.session_state.mutation_catalog = updated
        st.session_state.catalog_table_draft = None
        st.session_state.pending_delete_ids = []
        st.success(f"Deleted {len(ids)} mutation record(s).")
        st.rerun()
    if dc2.button("Cancel deletion"):
        st.session_state.pending_delete_ids = []
        st.rerun()

if st.session_state.pending_table_preview is not None:
    tpreview = st.session_state.pending_table_preview
    st.markdown("#### Pending table changes")
    tc1, tc2, tc3, tc4, tc5 = st.columns(5)
    tc1.metric("New", tpreview.new_rows)
    tc2.metric("Updated", tpreview.updated_rows)
    tc3.metric("Deleted", tpreview.deleted_rows)
    tc4.metric("Unchanged", tpreview.unchanged_rows)
    tc5.metric("Invalid", tpreview.invalid_rows)
    if tpreview.blocking_errors:
        st.error("; ".join(tpreview.blocking_errors))
    tc_apply, tc_cancel = st.columns(2)
    if tc_apply.button("Confirm save", type="primary", disabled=not tpreview.can_apply):
        try:
            updated = apply_table_preview(catalog_df, tpreview, source="editable_table")
            save_catalog(updated)
            st.session_state.mutation_catalog = updated
            st.session_state.catalog_table_draft = None
            st.session_state.pending_table_preview = None
            st.success("Mutation table saved.")
            st.rerun()
        except ValueError as exc:
            st.error(str(exc))
    if tc_cancel.button("Cancel save"):
        st.session_state.pending_table_preview = None
        st.rerun()

with st.expander("Audit history for selected record"):
    rec_options = active_df["mutation_record_id"].tolist()
    if rec_options:
        pick_rec = st.selectbox("mutation_record_id", rec_options)
        audit = load_audit_for_record(pick_rec)
        if audit:
            st.json(audit)
        else:
            st.caption("No audit entries for this record yet.")

with st.expander("Catalog management"):
    if st.button("Re-import legacy mutations"):
        from src.mutation_catalog import import_legacy_catalog
        st.session_state.mutation_catalog = import_legacy_catalog()
        save_catalog(st.session_state.mutation_catalog)
        st.session_state.catalog_table_draft = None
        st.rerun()
