"""Exon-skipping analysis — catalog-wide and target-based modeling."""

from __future__ import annotations

from typing import Any, Optional

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from src.coordinate_mapper import build_exon_records
from src.exon_skipping_catalog import (
    MUTATION_CLASS_FILTERS,
    build_catalog_options,
    catalog_row_to_map_row,
    format_exon_list,
    format_skip_candidate,
    frame_result_for_catalog_row,
    frame_status_label,
    mutation_description,
    overlay_from_row_and_candidate,
    rank_skip_candidates_for_row,
    validate_skip_target,
)
from src.mutation_catalog import active_catalog, load_catalog
from src.mutation_map_svg import (
    SKIP_COMBO_COLORS,
    SkipSchematicOverlay,
    build_interactive_map_html,
    build_skip_target_map_html,
    map_iframe_height,
)
from src.models import FrameStatus
from src.reporting import skip_candidates_dataframe
from src.streamlit_theme import apply_widescreen_theme, render_nav_bar
from src.ui_helpers import (
    check_reference_data,
    init_session_state,
    render_input_sidebar,
    render_research_warning,
    run_analysis_bundle,
)
from src.mutation_viz import sort_catalog_dataframe


def _current_target_analyzer():
    """Return the current target analyzer, refreshing a stale Streamlit import."""
    import importlib
    import inspect

    import src.exon_skipping_analysis as analysis_module
    import src.exon_skipping_catalog as catalog_module

    analyzer = catalog_module.analyze_mutations_for_skip_target
    parameters = inspect.signature(analyzer).parameters
    row_fields = getattr(
        catalog_module.RestoredMutationRow,
        "__dataclass_fields__",
        {},
    )
    required_fields = {
        "participant_id",
        "mutation_label",
        "skip_combination_rank",
    }

    if "target_mode" not in parameters or not required_fields.issubset(row_fields):
        importlib.reload(analysis_module)
        catalog_module = importlib.reload(catalog_module)
        analyzer = catalog_module.analyze_mutations_for_skip_target

    return analyzer


def _reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in df.columns if c != "mutation_record_id"]
    if "mutation_record_id" in df.columns:
        cols.append("mutation_record_id")
    return df[cols]


def _skip_combo_color_legend(candidate_count: int) -> None:
    """Show only skip-combo colors for candidates that exist."""
    if candidate_count <= 0:
        return
    parts = [
        f'<span style="color:{SKIP_COMBO_COLORS[i]};font-weight:700">'
        f"Skip combo {i + 1}</span>"
        for i in range(min(candidate_count, 3))
    ]
    st.markdown("&nbsp;&nbsp;|&nbsp;&nbsp;".join(parts), unsafe_allow_html=True)


FRAME_RESTORING_CAPTION = (
    "Only candidate strategies predicted to restore the reading frame are shown. "
    "Candidates are ranked by the fewest additional skipped exons, then by the "
    "least additional coding sequence removed."
)

NO_FRAME_CANDIDATES_MSG = (
    "No frame-restoring exon-skipping candidates were identified within the "
    "selected search limit."
)


def _colored_combo_label(combo_index: int) -> str:
    """HTML label for skip combo 1–3 with matching schematic color."""
    color = SKIP_COMBO_COLORS[combo_index]
    return (
        f'<span style="color:{color};font-weight:700">Skip combo {combo_index + 1}</span>'
    )


def _render_schematic_html(
    rows: list[dict[str, Any]],
    *,
    title: str,
    overlay: Optional[SkipSchematicOverlay] = None,
) -> None:
    html = build_interactive_map_html(
        rows,
        chart_width=2000,
        title=title,
        skip_overlay=overlay,
    )
    height = map_iframe_height(max(len(rows), 1))
    components.html(html, height=height, scrolling=True)


def _render_mutation_detail(
    row: dict[str, Any],
    *,
    exons: list[Any],
    max_additional_skips: int,
    advanced: bool,
    selected_candidate_idx: int = 0,
) -> None:
    """Shared detail view: metrics, candidates, schematic, assumptions."""
    frame, candidates = rank_skip_candidates_for_row(
        row,
        exons,
        max_additional_skips=max_additional_skips,
        advanced=advanced,
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("Original frame", frame_status_label(frame.status))
    resulting = "In Frame" if candidates else frame_status_label(frame.status)
    c2.metric("Resulting frame (best candidate)", resulting if candidates else "—")
    c3.metric("Candidates found", len(candidates))

    st.markdown(frame.explanation)

    if not candidates:
        if frame.status == FrameStatus.IN_FRAME:
            st.success("This mutation is already in frame; no skip strategy is required.")
        elif frame.status == FrameStatus.CANNOT_DETERMINE:
            st.warning("Frame status could not be determined for this mutation.")
        else:
            st.warning(NO_FRAME_CANDIDATES_MSG)
        return

    st.caption(FRAME_RESTORING_CAPTION)
    _skip_combo_color_legend(len(candidates))
    if len(candidates) > 1:
        st.caption(
            f"Skip combos 1–{len(candidates)} are outlined on the schematic "
            f"(red, green, blue)."
        )
    else:
        st.caption("Skip combo 1 is outlined on the schematic (red).")

    idx = st.selectbox(
        "Select candidate for junction emphasis",
        range(len(candidates)),
        index=min(selected_candidate_idx, len(candidates) - 1),
        format_func=lambda i: (
            f"Skip combo {i + 1}: {format_skip_candidate(candidates[i])} → "
            f"{candidates[i].junction_display()}"
        ),
        key=f"detail_candidate_{row.get('mutation_record_id', 'x')}",
    )
    selected = candidates[idx]

    header_cols = st.columns(min(len(candidates), 3))
    for i in range(len(header_cols)):
        with header_cols[i]:
            st.markdown(
                f'{_colored_combo_label(i)}: {format_skip_candidate(candidates[i])}',
                unsafe_allow_html=True,
            )

    m1, m2, m3 = st.columns(3)
    m1.metric("Additional exons skipped", len(selected.additional_skipped_exons))
    m2.metric("Additional bp removed", selected.additional_coding_bases_removed)
    m3.metric("Est. protein length (aa)", selected.estimated_protein_aa)

    cand_df = skip_candidates_dataframe(candidates).copy()
    cand_df.index = [f"Skip combo {i + 1}" for i in range(len(cand_df))]
    st.dataframe(cand_df, use_container_width=True)

    overlay = overlay_from_row_and_candidate(
        row, selected, all_candidates=candidates[:3],
    )
    _render_schematic_html(
        [catalog_row_to_map_row(row)],
        title=f"Skip strategy — {mutation_description(row)}",
        overlay=overlay,
    )

    with st.expander("Calculation explanation and assumptions", expanded=False):
        for item in selected.assumptions:
            st.markdown(f"- {item}")
        st.markdown(frame.explanation)
    if selected.affected_domains:
        st.markdown("**Affected protein domains (broad):** " + ", ".join(selected.affected_domains))


def _tab_model_by_target(catalog_df: pd.DataFrame, exons: list[Any]) -> None:
    st.markdown(
        "Select a proposed exon-skipping **target** and see which catalog mutations "
        "would be restored to in-frame if that target were skipped."
    )

    mode_label = st.radio(
        "Skipping strategy",
        [
            "All frame-restoring candidates",
            "Single-exon skipping",
            "Multi-exon skipping",
        ],
        horizontal=True,
        key="target_mode",
    )

    mode_lookup = {
        "All frame-restoring candidates": "all",
        "Single-exon skipping": "single",
        "Multi-exon skipping": "multi",
    }
    target_mode = mode_lookup[mode_label]

    n_exons = len(exons)
    target_exon_options = [None] + list(range(1, n_exons + 1))
    selected_target = st.selectbox(
        "Target exon",
        options=target_exon_options,
        index=0,
        format_func=lambda value: (
            "Select an exon"
            if value is None
            else f"Exon {value}"
        ),
        key="target_single_exon",
    )

    if selected_target is None:
        st.info(
            "Select an exon to identify mutations for which targeting that "
            "exon is part of a mathematically frame-restoring strategy."
        )
        return

    target_start = int(selected_target)
    target_end = target_start
    target_exons, err = validate_skip_target(target_start, target_end)
    if err:
        st.error(err)
        return

    st.subheader("Selected skip target")
    st.caption(format_exon_list(target_exons))
    components.html(
        build_skip_target_map_html(target_exons, title="Proposed skip target"),
        height=map_iframe_height(0),
        scrolling=True,
    )

    fc1, fc2 = st.columns([2, 1])
    with fc1:
        class_filter = st.multiselect(
            "Filter by mutation class",
            list(MUTATION_CLASS_FILTERS.keys()),
            key="target_class_filter",
        )
    with fc2:
        include_in_frame = st.checkbox(
            "Include mutations already in frame",
            value=False,
            key="target_include_in_frame",
        )

    restored = _current_target_analyzer()(
        catalog_df,
        target_exons,
        exons,
        include_in_frame=include_in_frame,
        class_filter=set(class_filter) if class_filter else None,
        target_mode=target_mode,
    )

    st.subheader("Mutations restored by this target")
    if not restored:
        st.warning(
            f"No {target_mode} frame-restoring candidates include "
            f"{format_exon_list(target_exons)}."
        )
        return

    table_rows = []
    for item in restored:
        table_rows.append({
            "ID": item.participant_id,
            "Mutation": item.mutation_label,
            "Class": item.mutation_class,
            "Original frame": item.original_frame,
            "Resulting frame": item.resulting_frame,
            "Skip combination": f"Skip combination {item.skip_combination_rank}",
            "Exons skipped": format_exon_list(item.skipped_exons),
            "mutation_record_id": item.mutation_record_id,
        })
    result_df = _reorder_columns(pd.DataFrame(table_rows))

    event = st.dataframe(
        result_df.drop(columns=["mutation_record_id"]),
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        key="target_results_table",
    )

    selected_idx = None
    if event and event.selection and event.selection.rows:
        selected_idx = event.selection.rows[0]

    if selected_idx is not None:
        item = restored[selected_idx]
        st.subheader("Mutation detail")
        frame, _ = frame_result_for_catalog_row(item.row, exons)
        c1, c2 = st.columns(2)
        c1.metric("Original frame", item.original_frame)
        c2.metric("Resulting frame", item.resulting_frame)
        overlay = overlay_from_row_and_candidate(
            item.row, item.candidate, all_candidates=[item.candidate],
        )
        _render_schematic_html(
            [catalog_row_to_map_row(item.row)],
            title=f"Restored by skipping {format_exon_list(item.skipped_exons)}",
            overlay=overlay,
        )
        with st.expander("Calculation explanation and assumptions", expanded=False):
            for note in item.candidate.assumptions:
                st.markdown(f"- {note}")
            st.markdown(frame.explanation)
        if item.candidate.affected_domains:
            st.markdown(
                "**Affected protein domains (broad):** "
                + ", ".join(item.candidate.affected_domains)
            )


def _tab_model_by_mutation(
    catalog_df: pd.DataFrame,
    exons: list[Any],
    settings: dict[str, Any],
) -> None:
    st.markdown(
        "Start from a mutation and review ranked skip-candidate strategies, "
        "or drill into a single mutation's detail view."
    )

    source = st.radio(
        "Mutation source",
        ["Current Mutation Explorer entry", "Mutation catalog"],
        horizontal=True,
        key="mutation_source",
    )

    detail_row: Optional[dict[str, Any]] = None

    if source == "Current Mutation Explorer entry":
        bundle = run_analysis_bundle(
            settings["variant_text"],
            settings["input_mode"],
            selected_exons=settings["selected_exons"],
            max_additional_skips=settings["max_additional_skips"],
            advanced_skip_search=settings["advanced_skip_search"],
        )
        variant = bundle.variant
        if not variant.is_valid or not variant.exon_range:
            st.info("Enter a valid whole-exon variant in the sidebar and click **Run analysis**.")
            return
        msub = "Deletion" if variant.variant_type.value == "deletion" else "Duplication"
        detail_row = {
            "id": variant.raw_input,
            "mutation_class": "Exonic",
            "mutation_subclass": msub,
            "start_region": f"e{variant.first_exon}",
            "stop_region": f"e{variant.last_exon}",
            "mutation_record_id": "__explorer_entry__",
        }
    else:
        options = build_catalog_options(catalog_df)
        if not options:
            st.warning("No active mutations in the catalog.")
            return
        labels = [o[0] for o in options]
        pick = st.selectbox("Select catalog mutation", labels, key="catalog_mutation_pick")
        rec_id = next(rid for lbl, rid in options if lbl == pick)
        active = active_catalog(catalog_df)
        match = active[active["mutation_record_id"] == rec_id]
        if not match.empty:
            detail_row = match.iloc[0].to_dict()

    st.subheader("Top skip candidates per mutation")
    st.caption(FRAME_RESTORING_CAPTION)

    active = sort_catalog_dataframe(active_catalog(catalog_df))
    summary_rows: list[dict[str, Any]] = []
    for _, raw in active.iterrows():
        row = raw.to_dict()
        frame, candidates = rank_skip_candidates_for_row(
            row,
            exons,
            max_additional_skips=settings["max_additional_skips"],
            advanced=settings["advanced_skip_search"],
        )
        entry = {
            "Participant": row.get("id", ""),
            "Mutation": mutation_description(row),
            "Class": row.get("mutation_subclass", ""),
            "Original frame": frame_status_label(frame.status),
            "Skip candidate #1": format_skip_candidate(candidates[0]) if len(candidates) > 0 else "—",
            "Skip candidate #2": format_skip_candidate(candidates[1]) if len(candidates) > 1 else "—",
            "Skip candidate #3": format_skip_candidate(candidates[2]) if len(candidates) > 2 else "—",
            "mutation_record_id": row.get("mutation_record_id", ""),
        }
        summary_rows.append(entry)

    summary_df = _reorder_columns(pd.DataFrame(summary_rows))
    display_df = summary_df.drop(columns=["mutation_record_id"]).copy()
    display_df.columns = [
        "Participant",
        "Mutation",
        "Class",
        "Original frame",
        f"Skip combo 1",
        f"Skip combo 2",
        f"Skip combo 3",
    ]
    event = st.dataframe(
        display_df,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        key="mutation_summary_table",
    )

    selected_summary_row: Optional[dict[str, Any]] = None
    if event and event.selection and event.selection.rows:
        idx = event.selection.rows[0]
        rec_id = summary_rows[idx]["mutation_record_id"]
        match = active[active["mutation_record_id"] == rec_id]
        if not match.empty:
            selected_summary_row = match.iloc[0].to_dict()

    st.subheader("Mutation detail")
    detail_source = selected_summary_row or detail_row
    if detail_source:
        _render_mutation_detail(
            detail_source,
            exons=exons,
            max_additional_skips=settings["max_additional_skips"],
            advanced=settings["advanced_skip_search"],
        )
    else:
        st.caption("Select a row in the table above or choose a catalog mutation.")


def main() -> None:
    init_session_state()
    apply_widescreen_theme()
    render_nav_bar(page_title="Exon Skipping Analysis")
    render_research_warning()

    if not check_reference_data():
        st.stop()

    if "mutation_catalog" not in st.session_state:
        st.session_state.mutation_catalog = load_catalog()

    catalog_df = st.session_state.mutation_catalog
    exons = build_exon_records()

    with st.sidebar:
        settings = render_input_sidebar()

    tab_target, tab_mutation = st.tabs(["Model by exon target", "Model by mutation"])

    with tab_target:
        _tab_model_by_target(catalog_df, exons)

    with tab_mutation:
        _tab_model_by_mutation(catalog_df, exons, settings)

    st.caption(
        "Evidence class: computational frame restoration only. "
        "Not an approved or experimentally validated therapy."
    )


if __name__ == "__main__":
    main()
