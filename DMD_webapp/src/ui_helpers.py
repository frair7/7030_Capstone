"""
Shared Streamlit UI helpers and session-state management.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import streamlit as st

from src.config import REFERENCE, reference_summary
from src.coordinate_mapper import build_exon_records, map_parsed_variant
from src.exon_skipping import find_skip_candidates
from src.frame_analysis import assess_variant_frame
from src.models import FrameStatus, InputMode, ParsedVariant, SkipCandidate, VariantType
from src.reference_data import EXONS_CSV, MAP_STYLES_CSV, PROTEIN_DOMAINS_CSV
from src.variant_parser import parse_variant
from src.visualization import ViewMode


INPUT_MODE_LABELS = {
    "Exon deletion": InputMode.EXON_DELETION,
    "Exon duplication": InputMode.EXON_DUPLICATION,
    "HGVS coding": InputMode.HGVS_CODING,
    "Genomic coordinate (GRCh38)": InputMode.GENOMIC,
    "Direct exon selection": InputMode.EXON_SELECTION,
}


@dataclass
class AnalysisBundle:
    """All analysis outputs for the current variant."""

    variant: ParsedVariant
    mapping: Any
    frame_result: Any
    skip_candidates: list[SkipCandidate]
    exons: list[Any]


def init_session_state() -> None:
    """Initialize shared session keys."""
    defaults = {
        "view_reset": 0,
        "variant_text": "del45-50",
        "input_mode_label": "Exon deletion",
        "max_additional_skips": 3,
        "view_mode_label": "Exon-order schematic",
        "selected_exons": [],
        "advanced_skip_search": False,
        "analysis_ran": False,
        "selected_candidate_idx": -1,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def plotly_config() -> dict:
    return {
        "scrollZoom": True,
        "displayModeBar": True,
        "modeBarButtonsToAdd": ["zoomIn2d", "zoomOut2d", "autoScale2d"],
    }


def render_research_warning() -> None:
    st.error(
        "**Research and educational use only.** This tool is **not** a clinical "
        "diagnostic device. Predicted exon-skipping strategies are theoretical "
        "calculations only — not evidence of therapeutic feasibility, approval, "
        "or clinical appropriateness."
    )


def render_reference_sidebar() -> None:
    st.header("Reference")
    for label, value in reference_summary().items():
        st.markdown(f"**{label}:** {value}")
    st.caption(
        f"Transcript direction on {REFERENCE.strand} strand runs opposite "
        "increasing genomic coordinates."
    )


def render_input_sidebar() -> dict[str, Any]:
    """Render shared sidebar controls; return current settings."""
    render_reference_sidebar()
    st.divider()
    st.header("Input")

    st.session_state.input_mode_label = st.selectbox(
        "Input mode",
        list(INPUT_MODE_LABELS.keys()),
        index=list(INPUT_MODE_LABELS.keys()).index(
            st.session_state.input_mode_label
        )
        if st.session_state.input_mode_label in INPUT_MODE_LABELS
        else 0,
    )
    input_mode = INPUT_MODE_LABELS[st.session_state.input_mode_label]

    selected_exons: list[int] = []
    if input_mode == InputMode.EXON_SELECTION:
        selected_exons = st.multiselect(
            "Select exons (1–79)",
            options=list(range(1, REFERENCE.coding_exon_count + 1)),
            default=st.session_state.selected_exons or [45, 46, 47, 48, 49, 50],
        )
        st.session_state.selected_exons = selected_exons
        variant_text = ""
    else:
        variant_text = st.text_input(
            "Variant input",
            value=st.session_state.variant_text,
            help="Examples: del45-50, dup2, c.5287C>T, X:31140000",
        )
        st.session_state.variant_text = variant_text

    st.header("Analysis options")
    st.session_state.max_additional_skips = st.slider(
        "Max additional skipped exons",
        min_value=0,
        max_value=5,
        value=st.session_state.max_additional_skips,
    )
    st.session_state.advanced_skip_search = st.checkbox(
        "Advanced skip search (non-contiguous)",
        value=st.session_state.advanced_skip_search,
    )
    st.session_state.view_mode_label = st.radio(
        "Map display mode",
        ["Exon-order schematic", "Transcript-scale (genomic)"],
        index=0 if st.session_state.view_mode_label.startswith("Exon-order") else 1,
    )

    col_run, col_reset = st.columns(2)
    run_clicked = col_run.button("Run analysis", type="primary", use_container_width=True)
    reset_view = col_reset.button("Reset view", use_container_width=True)
    if reset_view:
        st.session_state.view_reset += 1

    return {
        "input_mode": input_mode,
        "variant_text": variant_text,
        "selected_exons": selected_exons,
        "max_additional_skips": st.session_state.max_additional_skips,
        "advanced_skip_search": st.session_state.advanced_skip_search,
        "view_mode": (
            ViewMode.SCHEMATIC
            if st.session_state.view_mode_label.startswith("Exon-order")
            else ViewMode.GENOMIC
        ),
        "run_clicked": run_clicked,
    }


def get_view_mode() -> ViewMode:
    return (
        ViewMode.SCHEMATIC
        if st.session_state.view_mode_label.startswith("Exon-order")
        else ViewMode.GENOMIC
    )


def run_analysis_bundle(
    variant_text: str,
    input_mode: InputMode,
    *,
    selected_exons: Optional[list[int]] = None,
    max_additional_skips: int = 3,
    advanced_skip_search: bool = False,
) -> AnalysisBundle:
    """Parse input and run all analysis modules."""
    if not EXONS_CSV.exists():
        raise FileNotFoundError(
            "Reference exon table missing. Run scripts/fetch_reference_data.py"
        )

    exons = build_exon_records()
    if input_mode == InputMode.EXON_SELECTION:
        variant = parse_variant("", input_mode, selected_exons=selected_exons or [])
    else:
        variant = parse_variant(variant_text, input_mode)

    mapping = map_parsed_variant(variant, exons)
    frame_result = assess_variant_frame(variant, exons)

    skip_candidates: list[SkipCandidate] = []
    if (
        variant.is_valid
        and variant.variant_type == VariantType.DELETION
        and variant.exon_range
        and frame_result.status == FrameStatus.OUT_OF_FRAME
    ):
        first, last = variant.exon_range
        skip_candidates = find_skip_candidates(
            first,
            last,
            exons,
            max_additional_skips=max_additional_skips,
            mutation_label=variant.raw_input,
            advanced=advanced_skip_search,
        )

    return AnalysisBundle(
        variant=variant,
        mapping=mapping,
        frame_result=frame_result,
        skip_candidates=skip_candidates,
        exons=exons,
    )


def ensure_analysis(settings: dict[str, Any]) -> AnalysisBundle:
    """Run analysis when requested or return cached bundle."""
    if settings["run_clicked"] or not st.session_state.get("analysis_ran"):
        bundle = run_analysis_bundle(
            settings["variant_text"],
            settings["input_mode"],
            selected_exons=settings["selected_exons"],
            max_additional_skips=settings["max_additional_skips"],
            advanced_skip_search=settings["advanced_skip_search"],
        )
        st.session_state.analysis_bundle = bundle
        st.session_state.analysis_ran = True
        st.session_state.selected_candidate_idx = -1
    return st.session_state.analysis_bundle


def check_reference_data() -> bool:
    missing = [
        path.name
        for path in (EXONS_CSV, MAP_STYLES_CSV, PROTEIN_DOMAINS_CSV)
        if not path.exists()
    ]
    if missing:
        st.warning(
            "Reference tables not found in `reference_tables/`:\n\n"
            + ", ".join(f"`{name}`" for name in missing)
            + "\n\nRun `python scripts/fetch_reference_data.py` or see "
            "`reference_tables/README.md`."
        )
        return False
    return True
