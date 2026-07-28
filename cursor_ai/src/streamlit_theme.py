"""Widescreen Streamlit layout helpers for the DMD Explorer."""

from __future__ import annotations

import streamlit as st


def apply_widescreen_theme(*, compact_header: bool = True) -> None:
    """Inject CSS for full-width layout and readable typography."""
    padding_top = "0.75rem" if compact_header else "1.25rem"
    st.markdown(
        f"""
        <style>
          /* Full-width main content */
          .block-container {{
            padding-top: {padding_top};
            padding-bottom: 1rem;
            max-width: 98vw !important;
            padding-left: 1.5rem;
            padding-right: 1.5rem;
          }}
          /* Slightly narrower collapsed sidebar */
          [data-testid="stSidebar"][aria-expanded="false"] {{
            min-width: 3.5rem;
          }}
          [data-testid="stSidebar"] {{
            background-color: #f7f9fc;
          }}
          /* Map container */
          .dmd-map-shell {{
            width: 100%;
            overflow-x: auto;
            border: 1px solid #d8dee8;
            border-radius: 8px;
            background: #fff;
            margin-bottom: 0.5rem;
          }}
          .dmd-map-shell svg {{
            display: block;
            width: 100%;
            height: auto;
            font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
          }}
          .dmd-map-tip {{
            position: fixed;
            pointer-events: none;
            background: #1a2332;
            color: #fff;
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 13px;
            line-height: 1.45;
            z-index: 9999;
            max-width: 280px;
            box-shadow: 0 4px 14px rgba(0,0,0,0.22);
            display: none;
          }}
          /* GeneCards-inspired section headers */
          h1, h2, h3 {{
            letter-spacing: -0.01em;
          }}
          [data-testid="stMetricValue"] {{
            font-size: 1.35rem;
          }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_nav_bar(*, page_title: str) -> None:
    """Top bar with home toggle and page title."""
    home_col, title_col, _ = st.columns([1, 8, 1])
    with home_col:
        if st.button("Home", help="Return to home page", use_container_width=True):
            st.switch_page("app.py")
    with title_col:
        st.markdown(f"### {page_title}")
