"""Widescreen Streamlit layout helpers for the DMD Explorer."""

from __future__ import annotations

import streamlit as st


def apply_widescreen_theme(*, compact_header: bool = True) -> None:
    """Inject CSS for full-width layout, compact top spacing, and readable typography."""
    padding_top = "0.5rem" if compact_header else "1rem"
    st.markdown(
        f"""
        <style>
          /* Keep page content inside the visible viewport (no clipped titles). */
          header[data-testid="stHeader"] {{
            display: none !important;
            height: 0 !important;
            min-height: 0 !important;
          }}
          [data-testid="stToolbar"] {{
            display: none !important;
          }}
          [data-testid="stAppViewContainer"] {{
            margin-top: 0 !important;
            padding-top: 0 !important;
          }}
          section[data-testid="stMain"] > div {{
            padding-top: 0 !important;
            gap: 0 !important;
          }}
          .block-container {{
            padding-top: {padding_top} !important;
            padding-bottom: 0.75rem !important;
            max-width: 98vw !important;
            padding-left: 1.5rem !important;
            padding-right: 1.5rem !important;
            margin-top: 0 !important;
          }}
          /* Tighten vertical spacing between Streamlit blocks. */
          [data-testid="stVerticalBlockBorderWrapper"] {{
            padding-top: 0.15rem !important;
            padding-bottom: 0.15rem !important;
          }}
          /* Page and section titles — pull into view, minimal top margin. */
          h1, h2, h3, h4,
          [data-testid="stMarkdownContainer"] h1,
          [data-testid="stMarkdownContainer"] h2,
          [data-testid="stMarkdownContainer"] h3 {{
            margin-top: 0.15rem !important;
            margin-bottom: 0.35rem !important;
            padding-top: 0 !important;
            line-height: 1.25 !important;
          }}
          .dmd-page-title {{
            margin: 0 !important;
            padding: 0 !important;
            font-size: 1.35rem !important;
            font-weight: 700 !important;
            line-height: 1.2 !important;
            letter-spacing: -0.01em;
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
          [data-testid="stMetricValue"] {{
            font-size: 1.35rem;
          }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_nav_bar(*, page_title: str) -> None:
    """Compact top bar with home toggle and page title."""
    home_col, title_col, _ = st.columns([1, 8, 1])
    with home_col:
        if st.button("Home", help="Return to home page", use_container_width=True):
            st.switch_page("app.py")
    with title_col:
        st.markdown(
            f'<p class="dmd-page-title">{page_title}</p>',
            unsafe_allow_html=True,
        )
