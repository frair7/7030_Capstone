# Streamlit pages — `pages/`

Multi-page Streamlit modules. Navigation order follows numeric prefixes.

| File | Page |
|------|------|
| `1_Mutation_Explorer.py` | **Primary** — catalog intake, filters, interactive cohort map |
| `2_Exon_Skipping_Analysis.py` | Reading-frame skip candidate search |
| `3_Reference_Map.py` | Full exon table + Plotly reference visualization |
| `4_Methods_and_Limitations.py` | Methods, assumptions, disclaimers |

Entry point: `../app.py` (home page).

## Mutation Explorer flow

1. Map renders at top via `build_interactive_map_html()` from `src/mutation_map_svg.py`
2. Iframe height from `map_iframe_height()`; `scrolling=True` for wide map
3. User selects rows in catalog table → **Plot selected on map**
4. Optional intake **Preview on map** for draft row
5. PNG download in expander via `create_cohort_mutation_map()`

## Running

```bash
cd cursor_ai
python -m streamlit run app.py --server.port 8503
```

Open `/Mutation_Explorer` for the main map page.
