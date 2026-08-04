# Streamlit pages — `pages/`

Multi-page Streamlit modules. Navigation order follows numeric prefixes.

| File | Page |
|------|------|
| `1_Mutation_Explorer.py` | Interactive cohort mutation map and plot selection |
| `2_Mutation_Catalog.py` | Data entry, editing, CSV import/export, audit |
| `3_Exon_Skipping_Analysis.py` | Reading-frame skip candidate search |
| `4_Reference_Map.py` | Full exon table + reference visualization |
| `5_Methods_and_Limitations.py` | Methods, assumptions, disclaimers |

Entry point: `../app.py` (home page).

## Mutation Explorer flow

1. Map renders at top via `build_interactive_map_html()` from `src/mutation_map_svg.py`
2. User filters catalog rows and checks **Plot** in the selection table
3. **Plot selected on map** updates the interactive SVG map
4. Draft preview from Mutation Catalog appears when **Preview on map** was used there
5. PNG download in expander via `create_cohort_mutation_map()`

## Mutation Catalog flow

1. Enter mutations via quick entry or full intake form
2. Edit the catalog table inline; save, delete, or reset changes
3. Import/export CSV with validation preview and audit trail

## Running

```bash
cd cursor_ai
python -m streamlit run app.py --server.port 8503
```

Open `/Mutation_Explorer` for the map or `/Mutation_Catalog` for data management.
