# AI prompt — `pages/`

Folder-specific context for `cursor_ai/pages/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## Mutation Explorer (`1_Mutation_Explorer.py`)

Main UI page. Map changes usually happen in `../src/mutation_map_svg.py`, not here.

This page:

- Builds `map_rows` from `st.session_state.map_selection_rows` (+ optional draft preview)
- Calls `components.html(map_html, height=iframe_h, scrolling=True)`
- Catalog table uses `st.data_editor` with `selected` checkbox column

## Common user expectation

The transcript + domain rows always render. Patient deletion bars appear only after **Plot selected on map**.

## Do not

- Duplicate map rendering logic in the page — keep it in `src/`
- Reduce iframe height below `map_iframe_height()` (clips the SVG)

## After map changes

Restart Streamlit and hard-refresh the browser (Ctrl+Shift+R).
