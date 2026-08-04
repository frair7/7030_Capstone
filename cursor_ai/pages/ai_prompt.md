# AI prompt — `pages/`

Folder-specific context for `cursor_ai/pages/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## Mutation Explorer (`1_Mutation_Explorer.py`)

Main UI page. Map changes usually happen in `../src/mutation_map_svg.py`, not here.

This page:

- Builds `map_rows` from `st.session_state.map_selection_rows` (+ optional draft preview)
- Calls `components.html(map_html, height=iframe_h, scrolling=True)`
- Catalog table uses `st.data_editor` with `selected` checkbox column

## Exon Skipping Analysis (`3_Exon_Skipping_Analysis.py`)

Two tabs:

### Model by exon target (`_tab_model_by_target`)

- **Skipping strategy** radio: All frame-restoring / Single-exon / Multi-exon
- **Target exon** selectbox: defaults to "Select an exon" (no silent default to exon 45)
- Calls `analyze_mutations_for_skip_target()` with `target_mode` parameter
- Target results keep the catalog `ID` in its own column; `Mutation` excludes the ID
- Each result preserves its rank from the complete mutation-specific candidate list
- Table labels that rank as `Skip combination N` and shows the full set in `Exons skipped`
- Analysis runs automatically when an exon is selected (no Analyze button)
- Empty state: warning when no frame-restoring candidates include the selected exon
- `_current_target_analyzer()` conditionally reloads stale analysis modules only
  when the in-memory Streamlit API lacks current parameters/dataclass fields

**Correct question:** "For this mutation, does skipping this exon (or combination) restore the coding frame?"

**Wrong question:** "Which mutations overlap or are associated with exon N?"

### Model by mutation (`_tab_model_by_mutation`)

- Uses `rank_skip_candidates_for_row()` → mutation-specific `find_frame_restoring_candidates()`
- Shows top 3 ranked candidates per catalog mutation
- Mutation detail with schematic overlay from `overlay_from_row_and_candidate()`
- Deletions and duplications must never share a generic flanking-exon lookup

## Common user expectation

The transcript + domain rows always render. Patient deletion bars appear only after **Plot selected on map**.

## Do not

- Duplicate map rendering logic in the page — keep it in `src/`
- Reduce iframe height below `map_iframe_height()` (clips the SVG)
- Reintroduce exon 45 as default target exon without explicit user request
- Filter mutations by target exon overlap in the page — that logic belongs in `exon_skipping_analysis.py` via `analyze_mutations_for_skip_target()`

## After map changes

Restart Streamlit and hard-refresh the browser (Ctrl+Shift+R).
