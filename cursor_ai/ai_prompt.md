# AI context prompt — DMD Mutation and Exon-Skipping Explorer

Use this file to resume AI-assisted work on the capstone without re-explaining the project.

## Project summary

**BSGP 7030 Capstone** — locally hosted **Streamlit** web app for research/educational exploration of *DMD* / Dp427m mutations (NM_004006.3, GRCh38, 79 coding exons). **Not for clinical use.**

Active codebase: **`cursor_ai/`** subfolder only. Parent repo (`7030_Capstone/`) has older scripts and `data/input/DMD_Mutations_clean.csv`.

## Environment (OSC)

```bash
module load miniconda3/24.1.2-py310
conda activate 7030_capstone
conda env update --name 7030_capstone --file cursor_ai/environment.yml   # do NOT use --prune
cd ~/7030_Capstone/cursor_ai
export MPLCONFIGDIR=~/7030_Capstone/cursor_ai/.mplconfig
python -m streamlit run app.py --server.port 8503 --server.headless true
```

**SSH tunnel from laptop:**

```bash
ssh -N -L 8503:localhost:8503 frair7@ascend.osc.edu
```

Then open **http://localhost:8503/Mutation_Explorer**

Dependencies: streamlit, pandas, numpy, plotly, matplotlib, biopython, pytest, pydantic.

## Development phases

| Phase | Status | Notes |
|-------|--------|-------|
| 1 Scaffold | Done | `app.py`, `environment.yml`, `src/config.py` |
| 2 Reference data | Done | `data/dmd_exons_grch38.csv`, `data/dmd_domains.csv` |
| 3 Core analysis | Done | variant parser, coordinate mapper, frame analysis, exon skipping |
| 4 Visualization | Done | Plotly transcript map in `src/visualization.py` |
| 5 Multi-page app | Done | 4 Streamlit pages, reporting downloads |
| 6 Validation / polish | In progress | Cohort map styling complete; student review pending |

## App pages

1. **Mutation Explorer** (`pages/1_Mutation_Explorer.py`) — primary focus; catalog intake, filters, interactive SVG cohort map
2. **Exon Skipping Analysis** (`pages/3_Exon_Skipping_Analysis.py`) — two tabs: Model by exon target, Model by mutation
3. Reference Map
4. Methods and Limitations

## Mutation Explorer (current design)

### Catalog

- Legacy mutations in `data/mutation_catalog.csv`
- 15 catalog columns + `phenotype` (Healthy control, Asymptomatic, Paucisymptomatic, IMD, BMD, DMD)
- Intake form with quick entry autofill (`e45-e55` style)
- Filterable table with row checkboxes → **Plot selected on map**
- Optional map fields: `group`, `pct_dys_wb`

### Interactive cohort map (`src/mutation_map_svg.py`)

Rendered via `components.html` with fixed pixel SVG dimensions and horizontal scroll.

**Transcript row (3 layers, shared x-coordinates from coding position):**

1. Continuous domain color band (`build_part_segments`)
2. Puzzle-clipped fills into bump regions (`segments_for_exon` — no neighbor bleed)
3. Per-segment puzzle outlines using each part's `border` color

**Domain row (protein icons):**

- Rounded bubbles from `resolve_domains()` in `src/domain_alignment.py`
- Fill and border strips use the **same x-positions** as the transcript row
- Domain colors (fill / border):
  - ABD `#FFCCCC` / `#EB6F6F`
  - Hinges H1–H4 `#AEDCCA` / `#4BAF89`
  - Repeats R1–R24 `#FFFFCC` / `#F0C266`
  - CR Domain `#DAE9F8` / `#4D93D9`
  - C-terminal `#F0E3F4` / `#D86DCD`

**Layout:**

- Left metadata panel ends at `TABLE_RIGHT`; track starts at `TRACK_LEFT` (205 px)
- Table columns left of track: Group | Participant | MW | %Dys(WB)
- 12 px white gap between transcript and domain rows
- Mutation bars: white margin on all sides (no gray flanking boxes)
- 8 dashed hinge guides (start + end of H1–H4), drawn downward from domain row only

**Static PNG export:** `src/mutation_map_figure.py` via `create_cohort_mutation_map()`

### Authoritative map data

| File | Role |
|------|------|
| `scripts/reference_make_figure.py` | Authoritative `EXON_TABLE` + `DOMAIN_MAP` |
| `src/dp427m_exon_data.py` | Loads constants from reference script (`@lru_cache` on mtime) |
| `src/part_segments.py` | `PartSegment`, `build_part_segments`, `find_part_segments`, `segments_for_exon` |
| `src/puzzle_geometry.py` | Puzzle paths; junction edges use `M` (move), not `L`, when shared |
| `src/domain_alignment.py` | `resolve_domains`, `hinge_guide_x_positions` |
| `src/exon_shape_renderer.py` | Matplotlib transcript/domain row rendering |
| `src/mutation_map_svg.py` | Interactive HTML/SVG for Streamlit |
| `src/mutation_map_figure.py` | Matplotlib cohort figure for PNG export |

**Note:** `reference_tables/dmd_exon_map_styles.csv` (via `src/reference_data.py`) is a parallel CSV source used by Reference Map — **not** what feeds the cohort map. The cohort map uses `scripts/reference_make_figure.py`.

### Key DOMAIN_MAP fix (exon 17)

- R3: exons 15–17, but only the 5′ yellow sliver in exon 17 (via `find_part_segments`)
- H2: exon 17 only (interior green segment)
- R4: exons 18–20

### Multi-part junction exons

| Exon | Split |
|------|-------|
| 8 | ABD → H1 |
| 10 | H1 → R1 |
| 17 | R3 → H2 → R3 |
| 50 | R19 → H3 |
| 51 | H3 → R20 |
| 61 | R24 → H4 |
| 64 | H4 → CR Domain |

## Exon skipping analysis (critical architecture)

### Model by mutation tab

Uses the mutation-specific engine in `src/exon_skipping_analysis.py`:

- `rank_skip_candidates_for_row()` builds the correct deletion or duplication model
- It calls `find_frame_restoring_candidates()` and returns the top three
- `_legacy_candidate_from_analysis()` adapts results for the existing schematic UI
- Never send duplications through deletion-oriented `find_skip_candidates()`

### Model by exon target tab (mutation-specific simulation)

Uses **`src/exon_skipping_analysis.py`** — authoritative for candidate generation and target-exon filtering:

```
Mutation → generate adjacent candidate sets → evaluate coding delta →
discard non-restoring → rank → apply single/multi filter → apply target-exon filter → display
```

**Core rule:** The selected target exon is a **display filter only**. Never use it to decide whether a mutation is frame-restoring.

**Frame math (coding bp only, from `dmd_exons_grch38.csv`):**

```
Deletion:  mutation_delta = -(deleted coding bases)
Duplication: mutation_delta = +(duplicated coding bases)
For each skip set: final_delta = mutation_delta - sum(skipped coding bases)
Frame restored when: final_delta % 3 == 0
```

**Candidate generation rules:**

| Mutation class | Rules |
|----------------|-------|
| Deletion | Deleted exons cannot be skip targets. Generate contiguous upstream/downstream blocks bordering the deletion. Nonconsecutive exon numbers allowed across a deletion (e.g. del7 → skip 6,8). |
| Duplication | Duplicated interval is the base skip target. Extend with contiguous upstream/downstream blocks. |

`candidate_is_allowed()` enforces these rules even for directly evaluated
skip sets. In particular, duplication e3–e7 must reject skip 8–9 because the
candidate omits exons 3–7. Deletion e3–e7 may retain skip 8–9 only after its
coding delta passes modulo 3.

**Catalog bridge:** `analyze_mutations_for_skip_target()` in `src/exon_skipping_catalog.py` calls `find_frame_restoring_candidates()` per mutation, then `filter_candidates_by_target_exons()`.

**Target-results table:** Keep the participant `ID` separate from the mutation
label. Preserve each candidate's rank from the complete pre-filter candidate
list and display it as `Skip combination N`, alongside the full `Exons skipped`
set. Target filtering must not renumber the underlying mutation-specific rank.

**Validated examples (see `tests/test_exon_skipping_analysis.py`):**

| Mutation | Expected candidates |
|----------|---------------------|
| Dup exon 2 | Skip exon 2 (single) |
| Dup exons 3–4 | Skip 3,4 / 3,4,5 / 3,4,5,6,7,8 (no single-exon) |
| Del exon 7 | Skip 6,8 / 6,8,9 / 5,6,8 / 2,3,4,5,6 (no exon 7 in targets) |
| Del exon 45 | Skip exon 44 (single) |
| Del exons 3–7 + target exon 45 | Must NOT appear |

**Do not:**

- Filter mutations by whether they overlap a target exon range
- Default silently to exon 45 in the target tab UI
- Use `reconstruct_deletion_skip()` with user target exons as the skip set in `analyze_mutations_for_skip_target()`
- Use deletion-oriented candidate generation for a duplication
- Patch target-exon filtering without replacing the underlying mutation-specific candidate generation

## Testing

```bash
cd cursor_ai
module load miniconda3/24.1.2-py310 && conda activate 7030_capstone
export MPLCONFIGDIR=cursor_ai/.mplconfig
PYTHONPATH=. pytest -q
```

Last known: **154 tests passing** (includes `tests/test_exon_skipping_analysis.py`).

Exon-skipping regression tests:

```bash
PYTHONPATH=. pytest tests/test_exon_skipping_analysis.py tests/test_exon_skipping_catalog.py tests/test_skip_candidate_filtering.py -q
```

## Git / workflow notes

- Branch: `DMD-Mutation-Test`
- Only commit when user explicitly asks
- Do not use `conda env update --prune` on shared `7030_capstone` env
- After code changes, restart Streamlit so the browser picks up updates

## Folder-level AI prompts

| Folder | File |
|--------|------|
| `data/` | `data/ai_prompt.md` |
| `src/` | `src/ai_prompt.md` |
| `pages/` | `pages/ai_prompt.md` |
| `scripts/` | `scripts/ai_prompt.md` |
| `tests/` | `tests/ai_prompt.md` |

## Quick verification

```bash
# Generate standalone HTML preview (no Streamlit)
PYTHONPATH=. python -c "
from src.mutation_map_svg import build_interactive_map_html
from src.mutation_catalog import load_catalog
rows = load_catalog().head(3).to_dict('records')
open('outputs/map_preview.html','w').write(build_interactive_map_html(rows))
print('Wrote outputs/map_preview.html')
"
```

## Cursor handoff prompt (exon skipping — Model by Exon Target)

Paste when fixing or extending target-exon analysis:

```
Inspect: pages/3_Exon_Skipping_Analysis.py, src/exon_skipping_analysis.py,
src/exon_skipping_catalog.py (analyze_mutations_for_skip_target), data/dmd_exons_grch38.csv.

CORE RULE: Selected target exon is a DISPLAY FILTER only. For each mutation independently:
  mutation_delta = ± sum(coding lengths of affected exons)
  final_delta = mutation_delta - sum(coding lengths of skipped exons)
  Frame restored when final_delta % 3 == 0

Deletion candidates: contiguous upstream/downstream blocks; deleted exons cannot be targets.
Duplication candidates: duplicated interval required; extend upstream/downstream.
Rank by fewest targets, then fewest additional coding bases removed.

Do NOT filter mutations by overlap with target exon. Do NOT default to exon 45.
Run: pytest tests/test_exon_skipping_analysis.py -v
Validate: dup2→skip2, dup3-4→no single, del7→6,8 etc, del45→skip44, del3-7≠target45.
```
