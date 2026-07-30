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
2. Exon Skipping Analysis
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

## Testing

```bash
cd cursor_ai
module load miniconda3/24.1.2-py310 && conda activate 7030_capstone
export MPLCONFIGDIR=cursor_ai/.mplconfig
PYTHONPATH=. pytest -q
```

Last known: **82 tests passing**.

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
