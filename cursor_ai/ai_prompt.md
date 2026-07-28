# AI context prompt — DMD Mutation and Exon-Skipping Explorer

Use this file to resume AI-assisted work on the capstone without re-explaining the project.

## Project summary

**BSGP 7030 Capstone** — locally hosted **Streamlit** web app for research/educational exploration of *DMD* / Dp427m mutations (NM_004006.3, GRCh38, 79 coding exons). **Not for clinical use.**

Active codebase: `cursor_ai/` subfolder. Parent repo (`7030_Capstone/`) has older scripts and `data/input/DMD_Mutations_clean.csv`.

## Environment (OSC)

```bash
module load miniconda3/24.1.2-py310
conda activate 7030_capstone
conda env update --name 7030_capstone --file cursor_ai/environment.yml   # do NOT use --prune
cd cursor_ai
python -m streamlit run app.py --server.port 20001
```

User presents from OSC browser only (not Windows SSH tunnel). Port 8501 may be taken by another user.

Dependencies include: streamlit, pandas, numpy, plotly, matplotlib, biopython, pytest, pydantic.

## Development phases

| Phase | Status | Notes |
|-------|--------|-------|
| 1 Scaffold | Done | `app.py`, `environment.yml`, `src/config.py` |
| 2 Reference data | Done | `data/dmd_exons_grch38.csv` (79 exons, Ensembl), `data/dmd_domains.csv` |
| 3 Core analysis | Done | variant parser, coordinate mapper, frame analysis, exon skipping; 48+ tests |
| 4 Visualization | Done | Plotly transcript map in `src/visualization.py` |
| 5 Multi-page app | Done | 4 Streamlit pages, reporting downloads |
| 6 Validation / polish | Pending | README status table, final student review |

## App pages

1. **Mutation Explorer** (`pages/1_Mutation_Explorer.py`) — primary focus; mutation catalog intake, filters, cohort map
2. Exon Skipping Analysis
3. Reference Map
4. Methods and Limitations

## Mutation Explorer (current design)

### Catalog

- 101 legacy mutations imported to `data/mutation_catalog.csv`
- 15 catalog columns + `phenotype` (Healthy control, Asymptomatic, Paucisymptomatic, IMD, BMD, DMD)
- Intake form with quick entry autofill (`e45-e55` style)
- Filterable table with row checkboxes → **Plot selected on map**

Optional map sidebar fields (not yet in catalog CSV): `group`, `pct_dys_wb`

### Cohort mutation map (`src/mutation_map_figure.py`)

Uses authoritative data from `scripts/reference_make_figure.py` (user's Claude reference):

- **EXON_TABLE**: per-exon CDS bp, puzzle-piece shapes, fractional domain sub-segments
  (e.g. exon 8 ABD→H1 ~39%/61%, exon 17 R3→H2→R4 three-way, exon 64 H4→CR)
- **DOMAIN_MAP**: from user's domain/feature nucleotide table colors (31 segments, zero gaps)
- **Proportional exon widths** from true CDS bp with `MIN_EXON_W = 0.35` floor (not equal-width)
- **Multi-part clipped fills** via `src/exon_shape_renderer.py` — 87 colored segments across 79 exons
- Patient table: Group | Participant | MW | %Dys(WB) heatmap
- Phenotype boxes on right (groups D/E)

Catalog optional fields: `group`, `pct_dys_wb`

### Key map modules

| File | Role |
|------|------|
| `scripts/reference_make_figure.py` | Authoritative EXON_TABLE + DOMAIN_MAP (from user) |
| `src/dp427m_exon_data.py` | Loads constants from reference script |
| `src/mutation_map_figure.py` | Main cohort figure builder |
| `src/exon_shape_renderer.py` | Multi-part clipped puzzle-piece exons |

## User-provided reference materials

Located in `~/Downloads/` (not in repo):

- `make_figure.py` — Claude matplotlib reproduction (now superseded by `scripts/reference_make_figure.py`)
- `dmd_mutation_map.html` — ChatGPT interactive SVG (proportional bp widths, puzzle-piece shapes)
- Reference figure image (cohort deletion map with P-1…P-29)

**Important:** Both reference implementations now use **bp-proportional exon widths** (true CDS coding length, not mRNA/UTR length) with a minimum-width floor for short exons (71–79). The authoritative per-exon data lives in `scripts/reference_make_figure.py` as `EXON_TABLE`.

## Testing

```bash
cd cursor_ai
module load miniconda3/24.1.2-py310 && conda activate 7030_capstone
export MPLCONFIGDIR=cursor_ai/.mplconfig   # if matplotlib cache permission issues
PYTHONPATH=. pytest -q
```

Last known: **50 tests passing**.

## Git / workflow notes

- Branch: `DMD-Mutation-Test`
- User requested **do not commit** unless explicitly asked
- Do not use `conda env update --prune` on shared `7030_capstone` env

## Known follow-ups

- Add `group` and `pct_dys_wb` columns to catalog intake/table if user wants full reference sidebar
- Phase 6: README phase table still shows outdated status
- `figures/` folder in parent repo is empty
- Student manual review of frame logic and exon coordinates still pending per AI_USE.md

## Prompts that shaped this work

1. Build DMD Mutation Explorer capstone app in phases (scaffold → data → analysis → viz → UI)
2. Redesign Mutation Explorer as mutation catalog with 15 columns + phenotype
3. Rebuild mutation map to match published cohort figure (transcript + domain + deletion bars)
4. Incorporate Claude `make_figure.py` and ChatGPT `dmd_mutation_map.html` reference code
5. Fix map: wider figure, larger transcript/domain fonts, equal-width exons, puzzle pieces
6. Reconcile domain boundaries with Ensembl exon coordinates for multi-segment rendering (see below)

### Claude: Reconciled domain boundaries with exon coordinates for precise multi-segment rendering

Authoritative Ensembl-sourced exon table with sub-exon domain boundaries. Some exons span two or three domains within a single exon; those are split into correctly-colored sub-segments rather than treating each exon as one domain.

**Source data:** GRCh38, ENST00000357033.9 (Dp427m). Per-exon CDS coding length (bp), 5′/3′ jigsaw edge shapes (`Flat` / `Round` / `Point` / `N/A`), and 1–3 fractional `parts` per exon with exact domain fill/border hex colors from the user's domain/feature nucleotide coordinate table.

**Junction exons (multi-part):** 8 exons straddle a domain boundary and render as proportionally-split sub-segments inside a single jigsaw shape:

| Exon | Split | Approx. fractions |
|------|-------|-------------------|
| 8 | ABD → H1 | ~39% / ~61% |
| 10 | H1 → R1 | ~25% / ~75% |
| 17 | R3 → H2 → R4 | 3 parts |
| 50 | R19 → H3 | 2 parts |
| 51 | H3 → R20 | 2 parts |
| 61 | R24 → H4 | 2 parts |
| 64 | H4 → CR Domain | 2 parts |

**Verification (all passed):**
- All 79 exons match the given "# of Parts" and color columns exactly
- Zero gaps, zero overlaps across full CDS coverage after closing small linker gaps at domain boundaries
- 87 colored segments total (79 exons + 8 extra from junction splits)
- All 5 domain colors present end-to-end on the transcript row
- Pixel-level inspection confirmed transitions (e.g. ABD pink `#FFCCCC` through exons 1–7, H1 green `#AEDCCA` from exon 8's internal split onward)

**Domain colors (fill / border):**
- ABD `#FFCCCC` / `#EB6F6F`
- Hinges H1–H4 `#AEDCCA` / `#4BAF89`
- Repeats R1–R24 `#FFFFCC` / `#F0C266`
- CR Domain fill `#DAE9F8`, border `#4D93D9`
- C-terminal `#F0E3F4` / `#D86DCD`

**Implementation:**
- `scripts/reference_make_figure.py` — authoritative `EXON_TABLE` + `DOMAIN_MAP`; standalone matplotlib figure script
- `src/dp427m_exon_data.py` — loads constants from reference script via `exec()`
- `src/exon_shape_renderer.py` — `add_exon_patch()` draws clipped multi-part fills inside puzzle-piece outlines
- `src/mutation_map_figure.py` — Streamlit cohort map uses same table and proportional-width algorithm (min-width floor `MIN_EXON_W = 0.35`)
- Companion HTML (`dmd_mutation_map.html`, user Downloads) — interactive version with hover tooltips; uses same `EXON_TABLE` and clipPath-based multi-part rendering

**Layout algorithm (shared HTML + matplotlib):**
1. `raw_w = bp / total_bp * TOTAL_TRACK_WIDTH` per exon
2. Exons below `MIN_EXON_W` get fixed floor width; remaining width redistributed proportionally among larger exons
3. Each exon's `parts` are `[f0, f1, fill, border]` fractions of that exon's own width, clipped to jigsaw outline

**Static previews:**
- `data/preview_mutation_map.png` — Streamlit cohort map (single patient, regenerate via Quick verification below)
- `data/dp427m_figure.png` — full reference figure with all 29 patients (`python scripts/reference_make_figure.py`)

## Quick verification

```bash
# Generate static preview
python -c "
from src.coordinate_mapper import build_exon_records
from src.mutation_map_figure import create_cohort_mutation_map
rows = [{'id':'P-6','group':'B','start_region':'e3','stop_region':'e29',
         'expected_protein_size_kda':'273.4','pct_dys_wb':'49.16','phenotype':'DMD'}]
fig = create_cohort_mutation_map(rows, build_exon_records())
fig.savefig('data/preview_mutation_map.png', dpi=150, bbox_inches='tight')
"
```
