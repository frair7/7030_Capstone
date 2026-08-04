# AI-assisted development log

## Tools used

- **Cursor IDE** with AI coding assistance (Composer / agent mode)
- Reference implementations reviewed: Claude `make_figure.py`, ChatGPT `dmd_mutation_map.html`

## What AI assisted with

### Phase 1 — Scaffold
- `cursor_ai/` project scaffold and directory layout
- Lean `environment.yml` and `requirements.txt` for Streamlit stack
- Minimal Streamlit application shell (`app.py`)
- Centralized transcript configuration (`src/config.py`)

### Phase 2 — Reference data
- Ensembl REST API fetch script (`scripts/fetch_reference_data.py`)
- `data/dmd_exons_grch38.csv` — 79 exons with splice phases
- `data/dmd_domains.csv` — UniProt P11532 domain regions
- Validation script (`scripts/validate_reference_data.py`)

### Phase 3 — Core analysis
- `src/models.py`, `variant_parser.py`, `coordinate_mapper.py`
- `frame_analysis.py`, `exon_skipping.py` — frame logic from splice phases
- `exon_skipping_analysis.py` — mutation-specific candidate generation for Model by Exon Target tab
- pytest suite (154 tests)

### Phase 4 — Visualization
- Plotly combined transcript + domain figure (`src/visualization.py`)

### Phase 5 — Multi-page Streamlit app
- Pages: Mutation Explorer, Exon Skipping, Reference Map, Methods
- `src/ui_helpers.py`, `src/reporting.py` — CSV/text/HTML downloads

### Mutation Explorer cohort map (complete)
- Mutation catalog (`data/mutation_catalog.csv`)
- Intake form, autofill, filters, checkbox table
- **Interactive SVG cohort map** (`src/mutation_map_svg.py`):
  - Proportional exon widths from CDS bp
  - Puzzle-piece shapes with domain-aligned colors and per-segment borders
  - Domain row bubbles with matching fill/border x-coordinates
  - 8 hinge dashed guides, 12 px transcript/domain gap
  - White metadata panel; mutation bars without gray flanking boxes
- **Static PNG export** (`src/mutation_map_figure.py`)
- Authoritative data: `scripts/reference_make_figure.py` (`EXON_TABLE`, `DOMAIN_MAP`)

## Context files for future sessions

- **[ai_prompt.md](ai_prompt.md)** — full project handoff
- Per-folder `ai_prompt.md` in `data/`, `src/`, `pages/`, `scripts/`, `tests/`

## Major prompts

1. Full capstone specification: DMD Mutation and Exon-Skipping Explorer (phased build).
2. Mutation Explorer as catalog with 15 columns + phenotype and cohort map.
3. Rebuild map to match reference figure; incorporate Claude/ChatGPT reference code.
4. Transcript puzzle shapes over continuous color band; borders correlate with underlying domain hue.
5. Domain row alignment, hinge guides, left panel layout, R3/H2 exon 17 boundary fix.
6. Replace Model by Exon Target logic: mutation-specific frame restoration via `exon_skipping_analysis.py`; target exon as display filter only; validated dup2, dup3-4, del7, del45, del3-7≠exon45 cases.

## Manual review checklist

| Output | Reviewed by student | Notes |
|--------|---------------------|-------|
| Transcript IDs / assembly in `src/config.py` | ☐ | Cross-check against NCBI NM_004006.3 |
| `environment.yml` package list | ☐ | Confirm no `--prune` on shared env |
| Research-use disclaimer text | ☐ | |
| No fabricated exon coordinates | ✓ | Fetched from Ensembl ENST00000357033.9 |
| Exon table / splice phases | ☐ | Run `validate_reference_data.py` |
| Mutation map vs reference figure | ☐ | Compare to published cohort plot |
| Frame logic (e.g. del45-50) | ☐ | Independent verification |

## Biological correctness checks (ongoing)

- Phase 2: Validate all 79 exons against NCBI/Ensembl ✓
- Phase 3: Independently verify frame calculations for known cases
- Phase 6: Student sign-off before presentation

## Student responsibilities

- Approve git commits before pushing
- Run `conda env update` on OSC (no `--prune`)
- Verify authoritative coordinates before trusting frame outputs
- Final review of all AI-generated analysis logic and mutation map
