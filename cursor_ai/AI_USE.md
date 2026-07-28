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
- pytest suite (48+ tests)

### Phase 4 — Visualization
- Plotly combined transcript + domain figure (`src/visualization.py`)

### Phase 5 — Multi-page Streamlit app
- Pages: Mutation Explorer, Exon Skipping, Reference Map, Methods
- `src/ui_helpers.py`, `src/reporting.py` — CSV/text/HTML downloads

### Mutation Explorer redesign (in progress)
- Mutation catalog (`data/mutation_catalog.csv`, 101 rows from legacy CSV)
- Intake form, autofill, filters, checkbox table
- **Cohort mutation map** rebuilt with matplotlib to match published figure:
  - Equal-width exon units, puzzle-piece shapes (`src/exon_shape_renderer.py`)
  - `data/dp427m_exon_edge_shapes.csv` — curated 5′/3′ edge shapes
  - `src/mutation_map_figure.py` — transcript + domain + patient deletion bars
  - Aligned to user's Claude `make_figure.py` reference code

## Context file for future sessions

See **[ai_prompt.md](ai_prompt.md)** for full project handoff context, file map, and resume instructions.

## Major prompts

1. Full capstone specification: DMD Mutation and Exon-Skipping Explorer (phased build).
2. Phase 1-only: scaffold in `cursor_ai/` without touching parent repo scripts.
3. Mutation Explorer as catalog with 15 columns + phenotype and cohort map.
4. Rebuild map to match reference figure; incorporate Claude/ChatGPT reference code.
5. Fix transcript map: equal-width exons, puzzle pieces, reference colors.

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
