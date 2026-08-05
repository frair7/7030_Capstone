# Application source — `src/`

Python modules for the DMD Mutation and Exon-Skipping Explorer.

## Cohort mutation map (Mutation Explorer)

| Module | Role |
|--------|------|
| `mutation_map_svg.py` | Interactive HTML/SVG for Streamlit (`build_interactive_map_html`) |
| `mutation_map_figure.py` | Matplotlib cohort figure for PNG export |
| `exon_shape_renderer.py` | Transcript + domain row drawing (matplotlib) |
| `puzzle_geometry.py` | Puzzle-piece paths, bumps, domain bubble SVG |
| `part_segments.py` | `PartSegment`, `build_part_segments`, `find_part_segments`, `segments_for_exon` |
| `domain_alignment.py` | `resolve_domains`, hinge guide positions |
| `dp427m_exon_data.py` | Loads `EXON_TABLE` / `DOMAIN_MAP` from reference script |
| `exon_display_config.py` | Phenotype options, %Dys heatmap colors |
| `mutation_viz.py` | Exon range parsing from catalog rows |

### Rendering architecture

**Transcript row (SVG and matplotlib):**

1. Domain color band at coding-position x-coordinates
2. Puzzle-clipped fills (use `segments_for_exon`, not x-range overlap)
3. Per-segment stroke outlines using each part's `border` color

**Domain row:**

- Rounded bubbles clipped to segment-aligned fill/border strips
- Same x-coordinates as transcript row

## Mutation catalog

| Module | Role |
|--------|------|
| `mutation_catalog.py` | Load/save/filter catalog CSV |
| `mutation_autofill.py` | Computed fields from region input |
| `region_parser.py` | Quick entry parsing (`e45-e55`) |

## Core analysis

| Module | Role |
|--------|------|
| `variant_parser.py` | HGVS-style variant parsing |
| `coordinate_mapper.py` | Exon records from reference CSV |
| `frame_analysis.py` | Reading-frame consequence |
| `exon_skipping.py` | Skip candidate search |
| `models.py` | Pydantic data models |

## Visualization & UI

| Module | Role |
|--------|------|
| `visualization.py` | Plotly transcript + domain figure |
| `streamlit_theme.py` | Widescreen layout, nav bar |
| `ui_helpers.py` | Session state, research warning |
| `reporting.py` | CSV/text/HTML download helpers |
| `reference_data.py` | CSV loaders for Reference Map page |
| `config.py` | Transcript constants |

## Tests

Mirror modules under `../tests/`. Run:

```bash
PYTHONPATH=. pytest tests/ -q
```
