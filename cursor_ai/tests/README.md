# Tests — `tests/`

pytest suite for the DMD Mutation and Exon-Skipping Explorer.

## Run all tests

```bash
cd cursor_ai
export MPLCONFIGDIR=.mplconfig
PYTHONPATH=. pytest -q
```

Last known: **82 tests passing**.

## Map-related tests

| File | Covers |
|------|--------|
| `test_mutation_map_svg.py` | Interactive HTML/SVG output |
| `test_mutation_map_figure.py` | Matplotlib cohort figure |
| `test_puzzle_geometry.py` | Puzzle paths, junction edge behavior |
| `test_domain_alignment.py` | H1/R-repeat boundaries, R3/H2 exon 17 |

## Other modules

| File | Covers |
|------|--------|
| `test_variant_parser.py` | HGVS parsing |
| `test_coordinate_mapper.py` | Exon records |
| `test_frame_analysis.py` | Reading frame |
| `test_exon_skipping.py` | Skip candidates |
| `test_mutation_catalog.py` | Catalog I/O |
| `test_reference_data.py` | CSV loaders |
| `test_visualization.py` | Plotly figures |
| `test_reporting.py` | Download helpers |

## Convenience script

```bash
bash scripts/run_tests.sh
```
