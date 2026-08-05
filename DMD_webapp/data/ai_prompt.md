# AI prompt — `data/`

Folder-specific context for the `DMD_webapp/data/` directory. Full project context: [../ai_prompt.md](../ai_prompt.md).

## Purpose

Static CSV/JSON files used by the Streamlit app, tests, and validation scripts.

## What feeds the cohort mutation map

**Not** the CSVs in this folder (except `mutation_catalog.csv` for patient rows).

Cohort map exon shapes, domain colors, and fractional splits come from:

- `../scripts/reference_make_figure.py` → `EXON_TABLE`, `DOMAIN_MAP`
- Loaded via `../src/dp427m_exon_data.py`

## Key files

| File | Used by |
|------|---------|
| `dmd_exons_grch38.csv` | Frame analysis, coordinate mapper, Reference Map |
| `dmd_domains.csv` | Reference Map, broad domain summary |
| `mutation_catalog.csv` | Mutation Explorer catalog intake/table |
| `example_variants.csv` | Tests and demos |

## Do not edit casually

- `dmd_exons_grch38.csv` — regenerate via `scripts/fetch_reference_data.py`
- `mutation_catalog.csv` — user data; backup before bulk changes

## Regenerate validation

```bash
python scripts/validate_reference_data.py
```
