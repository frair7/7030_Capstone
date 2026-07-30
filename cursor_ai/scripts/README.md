# Scripts — `scripts/`

Utility and reference scripts for data and figure generation.

| Script | Purpose |
|--------|---------|
| `reference_make_figure.py` | **Authoritative** `EXON_TABLE` + `DOMAIN_MAP`; standalone matplotlib reference figure |
| `fetch_reference_data.py` | Download Ensembl exon table → `data/dmd_exons_grch38.csv` |
| `validate_reference_data.py` | Validate reference CSV integrity |
| `run_tests.sh` | Convenience wrapper for pytest |

## `reference_make_figure.py`

This is the **source of truth** for cohort map styling:

- 79 exons with CDS bp, `five_prime`/`three_prime` puzzle shapes
- Fractional `parts` per exon: `[f0, f1, fill, border]`
- `DOMAIN_MAP` — 31 domain segments with exon spans

Loaded at runtime by `src/dp427m_exon_data.py` via `exec()`.

### Regenerate static reference figure

```bash
cd cursor_ai
python scripts/reference_make_figure.py
# writes data/dp427m_figure.png
```

## Data fetch (network required)

```bash
python scripts/fetch_reference_data.py
python scripts/validate_reference_data.py
```

## Editing domain boundaries

When an exon spans multiple domains (e.g. exon 17: R3 → H2 → R3), update **both**:

1. `parts` fractions in `EXON_TABLE`
2. `DOMAIN_MAP` exon spans (e.g. H2 at exon 17, R4 from exon 18)

Then run `pytest tests/test_domain_alignment.py`.
