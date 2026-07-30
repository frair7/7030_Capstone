# Reference data — DMD Dp427m (GRCh38)

Validated reference tables and catalogs for the DMD Mutation and Exon-Skipping Explorer.

## Files

| File | Description |
|------|-------------|
| `dmd_exons_grch38.csv` | Coding exon coordinates, transcript/CDS positions, splice phases (79 exons) |
| `dmd_domains.csv` | Broad dystrophin protein domain regions (amino-acid coordinates) |
| `dp427m_exon_edge_shapes.csv` | Curated 5′/3′ puzzle edge shapes (legacy; cohort map uses `reference_make_figure.py`) |
| `mutation_catalog.csv` | Mutation catalog for Mutation Explorer intake/table |
| `example_variants.csv` | Curated example inputs for testing and demo |
| `cache/ensembl_ENST00000357033.json` | Raw Ensembl API response (reproducibility) |
| `preview_mutation_map.png` | Optional static map preview (regenerate via scripts) |
| `dp427m_figure.png` | Full reference figure from `scripts/reference_make_figure.py` |

## Canonical reference

| Field | Value |
|-------|-------|
| Gene | DMD |
| Isoform | Dp427m |
| RefSeq transcript | NM_004006.3 |
| Ensembl transcript | ENST00000357033.9 |
| Genome assembly | GRCh38 |
| Chromosome | X |
| Strand | reverse (minus) |
| Coding exons | 79 |
| Total CDS | 11,058 bp |

## Cohort map data (authoritative for Mutation Explorer)

The interactive and static cohort maps do **not** read edge shapes or domain colors from these CSVs directly. They load:

- `EXON_TABLE` — per-exon CDS bp, puzzle shapes, fractional `parts` with fill/border hex
- `DOMAIN_MAP` — 31 domain segments (ABD, H1–H4, R1–R24, CR, C-term)

Both are defined in `../scripts/reference_make_figure.py` and loaded by `../src/dp427m_exon_data.py`.

## Minus-strand note

DMD lies on the **minus (reverse) strand** of chromosome X.

- Genomic coordinates are stored in ascending order.
- Exon numbering follows **transcript order** (exon 1 = 5′ mRNA end).
- Exon 1 has the **highest** genomic coordinates.

## Regenerating Ensembl/UniProt tables

From `cursor_ai/`:

```bash
python scripts/fetch_reference_data.py
python scripts/validate_reference_data.py
```

Requires network access for the fetch step.

## Parent repository data (not primary reference)

- `../../data/input/DMD_Mutations_clean.csv` — legacy mutation source (imported into `mutation_catalog.csv`)
- `../../data/input/DMD_light.csv` — precomputed frame lookup
