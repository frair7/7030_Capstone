# Reference data — DMD Dp427m (GRCh38)

Validated reference tables for the DMD Mutation and Exon-Skipping Explorer.

## Files

| File | Description | Rows |
|------|-------------|------|
| `dmd_exons_grch38.csv` | Coding exon coordinates, transcript/CDS positions, splice phases | 79 |
| `dmd_domains.csv` | Broad dystrophin protein domain regions (amino-acid coordinates) | 5 |
| `example_variants.csv` | Curated example inputs for testing and demo | 5 |
| `cache/ensembl_ENST00000357033.json` | Raw Ensembl API response (reproducibility) | — |

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
| Total CDS | 11,058 bp (3,685 aa + stop codon) |

## Minus-strand note

DMD lies on the **minus (reverse) strand** of chromosome X.

- **Genomic coordinates** are stored in standard ascending order (`genomic_start_grch38` ≤ `genomic_end_grch38`).
- **Exon numbering** follows **transcript order** (exon 1 = 5′ end of mRNA, exon 79 = 3′ end).
- Transcript direction runs **opposite** to increasing genomic coordinates: exon 1 has the **highest** genomic coordinates.

## Provenance

### Exon table (`dmd_exons_grch38.csv`)

- **Source:** [Ensembl REST API](https://rest.ensembl.org/)
- **Transcript:** ENST00000357033.9 (GENCODE/Havana, matches Dp427m)
- **Assembly:** GRCh38
- **Cross-reference:** NCBI RefSeq NM_004006.3 / NP_003997.2
- **Splice phases:** Computed from cumulative coding length mod 3 (not manually entered)
- **Date accessed:** See `date_accessed` column in CSV (set when `fetch_reference_data.py` runs)

### Domain table (`dmd_domains.csv`)

- **Source:** [UniProt P11532](https://www.uniprot.org/uniprotkb/P11532) (dystrophin)
- **Coordinate system:** Amino-acid and CDS positions only — **no genomic coordinates**
- **Note:** Hinge regions overlap the rod domain by design; schema supports finer sub-features later

## Regenerating data

From `cursor_ai/`:

```bash
python scripts/fetch_reference_data.py
python scripts/validate_reference_data.py
```

Requires network access for the fetch step. After setup, the app can run offline using the cached CSV files.

## Validation checks

`scripts/validate_reference_data.py` verifies:

- Exactly 79 exons
- Unique, sequential exon numbers (1–79)
- Valid genomic and transcript coordinate ranges
- Contiguous transcript ordering
- Non-negative coding lengths
- Cumulative CDS consistency
- Splice phase consistency
- Minus-strand transcript ordering (exon 1 at higher genomic coordinates)
- Source and version metadata present

## Parent repository data (not used as primary reference)

- `../../data/input/DMD_Mutations_clean.csv` — mutation catalog
- `../../data/input/DMD_light.csv` — precomputed frame lookup table

These may inform Phase 3 validation but do **not** replace authoritative exon coordinates.
