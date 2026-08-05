# DMD Reference Tables

Validated reference tables for **DMD (dystrophin) Dp427m** map modeling, exon-skipping analysis, and mutation visualization. These tables are the authoritative inputs for building the DMD exon map figure and linking genomic coordinates to protein domains.

## Canonical reference

| Field | Value |
|-------|-------|
| Gene | DMD |
| Isoform | Dp427m |
| RefSeq transcript | NM_004006.3 |
| Ensembl transcript | ENST00000357033.9 |
| UniProt accession | P11532 |
| Genome assembly | GRCh38 |
| Chromosome | X |
| Strand | reverse (minus) |
| Coding exons | 79 |
| Total CDS | 11,058 bp |

## Files

| File | Rows | Description | Primary use |
|------|------|-------------|-------------|
| [`dmd_exons_grch38.csv`](dmd_exons_grch38.csv) | 79 | Genomic, transcript, and CDS coordinates for each coding exon, plus 5′/3′ puzzle-piece edge shapes | Map variants to exons; compute reading frame and splice phases; anchor exon blocks on the genomic map |
| [`dmd_exon_map_styles.csv`](dmd_exon_map_styles.csv) | 79 | Per-exon visualization styling: segment count, edge shapes/directions, and up to three border/fill color pairs | Render the DMD exon map figure with domain-colored, multi-segment exon tiles |
| [`dmd_protein_domains_dp427m.csv`](dmd_protein_domains_dp427m.csv) | 48 | Dp427m protein regions and functional features with amino-acid coordinates and map colors | Color exon segments by dystrophin domain; annotate UTRs, rod repeats, hinges, and C-terminal motifs |

## Data sources

### 1. Ensembl REST API — exon coordinates

- **URL:** https://rest.ensembl.org/
- **Transcript:** ENST00000357033.9 (GENCODE/Havana; matches Dp427m)
- **Assembly:** GRCh38
- **Cross-reference:** NCBI RefSeq NM_004006.3 / NP_003997.2
- **Used in:** `dmd_exons_grch38.csv` (`source`, `source_version`, `date_accessed` columns)
- **Provides:** Chromosome, genomic start/end, strand, Ensembl exon IDs, transcript positions, CDS boundaries, coding length, splice phases, cumulative CDS offsets

Splice phases are derived from cumulative coding length mod 3.

### 2. DMD map modeling reference — exon shapes and colors

- **Description:** Curated puzzle-piece edge shapes (Flat, Round, Point, N/A) and per-exon styling for the DMD map figure
- **Used in:** `dmd_exons_grch38.csv` (`five_prime_shape`, `three_prime_shape`) and `dmd_exon_map_styles.csv`
- **Provides:** Number of visual segments per exon, 5′/3′ edge geometry, shape direction glyphs (`|`, `)`, `>`), and primary/secondary/tertiary border and fill colors

Edge shapes encode splice junction geometry for the linear exon map. Exons that span multiple protein domains use `num_parts` > 1 with distinct color pairs per segment.

### 3. Leiden DMD database (dmd.nl) — Dp427m protein features

- **URL:** https://www.dmd.nl/seqs/Dp427m.html
- **Used in:** `dmd_protein_domains_dp427m.csv`
- **Provides:** Amino-acid coordinates for UTRs, actin-binding domain (ABD), central rod domain, 24 spectrin-repeat regions (R1–R24), four hinge regions (H1–H4), WW-domain, cysteine-rich domain, EF-hands, ZZ-domain, syntrophin binding sites, C-terminal region, and polyadenylation signals

Regions may overlap by design (e.g., hinge 4 overlaps repeat region 24; WW-domain sits within the rod domain).

### 4. UniProt — dystrophin annotation and variant context

- **URL:** https://www.uniprot.org/uniprotkb/P11532/variant-viewer?loadVariants=true
- **Accession:** P11532 (dystrophin)
- **Used for:** Cross-checking domain boundaries, functional motifs, and published variant positions against the Dp427m feature map

### 5. GeneCards — gene-level summary

- **URL:** https://www.genecards.org/search/results?q=DMD
- **Used for:** Gene nomenclature, disease associations, and transcript isoform overview

## Minus-strand note

DMD lies on the **minus (reverse) strand** of chromosome X.

- **Genomic coordinates** are stored in ascending order (`genomic_start_grch38` ≤ `genomic_end_grch38`).
- **Exon numbering** follows transcript order: exon 1 is the 5′ end of the mRNA, exon 79 is the 3′ end.
- Transcript direction runs opposite to increasing genomic coordinates; exon 1 has the highest genomic coordinates.

## Column reference

### `dmd_exons_grch38.csv`

| Column | Description |
|--------|-------------|
| `exon_number` | Transcript-order exon index (1–79) |
| `five_prime_shape` / `three_prime_shape` | Puzzle-piece edge type at 5′ and 3′ splice junctions |
| `chromosome`, `genomic_start_grch38`, `genomic_end_grch38`, `strand` | GRCh38 genomic location |
| `transcript_id`, `ensembl_exon_id` | Ensembl identifiers |
| `transcript_exon_start`, `transcript_exon_end` | Positions along the spliced transcript (nt) |
| `cds_start`, `cds_end`, `coding_length_bp` | CDS coordinates and exon coding length |
| `splice_phase_5prime`, `splice_phase_3prime` | Reading frame at each splice junction (0, 1, or 2) |
| `cumulative_cds_start`, `cumulative_cds_end` | Running CDS offsets across the transcript |
| `source`, `source_version`, `date_accessed` | Provenance metadata |

### `dmd_exon_map_styles.csv`

| Column | Description |
|--------|-------------|
| `exon_number` | Join key to `dmd_exons_grch38.csv` |
| `num_parts` | Number of colored segments within the exon tile (1–3) |
| `five_prime_shape`, `five_prime_shape_direction` | 5′ edge geometry and direction glyph |
| `primary_border_color`, `primary_fill_color` | First segment colors (hex) |
| `second_border_color`, `second_fill_color` | Second segment colors (when `num_parts` ≥ 2) |
| `third_border_color`, `third_fill_color` | Third segment colors (when `num_parts` = 3) |
| `three_prime_shape`, `three_prime_shape_direction` | 3′ edge geometry and direction glyph |

### `dmd_protein_domains_dp427m.csv`

| Column | Description |
|--------|-------------|
| `region_name` | Full feature name (e.g., "repeat region 12", "Cysteine-rich domain") |
| `feature` | Original coordinate notation from source (e.g., "4949..5272") |
| `first_amino_acid`, `last_amino_acid` | Inclusive amino-acid positions on the Dp427m protein |
| `short_code` | Abbreviated label (e.g., ABD, R12, H3, CR Domain) |
| `border_color`, `fill_color` | Hex colors for map rendering (blank when not color-coded) |
| `source`, `source_version`, `date_accessed` | Provenance metadata |

## Domain color legend

| Color (border / fill) | Domain class |
|-----------------------|--------------|
| `#EB6F6F` / `#FFCCCC` | Actin-binding domain (ABD) |
| `#4BAF89` / `#AEDCCA` | Hinge regions (H1–H4) |
| `#F0C266` / `#FFFFCC` | Rod spectrin-repeat regions (R1–R24) and WW-domain |
| `#4D93D9` / `#4D93D9` | Cysteine-rich domain |
| `#D86DCD` / `#F0E3F4` | Carboxy-terminal region |

## How the tables work together

1. **`dmd_exons_grch38.csv`** — Look up a variant's genomic position to find the affected exon, CDS offset, and reading frame.
2. **`dmd_protein_domains_dp427m.csv`** — Map a CDS/amino-acid position to a dystrophin functional region for clinical interpretation.
3. **`dmd_exon_map_styles.csv`** — Render each exon as a styled tile on the linear DMD map, splitting exons that cross domain boundaries into multiple colored segments.

## Date accessed

All tables were compiled on **2026-07-28**.

## Related project data

The `cursor_ai/` Streamlit app loads these tables via `src/reference_data.py`.
The summary domain file (`cursor_ai/data/dmd_domains.csv`) remains in the app
data folder for the protein-domain visualization track.
