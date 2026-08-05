# AI prompt — `src/`

Folder-specific context for `DMD_webapp/src/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## When editing the cohort map

Primary files (in dependency order):

1. `dp427m_exon_data.py` — data loader (exec `reference_make_figure.py`)
2. `part_segments.py` — segment x-ranges, domain resolution helpers
3. `puzzle_geometry.py` — exon/junction paths; junction-owned edges use `M` not `L`
4. `domain_alignment.py` — domain bubble extents, hinge guides
5. `exon_shape_renderer.py` — matplotlib rows
6. `mutation_map_svg.py` — Streamlit interactive map
7. `mutation_map_figure.py` — PNG export

## Reference atlas (`reference_atlas.py`)

Authoritative transforms and Plotly viewers for `pages/4_Reference_Map.py`:

- `build_exon_reference_frame()` — typed genomic, full transcript, CDS, UTR,
  amino-acid, and derived coding-genomic fields
- `create_genomic_viewer()` — standard chromosome orientation with reverse-strand arrow
- `create_transcript_viewer()` — equal-width schematic or complete spliced-transcript scale
- `build_detailed_protein_features()` — explicit Leiden transcript-nt → AA conversion
- `create_protein_overview()` / `create_protein_feature_viewer()` — AA-only protein axes
- `resource_table()` — validated and grouped external links

Critical coordinate facts:

- ENST00000357033.9 local transcript length: 13,992 bp
- Local CDS length: 11,058 bp including stop codon
- Derived UTR lengths from Ensembl table: 237 bp 5′ and 2,697 bp 3′
- Detailed feature CSV values 1–13,993 are transcript nucleotides despite legacy
  headers `first_amino_acid` / `last_amino_acid`; never plot those raw values as AA
- Genomic coding bounds are locally derived from full exon bounds, strand, and
  validated coding length and must remain inside the GRCh38 full-exon interval

## Exon skipping modules

### `exon_skipping_analysis.py` — authoritative mutation-specific frame math

- `ExonMutation`, `SkipCandidate` (analysis dataclasses — not `models.SkipCandidate`)
- `find_frame_restoring_candidates()` — generates adjacent skip sets, filters to `final_delta % 3 == 0`
- `filter_candidates_by_target_exons()` — display filter only; never determines frame restoration
- `evaluate_skip_set()` — single skip-set frame check
- `candidate_is_allowed()` — deleted exons cannot be targets; duplication candidates must include the full duplicated interval
- Uses coding lengths from `ExonRecord.coding_length_bp` / `dmd_exons_grch38.csv`

### `exon_skipping.py` + `transcript_reconstruction.py` — legacy junction-aware utilities

- `find_skip_candidates()` — boundary extensions + optional advanced combos
- `reconstruct_deletion_skip()` — builds retained exons, checks junction splice phases
- `evaluate_candidate()` — returns `models.SkipCandidate` only when frame restored
- `_candidate_from_transcript()` — builds legacy `SkipCandidate` for UI/schematic
- Do not call this deletion-oriented path for duplication recommendations

### `exon_skipping_catalog.py` — catalog page bridge

- `rank_skip_candidates_for_row()` — Model by Mutation ranking; calls `find_frame_restoring_candidates()` for both deletion and duplication
- `analyze_mutations_for_skip_target()` — Model by Exon Target tab; **must call `exon_skipping_analysis`**, not pass user targets directly to `reconstruct_deletion_skip()`
- `_legacy_candidate_from_analysis()` — converts analysis output to `models.SkipCandidate` for schematic overlays

## Critical rules (learned from debugging)

- Use `segments_for_exon(segments, exon_n)` for per-exon fills/outlines — **not** `segments_in_xrange` with bump padding (neighbor color bleed).
- Transcript colors follow **coding position** (`parts` fractions), not whole-exon domain labels.
- Domain row fill/border must use the **same segment x-coordinates** as the transcript row.
- Hinge guides: 8 positions (start + end of H1–H4), drawn from domain row downward only.
- `TRACK_LEFT` separates metadata table from exon track; mutation bars use white margins.
- **Target exon must never appear near the start of skip analysis** — calculate per-mutation candidates first, filter by target last.
- Duplication candidates must always contain every exon in the duplicated interval.

## Layout constants (`mutation_map_svg.py`)

- `TRANSCRIPT_DOMAIN_GAP = 12` (white break between rows)
- `TRACK_LEFT = 205` (exon track origin)
- SVG uses fixed pixel `width`/`height` + `.map-scroll` wrapper (not `width="100%"`)

## After changes

```bash
PYTHONPATH=. pytest tests/test_mutation_map_svg.py tests/test_puzzle_geometry.py tests/test_domain_alignment.py -q
PYTHONPATH=. pytest tests/test_exon_skipping_analysis.py tests/test_exon_skipping_catalog.py -q
```

Restart Streamlit to see updates in the browser.
