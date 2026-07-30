# AI prompt — `src/`

Folder-specific context for `cursor_ai/src/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## When editing the cohort map

Primary files (in dependency order):

1. `dp427m_exon_data.py` — data loader (exec `reference_make_figure.py`)
2. `part_segments.py` — segment x-ranges, domain resolution helpers
3. `puzzle_geometry.py` — exon/junction paths; junction-owned edges use `M` not `L`
4. `domain_alignment.py` — domain bubble extents, hinge guides
5. `exon_shape_renderer.py` — matplotlib rows
6. `mutation_map_svg.py` — Streamlit interactive map
7. `mutation_map_figure.py` — PNG export

## Critical rules (learned from debugging)

- Use `segments_for_exon(segments, exon_n)` for per-exon fills/outlines — **not** `segments_in_xrange` with bump padding (neighbor color bleed).
- Transcript colors follow **coding position** (`parts` fractions), not whole-exon domain labels.
- Domain row fill/border must use the **same segment x-coordinates** as the transcript row.
- Hinge guides: 8 positions (start + end of H1–H4), drawn from domain row downward only.
- `TRACK_LEFT` separates metadata table from exon track; mutation bars use white margins.

## Layout constants (`mutation_map_svg.py`)

- `TRANSCRIPT_DOMAIN_GAP = 12` (white break between rows)
- `TRACK_LEFT = 205` (exon track origin)
- SVG uses fixed pixel `width`/`height` + `.map-scroll` wrapper (not `width="100%"`)

## After changes

```bash
PYTHONPATH=. pytest tests/test_mutation_map_svg.py tests/test_puzzle_geometry.py tests/test_domain_alignment.py -q
```

Restart Streamlit to see updates in the browser.
