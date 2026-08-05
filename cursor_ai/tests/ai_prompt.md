# AI prompt — `tests/`

Folder-specific context for `cursor_ai/tests/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## Before pushing map changes

```bash
PYTHONPATH=. pytest tests/test_mutation_map_svg.py tests/test_puzzle_geometry.py tests/test_domain_alignment.py tests/test_mutation_map_figure.py -q
```

## Before pushing exon-skipping changes

```bash
PYTHONPATH=. pytest tests/test_exon_skipping_analysis.py tests/test_exon_skipping_catalog.py tests/test_skip_candidate_filtering.py tests/test_transcript_reconstruction.py -q
```

## Before pushing Reference Map changes

```bash
PYTHONPATH=. pytest tests/test_reference_atlas.py tests/test_reference_data.py tests/test_visualization.py -q
```

`test_reference_atlas.py` verifies:

- continuous section order with no `st.tabs()`
- GRCh38 coding genomic bounds and numeric CSV fields
- complete 13,992-bp spliced-transcript scale with explicit UTRs
- Leiden transcript-nucleotide feature conversion to Dp427m amino acids
- level-specific genomic, transcript, and protein axis labels
- external links built from the validated local DMD locus

## Exon skipping tests (`test_exon_skipping_analysis.py`)

Uses canonical Dp427m coding lengths via `conftest.py` fixture `exon_coding_lengths` (from `build_exon_records()` — do not hardcode exon lengths).

**Biological validation cases (must pass):**

| Test | Mutation | Expected |
|------|----------|----------|
| `test_deletion_45_single_exon_44_restores_frame` | Del 45 | Skip 44 restores frame |
| `test_deletion_7_expected_multi_exon_candidates` | Del 7 | 6,8 / 6,8,9 / 5,6,8 / 2,3,4,5,6; no single-exon; exon 7 never in skip targets |
| `test_duplication_3_4_expected_candidates` | Dup 3–4 | 3,4 / 3,4,5 / 3,4,5,6,7,8; no single-exon |
| `test_duplication_2_is_corrected_by_targeting_exon_2` | Dup 2 | Skip 2 restores frame |
| `test_unrelated_deletion_does_not_appear_for_exon_45` | Del 3–7 | No candidates include exon 45 |
| `test_duplication_rejects_downstream_only_skip_set` | Dup 3–7 | Reject skip 8–9 because it omits duplicated interval |
| `test_deletion_and_duplication_3_7_use_different_candidates` | S42 dup 3–7 / S79 del 3–7 | Dup ranks 3–7 first; deletion independently validates 8–9 |

**Invariant tests:**

- Every candidate: `final_delta % 3 == 0`
- No deleted exon in skip targets
- Every duplication candidate includes duplicated interval
- Target-exon filter never creates candidates (only removes)
- Candidate order is deterministic

Do not weaken these tests to make incorrect logic pass.

## Key assertions (map)

- `test_puzzle_geometry.py` — junction-owned edges use move (`M`), not line (`L`)
- `test_mutation_map_svg.py` — SVG contains segment clips, domain border colors, `hinge-guide`
- `test_domain_alignment.py` — H1 at exon 8/10 splits; R3/H2 boundary in exon 17

## Matplotlib cache

If font/cache errors occur on OSC:

```bash
export MPLCONFIGDIR=../.mplconfig
```

## Adding tests

Only add tests that cover real behavior (not trivial asserts). Match existing pytest style in sibling files.
