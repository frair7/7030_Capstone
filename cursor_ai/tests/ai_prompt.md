# AI prompt — `tests/`

Folder-specific context for `cursor_ai/tests/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## Before pushing map changes

```bash
PYTHONPATH=. pytest tests/test_mutation_map_svg.py tests/test_puzzle_geometry.py tests/test_domain_alignment.py tests/test_mutation_map_figure.py -q
```

## Key assertions

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
