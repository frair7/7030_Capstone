# AI prompt — `scripts/`

Folder-specific context for `cursor_ai/scripts/`. Full project context: [../ai_prompt.md](../ai_prompt.md).

## `reference_make_figure.py` is authoritative

`src/dp427m_exon_data.py` execs `EXON_TABLE` and `DOMAIN_MAP` from this file. Any change here affects:

- Mutation Explorer interactive map
- Mutation Explorer PNG export
- Domain alignment tests

## DOMAIN_MAP vs EXON_TABLE.parts

- `parts` = fine-grained color splits within an exon (coding-position fractions)
- `DOMAIN_MAP` = coarse domain labels for bubble boundaries and hinge guides
- They must agree at boundaries (e.g. R3 ends where H2 starts in exon 17)

## Known multi-part exons

8, 10, 17, 50, 51, 61, 64 — see root `ai_prompt.md` for split table.

## After editing

```bash
PYTHONPATH=. pytest tests/test_domain_alignment.py -q
```

Restart Streamlit if the app is running.

## Not the same as

`reference_tables/dmd_exon_map_styles.csv` (parent repo CSV pipeline) — used by Reference Map via `src/reference_data.py`, not the cohort map.
