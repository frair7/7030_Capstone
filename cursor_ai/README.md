# DMD Mutation and Exon-Skipping Explorer

**BSGP 7030 Capstone — research and educational web application**

Interactive Streamlit tool for mapping mutations in the human *DMD* gene relative to the Dp427m transcript, evaluating reading-frame consequences, and exploring computational exon-skipping strategies.

> **Not for clinical use.** Frame-restoration candidates are theoretical calculations only.

## Reference transcript

| Field | Value |
|-------|-------|
| Gene | DMD |
| Isoform | Dp427m |
| RefSeq transcript | NM_004006.3 |
| RefSeq protein | NP_003997.2 |
| Ensembl transcript | ENST00000357033.9 |
| Assembly | GRCh38 |
| Chromosome / strand | X / reverse |
| Coding exons | 79 |
| Protein length | 3,685 aa |

## Repository layout

```text
cursor_ai/
├── app.py                      # Streamlit entry point
├── ai_prompt.md                # AI session handoff (full project context)
├── environment.yml             # Conda environment
├── requirements.txt
├── README.md
├── DEMO.md
├── AI_USE.md
├── LICENSE
├── .streamlit/config.toml
├── data/                       # Reference tables + mutation catalog
├── src/                        # Application logic
├── pages/                      # Streamlit multi-page modules
├── scripts/                    # Reference figure + data utilities
├── tests/                      # pytest suite (82 tests)
└── outputs/                    # Generated previews (gitignored content)
```

Parent repository (`7030_Capstone/`) contains earlier batch scripts under `data/input/`.

## Environment setup

```bash
module load miniconda3/24.1.2-py310
conda activate 7030_capstone
conda env update --name 7030_capstone --file cursor_ai/environment.yml
```

> **Do not use `--prune`** — the environment is shared with prior coursework.

Verify:

```bash
python -c "import streamlit, plotly, Bio; print('OK')"
```

## Running the app

From `cursor_ai/`:

```bash
export MPLCONFIGDIR=~/7030_Capstone/cursor_ai/.mplconfig
python -m streamlit run app.py --server.port 8503 --server.headless true
```

**OSC + laptop tunnel:**

```bash
ssh -N -L 8503:localhost:8503 frair7@ascend.osc.edu
```

Open **http://localhost:8503/Mutation_Explorer**

Select catalog rows and click **Plot selected on map** to show patient deletion bars.

## Mutation Explorer map (highlights)

- **Transcript row:** puzzle-piece exons with domain colors aligned to coding position; per-segment border hues
- **Domain row:** rounded protein bubbles with matching fill/border x-coordinates
- **Hinge guides:** 8 dashed vertical lines (start + end of H1–H4)
- **Left panel:** Group, Participant, MW, %Dys(WB) left of the track boundary
- **Data source:** `scripts/reference_make_figure.py` (`EXON_TABLE`, `DOMAIN_MAP`)

Interactive SVG: `src/mutation_map_svg.py` · Static PNG: `src/mutation_map_figure.py`

## Development status

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Structure, environment, minimal app | Complete |
| 2 | Authoritative exon reference data | Complete |
| 3 | Parsing, frame analysis, skip search | Complete |
| 4 | Plotly visualization | Complete |
| 5 | Full Streamlit UI + cohort map | Complete |
| 6 | Tests, validation, documentation | In progress |

## Testing

```bash
cd cursor_ai
export MPLCONFIGDIR=.mplconfig
PYTHONPATH=. pytest -q
```

## Data provenance

- Exon coordinates: Ensembl ENST00000357033.9 → `data/dmd_exons_grch38.csv`
- Cohort map styling: `scripts/reference_make_figure.py` (user-verified fractional domain splits)
- Mutation catalog: `data/mutation_catalog.csv`

See [data/README.md](data/README.md) for details.

## AI use statement

See [AI_USE.md](AI_USE.md). For session handoff context, see [ai_prompt.md](ai_prompt.md).

## License

MIT — see [LICENSE](LICENSE).
