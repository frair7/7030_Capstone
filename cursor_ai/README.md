# DMD Mutation and Exon-Skipping Explorer

**BSGP 7030 Capstone — research and educational web application**

Interactive Streamlit tool for mapping mutations in the human *DMD* gene relative to the Dp427m transcript, evaluating reading-frame consequences, and exploring computational exon-skipping strategies.

> **Not for clinical use.** Frame-restoration candidates are theoretical calculations only.

![Screenshot placeholder](docs/screenshot_placeholder.png)

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

## Repository layout (this subfolder)

```text
cursor_ai/
├── app.py                  # Streamlit entry point
├── environment.yml         # Conda environment (adds web-app packages)
├── requirements.txt        # Pip-style direct dependencies
├── README.md
├── DEMO.md
├── LICENSE
├── .gitignore
├── .streamlit/config.toml
├── data/                   # Reference tables (Phase 2)
├── src/                    # Application logic
├── pages/                  # Streamlit multi-page modules (Phase 5)
├── scripts/                # Data validation & utilities
└── tests/                  # pytest suite (Phase 3+)
```

Parent repository (`7030_Capstone/`) contains earlier batch scripts and mutation tables under `data/input/` and `scripts/`.

## Environment setup

The shared Conda environment **`7030_capstone`** already exists on this system. It currently includes Python 3.10, pandas, numpy, Jupyter, R, and scikit-learn, but **not yet** Streamlit or Plotly.

**Do not use `--prune`** when updating — the environment is shared with prior coursework (R, JupyterLab, etc.) and pruning would remove those packages.

```bash
module load miniconda3/24.1.2-py310
conda activate 7030_capstone
conda env update --name 7030_capstone --file cursor_ai/environment.yml
```

Verify:

```bash
python -c "import streamlit, plotly, Bio; print('OK')"
```

## Running locally

From the `cursor_ai/` directory:

```bash
cd cursor_ai
python -m streamlit run app.py
```

Browser: `http://localhost:8501`

## Running on OSC (remote Linux server)

On the server:

```bash
cd cursor_ai
python -m streamlit run app.py \
  --server.address 0.0.0.0 \
  --server.port 8501 \
  --server.headless true
```

From your laptop, open an SSH tunnel (replace placeholders):

```bash
ssh -N -L 8501:localhost:8501 USERNAME@SERVER_HOSTNAME
```

Then browse to `http://localhost:8501`.

> Exact OSC login node, compute node, and port-forwarding steps depend on your account and job type. Consult OSC documentation for your cluster.

## Development status

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Structure, environment, minimal app | **Current** |
| 2 | Authoritative exon reference data | Pending |
| 3 | Parsing, frame analysis, skip search | **Complete** |
| 4 | Plotly visualization | **Complete** |
| 5 | Full Streamlit UI | **Complete** |
| 6 | Tests, validation, documentation | Pending |

## Testing (Phase 3+)

```bash
cd cursor_ai
pytest -q
```

## Data provenance

Reference exon coordinates will be retrieved from NCBI/Ensembl in Phase 2 — not manually guessed. Existing parent-repo files:

- `../data/input/DMD_Mutations_clean.csv`
- `../data/input/DMD_light.csv`

## AI use statement

See [AI_USE.md](AI_USE.md). For session handoff context, see [ai_prompt.md](ai_prompt.md).

## License

MIT — see [LICENSE](LICENSE).
