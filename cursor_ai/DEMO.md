# Demo script — DMD Mutation Explorer

Use this outline during your capstone presentation (~8 minutes).

## Before the demo

```bash
module load miniconda3/24.1.2-py310
conda activate 7030_capstone
cd ~/7030_Capstone/cursor_ai
export MPLCONFIGDIR=~/7030_Capstone/cursor_ai/.mplconfig
python -m streamlit run app.py --server.port 8503 --server.headless true
```

Open **http://localhost:8503/Mutation_Explorer** (OSC browser or SSH tunnel).

## Demo flow

### 1. Home page (30 sec)
- Project title and **research-use-only** warning
- Reference transcript: NM_004006.3 / GRCh38 / 79 exons

### 2. Mutation Explorer (4 min)
- Show **transcript row**: puzzle-piece exons with domain-colored fills and matching border hues
- Show **domain row**: protein bubbles aligned to the same x-coordinates
- Point out **hinge dashed guides** (8 lines at H1–H4 boundaries)
- Scroll catalog table → select 2–3 mutations → **Plot selected on map**
- Show deletion bars with Group / Participant / MW / %Dys(WB) panel on the left
- Hover an exon for CDS length; hover a bar for participant details
- Optional: expand **Download static map (PNG)**

### 3. Exon Skipping Analysis (2 min)
- Show ranked candidates (computational only)
- Read assumptions aloud — not a treatment recommendation

### 4. Reference Map (1 min)
- Full 79-exon table and domain track
- Mention Ensembl provenance

### 5. Methods & Limitations (1 min)
- Frame calculation from splice phases
- Known limitations and non-clinical disclaimer

## Q&A talking points
- Minus strand: exon 1 = highest genomic coordinate
- Transcript colors follow **coding position**, not whole-exon domain labels
- Multi-part exons (e.g. exon 8 ABD→H1, exon 17 R3→H2) split proportionally inside one puzzle shape
- Candidates ranked by fewest extra exons skipped
