# Demo script — DMD Mutation Explorer

Use this outline during your capstone presentation (~8 minutes).

## Before the demo

```bash
module load miniconda3/24.1.2-py310
conda activate 7030_capstone
cd ~/7030_Capstone/cursor_ai
python -m streamlit run app.py --server.port 20001
```

Open `http://localhost:20001` in your OSC browser.

## Demo flow

### 1. Home page (30 sec)
- Project title and **research-use-only** warning
- Reference transcript: NM_004006.3 / GRCh38 / 79 exons

### 2. Mutation Explorer (3 min)
- Input: `del45-50` → **Run analysis**
- Show parsed summary: exons 45–50, out of frame
- Point out transcript map: deleted exons in red
- Junction: exon 44 | 51
- Toggle **Exon-order schematic** ↔ **Transcript-scale**
- Download CSV / text report

### 3. Exon Skipping Analysis (2 min)
- Show ranked candidates (computational only)
- Select a candidate → updated map with purple skip highlights
- Read assumptions aloud — not a treatment recommendation

### 4. Reference Map (1 min)
- Full 79-exon table and domain track
- Mention Ensembl provenance

### 5. Methods & Limitations (1 min)
- Frame calculation from splice phases, not hard-coded tables
- Known limitations and non-clinical disclaimer

## Q&A talking points
- Minus strand: exon 1 = highest genomic coordinate
- Frame = coding bases removed mod 3 + phase cross-check
- Candidates ranked by fewest extra exons skipped
