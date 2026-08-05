"""Reference-atlas data transforms and Plotly viewers for DMD/Dp427m."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import pandas as pd
import plotly.graph_objects as go

from src.config import REFERENCE


UTR_FILL = "#E5E7EB"
UTR_BORDER = "#9CA3AF"
DEFAULT_FILL = "#D9E2EC"
DEFAULT_BORDER = "#64748B"

CATEGORY_COLORS = {
    "Major regions": ("#FFCCCC", "#EB6F6F"),
    "Spectrin repeats": ("#FFFFCC", "#F0C266"),
    "Hinges": ("#AEDCCA", "#4BAF89"),
    "Binding sites": ("#E8D5F2", "#8B5AA5"),
    "Functional motifs": ("#DAE9F8", "#4D93D9"),
}

STYLE_CLASSES = [
    {
        "Style Category": "Actin-binding domain",
        "Biological Meaning": "N-terminal actin-binding region",
        "Fill Color": "#FFCCCC",
        "Border Color": "#EB6F6F",
        "Opacity": 1.0,
        "Label Color": "#6B1F1F",
        "Notes": "Approved DMD project palette.",
    },
    {
        "Style Category": "Spectrin repeats and WW region",
        "Biological Meaning": "Central rod repeats and WW-domain",
        "Fill Color": "#FFFFCC",
        "Border Color": "#F0C266",
        "Opacity": 1.0,
        "Label Color": "#6B5200",
        "Notes": "Repeated regions may overlap other feature tracks.",
    },
    {
        "Style Category": "Hinge regions",
        "Biological Meaning": "H1–H4 flexible hinge regions",
        "Fill Color": "#AEDCCA",
        "Border Color": "#4BAF89",
        "Opacity": 1.0,
        "Label Color": "#174D39",
        "Notes": "Shown separately from spectrin repeats.",
    },
    {
        "Style Category": "Cysteine-rich region",
        "Biological Meaning": "Cysteine-rich, EF-hand, and ZZ region",
        "Fill Color": "#DAE9F8",
        "Border Color": "#4D93D9",
        "Opacity": 1.0,
        "Label Color": "#1E4F82",
        "Notes": "Approved DMD project palette.",
    },
    {
        "Style Category": "Carboxy-terminal region",
        "Biological Meaning": "Dp427m C-terminal region",
        "Fill Color": "#F0E3F4",
        "Border Color": "#D86DCD",
        "Opacity": 1.0,
        "Label Color": "#78306F",
        "Notes": "Approved DMD project palette.",
    },
    {
        "Style Category": "Untranslated region",
        "Biological Meaning": "5′ or 3′ UTR sequence",
        "Fill Color": UTR_FILL,
        "Border Color": UTR_BORDER,
        "Opacity": 1.0,
        "Label Color": "#4B5563",
        "Notes": "Neutral styling; never inferred from whitespace.",
    },
]


def _numeric_frame(
    rows: Iterable[dict[str, Any]],
    columns: Iterable[str],
) -> pd.DataFrame:
    frame = pd.DataFrame(list(rows))
    for column in columns:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce").astype(
                "Int64"
            )
    return frame


def build_exon_reference_frame(
    exon_rows: Iterable[dict[str, Any]],
    broad_domains: Iterable[dict[str, Any]] | None = None,
) -> pd.DataFrame:
    """Build typed genomic/transcript/CDS reference fields from local data."""
    numeric_columns = [
        "exon_number",
        "genomic_start_grch38",
        "genomic_end_grch38",
        "transcript_exon_start",
        "transcript_exon_end",
        "cds_start",
        "cds_end",
        "coding_length_bp",
        "splice_phase_5prime",
        "splice_phase_3prime",
        "cumulative_cds_start",
        "cumulative_cds_end",
    ]
    frame = _numeric_frame(exon_rows, numeric_columns).sort_values("exon_number")
    if frame.empty:
        return frame

    frame["full_exon_genomic_length_bp"] = (
        frame["genomic_end_grch38"] - frame["genomic_start_grch38"] + 1
    )
    frame["transcript_exon_length_bp"] = (
        frame["transcript_exon_end"] - frame["transcript_exon_start"] + 1
    )
    frame["noncoding_length_bp"] = (
        frame["transcript_exon_length_bp"] - frame["coding_length_bp"]
    )
    frame["utr_5prime_length_bp"] = 0
    frame["utr_3prime_length_bp"] = 0

    coding_rows = frame.loc[frame["coding_length_bp"] > 0]
    first_index = coding_rows.index[0]
    last_index = coding_rows.index[-1]
    frame.loc[first_index, "utr_5prime_length_bp"] = int(
        frame.loc[first_index, "noncoding_length_bp"]
    )
    frame.loc[last_index, "utr_3prime_length_bp"] = int(
        frame.loc[last_index, "noncoding_length_bp"]
    )

    frame["coding_transcript_start"] = (
        frame["transcript_exon_start"] + frame["utr_5prime_length_bp"]
    )
    frame["coding_transcript_end"] = (
        frame["coding_transcript_start"] + frame["coding_length_bp"] - 1
    )

    reverse = frame["strand"].str.lower().isin({"reverse", "-", "-1"})
    frame["coding_genomic_start"] = frame["genomic_start_grch38"]
    frame["coding_genomic_end"] = frame["genomic_end_grch38"]
    first_mask = frame.index == first_index
    last_mask = frame.index == last_index
    frame.loc[reverse & first_mask, "coding_genomic_end"] = (
        frame.loc[reverse & first_mask, "genomic_start_grch38"]
        + frame.loc[reverse & first_mask, "coding_length_bp"]
        - 1
    )
    frame.loc[reverse & last_mask, "coding_genomic_start"] = (
        frame.loc[reverse & last_mask, "genomic_end_grch38"]
        - frame.loc[reverse & last_mask, "coding_length_bp"]
        + 1
    )
    frame.loc[~reverse & first_mask, "coding_genomic_start"] = (
        frame.loc[~reverse & first_mask, "genomic_end_grch38"]
        - frame.loc[~reverse & first_mask, "coding_length_bp"]
        + 1
    )
    frame.loc[~reverse & last_mask, "coding_genomic_end"] = (
        frame.loc[~reverse & last_mask, "genomic_start_grch38"]
        + frame.loc[~reverse & last_mask, "coding_length_bp"]
        - 1
    )
    frame["coding_genomic_length_bp"] = (
        frame["coding_genomic_end"] - frame["coding_genomic_start"] + 1
    )

    frame["encoded_amino_acid_start"] = (
        (frame["cumulative_cds_start"] - 1) // 3 + 1
    ).clip(upper=REFERENCE.protein_length_aa)
    frame["encoded_amino_acid_end"] = (
        (frame["cumulative_cds_end"] + 2) // 3
    ).clip(upper=REFERENCE.protein_length_aa)
    frame["coding_frame"] = frame["splice_phase_5prime"].map(
        lambda value: (
            f"Frame {int(value)}"
            if pd.notna(value)
            else "Not available"
        )
    )

    frame["primary_protein_region"] = "Unassigned"
    if broad_domains:
        domains = list(broad_domains)
        for index, row in frame.iterrows():
            aa_start = int(row["encoded_amino_acid_start"])
            aa_end = int(row["encoded_amino_acid_end"])
            matches = [
                domain["domain_name"]
                for domain in domains
                if aa_start <= int(domain["protein_end_aa"])
                and aa_end >= int(domain["protein_start_aa"])
            ]
            frame.at[index, "primary_protein_region"] = (
                "; ".join(matches) if matches else "Unassigned"
            )

    return frame.reset_index(drop=True)


def genomic_reference_table(exons: pd.DataFrame) -> pd.DataFrame:
    """Return genomic coding-region columns with numeric coordinate dtypes."""
    table = pd.DataFrame(
        {
            "Exon": exons["exon_number"],
            "Chromosome": exons["chromosome"],
            "Coding Genomic Start": exons["coding_genomic_start"],
            "Coding Genomic End": exons["coding_genomic_end"],
            "Coding Genomic Length": exons["coding_genomic_length_bp"],
            "Full Exon Genomic Start": exons["genomic_start_grch38"],
            "Full Exon Genomic End": exons["genomic_end_grch38"],
            "Strand": exons["strand"],
            "Genome Assembly": REFERENCE.genome_assembly,
            "Transcript Accession": exons["transcript_id"],
            "Reference Sequence Accession": REFERENCE.refseq_transcript,
        }
    )
    return table


def transcript_reference_table(exons: pd.DataFrame) -> pd.DataFrame:
    """Return transcript/CDS fields supported by the local exon table."""
    return pd.DataFrame(
        {
            "Exon": exons["exon_number"],
            "Transcript Start": exons["transcript_exon_start"],
            "Transcript End": exons["transcript_exon_end"],
            "Total Transcript Length": exons["transcript_exon_length_bp"],
            "CDS Start": exons["cds_start"],
            "CDS End": exons["cds_end"],
            "Coding Length": exons["coding_length_bp"],
            "5′ UTR Length": exons["utr_5prime_length_bp"],
            "3′ UTR Length": exons["utr_3prime_length_bp"],
            "Coding Frame": exons["coding_frame"],
            "5′ Phase": exons["splice_phase_5prime"],
            "3′ Phase": exons["splice_phase_3prime"],
            "Encoded Amino-Acid Start": exons["encoded_amino_acid_start"],
            "Encoded Amino-Acid End": exons["encoded_amino_acid_end"],
            "Primary Protein Region": exons["primary_protein_region"],
        }
    )


def _feature_category(name: str) -> str:
    lowered = name.lower()
    if "actin binding domain" in lowered:
        return "Major regions"
    if "repeat region" in lowered:
        return "Spectrin repeats"
    if "hinge region" in lowered:
        return "Hinges"
    if "binding" in lowered:
        return "Binding sites"
    if any(token in lowered for token in ("ef-hand", "ww-domain", "zz-domain", "heptad")):
        return "Functional motifs"
    return "Major regions"


def _feature_colors(
    name: str,
    category: str,
) -> tuple[str, str]:
    lowered = name.lower()
    if "actin binding domain" in lowered or "unique n-terminus" in lowered:
        return "#FFCCCC", "#EB6F6F"
    if "central rod" in lowered:
        return "#FFFFCC", "#F0C266"
    if "cysteine-rich" in lowered:
        return "#DAE9F8", "#4D93D9"
    if "carboxy-terminal" in lowered:
        return "#F0E3F4", "#D86DCD"
    return CATEGORY_COLORS[category]


def _encoded_exon_range(
    aa_start: int,
    aa_end: int,
    exons: pd.DataFrame,
) -> tuple[int | None, int | None]:
    cds_start = (aa_start - 1) * 3 + 1
    cds_end = aa_end * 3
    overlap = exons.loc[
        (exons["cumulative_cds_start"] <= cds_end)
        & (exons["cumulative_cds_end"] >= cds_start)
    ]
    if overlap.empty:
        return None, None
    return int(overlap["exon_number"].min()), int(overlap["exon_number"].max())


def build_detailed_protein_features(
    feature_rows: Iterable[dict[str, Any]],
    exons: pd.DataFrame,
) -> pd.DataFrame:
    """Convert the Leiden transcript-nt feature map to explicit protein AA fields.

    The source file's numeric feature positions use Dp427m transcript
    nucleotides (5′ UTR 1–244; coding protein 245–11299), despite legacy
    headers naming the values amino-acid positions.
    """
    rows: list[dict[str, Any]] = []
    coding_nt_start = 245
    coding_nt_end = 11299

    for source_row in feature_rows:
        name = str(source_row["region_name"]).strip()
        source_start = int(source_row["first_amino_acid"])
        source_end = int(source_row["last_amino_acid"])
        clipped_start = max(source_start, coding_nt_start)
        clipped_end = min(source_end, coding_nt_end)
        if clipped_start > clipped_end:
            continue
        if name.lower() in {"codon start", "dp427m", "cds"}:
            continue

        aa_start = (clipped_start - coding_nt_start) // 3 + 1
        aa_end = (clipped_end - coding_nt_start) // 3 + 1
        aa_start = max(1, min(aa_start, REFERENCE.protein_length_aa))
        aa_end = max(aa_start, min(aa_end, REFERENCE.protein_length_aa))
        exon_start, exon_end = _encoded_exon_range(aa_start, aa_end, exons)
        category = _feature_category(name)
        default_fill, default_border = _feature_colors(name, category)
        fill = str(source_row.get("fill_color", "")).strip() or default_fill
        border = str(source_row.get("border_color", "")).strip() or default_border

        rows.append(
            {
                "Feature": name,
                "Feature Category": category,
                "Amino-Acid Start": aa_start,
                "Amino-Acid End": aa_end,
                "Length": aa_end - aa_start + 1,
                "Encoded Exon Start": exon_start,
                "Encoded Exon End": exon_end,
                "Fill Color": fill,
                "Border Color": border,
                "Source Transcript-nt Start": source_start,
                "Source Transcript-nt End": source_end,
                "Source": source_row.get("source", ""),
                "Source Version": source_row.get("source_version", ""),
                "Date Accessed": source_row.get("date_accessed", ""),
                "Notes": "AA coordinates derived from source transcript-nt feature positions.",
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["Amino-Acid Start", "Amino-Acid End", "Feature"],
        kind="stable",
    ).reset_index(drop=True)


def _base_layout(
    figure: go.Figure,
    *,
    height: int,
    x_title: str,
    x_range: list[int | float],
    y_ticks: tuple[list[float], list[str]],
    range_slider: bool = False,
) -> go.Figure:
    figure.update_layout(
        height=height,
        margin=dict(l=155, r=25, t=35, b=55),
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=False,
        dragmode="pan",
        hovermode="closest",
    )
    figure.update_xaxes(
        title=x_title,
        range=x_range,
        showgrid=True,
        gridcolor="#EEF2F7",
        zeroline=False,
        rangeslider=dict(visible=range_slider, thickness=0.08),
    )
    figure.update_yaxes(
        tickmode="array",
        tickvals=y_ticks[0],
        ticktext=y_ticks[1],
        range=[0.2, max(y_ticks[0]) + 0.7],
        showgrid=False,
        zeroline=False,
    )
    return figure


def create_genomic_viewer(exons: pd.DataFrame) -> go.Figure:
    """Create a browser-like genomic viewer in standard chromosome orientation."""
    figure = go.Figure()
    gene_start = int(exons["genomic_start_grch38"].min())
    gene_end = int(exons["genomic_end_grch38"].max())

    figure.add_trace(
        go.Scatter(
            x=[gene_start, gene_end],
            y=[3, 3],
            mode="lines",
            line=dict(color="#475569", width=2),
            hoverinfo="skip",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[gene_start, gene_end],
            y=[2, 2],
            mode="lines",
            line=dict(color="#6B7280", width=6),
            hovertemplate=(
                f"<b>DMD gene span</b><br>chrX:{gene_start:,}–{gene_end:,}"
                f"<br>{REFERENCE.genome_assembly}<br>Reverse strand<extra></extra>"
            ),
        )
    )
    figure.add_annotation(
        x=gene_start + (gene_end - gene_start) * 0.28,
        y=2.23,
        text="Transcription direction  ←  (reverse strand)",
        showarrow=False,
        font=dict(color="#7C2D12", size=12),
    )
    figure.add_trace(
        go.Scatter(
            x=[gene_start, gene_end],
            y=[1, 1],
            mode="lines",
            line=dict(color="#CBD5E1", width=2),
            hoverinfo="skip",
        )
    )

    hover_x: list[float] = []
    hover_text: list[str] = []
    for _, exon in exons.iterrows():
        start = int(exon["coding_genomic_start"])
        end = int(exon["coding_genomic_end"])
        figure.add_shape(
            type="rect",
            x0=start,
            x1=end,
            y0=0.72,
            y1=1.28,
            fillcolor="#4C78A8",
            line=dict(color="#1F4E79", width=1),
        )
        hover_x.append((start + end) / 2)
        hover_text.append(
            "<br>".join(
                [
                    f"<b>Exon {int(exon['exon_number'])}</b>",
                    f"Chromosome: {exon['chromosome']}",
                    f"Coding genomic: {start:,}–{end:,}",
                    f"Coding genomic length: {int(exon['coding_genomic_length_bp']):,} bp",
                    (
                        "Full exon genomic: "
                        f"{int(exon['genomic_start_grch38']):,}–"
                        f"{int(exon['genomic_end_grch38']):,}"
                    ),
                    f"Strand: {exon['strand']}",
                    f"Assembly: {REFERENCE.genome_assembly}",
                    f"Transcript: {exon['transcript_id']}",
                    (
                        "Transcript exon: "
                        f"{int(exon['transcript_exon_start']):,}–"
                        f"{int(exon['transcript_exon_end']):,}"
                    ),
                    f"CDS offset: {int(exon['cds_start']):,}–{int(exon['cds_end']):,}",
                ]
            )
        )
    figure.add_trace(
        go.Scatter(
            x=hover_x,
            y=[1] * len(hover_x),
            mode="markers",
            marker=dict(size=10, opacity=0),
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
        )
    )
    figure.add_annotation(
        x=gene_end,
        y=1.42,
        text="Exon 1 (5′ transcript end)",
        showarrow=True,
        arrowhead=2,
        ax=-55,
        ay=-25,
        font=dict(size=10),
    )
    return _base_layout(
        figure,
        height=390,
        x_title="Chromosome X genomic position (bp)",
        x_range=[gene_start, gene_end],
        y_ticks=([1, 2, 3], ["Dp427m coding exons", "DMD gene span", "Chromosome coordinates"]),
        range_slider=True,
    )


def _domain_parts_by_exon(
    exon_domain_rows: Iterable[dict[str, Any]],
) -> dict[int, list[tuple[float, float, str, str]]]:
    return {
        int(row["n"]): list(row["parts"])
        for row in exon_domain_rows
    }


def create_transcript_viewer(
    exons: pd.DataFrame,
    exon_domain_rows: Iterable[dict[str, Any]],
    *,
    schematic: bool,
) -> go.Figure:
    """Create a rectangular exon viewer in full mature-transcript coordinates."""
    figure = go.Figure()
    parts_by_exon = _domain_parts_by_exon(exon_domain_rows)
    hover_x: list[float] = []
    hover_text: list[str] = []

    for _, exon in exons.iterrows():
        number = int(exon["exon_number"])
        transcript_start = int(exon["transcript_exon_start"])
        transcript_end = int(exon["transcript_exon_end"])
        transcript_length = int(exon["transcript_exon_length_bp"])
        coding_start = int(exon["coding_transcript_start"])
        coding_end = int(exon["coding_transcript_end"])

        if schematic:
            exon_x0 = float(number - 1)
            exon_x1 = float(number)

            def map_position(position: int) -> float:
                return exon_x0 + (position - transcript_start) / transcript_length
        else:
            exon_x0 = float(transcript_start)
            exon_x1 = float(transcript_end + 1)

            def map_position(position: int) -> float:
                return float(position)

        if int(exon["utr_5prime_length_bp"]) > 0:
            figure.add_shape(
                type="rect",
                x0=exon_x0,
                x1=map_position(coding_start),
                y0=0.68,
                y1=1.32,
                fillcolor=UTR_FILL,
                line=dict(color=UTR_BORDER, width=1),
            )

        coding_x0 = map_position(coding_start)
        coding_x1 = map_position(coding_end + 1)
        coding_width = coding_x1 - coding_x0
        parts = parts_by_exon.get(number, [(0.0, 1.0, DEFAULT_FILL, DEFAULT_BORDER)])
        for fraction_start, fraction_end, fill, border in parts:
            figure.add_shape(
                type="rect",
                x0=coding_x0 + float(fraction_start) * coding_width,
                x1=coding_x0 + float(fraction_end) * coding_width,
                y0=0.68,
                y1=1.32,
                fillcolor=fill,
                line=dict(color=border, width=1),
            )

        if int(exon["utr_3prime_length_bp"]) > 0:
            figure.add_shape(
                type="rect",
                x0=map_position(coding_end + 1),
                x1=exon_x1,
                y0=0.68,
                y1=1.32,
                fillcolor=UTR_FILL,
                line=dict(color=UTR_BORDER, width=1),
            )

        figure.add_shape(
            type="line",
            x0=exon_x0,
            x1=exon_x0,
            y0=0.58,
            y1=1.42,
            line=dict(color="#475569", width=0.8),
        )
        center = (exon_x0 + exon_x1) / 2
        hover_x.append(center)
        hover_text.append(
            "<br>".join(
                [
                    f"<b>Exon {number}</b>",
                    f"Transcript: {transcript_start:,}–{transcript_end:,}",
                    f"Exon transcript length: {transcript_length:,} bp",
                    f"CDS offset: {int(exon['cds_start']):,}–{int(exon['cds_end']):,}",
                    f"Coding length: {int(exon['coding_length_bp']):,} bp",
                    f"5′ UTR: {int(exon['utr_5prime_length_bp']):,} bp",
                    f"3′ UTR: {int(exon['utr_3prime_length_bp']):,} bp",
                    (
                        "Genomic: "
                        f"{int(exon['genomic_start_grch38']):,}–"
                        f"{int(exon['genomic_end_grch38']):,}"
                    ),
                    f"Assembly / strand: {REFERENCE.genome_assembly} / {exon['strand']}",
                    (
                        "Encoded amino acids: "
                        f"{int(exon['encoded_amino_acid_start'])}–"
                        f"{int(exon['encoded_amino_acid_end'])}"
                    ),
                    f"Protein region: {exon['primary_protein_region']}",
                    f"Transcript accession: {exon['transcript_id']}",
                ]
            )
        )
        if schematic and (number == 1 or number == 79 or number % 5 == 0):
            figure.add_annotation(
                x=center,
                y=1,
                text=str(number),
                showarrow=False,
                font=dict(size=8, color="#1F2937"),
            )

    figure.add_trace(
        go.Scatter(
            x=hover_x,
            y=[1] * len(hover_x),
            mode="markers",
            marker=dict(size=12, opacity=0),
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
        )
    )
    figure.add_shape(
        type="line",
        x0=(0 if schematic else int(exons["transcript_exon_start"].min())),
        x1=(79 if schematic else int(exons["transcript_exon_end"].max()) + 1),
        y0=1,
        y1=1,
        line=dict(color="#CBD5E1", width=1),
        layer="below",
    )
    x_title = (
        "Exon order (5′ → 3′, schematic)"
        if schematic
        else "Dp427m mature transcript position (bp)"
    )
    x_range = (
        [0, 79]
        if schematic
        else [
            int(exons["transcript_exon_start"].min()),
            int(exons["transcript_exon_end"].max()) + 1,
        ]
    )
    return _base_layout(
        figure,
        height=260,
        x_title=x_title,
        x_range=x_range,
        y_ticks=([1], ["Spliced Dp427m transcript"]),
        range_slider=not schematic,
    )


def create_protein_overview(
    broad_domains: Iterable[dict[str, Any]],
) -> go.Figure:
    """Create a full-length protein orientation viewer in amino acids."""
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=[1, REFERENCE.protein_length_aa],
            y=[1, 1],
            mode="lines",
            line=dict(color="#CBD5E1", width=12),
            hoverinfo="skip",
        )
    )
    major_rows = [
        row
        for row in broad_domains
        if row["domain_name"] != "Hinge regions"
    ]
    for row in major_rows:
        name = row["domain_name"]
        start = int(row["protein_start_aa"])
        end = int(row["protein_end_aa"])
        default_fill, default_border = _feature_colors(name, "Major regions")
        fill = row.get("display_fill_hex") or default_fill
        border = row.get("display_border_hex") or default_border
        figure.add_trace(
            go.Scatter(
                x=[start, end],
                y=[1, 1],
                mode="lines",
                line=dict(color=fill, width=20),
                text=[
                    f"<b>{name}</b><br>AA {start:,}–{end:,}<br>Length: {end-start+1:,} aa"
                ] * 2,
                hovertemplate="%{text}<extra></extra>",
            )
        )
        figure.add_shape(
            type="line",
            x0=start,
            x1=end,
            y0=0.79,
            y1=0.79,
            line=dict(color=border, width=2),
        )
    return _base_layout(
        figure,
        height=230,
        x_title="Dp427m protein position (amino acids)",
        x_range=[1, REFERENCE.protein_length_aa],
        y_ticks=([1], ["Full-length Dp427m"]),
        range_slider=True,
    )


def create_protein_feature_viewer(features: pd.DataFrame) -> go.Figure:
    """Create separate amino-acid tracks for overlapping feature categories."""
    categories = [
        "Major regions",
        "Spectrin repeats",
        "Hinges",
        "Binding sites",
        "Functional motifs",
    ]
    y_by_category = {category: len(categories) - index for index, category in enumerate(categories)}
    figure = go.Figure()

    for _, feature in features.iterrows():
        category = str(feature["Feature Category"])
        y = y_by_category[category]
        start = int(feature["Amino-Acid Start"])
        end = int(feature["Amino-Acid End"])
        figure.add_shape(
            type="rect",
            x0=start,
            x1=end,
            y0=y - 0.28,
            y1=y + 0.28,
            fillcolor=feature["Fill Color"],
            line=dict(color=feature["Border Color"], width=1.2),
        )
        figure.add_trace(
            go.Scatter(
                x=[(start + end) / 2],
                y=[y],
                mode="markers",
                marker=dict(size=12, opacity=0),
                text=[
                    "<br>".join(
                        [
                            f"<b>{feature['Feature']}</b>",
                            f"Category: {category}",
                            f"Amino acids: {start:,}–{end:,}",
                            f"Length: {int(feature['Length']):,} aa",
                            (
                                "Encoded exons: "
                                f"{feature['Encoded Exon Start']}–"
                                f"{feature['Encoded Exon End']}"
                            ),
                            f"Source: {feature['Source']}",
                            f"Notes: {feature['Notes']}",
                        ]
                    )
                ],
                hovertemplate="%{text}<extra></extra>",
            )
        )

    return _base_layout(
        figure,
        height=420,
        x_title="Dp427m protein position (amino acids)",
        x_range=[1, REFERENCE.protein_length_aa],
        y_ticks=(
            [y_by_category[category] for category in categories],
            categories,
        ),
        range_slider=True,
    )


def create_style_summary_figure() -> go.Figure:
    """Create a compact visual legend for project styling classes."""
    styles = pd.DataFrame(STYLE_CLASSES)
    figure = go.Figure()
    for index, row in styles.iterrows():
        y = len(styles) - index
        figure.add_shape(
            type="rect",
            x0=0,
            x1=1,
            y0=y - 0.3,
            y1=y + 0.3,
            fillcolor=row["Fill Color"],
            line=dict(color=row["Border Color"], width=2),
        )
        figure.add_annotation(
            x=1.15,
            y=y,
            text=f"<b>{row['Style Category']}</b> — {row['Biological Meaning']}",
            showarrow=False,
            xanchor="left",
            font=dict(size=11),
        )
    figure.update_layout(
        height=330,
        margin=dict(l=25, r=20, t=15, b=15),
        xaxis=dict(visible=False, range=[0, 6]),
        yaxis=dict(visible=False, range=[0.3, len(styles) + 0.7]),
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=False,
    )
    return figure


def style_reference_table() -> pd.DataFrame:
    """Return one row per meaningful visual style class."""
    return pd.DataFrame(STYLE_CLASSES)


def resource_table(
    genomic_start: int,
    genomic_end: int,
) -> pd.DataFrame:
    """Return grouped, labeled external research resources."""
    ucsc_region = (
        "https://genome.ucsc.edu/cgi-bin/hgTracks?"
        f"db=hg38&position=chrX%3A{genomic_start}-{genomic_end}"
    )
    ensembl_region = (
        "https://www.ensembl.org/Homo_sapiens/Location/View?"
        f"r=X%3A{genomic_start}-{genomic_end}"
    )
    resources = [
        ("DMD-Specific Resources", "Leiden DMD Resources", "Leiden", "Gene / transcript", "DMD-specific transcript, exon, isoform, and reading-frame reference", "https://www.dmd.nl/DMD_home.html"),
        ("DMD-Specific Resources", "DMD Open-access Variant Explorer", "Leiden", "Transcript / variant", "DMD transcript-based variant interpretation and reading-frame context", "https://www.dmd.nl/DOVE/"),
        ("Genome and Coordinate Browsers", "Open DMD in UCSC (GRCh38)", "UCSC", "Genome", "Inspect the validated local DMD GRCh38 locus", ucsc_region),
        ("Genome and Coordinate Browsers", "Open UCSC Table Browser", "UCSC", "Genome", "Export and compare assembly-specific annotation tables", "https://genome.ucsc.edu/cgi-bin/hgTables"),
        ("Genome and Coordinate Browsers", "Open UCSC LiftOver", "UCSC", "Genome", "Convert coordinates explicitly between genome assemblies", "https://genome.ucsc.edu/cgi-bin/hgLiftOver"),
        ("Genome and Coordinate Browsers", "Open Ensembl DMD Gene", "Ensembl", "Genome / transcript", "Inspect DMD gene and transcript models", "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000198947"),
        ("Genome and Coordinate Browsers", "Open Ensembl Region", "Ensembl", "Genome", "Inspect the validated local GRCh38 DMD region", ensembl_region),
        ("NCBI Sequence and Variant Resources", "Open NCBI DMD Gene", "NCBI", "Gene", "Authoritative NCBI gene record and genomic context", "https://www.ncbi.nlm.nih.gov/gene/1756"),
        ("NCBI Sequence and Variant Resources", "Open Dp427m RefSeq", "NCBI", "Transcript", f"RefSeq mRNA record used by project configuration ({REFERENCE.refseq_transcript})", f"https://www.ncbi.nlm.nih.gov/nuccore/{REFERENCE.refseq_transcript}"),
        ("NCBI Sequence and Variant Resources", "Open NCBI Genome Data Viewer", "NCBI", "Genome", "Assembly-specific RefSeq genomic browsing", "https://www.ncbi.nlm.nih.gov/genome/gdv/browser/gene/?id=1756"),
        ("NCBI Sequence and Variant Resources", "Open NCBI BLAST", "NCBI", "Sequence validation", "Compare nucleotide or protein sequences; not a coordinate-table source", "https://blast.ncbi.nlm.nih.gov/Blast.cgi"),
        ("NCBI Sequence and Variant Resources", "Open Nucleotide BLAST", "NCBI", "Sequence validation", "Align nucleotide sequence against nucleotide databases", "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch"),
        ("NCBI Sequence and Variant Resources", "Open Protein BLAST", "NCBI", "Sequence validation", "Align dystrophin protein sequence or fragments", "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch"),
        ("NCBI Sequence and Variant Resources", "Open ClinVar DMD Search", "NCBI", "Clinical variants", "Clinical submissions and review status; not a normal exon-coordinate source", "https://www.ncbi.nlm.nih.gov/clinvar/?term=DMD%5Bgene%5D"),
        ("Protein Resources", "Open UniProt Dystrophin", "UniProt", "Protein", "Curated dystrophin sequence, domains, motifs, and binding sites", "https://www.uniprot.org/uniprotkb/P11532/entry"),
        ("Protein Resources", "Open InterPro Dystrophin", "EMBL-EBI", "Protein", "Domain signatures across multiple protein databases", "https://www.ebi.ac.uk/interpro/protein/UniProt/P11532/"),
        ("Protein Resources", "Open NCBI Dp427m Protein", "NCBI", "Protein", f"RefSeq protein record ({REFERENCE.refseq_protein})", f"https://www.ncbi.nlm.nih.gov/protein/{REFERENCE.refseq_protein}"),
        ("Aggregated Portals", "Open GeneCards DMD", "GeneCards", "Aggregated", "Aggregated gene portal; not the primary project coordinate source", "https://www.genecards.org/cgi-bin/carddisp.pl?gene=DMD"),
    ]
    return pd.DataFrame(
        resources,
        columns=[
            "Group",
            "Resource",
            "Organization",
            "Biological Level",
            "Purpose",
            "Link",
        ],
    )
