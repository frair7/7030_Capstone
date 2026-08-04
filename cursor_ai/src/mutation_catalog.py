"""
Mutation catalog — load, save, filter, and manage cataloged mutations.
"""

from __future__ import annotations

import csv
import re
import uuid
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from src.config import REFERENCE
from src.mutation_viz import sort_catalog_dataframe

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
CATALOG_CSV = DATA_DIR / "mutation_catalog.csv"
PARENT_CATALOG = (
    Path(__file__).resolve().parents[2] / "data" / "input" / "DMD_Mutations_clean.csv"
)

CATALOG_COLUMNS = [
    "mutation_record_id",
    "id",
    "mutation_class",
    "mutation_subclass",
    "start_region",
    "stop_region",
    "frame",
    "phenotype",
    "group",
    "pct_dys_wb",
    "cdna",
    "rna",
    "protein",
    "domains_affected",
    "expected_protein_size_kda",
    "molecular_consequence",
    "experimental_research_category",
    "tissue_category_comments",
    "general_comments",
    "is_deleted",
    "deleted_at",
    "deleted_by",
    "deletion_reason",
]

DISPLAY_COLUMNS = (
    ["selected", "delete_row"]
    + [c for c in CATALOG_COLUMNS if c != "mutation_record_id"]
    + ["mutation_record_id"]
)

EXPLORER_DISPLAY_COLUMNS = [
    "selected",
    "id",
    "mutation_class",
    "mutation_subclass",
    "start_region",
    "stop_region",
    "frame",
    "phenotype",
    "group",
    "mutation_record_id",
]

CATALOG_EDITOR_COLUMNS = (
    ["delete_row"]
    + [c for c in CATALOG_COLUMNS if c != "mutation_record_id"]
    + ["mutation_record_id"]
)


def new_mutation_record_id() -> str:
    return str(uuid.uuid4())


def _ensure_schema(df: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    """Add missing columns and assign mutation_record_id where absent."""
    out = df.copy()
    changed = False
    for col in CATALOG_COLUMNS:
        if col not in out.columns:
            out[col] = ""
            changed = True
    out = out[CATALOG_COLUMNS].fillna("").astype(str)
    for idx in out.index:
        if not str(out.at[idx, "mutation_record_id"]).strip():
            out.at[idx, "mutation_record_id"] = new_mutation_record_id()
            changed = True
    return out, changed


def _normalize_frame(raw: str) -> str:
    val = (raw or "").strip()
    mapping = {
        "IF": "In-frame",
        "OOF": "Out-of-frame",
        "IN-FRAME": "In-frame",
        "OUT-OF-FRAME": "Out-of-frame",
        "UNK": "Cannot determine",
        "UNKNOWN": "Cannot determine",
        "N/A": "N/A",
        "(BLANK)": "",
        "BLANK": "",
    }
    return mapping.get(val.upper(), val)


def _split_legacy_region(raw: str) -> tuple[str, str]:
    """Split legacy Exon(intron) field into start/stop tokens."""
    text = (raw or "").strip()
    if not text:
        return "", ""

    junction = re.match(r"^\s*(e|i)(\d+)\s*[-–]\s*(e|i)(\d+)\s*$", text, re.I)
    if junction:
        p1, n1, p2, n2 = junction.groups()
        return f"{p1.lower()}{n1}", f"{p2.lower()}{n2}"

    rng = re.match(r"^\s*e?(\d+)\s*[-–]\s*e?(\d+)\s*$", text, re.I)
    if rng:
        return f"e{rng.group(1)}", f"e{rng.group(2)}"

    single = re.match(r"^\s*(e|i)(\d+)\s*$", text, re.I)
    if single:
        tok = f"{single.group(1).lower()}{single.group(2)}"
        return tok, tok

    return text, text


def import_legacy_catalog(path: Optional[Path] = None) -> pd.DataFrame:
    """Import from parent-repo DMD_Mutations_clean.csv."""
    src = path or PARENT_CATALOG
    rows: list[dict[str, str]] = []
    with src.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        reader.fieldnames = [h.strip() for h in (reader.fieldnames or [])]
        for row in reader:
            row = {k.strip(): (v.strip() if v else v) for k, v in row.items()}
            region = row.get("Exon(intron)", "") or row.get("Exon (intron)", "")
            start_r, stop_r = _split_legacy_region(region)
            rows.append(
                {
                    "mutation_record_id": new_mutation_record_id(),
                    "id": row.get("Psuedo", row.get("Pseudo", "")).strip(),
                    "mutation_class": (row.get("Mutation Class") or "").strip(),
                    "mutation_subclass": (row.get("Mutation Subclass") or "").strip(),
                    "start_region": start_r,
                    "stop_region": stop_r,
                    "frame": _normalize_frame(row.get("Frame", "")),
                    "phenotype": "",
                    "cdna": "",
                    "rna": "",
                    "protein": "",
                    "domains_affected": "",
                    "expected_protein_size_kda": "",
                    "molecular_consequence": "",
                    "experimental_research_category": "",
                    "tissue_category_comments": "",
                    "general_comments": "",
                    "is_deleted": "",
                    "deleted_at": "",
                    "deleted_by": "",
                    "deletion_reason": "",
                }
            )
    return pd.DataFrame(rows, columns=CATALOG_COLUMNS)


def load_catalog(path: Optional[Path] = None) -> pd.DataFrame:
    """Load mutation catalog; create from legacy import if missing."""
    csv_path = path or CATALOG_CSV
    if not csv_path.exists():
        df = import_legacy_catalog()
        save_catalog(df, csv_path)
        return df
    df = pd.read_csv(csv_path, dtype=str).fillna("")
    df, changed = _ensure_schema(df)
    if changed:
        save_catalog(df, csv_path)
    return df


def save_catalog(df: pd.DataFrame, path: Optional[Path] = None) -> None:
    """Persist catalog to CSV."""
    csv_path = path or CATALOG_CSV
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    out = df[CATALOG_COLUMNS].copy()
    out.to_csv(csv_path, index=False)


def next_id(df: pd.DataFrame) -> str:
    """Generate next subject ID like S102."""
    if df.empty:
        return "S1"
    nums = []
    for val in df["id"]:
        m = re.match(r"[Ss]?(\d+)", str(val))
        if m:
            nums.append(int(m.group(1)))
    n = max(nums, default=0) + 1
    return f"S{n}"


def ensure_record_ids(df: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    return _ensure_schema(df)


def active_catalog(df: pd.DataFrame) -> pd.DataFrame:
    """Return non-soft-deleted rows."""
    if "is_deleted" not in df.columns:
        return df.copy()
    deleted = df["is_deleted"].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    return df[~deleted].reset_index(drop=True)


def filter_catalog(
    df: pd.DataFrame,
    *,
    mutation_classes: Optional[list[str]] = None,
    mutation_subclasses: Optional[list[str]] = None,
    frames: Optional[list[str]] = None,
    search_text: str = "",
) -> pd.DataFrame:
    """Apply filters to the catalog table."""
    out = active_catalog(df.copy())
    if mutation_classes:
        out = out[out["mutation_class"].isin(mutation_classes)]
    if mutation_subclasses:
        out = out[out["mutation_subclass"].isin(mutation_subclasses)]
    if frames:
        out = out[out["frame"].isin(frames)]
    if search_text.strip():
        q = search_text.strip().lower()
        mask = out.apply(
            lambda row: q in " ".join(str(v).lower() for v in row.values),
            axis=1,
        )
        out = out[mask]
    return sort_catalog_dataframe(out)


def catalog_for_display(df: pd.DataFrame) -> pd.DataFrame:
    """Add selection columns for the data editor (legacy combined view)."""
    display = df.copy()
    display.insert(0, "selected", False)
    display.insert(1, "delete_row", False)
    return display[DISPLAY_COLUMNS]


def catalog_for_plot_selection(df: pd.DataFrame) -> pd.DataFrame:
    """Add Plot checkbox column for Mutation Explorer selection table."""
    display = df.copy()
    if "selected" not in display.columns:
        display.insert(0, "selected", False)
    return display[EXPLORER_DISPLAY_COLUMNS]


def catalog_for_editor(df: pd.DataFrame) -> pd.DataFrame:
    """Add delete checkbox column for Mutation Catalog editing table."""
    display = df.copy()
    if "delete_row" not in display.columns:
        display.insert(0, "delete_row", False)
    return display[CATALOG_EDITOR_COLUMNS]
