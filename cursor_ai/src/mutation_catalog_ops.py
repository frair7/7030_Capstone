"""
Mutation catalog data-management: validation, diff, CSV import, and audit trail.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

import pandas as pd

from src.exon_display_config import PHENOTYPE_OPTIONS
from src.mutation_catalog import (
    CATALOG_COLUMNS,
    CATALOG_CSV,
    active_catalog,
    ensure_record_ids,
    new_mutation_record_id,
)
from src.mutation_viz import exon_range_from_row

AUDIT_LOG = CATALOG_CSV.parent / "mutation_catalog_audit.jsonl"

ALLOWED_MUTATION_CLASSES = {"", "Exonic", "Subexonic", "Intronic"}
ALLOWED_FRAMES = {"", "In-frame", "Out-of-frame", "Cannot determine", "N/A"}
READ_ONLY_FIELDS = frozenset({"mutation_record_id", "is_deleted", "deleted_at", "deleted_by", "deletion_reason"})


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _norm(val: Any) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    return str(val).strip()


def is_active_row(row: dict[str, Any]) -> bool:
    return _norm(row.get("is_deleted")).lower() not in ("true", "1", "yes")


def normalize_upload_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map CSV headers to canonical catalog columns."""
    rename = {c: c.strip() for c in df.columns}
    out = df.rename(columns=rename).copy()
    for col in CATALOG_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    extra = [c for c in out.columns if c not in CATALOG_COLUMNS]
    if extra:
        out = out.drop(columns=extra)
    return out[CATALOG_COLUMNS].fillna("").astype(str)


def validate_row(row: dict[str, Any], *, existing_ids: set[str]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []

    rec_id = _norm(row.get("mutation_record_id"))
    if rec_id and rec_id in existing_ids:
        errors.append(f"Duplicate mutation_record_id in batch: {rec_id}")

    if not _norm(row.get("id")):
        errors.append("Participant id is required")

    mclass = _norm(row.get("mutation_class"))
    if not mclass:
        errors.append("mutation_class is required")
    elif mclass not in ALLOWED_MUTATION_CLASSES:
        errors.append(f"Invalid mutation_class: {mclass}")

    frame = _norm(row.get("frame"))
    if frame and frame not in ALLOWED_FRAMES:
        warnings.append(f"Unusual frame value: {frame}")

    pheno = _norm(row.get("phenotype"))
    if pheno and pheno not in PHENOTYPE_OPTIONS:
        warnings.append(f"Unusual phenotype: {pheno}")

    if _norm(row.get("start_region")) or _norm(row.get("stop_region")):
        if exon_range_from_row(row) is None:
            errors.append("Invalid or unparseable exon region")

    pct = _norm(row.get("pct_dys_wb"))
    if pct:
        try:
            float(pct)
        except ValueError:
            errors.append(f"Invalid pct_dys_wb: {pct}")

    kda = _norm(row.get("expected_protein_size_kda"))
    if kda:
        try:
            float(kda)
        except ValueError:
            warnings.append(f"Invalid expected_protein_size_kda: {kda}")

    if not _norm(row.get("protein")):
        warnings.append("Missing optional protein field")

    return errors, warnings


def _field_changes(old: dict[str, Any], new: dict[str, Any]) -> list[dict[str, str]]:
    changes: list[dict[str, str]] = []
    for col in CATALOG_COLUMNS:
        if col in READ_ONLY_FIELDS:
            continue
        ov, nv = _norm(old.get(col)), _norm(new.get(col))
        if ov != nv:
            changes.append({"field": col, "existing": ov, "uploaded": nv})
    return changes


@dataclass
class ClassifiedRow:
    index: int
    mutation_record_id: str
    status: str
    row: dict[str, Any]
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    field_changes: list[dict[str, str]] = field(default_factory=list)


@dataclass
class CatalogChangePreview:
    total_rows: int
    new_rows: int
    updated_rows: int
    unchanged_rows: int
    invalid_rows: int
    duplicate_ids: int
    unknown_ids: int
    deleted_rows: int
    rows: list[ClassifiedRow]
    blocking_errors: list[str]

    @property
    def can_apply(self) -> bool:
        return not self.blocking_errors and self.invalid_rows == 0 and self.duplicate_ids == 0 and self.unknown_ids == 0


def _stored_by_id(stored: pd.DataFrame) -> dict[str, dict[str, Any]]:
    active = active_catalog(stored)
    return {
        _norm(r["mutation_record_id"]): r.to_dict()
        for _, r in active.iterrows()
        if _norm(r["mutation_record_id"])
    }


def preview_csv_import(stored: pd.DataFrame, uploaded: pd.DataFrame) -> CatalogChangePreview:
    uploaded = normalize_upload_columns(uploaded)
    stored_map = _stored_by_id(stored)
    seen_ids: set[str] = set()
    rows: list[ClassifiedRow] = []
    blocking: list[str] = []
    counts = dict(new_rows=0, updated_rows=0, unchanged_rows=0, invalid_rows=0,
                  duplicate_ids=0, unknown_ids=0)

    for i, raw in uploaded.iterrows():
        row = {c: _norm(raw[c]) for c in CATALOG_COLUMNS}
        rec_id = row["mutation_record_id"]
        errs: list[str] = []
        warns: list[str] = []

        if rec_id:
            if rec_id in seen_ids:
                counts["duplicate_ids"] += 1
                errs.append(f"Duplicate mutation_record_id in upload: {rec_id}")
            seen_ids.add(rec_id)

        ve, vw = validate_row(row, existing_ids=seen_ids - {rec_id} if rec_id else seen_ids)
        errs.extend(ve)
        warns.extend(vw)

        if not rec_id:
            status = "new"
            counts["new_rows"] += 1
            changes: list[dict[str, str]] = []
        elif rec_id not in stored_map:
            counts["unknown_ids"] += 1
            status = "unknown_id"
            errs.append(f"Unknown mutation_record_id: {rec_id}")
            changes = []
        else:
            changes = _field_changes(stored_map[rec_id], row)
            if changes:
                status = "updated"
                counts["updated_rows"] += 1
            else:
                status = "unchanged"
                counts["unchanged_rows"] += 1

        if errs:
            status = "invalid" if status != "unknown_id" else status
            if status == "invalid":
                counts["invalid_rows"] += 1

        rows.append(ClassifiedRow(i, rec_id, status, row, errs, warns,
                                  changes if rec_id in stored_map else []))

    if counts["duplicate_ids"]:
        blocking.append(f"{counts['duplicate_ids']} duplicate mutation_record_id(s) in upload")
    if counts["unknown_ids"]:
        blocking.append(f"{counts['unknown_ids']} unknown mutation_record_id(s)")
    if counts["invalid_rows"]:
        blocking.append(f"{counts['invalid_rows']} invalid row(s)")

    return CatalogChangePreview(
        total_rows=len(uploaded),
        deleted_rows=0,
        rows=rows,
        blocking_errors=blocking,
        **counts,
    )


def preview_table_edit(stored: pd.DataFrame, edited: pd.DataFrame) -> CatalogChangePreview:
    """Diff active stored catalog against edited table (add/update/delete)."""
    stored = ensure_record_ids(stored.copy())[0]
    edited = normalize_upload_columns(edited)
    stored_map = _stored_by_id(stored)
    edited_map: dict[str, dict[str, Any]] = {}
    seen: set[str] = set()
    rows: list[ClassifiedRow] = []
    blocking: list[str] = []
    counts = dict(new_rows=0, updated_rows=0, unchanged_rows=0, invalid_rows=0,
                  duplicate_ids=0, unknown_ids=0, deleted_rows=0)

    for i, raw in edited.iterrows():
        row = {c: _norm(raw[c]) for c in CATALOG_COLUMNS}
        rec_id = row["mutation_record_id"]
        errs: list[str] = []
        warns: list[str] = []

        if not rec_id:
            row["mutation_record_id"] = new_mutation_record_id()
            rec_id = row["mutation_record_id"]
            status = "new"
            counts["new_rows"] += 1
        elif rec_id in seen:
            counts["duplicate_ids"] += 1
            errs.append(f"Duplicate mutation_record_id: {rec_id}")
            status = "invalid"
        elif rec_id not in stored_map:
            counts["unknown_ids"] += 1
            errs.append(f"Unknown mutation_record_id: {rec_id}")
            status = "unknown_id"
        else:
            changes = _field_changes(stored_map[rec_id], row)
            if changes:
                status = "updated"
                counts["updated_rows"] += 1
            else:
                status = "unchanged"
                counts["unchanged_rows"] += 1

        seen.add(rec_id)
        edited_map[rec_id] = row
        ve, vw = validate_row(row, existing_ids=seen - {rec_id})
        errs.extend(ve)
        warns.extend(vw)
        if errs and status not in ("unknown_id",):
            status = "invalid"
            counts["invalid_rows"] += 1

        rows.append(ClassifiedRow(i, rec_id, status, row, errs, warns,
                                  _field_changes(stored_map[rec_id], row) if rec_id in stored_map else []))

    for rec_id, old in stored_map.items():
        if rec_id not in edited_map:
            counts["deleted_rows"] += 1
            del_row = {**old, "is_deleted": "true", "deleted_at": _utc_now_iso()}
            rows.append(ClassifiedRow(-1, rec_id, "deleted", del_row))

    if counts["duplicate_ids"]:
        blocking.append(f"{counts['duplicate_ids']} duplicate mutation_record_id(s)")
    if counts["unknown_ids"]:
        blocking.append(f"{counts['unknown_ids']} unknown mutation_record_id(s)")
    if counts["invalid_rows"]:
        blocking.append(f"{counts['invalid_rows']} invalid row(s)")

    return CatalogChangePreview(
        total_rows=len(edited),
        rows=rows,
        blocking_errors=blocking,
        **counts,
    )


def apply_import_preview(stored: pd.DataFrame, preview: CatalogChangePreview, *, source: str, changed_by: str = "") -> pd.DataFrame:
    if not preview.can_apply:
        raise ValueError("Cannot apply import with blocking errors")

    stored, _ = ensure_record_ids(stored.copy())
    stored_map = {
        _norm(r["mutation_record_id"]): r.to_dict()
        for _, r in stored.iterrows()
        if _norm(r["mutation_record_id"])
    }

    for item in preview.rows:
        if item.status not in ("new", "updated", "unchanged"):
            continue
        row = dict(item.row)
        if item.status == "new" and not _norm(row.get("mutation_record_id")):
            row["mutation_record_id"] = new_mutation_record_id()
        rec_id = _norm(row["mutation_record_id"])
        old = stored_map.get(rec_id, {})
        stored_map[rec_id] = row
        if item.status in ("new", "updated"):
            _append_audit(item.status, row, old, source, changed_by)

    out = pd.DataFrame(list(stored_map.values()), columns=CATALOG_COLUMNS)
    out, _ = ensure_record_ids(out)
    return out


def apply_table_preview(stored: pd.DataFrame, preview: CatalogChangePreview, *, source: str, changed_by: str = "") -> pd.DataFrame:
    if not preview.can_apply:
        raise ValueError("Cannot apply table changes with blocking errors")

    stored, _ = ensure_record_ids(stored.copy())
    stored_map = {
        _norm(r["mutation_record_id"]): r.to_dict()
        for _, r in stored.iterrows()
        if _norm(r["mutation_record_id"])
    }

    for item in preview.rows:
        rec_id = _norm(item.mutation_record_id)
        if item.status == "deleted":
            old = stored_map.get(rec_id, {})
            new_row = {**old, **item.row}
            stored_map[rec_id] = new_row
            _append_audit("deleted", new_row, old, source, changed_by)
        elif item.status in ("new", "updated", "unchanged"):
            old = stored_map.get(rec_id, {})
            stored_map[rec_id] = item.row
            if item.status in ("new", "updated"):
                _append_audit(item.status, item.row, old, source, changed_by)

    out = pd.DataFrame(list(stored_map.values()), columns=CATALOG_COLUMNS)
    out, _ = ensure_record_ids(out)
    return out


def soft_delete_rows(stored: pd.DataFrame, record_ids: list[str], *, changed_by: str = "", reason: str = "") -> pd.DataFrame:
    stored, _ = ensure_record_ids(stored.copy())
    ids = {_norm(i) for i in record_ids}
    for idx in stored.index:
        rec_id = _norm(stored.at[idx, "mutation_record_id"])
        if rec_id in ids:
            old = stored.loc[idx].to_dict()
            stored.at[idx, "is_deleted"] = "true"
            stored.at[idx, "deleted_at"] = _utc_now_iso()
            stored.at[idx, "deleted_by"] = changed_by
            stored.at[idx, "deletion_reason"] = reason
            _append_audit("deleted", stored.loc[idx].to_dict(), old, "table", changed_by)
    return stored


def export_catalog_csv(df: pd.DataFrame) -> bytes:
    out, _ = ensure_record_ids(df.copy())
    out = active_catalog(out)
    from src.mutation_viz import sort_catalog_dataframe
    out = sort_catalog_dataframe(out)
    return out[CATALOG_COLUMNS].to_csv(index=False).encode("utf-8")


def validation_errors_csv(preview: CatalogChangePreview) -> bytes:
    lines = ["row_index,mutation_record_id,status,errors,warnings"]
    for item in preview.rows:
        if item.errors or item.warnings:
            lines.append(
                f'{item.index},{item.mutation_record_id},{item.status},'
                f'{"|".join(item.errors)},{"|".join(item.warnings)}'
            )
    return "\n".join(lines).encode("utf-8")


def _append_audit(action: str, new_row: dict[str, Any], old_row: dict[str, Any], source: str, changed_by: str) -> None:
    entry = {
        "mutation_record_id": _norm(new_row.get("mutation_record_id")),
        "action": action,
        "changed_at": _utc_now_iso(),
        "changed_by": changed_by,
        "source": source,
        "previous": {k: _norm(old_row.get(k)) for k in CATALOG_COLUMNS if _norm(old_row.get(k)) != _norm(new_row.get(k))},
        "new": {k: _norm(new_row.get(k)) for k in CATALOG_COLUMNS if _norm(old_row.get(k)) != _norm(new_row.get(k))},
    }
    AUDIT_LOG.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT_LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry) + "\n")


def load_audit_for_record(record_id: str, *, limit: int = 50) -> list[dict[str, Any]]:
    if not AUDIT_LOG.exists():
        return []
    rid = _norm(record_id)
    entries: list[dict[str, Any]] = []
    with AUDIT_LOG.open(encoding="utf-8") as handle:
        for line in handle:
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue
            if _norm(entry.get("mutation_record_id")) == rid:
                entries.append(entry)
    return entries[-limit:]
