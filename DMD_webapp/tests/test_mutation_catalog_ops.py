"""Tests for mutation catalog data-management operations."""

from __future__ import annotations

import uuid

import pandas as pd
import pytest

from src.mutation_catalog import CATALOG_COLUMNS, ensure_record_ids
from src.mutation_catalog_ops import (
    apply_import_preview,
    apply_table_preview,
    export_catalog_csv,
    preview_csv_import,
    preview_table_edit,
    soft_delete_rows,
)


def _row(
    *,
    record_id: str = "",
    participant: str = "P1",
    start: str = "e3",
    stop: str = "e3",
    subclass: str = "Deletion",
) -> dict[str, str]:
    row = {col: "" for col in CATALOG_COLUMNS}
    row.update({
        "mutation_record_id": record_id,
        "id": participant,
        "mutation_class": "Exonic",
        "mutation_subclass": subclass,
        "start_region": start,
        "stop_region": stop,
    })
    return row


def _df(*rows: dict[str, str]) -> pd.DataFrame:
    out, _ = ensure_record_ids(pd.DataFrame(list(rows), columns=CATALOG_COLUMNS))
    return out


def test_export_includes_permanent_ids() -> None:
    df = _df(_row(record_id="11111111-1111-1111-1111-111111111111"))
    csv = export_catalog_csv(df).decode("utf-8")
    assert "mutation_record_id" in csv.splitlines()[0]
    assert "11111111-1111-1111-1111-111111111111" in csv


def test_reimport_with_no_edits_is_unchanged() -> None:
    stored = _df(_row(participant="S1", start="e5", stop="e5"))
    preview = preview_csv_import(stored, stored[CATALOG_COLUMNS].copy())
    assert preview.new_rows == 0
    assert preview.updated_rows == 0
    assert preview.unchanged_rows == len(stored)
    assert preview.can_apply


def test_edit_one_field_updates_correct_record() -> None:
    rec_id = str(uuid.uuid4())
    stored = _df(_row(record_id=rec_id, participant="S1", start="e3", stop="e3"))
    edited = stored.copy()
    edited.at[0, "phenotype"] = "DMD"
    preview = preview_table_edit(stored, edited)
    assert preview.updated_rows == 1
    assert preview.rows[0].mutation_record_id == rec_id
    saved = apply_table_preview(stored, preview, source="test")
    assert saved.loc[saved["mutation_record_id"] == rec_id, "phenotype"].iloc[0] == "DMD"


def test_csv_row_reorder_does_not_change_matching() -> None:
    r1 = _row(participant="A", start="e3", stop="e3")
    r2 = _row(participant="B", start="e5", stop="e5")
    stored, _ = ensure_record_ids(_df(r1, r2))
    upload = stored.iloc[::-1][CATALOG_COLUMNS]
    preview = preview_csv_import(stored, upload)
    assert preview.updated_rows == 0
    assert preview.unchanged_rows == 2


def test_blank_id_creates_new_record() -> None:
    stored = _df(_row(participant="S1", start="e1", stop="e1"))
    new_row = _row(participant="S2", start="e2", stop="e2")
    upload = pd.concat([stored[CATALOG_COLUMNS], pd.DataFrame([new_row])], ignore_index=True)
    preview = preview_csv_import(stored, upload)
    assert preview.new_rows == 1
    saved = apply_import_preview(stored, preview, source="test")
    assert len(saved) == len(stored) + 1


def test_unknown_id_is_rejected() -> None:
    stored = _df(_row(participant="S1"))
    upload = stored[CATALOG_COLUMNS].copy()
    upload.at[0, "mutation_record_id"] = str(uuid.uuid4())
    preview = preview_csv_import(stored, upload)
    assert preview.unknown_ids == 1
    assert not preview.can_apply


def test_duplicate_ids_are_rejected() -> None:
    rec_id = str(uuid.uuid4())
    stored = _df(_row(record_id=rec_id, participant="S1"))
    upload = pd.concat([stored[CATALOG_COLUMNS], stored[CATALOG_COLUMNS]], ignore_index=True)
    preview = preview_csv_import(stored, upload)
    assert preview.duplicate_ids >= 1
    assert not preview.can_apply


def test_table_edit_persists_new_row_with_id() -> None:
    stored = _df(_row(participant="S1", start="e1", stop="e1"))
    edited = pd.concat(
        [stored[CATALOG_COLUMNS], pd.DataFrame([_row(participant="S2", start="e2", stop="e2")])],
        ignore_index=True,
    )
    preview = preview_table_edit(stored, edited)
    assert preview.new_rows == 1
    saved = apply_table_preview(stored, preview, source="test")
    assert len(saved) == 2
    assert all(saved["mutation_record_id"].str.strip())


def test_soft_delete_only_targets_record_id() -> None:
    r1 = _row(participant="S1", start="e1", stop="e1")
    r2 = _row(participant="S1", start="e2", stop="e2")
    stored, _ = ensure_record_ids(_df(r1, r2))
    target = stored.iloc[0]["mutation_record_id"]
    saved = soft_delete_rows(stored, [target])
    active = saved[saved["is_deleted"].str.lower() != "true"]
    assert len(active) == 1
    assert active.iloc[0]["mutation_record_id"] != target


def test_invalid_upload_does_not_apply() -> None:
    stored = _df(_row(participant="S1"))
    bad = stored[CATALOG_COLUMNS].copy()
    bad.at[0, "mutation_class"] = ""
    preview = preview_csv_import(stored, bad)
    assert not preview.can_apply
    with pytest.raises(ValueError):
        apply_import_preview(stored, preview, source="test")


def test_multiple_mutations_same_participant_remain_distinct() -> None:
    r1 = _row(participant="S1", start="e3", stop="e3")
    r2 = _row(participant="S1", start="e4", stop="e5")
    stored, _ = ensure_record_ids(_df(r1, r2))
    ids = set(stored["mutation_record_id"])
    assert len(ids) == 2


def test_export_after_edit_contains_saved_values(tmp_path) -> None:
    rec_id = str(uuid.uuid4())
    stored = _df(_row(record_id=rec_id, participant="S1", start="e3", stop="e3"))
    edited = stored.copy()
    edited.at[0, "phenotype"] = "BMD"
    preview = preview_table_edit(stored, edited)
    saved = apply_table_preview(stored, preview, source="test")
    csv = export_catalog_csv(saved).decode("utf-8")
    assert rec_id in csv
    assert "BMD" in csv


def test_new_mutation_appears_in_export() -> None:
    stored = _df(_row(participant="S1", start="e1", stop="e1"))
    new_row = _row(participant="S2", start="e2", stop="e2")
    upload = pd.concat([stored[CATALOG_COLUMNS], pd.DataFrame([new_row])], ignore_index=True)
    preview = preview_csv_import(stored, upload)
    saved = apply_import_preview(stored, preview, source="test")
    csv = export_catalog_csv(saved).decode("utf-8")
    assert csv.count("S2") >= 1
    assert len(saved) == len(stored) + 1


def test_sorted_table_rows_match_by_record_id_not_position() -> None:
    r1 = _row(participant="A", start="e10", stop="e10")
    r2 = _row(participant="B", start="e3", stop="e3")
    stored, _ = ensure_record_ids(_df(r1, r2))
    rec_b = stored.loc[stored["id"] == "B", "mutation_record_id"].iloc[0]
    edited = stored.iloc[::-1].copy()
    edited.at[edited.index[0], "phenotype"] = "DMD"
    preview = preview_table_edit(stored, edited)
    updated = [r for r in preview.rows if r.status == "updated"]
    assert len(updated) == 1
    assert updated[0].mutation_record_id == rec_b


def test_direct_table_edit_persists_after_reload() -> None:
    rec_id = str(uuid.uuid4())
    stored = _df(_row(record_id=rec_id, participant="S1", start="e5", stop="e5"))
    edited = stored.copy()
    edited.at[0, "mutation_subclass"] = "Duplication"
    preview = preview_table_edit(stored, edited)
    saved = apply_table_preview(stored, preview, source="test")
    reloaded, _ = ensure_record_ids(saved.copy())
    assert reloaded.loc[reloaded["mutation_record_id"] == rec_id, "mutation_subclass"].iloc[0] == "Duplication"
