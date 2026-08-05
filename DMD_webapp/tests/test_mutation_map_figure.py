"""Tests for cohort mutation map."""

from __future__ import annotations

from src.coordinate_mapper import build_exon_records
from src.dp427m_exon_data import bar_style_for_phenotype, get_domain_map, get_exon_table
from src.mutation_map_figure import create_cohort_mutation_map


def test_exon_table_loaded() -> None:
    table = get_exon_table()
    assert len(table) == 79
    assert table[0]["n"] == 1
    assert len(table[7]["parts"]) == 2  # exon 8 ABD/H1 split


def test_domain_map_loaded() -> None:
    domains = get_domain_map()
    assert len(domains) == 31
    assert domains[0][2] == "Actin Binding Domain"


def test_phenotype_bar_styles() -> None:
    assert bar_style_for_phenotype("DMD") == "black"
    assert bar_style_for_phenotype("BMD") == "grey"
    assert bar_style_for_phenotype("Asymptomatic") == "outline"


def test_cohort_map_builds() -> None:
    exons = build_exon_records()
    row = {
        "id": "P-6",
        "start_region": "e3",
        "stop_region": "e29",
        "expected_protein_size_kda": "273.4",
        "phenotype": "DMD",
        "group": "B",
        "pct_dys_wb": "49.16",
    }
    fig = create_cohort_mutation_map([row], exons)
    assert fig.get_figwidth() >= 20


def test_cr_domain_colors() -> None:
    table = get_exon_table()
    e65 = next(e for e in table if e["n"] == 65)
    assert e65["parts"][0][2] == "#DAE9F8"
    assert e65["parts"][0][3] == "#4D93D9"
    domains = get_domain_map()
    cr = next(d for d in domains if d[2] == "CR Domain")
    assert cr[3] == "#DAE9F8"


def test_multi_mutation_map() -> None:
    exons = build_exon_records()
    rows = [
        {"id": "P-5", "start_region": "e3", "stop_region": "e27", "phenotype": "BMD", "group": "B"},
        {"id": "P-6", "start_region": "e3", "stop_region": "e29", "phenotype": "DMD", "group": "B"},
    ]
    fig = create_cohort_mutation_map(rows, exons)
    assert len(fig.axes[0].patches) > 50
