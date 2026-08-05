"""Verify transcript part edges match domain bubble edges."""

from __future__ import annotations

from src.domain_alignment import (
    build_aligned_domain_bubbles,
    draw_coords,
    hinge_guide_x_positions,
    resolve_domains,
    validate_domain_alignment,
)
from src.dp427m_exon_data import get_domain_map, get_exon_table
from src.part_segments import build_part_segments, part_rects_for_exon
from src.mutation_map_figure import LEFT_MARGIN, MIN_EXON_W, TOTAL_TRACK_WIDTH


def _layout() -> tuple[list[float], dict[int, tuple[float, float]]]:
    exon_table = get_exon_table()
    total_bp = sum(e["bp"] for e in exon_table)
    raw_w = [e["bp"] / total_bp * TOTAL_TRACK_WIDTH for e in exon_table]
    is_fixed = [w < MIN_EXON_W for w in raw_w]
    fixed_total = sum(MIN_EXON_W for f in is_fixed if f)
    raw_large_total = sum(w for w, f in zip(raw_w, is_fixed) if not f)
    remaining = max(10.0, TOTAL_TRACK_WIDTH - fixed_total)
    widths = [
        MIN_EXON_W if f else (w / raw_large_total * remaining)
        for w, f in zip(raw_w, is_fixed)
    ]
    exon_x: dict[int, tuple[float, float]] = {}
    x = LEFT_MARGIN
    for e, w in zip(exon_table, widths):
        exon_x[e["n"]] = (x, x + w)
        x += w
    return widths, exon_x


def test_h1_edges_match_part_segments() -> None:
    _, exon_x = _layout()
    exon_table = get_exon_table()
    bubbles = build_aligned_domain_bubbles(get_domain_map(), exon_table, exon_x)
    h1 = next(b for b in bubbles if b.label == "H1")

    e8 = next(e for e in exon_table if e["n"] == 8)
    rects8 = part_rects_for_exon(8, exon_x, e8["parts"])
    h1_part8 = next(r for r in rects8 if r[2].upper() == "#AEDCCA")
    assert abs(h1.x0 - h1_part8[0]) < 1e-6

    e10 = next(e for e in exon_table if e["n"] == 10)
    rects10 = part_rects_for_exon(10, exon_x, e10["parts"])
    h1_part10 = next(r for r in rects10 if r[2].upper() == "#AEDCCA")
    assert abs(h1.x1 - h1_part10[1]) < 1e-6


def test_r3_ends_where_h2_begins() -> None:
    """R3 must not include the trailing R4 sliver in exon 17."""
    _, exon_x = _layout()
    exon_table = get_exon_table()
    dm = get_domain_map()
    r3 = next(b for b in build_aligned_domain_bubbles(dm, exon_table, exon_x) if b.label == "R3")
    h2 = next(b for b in build_aligned_domain_bubbles(dm, exon_table, exon_x) if b.label == "H2")
    assert abs(r3.x1 - h2.x0) < 1e-6

    e17 = next(e for e in exon_table if e["n"] == 17)
    xa, xb = exon_x[17]
    h2_start_x = xa + e17["parts"][1][0] * (xb - xa)
    assert abs(r3.x1 - h2_start_x) < 1e-6


def test_r4_includes_exon_17_trailing_sliver() -> None:
    """R4 begins at the 3′ yellow segment in exon 17 after H2."""
    _, exon_x = _layout()
    exon_table = get_exon_table()
    dm = get_domain_map()
    r4 = next(b for b in build_aligned_domain_bubbles(dm, exon_table, exon_x) if b.label == "R4")
    h2 = next(b for b in build_aligned_domain_bubbles(dm, exon_table, exon_x) if b.label == "H2")
    assert abs(r4.x0 - h2.x1) < 1e-6


def test_hinge_guides_include_start_and_end() -> None:
    _, exon_x = _layout()
    xs = hinge_guide_x_positions(get_domain_map(), get_exon_table(), exon_x)
    assert len(xs) == 8
    by_label = {
        b.label: b
        for b in build_aligned_domain_bubbles(get_domain_map(), get_exon_table(), exon_x)
        if b.label in ("H1", "H2", "H3", "H4")
    }
    for label in ("H1", "H2", "H3", "H4"):
        b = by_label[label]
        assert any(abs(x - b.x0) < 1e-6 for x in xs)
        assert any(abs(x - b.x1) < 1e-6 for x in xs)


def test_repeat_bubbles_do_not_overlap() -> None:
    _, exon_x = _layout()
    bubbles = build_aligned_domain_bubbles(get_domain_map(), get_exon_table(), exon_x)
    repeats = [b for b in bubbles if b.is_repeat]
    assert len(repeats) == 24
    for prev, curr in zip(repeats, repeats[1:]):
        assert prev.x1 <= curr.x0 + 1e-6


def test_r1_starts_after_h1_r2_starts_after_r1() -> None:
    _, exon_x = _layout()
    bubbles = build_aligned_domain_bubbles(get_domain_map(), get_exon_table(), exon_x)
    h1 = next(b for b in bubbles if b.label == "H1")
    r1 = next(b for b in bubbles if b.label == "R1")
    r2 = next(b for b in bubbles if b.label == "R2")
    assert r1.x0 >= h1.x0
    assert r2.x0 >= r1.x1 - 1e-6
    assert r2.x0 >= exon_x[13][0] - 1e-6


def test_validate_domain_alignment_passes() -> None:
    _, exon_x = _layout()
    validate_domain_alignment(get_domain_map(), get_exon_table(), exon_x)


def test_draw_coords_inset_only_when_wide_enough() -> None:
    x0, x1 = draw_coords(10.0, 20.0, gap=0.75)
    assert x0 == 10.375
    assert x1 == 19.625
    narrow_x0, narrow_x1 = draw_coords(10.0, 10.5, gap=0.75)
    assert narrow_x0 == 10.0
    assert narrow_x1 == 10.5


def test_resolved_domains_have_draw_inset() -> None:
    _, exon_x = _layout()
    domains = resolve_domains(get_domain_map(), get_exon_table(), exon_x, gap=0.75)
    wide = next(d for d in domains if d.label == "R1")
    assert wide.draw_x0 > wide.biological_x0
    assert wide.draw_x1 < wide.biological_x1


def test_adjacent_domain_bubbles_share_edges() -> None:
    """H1 ends where R1 begins at exon-10 junction."""
    _, exon_x = _layout()
    bubbles = build_aligned_domain_bubbles(get_domain_map(), get_exon_table(), exon_x)
    h1 = next(b for b in bubbles if b.label == "H1")
    r1 = next(b for b in bubbles if b.label == "R1")
    assert abs(h1.x1 - r1.x0) < 1e-6
