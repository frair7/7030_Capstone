#!/usr/bin/env python3
"""
find_exon_skips_expanded.py

Same logic as find_exon_skips.py, but lists every exon in a skip range
individually (e.g. "e4, e5, e6, e7") instead of using range notation
("e4-e7"). Useful when you need each skipped exon spelled out.

Usage:
    python find_exon_skips_expanded.py \
        --mutations data/DMD_Mutations_clean.csv \
        --light data/DMD_light.csv \
        --out data/exon_skip_results_expanded.csv
"""

import argparse
import csv
import re
from collections import defaultdict


def parse_exon_bounds(raw):
    """
    Extract the min/max exon numbers from strings like:
        'e5'          -> (5, 5, False)
        'e1-e7'       -> (1, 7, False)
        'e3-e5 '      -> (3, 5, False)
        'i11'         -> (11, 11, True)   <- intronic, flagged approximate
        'i7-e8'       -> (7, 8, True)     <- mixed, flagged approximate
        'e45-ei45'    -> (45, 45, True)   <- junctional, flagged approximate
    Returns (start, end, is_approx) or (None, None, True) if unparseable.
    """
    if not raw:
        return None, None, True

    raw = raw.strip()
    is_approx = "i" in raw.lower()  # any intron marker present

    # find all standalone numbers in the string
    numbers = [int(n) for n in re.findall(r"\d+", raw)]
    if not numbers:
        return None, None, True

    return min(numbers), max(numbers), is_approx


def load_light_table(path):
    """
    Load DMD_light.csv into two lists of rows (Deletion / Duplication),
    each row: dict with start, end, frame, early_stop, mw, junction_peptide.
    """
    table = defaultdict(list)
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        # normalize header names (strip whitespace)
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for row in reader:
            row = {k.strip(): (v.strip() if v else v) for k, v in row.items()}
            mtype = row.get("Deletion Duplication?", "")
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (ValueError, TypeError, KeyError):
                continue
            table[mtype].append({
                "start": start,
                "end": end,
                "frame": (row.get("frame") or "").lower(),
                "early_stop": row.get("early_stop"),
                "mw": row.get("MW"),
                "junction_peptide": row.get("junction_peptide"),
            })
    return table


def find_skip_candidates(mtype_table, mut_start, mut_end, top_n=3):
    """
    Find rows in mtype_table where [start, end] fully encloses
    [mut_start, mut_end] and frame == 'in'. Rank by:
      1. smallest total exon range (fewest extra exons skipped)
      2. closest fit to original mutation boundaries
    """
    candidates = []
    for row in mtype_table:
        if row["start"] <= mut_start and row["end"] >= mut_end and row["frame"] == "in":
            span = row["end"] - row["start"]
            closeness = (mut_start - row["start"]) + (row["end"] - mut_end)
            candidates.append((span, closeness, row))

    candidates.sort(key=lambda x: (x[0], x[1]))
    return [c[2] for c in candidates[:top_n]]


def format_exon_range(start, end):
    """
    Format an exon range as a comma-separated list of every exon in it,
    e.g. (4, 7) -> 'e4, e5, e6, e7' and (40, 40) -> 'e40'.
    """
    return ", ".join(f"e{i}" for i in range(start, end + 1))


def get_native_frame(mtype_table, mut_start, mut_end):
    """
    Look up the frame ('in'/'out') that DMD_light.csv assigns to the
    mutation's exact native exon range. Returns None if that exact
    range isn't present in the table.
    """
    for row in mtype_table:
        if row["start"] == mut_start and row["end"] == mut_end:
            return row["frame"]
    return None


def normalize_orig_frame(frame_orig, msub=None):
    """
    Determine the expected frame ('in'/'out'/None) for comparison against
    the light table's predicted frame.

    Some mutation subclasses have a definitional frame consequence that
    overrides whatever is literally written in the Frame column (since
    that column can contain data errors):
      - Nonsense (subexonic) mutations are, by definition, out-of-frame
        (premature stop / OOF).
      - Missense (subexonic) mutations are, by definition, in-frame
        (single amino acid substitution, no frame disruption).
    Any other subclass falls back to parsing the literal Frame column
    (IF/OOF/Unk/N/A/(blank)/etc.).
    """
    if msub:
        msub_l = msub.strip().lower()
        if msub_l == "nonsense":
            return "out"
        if msub_l == "missense":
            return "in"

    if not frame_orig:
        return None
    f = frame_orig.strip().lower()
    if f in ("if",):
        return "in"
    if f in ("oof",):
        return "out"
    return None  # Unk, N/A, (blank), or anything else - can't compare


def is_already_in_frame(mtype_table, mut_start, mut_end):
    """
    Check whether the mutation's own exact range is already in-frame
    according to the light table (no skip needed).
    """
    for row in mtype_table:
        if row["start"] == mut_start and row["end"] == mut_end and row["frame"] == "in":
            return True
    return False


def exon_range_is_modeled(mtype_table, mut_start, mut_end):
    """
    Check whether ANY row in mtype_table covers this mutation's range at all
    (regardless of frame). If not, the range likely isn't modeled in
    DMD_light.csv (e.g. exon 1 / start-codon-containing exon is often excluded).
    """
    for row in mtype_table:
        if row["start"] <= mut_start and row["end"] >= mut_end:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutations", default="data/DMD_Mutations_clean.csv")
    ap.add_argument("--light", default="data/DMD_light.csv")
    ap.add_argument("--out", default="data/exon_skip_results_expanded.csv")
    ap.add_argument("--top", type=int, default=3)
    args = ap.parse_args()

    light_table = load_light_table(args.light)

    results = []
    with open(args.mutations, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [h.strip() for h in reader.fieldnames]

        for row in reader:
            row = {k.strip(): (v.strip() if v else v) for k, v in row.items()}
            mut_id = row.get("Psuedo", "").strip()
            if not mut_id:
                continue

            mclass = row.get("Mutation Class", "")
            msub = row.get("Mutation Subclass", "")
            exon_raw = row.get("Exon(intron)", "")
            frame_orig = row.get("Frame", "")

            start, end, is_approx = parse_exon_bounds(exon_raw)

            note = ""
            if start is None:
                note = "Could not parse exon range - skipped"
                results.append([
                    mut_id, f"{mclass}/{msub} {exon_raw}",
                    "N/A", "N/A", *["" for _ in range(args.top * 2)], note
                ])
                continue
            if is_approx:
                note = "Intronic/junctional breakpoint - exon bounds are approximate"

            # Which DMD_light table to search: Duplication subclass -> Duplication table,
            # everything else (Deletion, Nonsense, Splice, Frameshift, Pseudoexon,
            # Junctional Deletion) is modeled as an effective deletion for skip purposes.
            lookup_type = "Duplication" if msub == "Duplication" else "Deletion"
            mtype_rows = light_table.get(lookup_type, [])

            msub_l = (msub or "").strip().lower()
            is_definitional = msub_l in ("nonsense", "missense")

            extra_notes = []

            if is_definitional:
                # Nonsense/Missense have a fixed frame consequence by definition
                # (premature stop = OOF, single-AA substitution = IF). This does
                # NOT depend on whether that exon's length happens to be a
                # multiple of 3, so we do not use the light table's native
                # single-exon-deletion lookup to determine it.
                definitional_frame = "out" if msub_l == "nonsense" else "in"
                predicted_frame = definitional_frame
                literal_frame = normalize_orig_frame(frame_orig, None)
                if literal_frame is None:
                    frame_match = "orig frame unclear"
                elif literal_frame == definitional_frame:
                    frame_match = "match"
                else:
                    frame_match = "MISMATCH"
                    extra_notes.append(
                        f"Frame column says '{frame_orig}' but {msub} mutations are "
                        f"'{definitional_frame}' by definition - used definitional value"
                    )

                if definitional_frame == "in":
                    # Missense: no frame disruption at all, so no skip is needed
                    # and no skip candidates should be searched or shown.
                    extra_notes.append(
                        "Missense mutations are in-frame by definition - no skip needed"
                    )
                    already_in_frame = True
                    candidates = []
                else:
                    # Nonsense: genuinely OOF - search the light table for
                    # real exon-skip rescue candidates around this exon.
                    already_in_frame = False
                    candidates = find_skip_candidates(mtype_rows, start, end, top_n=args.top)
            else:
                # Predicted frame from the light table for the mutation's native range
                predicted_frame = get_native_frame(mtype_rows, start, end)
                expected_frame = normalize_orig_frame(frame_orig, msub)
                if predicted_frame is None:
                    frame_match = "not modeled"
                elif expected_frame is None:
                    frame_match = "orig frame unclear"
                elif predicted_frame == expected_frame:
                    frame_match = "match"
                else:
                    frame_match = "MISMATCH"
                already_in_frame = is_already_in_frame(mtype_rows, start, end)
                candidates = find_skip_candidates(mtype_rows, start, end, top_n=args.top)

            skip_strs = [format_exon_range(c["start"], c["end"]) for c in candidates]
            skip_counts = [str(c["end"] - c["start"] + 1) for c in candidates]
            while len(skip_strs) < args.top:
                skip_strs.append("")
                skip_counts.append("")

            if already_in_frame and not is_definitional:
                extra_notes.append("Already in-frame at native boundaries - no skip needed")
            elif not already_in_frame and not candidates:
                if start == 1 or end == 1:
                    extra_notes.append(
                        "Not amenable to exon skipping - involves exon 1 "
                        "(contains start codon; no transcript remains to skip)"
                    )
                elif not exon_range_is_modeled(mtype_rows, start, end):
                    extra_notes.append(
                        "No candidates - exon range not present in DMD_light.csv - "
                        "data gap, review manually"
                    )
                else:
                    extra_notes.append("No in-frame rescue found within modeled range")

            if extra_notes:
                note = (note + "; " if note else "") + "; ".join(extra_notes)

            mutation_desc = f"{mclass}/{msub} {exon_raw} (orig frame: {frame_orig})"
            skip_cols = []
            for s, n in zip(skip_strs, skip_counts):
                skip_cols.extend([s, n])
            results.append([
                mut_id, mutation_desc,
                predicted_frame if predicted_frame else "N/A",
                frame_match,
                *skip_cols, note,
            ])

    header_skip_cols = []
    for i in range(args.top):
        header_skip_cols.extend([f"Skip_{i+1}", f"Skip_{i+1}_Exons_Skipped"])

    header = (
        ["ID", "Mutation", "Predicted_Frame", "Frame_Check"]
        + header_skip_cols
        + ["Notes"]
    )
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {args.out}")


if __name__ == "__main__":
    main()
