#!/usr/bin/env python3
"""
find_exon_skips.py

For each mutation in DMD_Mutations_clean.csv, find the first 3 exon-skipping
combinations (from DMD_light.csv) that would restore the reading frame.

Usage:
    python find_exon_skips.py \
        --mutations data/DMD_Mutations_clean.csv \
        --light data/DMD_light.csv \
        --out data/exon_skip_results.csv
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutations", default="data/DMD_Mutations_clean.csv")
    ap.add_argument("--light", default="data/DMD_light.csv")
    ap.add_argument("--out", default="data/exon_skip_results.csv")
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
                results.append([mut_id, f"{mclass}/{msub} {exon_raw}", "", "", "", note])
                continue
            if is_approx:
                note = "Intronic/junctional breakpoint - exon bounds are approximate"

            # Which DMD_light table to search: Duplication subclass -> Duplication table,
            # everything else (Deletion, Nonsense, Splice, Frameshift, Pseudoexon,
            # Junctional Deletion) is modeled as an effective deletion for skip purposes.
            lookup_type = "Duplication" if msub == "Duplication" else "Deletion"
            mtype_rows = light_table.get(lookup_type, [])

            candidates = find_skip_candidates(mtype_rows, start, end, top_n=args.top)
            skip_strs = [f"e{c['start']}-e{c['end']}" for c in candidates]
            while len(skip_strs) < args.top:
                skip_strs.append("")

            mutation_desc = f"{mclass}/{msub} {exon_raw} (orig frame: {frame_orig})"
            results.append([mut_id, mutation_desc, *skip_strs, note])

    header = ["ID", "Mutation"] + [f"Skip_{i+1}" for i in range(args.top)] + ["Notes"]
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {args.out}")


if __name__ == "__main__":
    main()
