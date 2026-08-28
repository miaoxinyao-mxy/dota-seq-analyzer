#!/usr/bin/env python3
"""Call SSR phase variation from reconstructed target sequences."""

# 2026-08-11: Restrict phase-variation analysis to SSR targets.
# Reason: SSR is the only supported phase-variation mode.

import argparse
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from helper_functions import get_target_modes


OUTPUT_COLUMNS = [
    "Barcode",
    "Target",
    "Mode",
    "Sequence_assignment",
    "Call",
    "Motif",
    "Repeat_count",
    "Start",
    "End",
    "Direct_similarity",
    "Evidence",
]


def find_tandem_repeats(sequence: str):
    """Return non-overlapping motifs repeated at least twice in tandem."""
    repeats = []
    index = 0
    while index < len(sequence):
        best = None
        # 2026-08-10: Search every motif length that can repeat at least twice.
        # Reason: no PV dataset is available to justify motif-length or repeat-count cutoffs.
        for motif_length in range(1, (len(sequence) - index) // 2 + 1):
            motif = sequence[index:index + motif_length]
            repeat_count = 1
            while sequence[
                index + repeat_count * motif_length:
                index + (repeat_count + 1) * motif_length
            ] == motif:
                repeat_count += 1

            if repeat_count >= 2:
                candidate = {
                    "motif": motif,
                    "repeat_count": repeat_count,
                    "start": index,
                    "end": index + repeat_count * motif_length,
                }
                if best is None or candidate["end"] - candidate["start"] > best["end"] - best["start"]:
                    best = candidate
        if best is None:
            index += 1
        else:
            repeats.append(best)
            index = best["end"]
    return repeats


def _target_for_assignment(assignment: str, targets) -> str | None:
    """Match an arbitrary reconstructed-sequence name back to its primer target."""
    for target in sorted(targets, key=len, reverse=True):
        if assignment == target or assignment.startswith(f"{target}_seq_"):
            return target
    return None


def load_sequence_records(sequence_list: str, target_modes):
    """Load reconstructed R1/R2 sequences for targets marked ssr."""
    dataframe = pd.read_csv(sequence_list, sep="\t")
    records = {}
    for row in dataframe.itertuples(index=False):
        assignment = str(row[0])
        target = _target_for_assignment(assignment, target_modes)
        if target is None or target_modes[target] != "ssr":
            continue
        sequences = str(row[2]).split("|", 1)
        if len(sequences) != 2:
            continue
        records[assignment] = {
            "target": target,
            "mode": target_modes[target],
            "cell_count": int(row[1]),
            "r1": sequences[0],
            "r2": sequences[1],
        }
    return records


def call_sequence_phase_variation(records):
    """Return one provisional phase-variation call for each reconstructed sequence."""
    calls = {}
    by_target = {}
    for assignment, record in records.items():
        by_target.setdefault(record["target"], []).append((assignment, record))

    for target, target_records in by_target.items():
        mode = target_records[0][1]["mode"]
        if mode == "ssr":
            for assignment, record in target_records:
                combined = record["r1"] + str(Seq(record["r2"]).reverse_complement())
                repeats = find_tandem_repeats(combined)
                best = max(repeats, key=lambda item: item["end"] - item["start"], default=None)
                calls[assignment] = {
                    "call": "ssr_detected" if best else "no_ssr_detected",
                    "motif": best["motif"] if best else None,
                    "repeat_count": best["repeat_count"] if best else None,
                    "start": best["start"] if best else None,
                    "end": best["end"] if best else None,
                    "direct_similarity": None,
                    "evidence": "tandem_repeat_scan",
                }
            continue

    return calls


def write_cell_calls(
    sequence_list: str, cell_matrix: str, primers_file: str, output_tsv: str,
) -> int:
    """Map sequence-level PV calls back to cell barcodes and write a report."""
    target_modes = {
        target: mode
        for target, mode in get_target_modes(primers_file).items()
        if mode == "ssr"
    }
    records = load_sequence_records(sequence_list, target_modes)
    sequence_calls = call_sequence_phase_variation(records)
    matrix = pd.read_csv(cell_matrix, sep="\t", index_col="Barcode")

    rows = []
    for barcode, cell in matrix.iterrows():
        for target, mode in target_modes.items():
            assignment = str(cell.get(target, "0"))
            if assignment not in sequence_calls:
                continue
            call = sequence_calls[assignment]
            rows.append({
                "Barcode": barcode,
                "Target": target,
                "Mode": mode,
                "Sequence_assignment": assignment,
                "Call": call["call"],
                "Motif": call["motif"],
                "Repeat_count": call["repeat_count"],
                "Start": call["start"],
                "End": call["end"],
                "Direct_similarity": call["direct_similarity"],
                "Evidence": call["evidence"],
            })

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(output_path, sep="\t", index=False)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Call cell-level SSR phase variation.")
    parser.add_argument("--sequence_list", required=True)
    parser.add_argument("--cell_matrix", required=True)
    parser.add_argument("--primers_file", required=True)
    parser.add_argument("--output_tsv", default="reports/cell_phase_variation.tsv")
    args = parser.parse_args()

    count = write_cell_calls(
        args.sequence_list, args.cell_matrix, args.primers_file, args.output_tsv)
    print(f"Wrote {count} cell-level phase-variation calls to {args.output_tsv}")


if __name__ == "__main__":
    main()
