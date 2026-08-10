#!/usr/bin/env python3
"""Call simple SSR and inversion phase variation from reconstructed target sequences."""

# 2026-08-10: Add one cell-aware phase-variation analysis module.
# Reason: primer-selected SSR and inversion targets should reuse reconstructed R1/R2 sequences.

import argparse
import csv
from difflib import SequenceMatcher
from pathlib import Path
import subprocess

import pandas as pd
from Bio import SeqIO
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
    "Inverted_similarity",
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
        if assignment == target or assignment.startswith(f"{target}_<"):
            return target
    return None


def load_sequence_records(sequence_list: str, target_modes):
    """Load reconstructed R1/R2 sequences for targets marked ssr or inv."""
    dataframe = pd.read_csv(sequence_list, sep="\t")
    records = {}
    for row in dataframe.itertuples(index=False):
        assignment = str(row[0])
        target = _target_for_assignment(assignment, target_modes)
        if target is None or target_modes[target] not in {"ssr", "inv"}:
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


def find_reference_covered_inv_targets(primers_file: str, reference: str, inv_targets):
    """Return INV targets whose two primer regions occur in one reference record."""
    reference_sequences = [
        str(record.seq).upper() for record in SeqIO.parse(reference, "fasta")
    ]
    covered = set()
    with open(primers_file, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        target_column = "Target" if reader.fieldnames and "Target" in reader.fieldnames else "Primer"
        for row in reader:
            target = (row.get(target_column) or "").strip()
            if target not in inv_targets:
                continue
            primer_f = (row.get("F") or "").strip().upper()
            primer_r = (row.get("R") or "").strip().upper()
            reverse_f = str(Seq(primer_f).reverse_complement())
            reverse_r = str(Seq(primer_r).reverse_complement())
            primer_pairs = [
                (primer_f, reverse_r),
                (primer_r, reverse_f),
                (reverse_f, primer_r),
                (reverse_r, primer_f),
            ]
            if any(
                left in sequence and right in sequence
                for sequence in reference_sequences
                for left, right in primer_pairs
            ):
                covered.add(target)
    return covered


def align_inversion_targets(records, reference: str, work_dir: str):
    """Align reconstructed INV sequences with BWA and return orientation evidence."""
    inv_records = {
        assignment: record
        for assignment, record in records.items()
        if record["mode"] == "inv"
    }
    if not inv_records:
        return {}

    # 2026-08-10: Build the BWA index under tmp instead of beside the user reference.
    # Reason: analysis must not modify a supplied reference FASTA or its directory.
    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)
    index_prefix = work_path / "reference"
    r1_fasta = work_path / "inv_R1.fasta"
    r2_fasta = work_path / "inv_R2.fasta"
    sam_path = work_path / "inv_alignments.sam"

    with r1_fasta.open("w", encoding="utf-8") as r1_handle, \
            r2_fasta.open("w", encoding="utf-8") as r2_handle:
        for assignment, record in inv_records.items():
            r1_handle.write(f">{assignment}\n{record['r1']}\n")
            r2_handle.write(f">{assignment}\n{record['r2']}\n")

    subprocess.run(
        ["bwa", "index", "-p", str(index_prefix), str(reference)], check=True)
    with sam_path.open("w", encoding="utf-8") as sam_handle:
        subprocess.run(
            ["bwa", "mem", str(index_prefix), str(r1_fasta), str(r2_fasta)],
            stdout=sam_handle,
            check=True,
        )

    alignments = {}
    with sam_path.open("r", encoding="utf-8") as sam_handle:
        for line in sam_handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            assignment = fields[0]
            flag = int(fields[1])
            if flag & 4:
                continue
            alignments.setdefault(assignment, []).append({
                "flag": flag,
                "reverse": bool(flag & 16),
                "supplementary": bool(flag & 2048) or any(
                    tag.startswith("SA:Z:") for tag in fields[11:]
                ),
            })

    evidence = {}
    for assignment in inv_records:
        hits = alignments.get(assignment, [])
        primary = [hit for hit in hits if not hit["flag"] & 256 and not hit["flag"] & 2048]
        supplementary = any(hit["supplementary"] for hit in hits)
        same_strand_pair = len(primary) >= 2 and len({hit["reverse"] for hit in primary[:2]}) == 1
        evidence[assignment] = {
            "inverted": supplementary or same_strand_pair,
            "description": (
                "supplementary_alignment" if supplementary
                else "same_strand_pair" if same_strand_pair
                else "reference_orientation" if primary
                else "unmapped"
            ),
        }
    return evidence


def call_sequence_phase_variation(records, inversion_evidence=None):
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
                    "inverted_similarity": None,
                    "evidence": "tandem_repeat_scan",
                }
            continue

        # 2026-08-10: Use the most abundant sequence as a provisional inversion orientation baseline.
        # Reason: real PV reference/SAM data are not yet available for calibrating breakpoint-aware calls.
        baseline_assignment, baseline_record = max(
            target_records, key=lambda item: item[1]["cell_count"])
        baseline = baseline_record["r1"] + str(Seq(baseline_record["r2"]).reverse_complement())
        for assignment, record in target_records:
            if inversion_evidence and assignment in inversion_evidence:
                evidence = inversion_evidence[assignment]
                calls[assignment] = {
                    "call": "inversion_detected" if evidence["inverted"] else "reference_orientation",
                    "motif": None,
                    "repeat_count": None,
                    "start": None,
                    "end": None,
                    "direct_similarity": None,
                    "inverted_similarity": None,
                    "evidence": evidence["description"],
                }
                continue

            candidate = record["r1"] + str(Seq(record["r2"]).reverse_complement())
            direct = SequenceMatcher(None, baseline, candidate, autojunk=False).ratio()
            inverted = SequenceMatcher(
                None, baseline, str(Seq(candidate).reverse_complement()), autojunk=False
            ).ratio()
            # 2026-08-10: Compare orientations directly without an unvalidated score cutoff.
            # Reason: real inversion data are required before setting a biological threshold.
            is_inverted = assignment != baseline_assignment and inverted > direct
            calls[assignment] = {
                "call": "inversion_detected" if is_inverted else "reference_orientation",
                "motif": None,
                "repeat_count": None,
                "start": None,
                "end": None,
                "direct_similarity": round(direct, 6),
                "inverted_similarity": round(inverted, 6),
                "evidence": "dominant_sequence_orientation",
            }
    return calls


def write_cell_calls(
    sequence_list: str, cell_matrix: str, primers_file: str, output_tsv: str,
    reference: str | None = None, work_dir: str = "tmp/phase_variation",
) -> int:
    """Map sequence-level PV calls back to cell barcodes and write a report."""
    target_modes = {
        target: mode
        for target, mode in get_target_modes(primers_file).items()
        if mode in {"ssr", "inv"}
    }
    records = load_sequence_records(sequence_list, target_modes)
    # 2026-08-10: Use BWA only when INV analysis and a reference are both available.
    # Reason: SSR does not require alignment and INV must still produce a provisional call without reference data.
    if reference:
        # 2026-08-10: Align only INV targets fully covered by one supplied reference record.
        # Reason: partial reference panels must not generate unsupported inversion calls.
        inv_targets = {target for target, mode in target_modes.items() if mode == "inv"}
        covered_inv_targets = find_reference_covered_inv_targets(
            primers_file, reference, inv_targets)
        skipped_targets = sorted(inv_targets - covered_inv_targets)
        if skipped_targets:
            print("Skipping INV targets without complete reference coverage: " + ", ".join(skipped_targets))
        records = {
            assignment: record
            for assignment, record in records.items()
            if record["mode"] != "inv" or record["target"] in covered_inv_targets
        }
        inversion_evidence = align_inversion_targets(records, reference, work_dir)
    else:
        inversion_evidence = {}
    sequence_calls = call_sequence_phase_variation(records, inversion_evidence)
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
                "Inverted_similarity": call["inverted_similarity"],
                "Evidence": call["evidence"],
            })

    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=OUTPUT_COLUMNS).to_csv(output_path, sep="\t", index=False)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Call cell-level SSR and inversion phase variation.")
    parser.add_argument("--sequence_list", required=True)
    parser.add_argument("--cell_matrix", required=True)
    parser.add_argument("--primers_file", required=True)
    parser.add_argument("--output_tsv", default="reports/cell_phase_variation.tsv")
    parser.add_argument("--reference")
    parser.add_argument("--work_dir", default="tmp/phase_variation")
    args = parser.parse_args()

    count = write_cell_calls(
        args.sequence_list, args.cell_matrix, args.primers_file, args.output_tsv,
        args.reference, args.work_dir)
    print(f"Wrote {count} cell-level phase-variation calls to {args.output_tsv}")


if __name__ == "__main__":
    main()
