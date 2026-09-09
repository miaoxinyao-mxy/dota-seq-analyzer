#!/usr/bin/env python3
"""Generate human-readable comparison reports from canonical v3 results."""

import argparse
from collections import Counter
import json
from pathlib import Path

import pandas as pd

TAXONOMY_RANKS = (("R1", "domain"), ("P", "phylum"), ("C", "class"),
                  ("O", "order"), ("F", "family"), ("G", "genus"),
                  ("S", "species"))
CELL_COLUMNS = (
    "Barcode", "Predicted taxonomy", "Confidence", "Contamination",
    "Total # of 16s reads", "Technical noise count", "Reads_used_for_ASV",
    "Raw unique_core_sequences", "Dominant_raw_read_count",
    "Coexisting_2bp_reads", "Unauthorized_secondary_reads",
    "Final_cell_asv_reads", "Max_internal_distance", "Assigned_core_asv",
    "Status",
)
PV_COLUMNS = (
    "Barcode", "Target", "Mode", "Sequence_assignment", "Call", "Motif",
    "Repeat_count", "Start", "End", "Direct_similarity", "Evidence",
)
REFERENCE_COLUMNS = (
    "sub-ARG_name", "subject", "read1_start", "read1_end",
    "read1_percent_identity", "read2_start", "read2_end",
    "read2_percent_identity", "read1_seq", "read2_seq",
)


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _lineage(ranks: dict) -> str:
    return " | ".join(
        f"{code} - {ranks[field] if ranks[field] is not None else 'None'}"
        for code, field in TAXONOMY_RANKS)


def _sequence_sort_key(record: dict, target_order: dict[str, int]):
    sequence_id = record["sequence_cluster_id"]
    _, _, suffix = sequence_id.rpartition("_seq_")
    number = int(suffix) if suffix.isdigit() else suffix
    return target_order[record["parent_target"]], number


def _load_canonical(canonical_dir: Path):
    run_summary = _read_json(canonical_dir / "run_summary.json")
    cells = _read_jsonl(canonical_dir / "cells.jsonl")
    asvs = _read_jsonl(canonical_dir / "asvs.jsonl")
    target_path = canonical_dir / "target_sequences.jsonl"
    target_sequences = _read_jsonl(target_path) if target_path.exists() else None

    # 2026-09-09: Reject mixed canonical inputs before deriving reports.
    # Reason: report rows from different runs must never be combined silently.
    run_id = run_summary["run_id"]
    groups = (("cell", cells), ("ASV", asvs),
              ("target sequence", target_sequences or []))
    for label, records in groups:
        if any(record.get("run_id") != run_id for record in records):
            raise ValueError(f"{label} run_id does not match run_summary.json")
    return run_summary, cells, asvs, target_sequences


def _cell_target_matrix(run_summary: dict, cells: list[dict]) -> pd.DataFrame:
    target_names = [item["target_name"] for item in run_summary["target_panel"]]
    rows = []
    for cell in cells:
        asv = cell["asv"]
        evidence = cell["taxonomy_evidence"]
        row = {
            "Barcode": cell["cell_barcode"],
            "Predicted taxonomy": _lineage(cell["taxonomy"]["ranks"]),
            "Confidence": cell["taxonomy"]["confidence"],
            "Contamination": cell["taxonomy"]["contamination"],
            "Total # of 16s reads": evidence["total_16s_reads"],
            "Technical noise count": evidence["technical_noise_reads"],
            "Reads_used_for_ASV": asv["reads_used"],
            "Raw unique_core_sequences": asv["raw_unique_core_sequences"],
            "Dominant_raw_read_count": asv["dominant_raw_read_count"],
            "Coexisting_2bp_reads": asv["coexisting_2bp_reads"],
            "Unauthorized_secondary_reads": asv["unauthorized_secondary_reads"],
            "Final_cell_asv_reads": asv["final_cell_asv_reads"],
            "Max_internal_distance": asv["max_internal_distance"],
            "Assigned_core_asv": asv["asv_id"],
            "Status": asv["status"],
        }
        observed_targets = {item["target_name"]: item for item in cell["targets"]}
        for target_name in target_names:
            target = observed_targets.get(target_name)
            if target is None or target["assignment_type"] == "filtered_out":
                value = 0
            elif target["assignment_type"] == "sequence_cluster":
                value = target["sequence_cluster_id"]
            elif target["assignment_type"] == "sequence_unresolved":
                value = f"{target_name}_parent"
            else:
                value = target["filtered_read_count"]
            row[target_name] = value
        rows.append(row)
    return pd.DataFrame(rows, columns=[*CELL_COLUMNS, *target_names])


def _asv_summary(asvs: list[dict]) -> pd.DataFrame:
    ordered = sorted(asvs, key=lambda item: int(item["asv_id"].split("_")[-1]))
    rows = [{
        "Core_ASV_ID": item["asv_id"],
        "cell_count": item["final_surviving_cell_count"],
        "core_sequence": item["core_sequence"]["r1"] + "|" + item["core_sequence"]["r2"],
    } for item in ordered]
    return pd.DataFrame(rows, columns=["Core_ASV_ID", "cell_count", "core_sequence"])


def _taxonomy_summary(cells: list[dict]) -> pd.DataFrame:
    counts = Counter(_lineage(cell["taxonomy"]["ranks"]) for cell in cells)
    return pd.DataFrame(
        counts.items(), columns=["MLE_Taxonomic_Classification", "Cell_Count"]
    ).sort_values("Cell_Count", ascending=False, kind="stable", ignore_index=True)


def _target_summary(run_summary: dict) -> pd.DataFrame:
    rows = [{
        "ARG": item["target_name"],
        "Original": item["original_positive_cells"],
        "Filtered Out": item["filtered_out_cells"],
        "Remaining": item["remaining_positive_cells"],
        "% Retention": item["retention_percent"],
    } for item in run_summary["target_filtering_summary"]]
    return pd.DataFrame(
        rows, columns=["ARG", "Original", "Filtered Out", "Remaining", "% Retention"])


def _phase_variation_report(
    run_summary: dict, cells: list[dict], target_sequences: list[dict],
) -> pd.DataFrame:
    sequences = {item["sequence_cluster_id"]: item for item in target_sequences}
    target_order = [item["target_name"] for item in run_summary["target_panel"]]
    rows = []
    for cell in cells:
        by_target = {item["target_name"]: item for item in cell["targets"]}
        for target_name in target_order:
            target = by_target.get(target_name)
            if target is None or target["assignment_type"] != "sequence_cluster":
                continue
            sequence_id = target["sequence_cluster_id"]
            phase = sequences[sequence_id].get("phase_variation")
            if phase is None:
                continue
            rows.append({
                "Barcode": cell["cell_barcode"], "Target": target_name,
                "Mode": phase["mode"], "Sequence_assignment": sequence_id,
                "Call": phase["call"], "Motif": phase["motif"],
                "Repeat_count": phase["repeat_count"], "Start": phase["start"],
                "End": phase["end"],
                "Direct_similarity": phase["direct_similarity"],
                "Evidence": phase["evidence"],
            })
    return pd.DataFrame(rows, columns=PV_COLUMNS)


def _reference_report(run_summary: dict, target_sequences: list[dict]) -> pd.DataFrame:
    target_order = {item["target_name"]: index
                    for index, item in enumerate(run_summary["target_panel"])}
    rows = []
    ordered = sorted(
        target_sequences,
        key=lambda item: _sequence_sort_key(item, target_order),
    )
    for sequence in ordered:
        for match in sequence.get("reference_matches", []):
            rows.append({
                "sub-ARG_name": sequence["sequence_cluster_id"],
                "subject": match["subject"],
                "read1_start": match["r1_start"], "read1_end": match["r1_end"],
                "read1_percent_identity": f'{match["r1_percent_identity"]:g}%',
                "read2_start": match["r2_start"], "read2_end": match["r2_end"],
                "read2_percent_identity": f'{match["r2_percent_identity"]:g}%',
                "read1_seq": sequence["core_sequence"]["r1"],
                "read2_seq": sequence["core_sequence"]["r2"],
            })
    return pd.DataFrame(rows, columns=REFERENCE_COLUMNS)


def _run_overview(run_summary: dict) -> pd.DataFrame:
    read_qc = run_summary["read_qc"]
    classification = run_summary["read_classification"]
    funnel = run_summary["barcode_funnel"]
    rows = [
        ("reads", "raw_read_pairs", read_qc["raw_read_pairs"]),
        ("reads", "passed_read_pairs", read_qc["passed_read_pairs"]),
        ("reads", "filtered_read_pairs", read_qc["filtered_read_pairs"]),
        ("reads", "accepted_16s_reads", classification["accepted_16s_reads"]),
        ("reads", "target_reads", classification["target_reads"]),
        ("reads", "unclassified_reads", classification["unclassified_reads"]),
    ]
    rows.extend(("cells", metric, value) for metric, value in funnel.items())
    rows.append(("asvs", "final_asv_count",
                 run_summary["asv_summary"]["final_asv_count"]))
    rows.extend(
        ("targets", f'{item["target_name"]}.positive_cells',
         item["remaining_positive_cells"])
        for item in run_summary["target_filtering_summary"])
    if "phase_variation_summary" in run_summary:
        phase = run_summary["phase_variation_summary"]
        rows.append(("phase_variation", "cell_calls", phase["cell_calls"]))
        rows.extend(("phase_variation", f"calls.{call}", count)
                    for call, count in phase["calls"].items())
    return pd.DataFrame(rows, columns=["section", "metric", "value"])


def generate_reports(canonical_dir: str | Path, output_dir: str | Path) -> list[str]:
    canonical_dir = Path(canonical_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    run_summary, cells, asvs, target_sequences = _load_canonical(canonical_dir)
    reports = {
        "summary.tsv": _run_overview(run_summary),
        "cell_target_matrix.tsv": _cell_target_matrix(run_summary, cells),
        "asv_summary.tsv": _asv_summary(asvs),
        "taxonomy_summary.tsv": _taxonomy_summary(cells),
        "target_summary.tsv": _target_summary(run_summary),
    }
    modes = {item["mode"] for item in run_summary["target_panel"]}
    if "ssr" in modes:
        # 2026-09-09: Require the canonical sequence source for optional PV reports.
        # Reason: a requested report must not be synthesized from a missing canonical file.
        if target_sequences is None:
            raise ValueError(
                "target_sequences.jsonl is required for phase-variation reports")
        reports["cell_phase_variation.tsv"] = _phase_variation_report(
            run_summary, cells, target_sequences)
    if run_summary["inputs"]["reference_fasta"] is not None:
        # 2026-09-09: Require the canonical sequence source for optional reference reports.
        # Reason: reference annotations belong to target-sequence records, not run metadata alone.
        if target_sequences is None:
            raise ValueError(
                "target_sequences.jsonl is required for reference-match reports")
        reports["reference_matches.tsv"] = _reference_report(
            run_summary, target_sequences)
    for name, dataframe in reports.items():
        dataframe.to_csv(output_dir / name, sep="	", index=False)
    return list(reports)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate comparison reports from canonical v3 results.")
    parser.add_argument("--canonical-dir", default=".")
    parser.add_argument("--output-dir", default=".phase2a_reports")
    args = parser.parse_args()
    reports = generate_reports(args.canonical_dir, args.output_dir)
    print(f"Generated {len(reports)} canonical-derived reports in {args.output_dir}")


if __name__ == "__main__":
    main()
