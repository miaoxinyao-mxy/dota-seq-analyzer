#!/usr/bin/env python3
"""Build and validate the version 3.0.0 canonical result set."""

import argparse
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
from typing import Any

import pandas as pd

SCHEMA_VERSION = "3.0.0"
TAXONOMY_RANKS = {
    "R1": "domain", "P": "phylum", "C": "class", "O": "order",
    "F": "family", "G": "genus", "S": "species",
}
FINAL_ASV_STATUSES = {"single_ASV", "corrupted_single_ASV"}
ASSIGNMENT_TYPES = {
    "filtered_out", "detected", "sequence_unresolved", "sequence_cluster",
}


def _native(value: Any):
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def _integer(value: Any) -> int:
    value = _native(value)
    if value is None:
        return 0
    return int(float(value))


def _number_or_none(value: Any):
    value = _native(value)
    return None if value is None else float(value)


def _parse_taxonomy(lineage: Any) -> dict[str, str | None]:
    ranks = {rank: None for rank in TAXONOMY_RANKS.values()}
    if pd.isna(lineage):
        return ranks
    for item in str(lineage).split(" | "):
        code, separator, taxon = item.partition(" - ")
        if separator and code in TAXONOMY_RANKS:
            ranks[TAXONOMY_RANKS[code]] = None if taxon == "None" else taxon
    return ranks


def _read_json(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _count_lines(path: str | Path) -> int:
    with open(path, encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def _read_stagger_counts(path: str | Path) -> Counter:
    counts = Counter()
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n")
        if header != "read_index\tread_id\tr1_primer_start":
            raise ValueError(f"Invalid 16S manifest header: {header!r}")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) != 3:
                raise ValueError(f"Invalid 16S manifest row: {line!r}")
            counts[int(values[2])] += 1
    return counts


def _load_sequence_catalog(path: str | Path | None) -> dict[str, dict]:
    if not path:
        return {}
    dataframe = pd.read_csv(path, sep="\t")
    catalog = {}
    for row in dataframe.itertuples(index=False, name=None):
        sequence_id = str(row[0])
        paired = str(row[2]).split("|", 1)
        if len(paired) != 2:
            raise ValueError(f"Invalid paired sequence for {sequence_id}")
        parent, separator, number = sequence_id.rpartition("_seq_")
        if not separator or not parent or not number.isdigit() or int(number) < 1:
            raise ValueError(f"Invalid reconstructed sequence ID: {sequence_id}")
        if sequence_id in catalog:
            raise ValueError(f"Duplicate reconstructed sequence ID: {sequence_id}")
        catalog[sequence_id] = {
            "parent_target": parent,
            "r1": paired[0],
            "r2": paired[1],
        }
    return catalog


def _load_phase_variation(path: str | Path | None) -> dict[str, dict]:
    if not path:
        return {}
    dataframe = pd.read_csv(path, sep="\t")
    by_sequence = {}
    for row in dataframe.to_dict(orient="records"):
        sequence_id = str(row["Sequence_assignment"])
        call = {
            "mode": str(row["Mode"]),
            "call": str(row["Call"]),
            "motif": _native(row["Motif"]),
            "repeat_count": _native(row["Repeat_count"]),
            "start": _native(row["Start"]),
            "end": _native(row["End"]),
            "direct_similarity": _native(row["Direct_similarity"]),
            "evidence": str(row["Evidence"]),
        }
        for key in ("repeat_count", "start", "end"):
            if call[key] is not None:
                call[key] = int(call[key])
        previous = by_sequence.setdefault(sequence_id, call)
        if previous != call:
            raise ValueError(f"Inconsistent phase-variation calls for {sequence_id}")
    return by_sequence


def _percent(value: Any) -> float:
    return float(str(value).rstrip("%"))


def _load_reference_matches(path: str | Path | None) -> dict[str, list[dict]]:
    if not path:
        return {}
    dataframe = pd.read_csv(path, sep="\t")
    matches = defaultdict(list)
    for row in dataframe.to_dict(orient="records"):
        matches[str(row["sub-ARG_name"])].append({
            "subject": str(row["subject"]),
            "r1_start": int(row["read1_start"]),
            "r1_end": int(row["read1_end"]),
            "r1_percent_identity": _percent(row["read1_percent_identity"]),
            "r2_start": int(row["read2_start"]),
            "r2_end": int(row["read2_end"]),
            "r2_percent_identity": _percent(row["read2_percent_identity"]),
        })
    return dict(matches)


def _load_sequence_qc(path: str | Path | None, target_names: list[str]) -> dict:
    if not path:
        return {}
    dataframe = pd.read_csv(path, sep="\t", index_col="Barcode")
    result = {}
    for barcode, row in dataframe.iterrows():
        per_target = {}
        for target in target_names:
            columns = {
                "confidence": f"{target}: Confidence",
                "contamination": f"{target}: Contamination",
                "total_reads": f"{target}: Total_#_of_ARG_reads",
                "technical_noise_reads": f"{target}: Technical_noise_count",
            }
            if columns["confidence"] not in dataframe.columns:
                continue
            values = {key: _native(row[column]) for key, column in columns.items()}
            if values["confidence"] in {None, "-"}:
                continue
            per_target[target] = {
                "confidence": float(values["confidence"]),
                "contamination": float(values["contamination"]),
                "total_reads": int(float(values["total_reads"])),
                "technical_noise_reads": int(float(values["technical_noise_reads"])),
            }
        if per_target:
            result[str(barcode)] = per_target
    return result


def _load_asv_sequences(path: str | Path) -> dict[str, dict]:
    dataframe = pd.read_csv(path, sep="\t")
    result = {}
    for row in dataframe.itertuples(index=False, name=None):
        asv_id = str(row[0])
        paired = str(row[2]).split("|", 1)
        if len(paired) != 2:
            raise ValueError(f"Invalid paired ASV sequence for {asv_id}")
        result[asv_id] = {"r1": paired[0], "r2": paired[1]}
    return result


def _target_assignment(
    target: str, raw_count: int, filtered_count: int, final_value: Any,
    sequence_catalog: dict[str, dict],
) -> tuple[str, str | None]:
    if filtered_count == 0:
        return "filtered_out", None
    value = str(_native(final_value))
    if value == f"{target}_parent":
        return "sequence_unresolved", None
    if value in sequence_catalog:
        if sequence_catalog[value]["parent_target"] != target:
            raise ValueError(f"Sequence {value} does not belong to target {target}")
        return "sequence_cluster", value
    return "detected", None


def build_cells(
    run_id: str, raw_table_path: str, filtered_table_path: str,
    final_table_path: str, target_names: list[str], sequence_catalog: dict,
    sequence_qc_path: str | None,
) -> list[dict]:
    raw = pd.read_csv(raw_table_path, sep="\t", index_col="Barcode")
    filtered = pd.read_csv(filtered_table_path, sep="\t", index_col="Barcode")
    final = pd.read_csv(final_table_path, sep="\t", index_col="Barcode")
    if not final.index.equals(raw.index) or not final.index.equals(filtered.index):
        raise ValueError("Final, raw-count, and filtered-count cell tables have different cells or order")
    missing = [target for target in target_names if target not in raw or target not in filtered or target not in final]
    if missing:
        raise ValueError("Cell tables are missing targets: " + ", ".join(missing))
    sequence_qc = _load_sequence_qc(sequence_qc_path, target_names)
    cells = []
    for barcode, row in final.iterrows():
        asv_id = str(row["Assigned_core_asv"])
        if asv_id in {"", "NA", "None", "nan"}:
            raise ValueError(f"Final cell {barcode} has no valid ASV assignment")
        status = str(row["Status"])
        if status not in FINAL_ASV_STATUSES:
            raise ValueError(f"Final cell {barcode} has non-final ASV status {status}")
        targets = []
        for target in target_names:
            raw_count = _integer(raw.at[barcode, target])
            if raw_count == 0:
                continue
            filtered_count = _integer(filtered.at[barcode, target])
            assignment_type, sequence_id = _target_assignment(
                target, raw_count, filtered_count, row[target], sequence_catalog)
            target_record = {
                "target_name": target,
                "raw_read_count": raw_count,
                "filtered_read_count": filtered_count,
                "assignment_type": assignment_type,
            }
            if sequence_id is not None:
                target_record["sequence_cluster_id"] = sequence_id
            call_qc = sequence_qc.get(str(barcode), {}).get(target)
            if call_qc is not None:
                target_record["sequence_call_qc"] = call_qc
            targets.append(target_record)
        cells.append({
            "schema_version": SCHEMA_VERSION,
            "record_type": "cell",
            "run_id": run_id,
            "cell_barcode": str(barcode),
            "taxonomy": {
                "ranks": _parse_taxonomy(row["Predicted taxonomy"]),
                "confidence": _number_or_none(row["Confidence"]),
                "contamination": _number_or_none(row["Contamination"]),
            },
            "taxonomy_evidence": {
                "total_16s_reads": _integer(row["Total # of 16s reads"]),
                "technical_noise_reads": _integer(row["Technical noise count"]),
            },
            "asv": {
                "asv_id": asv_id,
                "status": status,
                "reads_used": _integer(row["Reads_used_for_ASV"]),
                "raw_unique_core_sequences": _integer(row["Raw unique_core_sequences"]),
                "dominant_raw_read_count": _integer(row["Dominant_raw_read_count"]),
                "coexisting_2bp_reads": _integer(row["Coexisting_2bp_reads"]),
                "unauthorized_secondary_reads": _integer(row["Unauthorized_secondary_reads"]),
                "final_cell_asv_reads": _integer(row["Final_cell_asv_reads"]),
                "max_internal_distance": _integer(row["Max_internal_distance"]),
            },
            "targets": targets,
        })
    return cells


def build_asvs(run_id: str, cells: list[dict], global_asv_path: str) -> list[dict]:
    sequences = _load_asv_sequences(global_asv_path)
    counts = Counter(cell["asv"]["asv_id"] for cell in cells)
    records = []
    for asv_id, count in counts.items():
        if asv_id not in sequences:
            raise ValueError(f"Cell references missing ASV {asv_id}")
        records.append({
            "schema_version": SCHEMA_VERSION,
            "record_type": "asv",
            "run_id": run_id,
            "asv_id": asv_id,
            "core_sequence": sequences[asv_id],
            "final_surviving_cell_count": count,
        })
    return records


def build_target_sequences(
    run_id: str, cells: list[dict], sequence_catalog: dict,
    phase_variation_path: str | None, reference_matches_path: str | None,
) -> list[dict]:
    counts = Counter(
        target["sequence_cluster_id"]
        for cell in cells for target in cell["targets"]
        if target["assignment_type"] == "sequence_cluster"
    )
    pv = _load_phase_variation(phase_variation_path)
    references = _load_reference_matches(reference_matches_path)
    records = []
    for sequence_id, count in counts.items():
        if sequence_id not in sequence_catalog:
            raise ValueError(f"Cell references missing target sequence {sequence_id}")
        source = sequence_catalog[sequence_id]
        record = {
            "schema_version": SCHEMA_VERSION,
            "record_type": "target_sequence",
            "run_id": run_id,
            "sequence_cluster_id": sequence_id,
            "parent_target": source["parent_target"],
            "core_sequence": {"r1": source["r1"], "r2": source["r2"]},
            "final_surviving_cell_count": count,
        }
        if sequence_id in pv:
            record["phase_variation"] = pv[sequence_id]
        if reference_matches_path is not None:
            record["reference_matches"] = references.get(sequence_id, [])
        records.append(record)
    return records


def _target_filtering_summary(path: str, target_names: list[str]) -> list[dict]:
    dataframe = pd.read_csv(path, sep="\t", index_col=0)
    records = []
    for target in target_names:
        row = dataframe.loc[target]
        records.append({
            "target_name": target,
            "original_positive_cells": _integer(row["Original"]),
            "filtered_out_cells": _integer(row["Filtered Out"]),
            "remaining_positive_cells": _integer(row["Remaining"]),
            "retention_percent": float(row["% Retention"]),
        })
    return records


def build_run_summary(
    run_id: str, runtime_config: dict, cells: list[dict], asvs: list[dict],
    read_qc_path: str, barcode_cluster_path: str, barcode_filter_path: str,
    asv_stats_path: str, manifest_path: str, packet_paths: tuple[str, str, str],
    target_filter_stats_path: str, phase_variation_path: str | None,
) -> dict:
    read_qc = _read_json(read_qc_path)
    barcode_cluster = _read_json(barcode_cluster_path)
    barcode_filter = _read_json(barcode_filter_path)
    asv_stats = _read_json(asv_stats_path)
    classifications = {
        "accepted_16s_reads": _count_lines(packet_paths[0]),
        "target_reads": _count_lines(packet_paths[1]),
        "unclassified_reads": _count_lines(packet_paths[2]),
    }
    if read_qc["passed_pairs"] != sum(classifications.values()):
        raise ValueError(
            "Read classification counts do not equal passed read pairs: "
            f"{read_qc['passed_pairs']} != {sum(classifications.values())}")
    stagger = _read_stagger_counts(manifest_path)
    if sum(stagger.values()) != classifications["accepted_16s_reads"]:
        raise ValueError("16S stagger counts do not equal accepted 16S packet count")
    final_cells = len(cells)
    if final_cells != asv_stats["final_cells"]:
        raise ValueError("Final cell count differs between ASV stage and canonical cells")
    summary = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "run_summary",
        "run_id": run_id,
        "status": "completed",
        "software": runtime_config["software"],
        "inputs": runtime_config["inputs"],
        "target_panel": runtime_config["target_panel"],
        "parameters": runtime_config["parameters"],
        "read_qc": {
            "raw_read_pairs": read_qc["total_pairs"],
            "passed_read_pairs": read_qc["passed_pairs"],
            "filtered_read_pairs": read_qc["filtered_pairs"],
            "rejected": {
                key: read_qc[key] for key in (
                    "read_length_or_quality_length", "mean_phred",
                    "barcode_length", "barcode_q25")
            },
        },
        "read_classification": classifications,
        "16s_r1_primer_starts": {
            str(start): stagger.get(start, 0)
            for start in runtime_config["parameters"]["primer_matching"]["valid_16s_r1_starts"]
        },
        "barcode_funnel": {
            "raw_unique_barcodes": barcode_cluster["raw_unique_barcodes"],
            "clustered_barcodes": barcode_cluster["clustered_barcodes"],
            "after_stage1_taxonomy_filter": barcode_filter["after_stage1_taxonomy_filter"],
            "after_minimum_cells_per_taxon": barcode_filter["after_minimum_cells_per_taxon"],
            "before_asv_filter": asv_stats["before_asv_filter"],
            "after_asv_status_filter": asv_stats["after_asv_status_filter"],
            "asv_taxonomy_conflicts_removed": asv_stats["asv_taxonomy_conflicts_removed"],
            "final_cells": final_cells,
        },
        "asv_summary": {"final_asv_count": len(asvs)},
        "target_filtering_summary": _target_filtering_summary(
            target_filter_stats_path,
            [item["target_name"] for item in runtime_config["target_panel"]]),
    }
    if phase_variation_path:
        pv = pd.read_csv(phase_variation_path, sep="\t")
        summary["phase_variation_summary"] = {
            "cell_calls": len(pv),
            "calls": {str(k): int(v) for k, v in pv["Call"].value_counts().items()},
        }
    return summary


def _require(condition: bool, message: str):
    if not condition:
        raise ValueError(message)


def _require_keys(record: dict, required: set[str], optional: set[str], label: str):
    keys = set(record)
    missing = required - keys
    unexpected = keys - required - optional
    _require(not missing, f"{label} missing required fields: {sorted(missing)}")
    _require(not unexpected, f"{label} has unexpected fields: {sorted(unexpected)}")


def _nonnegative_integer(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _number_or_null(value: Any) -> bool:
    return value is None or (
        isinstance(value, (int, float)) and not isinstance(value, bool)
        and math.isfinite(value)
    )



def _require_nonnegative_fields(record: dict, fields: set[str], label: str):
    _require(
        all(_nonnegative_integer(record[field]) for field in fields),
        f"{label} contains an invalid count",
    )


def _validate_run_summary_schema(run_summary: dict) -> set[str]:
    """Validate the exact v3.0.0 run-summary field contract."""
    _require_keys(run_summary["software"], {"name", "version", "git_commit"}, set(), "software")
    _require(isinstance(run_summary["software"]["name"], str), "Invalid software name")
    _require(isinstance(run_summary["software"]["version"], str), "Invalid software version")
    _require(
        run_summary["software"]["git_commit"] is None
        or isinstance(run_summary["software"]["git_commit"], str),
        "Invalid Git commit",
    )

    input_fields = {"r1", "r2", "primers", "taxonomy_database", "reference_fasta"}
    _require_keys(run_summary["inputs"], input_fields, set(), "inputs")
    _require(
        all(isinstance(run_summary["inputs"][field], str) for field in input_fields - {"reference_fasta"}),
        "Invalid required input path",
    )
    _require(
        run_summary["inputs"]["reference_fasta"] is None
        or isinstance(run_summary["inputs"]["reference_fasta"], str),
        "Invalid reference FASTA path",
    )

    target_names = []
    for target in run_summary["target_panel"]:
        _require_keys(target, {"target_name", "mode"}, set(), "target-panel record")
        _require(isinstance(target["target_name"], str) and target["target_name"], "Invalid target name")
        _require(isinstance(target["mode"], str) and target["mode"], "Invalid target mode")
        target_names.append(target["target_name"])
    _require(len(target_names) == len(set(target_names)), "Duplicate target-panel names")

    parameters = run_summary["parameters"]
    parameter_shapes = {
        "read_qc": {"minimum_read_length", "minimum_mean_phred", "barcode_length", "minimum_barcode_q25_bases"},
        "primer_matching": {"maximum_shift", "maximum_mismatches", "r2_primer_start", "valid_16s_r1_starts"},
        "barcode_clustering": {"barcode_length", "maximum_shift"},
        "cell_taxonomy_filtering": {"minimum_16s_reads", "maximum_contamination", "minimum_cells_per_taxon"},
        "taxonomy_mle": {"p_match", "p_none", "p_error", "alpha_prior", "beta_prior", "minimum_confidence", "minimum_noise_reads", "noise_cutoff_ratio"},
        "asv": {"r1_start", "r1_end", "r2_start", "r2_end", "maximum_distance", "maximum_shift", "minimum_reads", "mixed_ratio_threshold", "filter_corrupted_single_asv", "taxonomy_conflict_minimum_cells", "taxonomy_conflict_dominant_phylum_fraction"},
        "target_background_filtering": {"alpha"},
        "target_sequence_reconstruction": {"performed", "maximum_shift", "maximum_mismatches", "alpha", "r1_start", "r1_end", "r2_start", "r2_end", "include_all_targets", "mle"},
    }
    _require_keys(parameters, {"analysis_workers", *parameter_shapes}, set(), "parameters")
    _require(
        isinstance(parameters["analysis_workers"], int)
        and not isinstance(parameters["analysis_workers"], bool)
        and parameters["analysis_workers"] >= 1,
        "Invalid analysis worker count",
    )
    for section, fields in parameter_shapes.items():
        _require(isinstance(parameters[section], dict), f"parameters.{section} must be an object")
        _require_keys(parameters[section], fields, set(), f"parameters.{section}")
    starts = parameters["primer_matching"]["valid_16s_r1_starts"]
    _require(
        isinstance(starts, list) and all(_nonnegative_integer(value) for value in starts),
        "Invalid 16S primer-start configuration",
    )
    _require(isinstance(parameters["asv"]["filter_corrupted_single_asv"], bool), "Invalid ASV filter setting")
    reconstruction = parameters["target_sequence_reconstruction"]
    _require(isinstance(reconstruction["performed"], bool), "Invalid reconstruction setting")
    _require(isinstance(reconstruction["include_all_targets"], bool), "Invalid reconstruction target setting")
    _require_keys(reconstruction["mle"], parameter_shapes["taxonomy_mle"], set(), "reconstruction MLE")

    read_qc = run_summary["read_qc"]
    _require_keys(read_qc, {"raw_read_pairs", "passed_read_pairs", "filtered_read_pairs", "rejected"}, set(), "read QC")
    read_count_fields = {"raw_read_pairs", "passed_read_pairs", "filtered_read_pairs"}
    _require_nonnegative_fields(read_qc, read_count_fields, "read QC")
    rejected_fields = {"read_length_or_quality_length", "mean_phred", "barcode_length", "barcode_q25"}
    _require_keys(read_qc["rejected"], rejected_fields, set(), "read-QC rejection counts")
    _require_nonnegative_fields(read_qc["rejected"], rejected_fields, "read-QC rejection counts")
    _require(read_qc["raw_read_pairs"] == read_qc["passed_read_pairs"] + read_qc["filtered_read_pairs"], "Raw-read QC counts do not balance")
    _require(read_qc["filtered_read_pairs"] == sum(read_qc["rejected"].values()), "Read-QC rejection counts do not balance")

    classification_fields = {"accepted_16s_reads", "target_reads", "unclassified_reads"}
    _require_keys(run_summary["read_classification"], classification_fields, set(), "read classification")
    _require_nonnegative_fields(run_summary["read_classification"], classification_fields, "read classification")
    expected_starts = {str(value) for value in starts}
    _require_keys(run_summary["16s_r1_primer_starts"], expected_starts, set(), "16S primer starts")
    _require_nonnegative_fields(run_summary["16s_r1_primer_starts"], expected_starts, "16S primer starts")

    funnel_fields = {"raw_unique_barcodes", "clustered_barcodes", "after_stage1_taxonomy_filter", "after_minimum_cells_per_taxon", "before_asv_filter", "after_asv_status_filter", "asv_taxonomy_conflicts_removed", "final_cells"}
    _require_keys(run_summary["barcode_funnel"], funnel_fields, set(), "barcode funnel")
    _require_nonnegative_fields(run_summary["barcode_funnel"], funnel_fields, "barcode funnel")
    _require_keys(run_summary["asv_summary"], {"final_asv_count"}, set(), "ASV summary")
    _require_nonnegative_fields(run_summary["asv_summary"], {"final_asv_count"}, "ASV summary")

    target_filter_fields = {"target_name", "original_positive_cells", "filtered_out_cells", "remaining_positive_cells", "retention_percent"}
    observed_targets = []
    for target in run_summary["target_filtering_summary"]:
        _require_keys(target, target_filter_fields, set(), "target-filtering summary")
        observed_targets.append(target["target_name"])
        _require_nonnegative_fields(target, {"original_positive_cells", "filtered_out_cells", "remaining_positive_cells"}, f"target-filtering summary for {target['target_name']}")
        _require(_number_or_null(target["retention_percent"]), "Invalid target retention percent")
    _require(observed_targets == target_names, "Target-filtering summary does not match target-panel order")

    if "phase_variation_summary" in run_summary:
        phase = run_summary["phase_variation_summary"]
        _require_keys(phase, {"cell_calls", "calls"}, set(), "phase-variation summary")
        _require(_nonnegative_integer(phase["cell_calls"]), "Invalid phase-variation call count")
        _require(isinstance(phase["calls"], dict), "Invalid phase-variation counts")
        _require(all(isinstance(key, str) and _nonnegative_integer(value) for key, value in phase["calls"].items()), "Invalid phase-variation counts")
    return set(target_names)

def validate_canonical(
    run_id: str, cells: list[dict], asvs: list[dict],
    target_sequences: list[dict] | None, run_summary: dict,
):
    import uuid
    parsed_run_id = uuid.UUID(run_id)
    _require(parsed_run_id.version == 4 and str(parsed_run_id) == run_id, "run_id is not a canonical UUIDv4")
    run_required = {
        "schema_version", "record_type", "run_id", "status", "software",
        "inputs", "target_panel", "parameters", "read_qc",
        "read_classification", "16s_r1_primer_starts", "barcode_funnel",
        "asv_summary", "target_filtering_summary",
    }
    _require_keys(run_summary, run_required, {"phase_variation_summary"}, "run_summary")
    _require(run_summary["schema_version"] == SCHEMA_VERSION, "Invalid run schema version")
    _require(run_summary["record_type"] == "run_summary", "Invalid run record type")
    _require(run_summary["run_id"] == run_id, "Invalid run ID")
    _require(run_summary["status"] == "completed", "Canonical run is not completed")
    for section in ("software", "inputs", "parameters", "read_qc", "read_classification", "barcode_funnel", "asv_summary"):
        _require(isinstance(run_summary[section], dict), f"run_summary.{section} must be an object")
    _require(isinstance(run_summary["target_panel"], list), "target_panel must be an array")
    _require(isinstance(run_summary["target_filtering_summary"], list), "target_filtering_summary must be an array")
    target_panel = _validate_run_summary_schema(run_summary)

    barcodes = set()
    asv_references = Counter()
    sequence_references = Counter()
    cell_required = {
        "schema_version", "record_type", "run_id", "cell_barcode",
        "taxonomy", "taxonomy_evidence", "asv", "targets",
    }
    asv_cell_fields = {
        "asv_id", "status", "reads_used", "raw_unique_core_sequences",
        "dominant_raw_read_count", "coexisting_2bp_reads",
        "unauthorized_secondary_reads", "final_cell_asv_reads",
        "max_internal_distance",
    }
    for cell in cells:
        _require_keys(cell, cell_required, set(), "cell")
        _require(cell["schema_version"] == SCHEMA_VERSION, "Invalid cell schema version")
        _require(cell["record_type"] == "cell", "Invalid cell record type")
        _require(cell["run_id"] == run_id, "Cell run ID mismatch")
        barcode = cell["cell_barcode"]
        _require(isinstance(barcode, str) and barcode and barcode not in barcodes, f"Duplicate/invalid barcode {barcode}")
        barcodes.add(barcode)
        _require_keys(cell["taxonomy"], {"ranks", "confidence", "contamination"}, set(), f"taxonomy for {barcode}")
        ranks = cell["taxonomy"]["ranks"]
        _require(set(ranks) == set(TAXONOMY_RANKS.values()), f"Invalid taxonomy ranks for {barcode}")
        _require(all(value is None or isinstance(value, str) for value in ranks.values()), f"Invalid taxonomy value for {barcode}")
        _require(_number_or_null(cell["taxonomy"]["confidence"]), f"Invalid taxonomy confidence for {barcode}")
        _require(_number_or_null(cell["taxonomy"]["contamination"]), f"Invalid taxonomy contamination for {barcode}")
        _require_keys(cell["taxonomy_evidence"], {"total_16s_reads", "technical_noise_reads"}, set(), f"taxonomy evidence for {barcode}")
        _require(all(_nonnegative_integer(value) for value in cell["taxonomy_evidence"].values()), f"Invalid taxonomy evidence for {barcode}")
        asv = cell["asv"]
        _require_keys(asv, asv_cell_fields, set(), f"ASV assignment for {barcode}")
        _require(isinstance(asv["asv_id"], str) and asv["asv_id"], f"Invalid ASV ID for {barcode}")
        _require(asv["status"] in FINAL_ASV_STATUSES, f"Invalid final ASV status for {barcode}")
        _require(all(_nonnegative_integer(asv[field]) for field in asv_cell_fields - {"asv_id", "status"}), f"Invalid ASV counts for {barcode}")
        asv_references[asv["asv_id"]] += 1
        _require(isinstance(cell["targets"], list), f"Targets must be an array for {barcode}")
        seen_targets = set()
        for target in cell["targets"]:
            _require_keys(
                target,
                {"target_name", "raw_read_count", "filtered_read_count", "assignment_type"},
                {"sequence_cluster_id", "sequence_call_qc"},
                f"target for {barcode}")
            name = target["target_name"]
            _require(name in target_panel and name not in seen_targets, f"Invalid/duplicate target {name} for {barcode}")
            seen_targets.add(name)
            _require(_nonnegative_integer(target["raw_read_count"]) and target["raw_read_count"] > 0, f"Sparse target {name} has no raw reads")
            _require(_nonnegative_integer(target["filtered_read_count"]), f"Invalid filtered count for {name}")
            assignment = target["assignment_type"]
            _require(assignment in ASSIGNMENT_TYPES, f"Invalid target assignment {assignment}")
            if assignment == "filtered_out":
                _require(target["filtered_read_count"] == 0, "filtered_out target has nonzero count")
            else:
                _require(target["filtered_read_count"] > 0, "positive target has zero count")
            if assignment == "sequence_cluster":
                sequence_id = target.get("sequence_cluster_id")
                _require(isinstance(sequence_id, str) and sequence_id, "sequence_cluster assignment has no ID")
                sequence_references[sequence_id] += 1
            else:
                _require("sequence_cluster_id" not in target, "Non-sequence target has a sequence ID")
            if "sequence_call_qc" in target:
                qc = target["sequence_call_qc"]
                _require_keys(qc, {"confidence", "contamination", "total_reads", "technical_noise_reads"}, set(), f"sequence call QC for {barcode}/{name}")
                _require(_number_or_null(qc["confidence"]) and _number_or_null(qc["contamination"]), f"Invalid sequence call QC probabilities for {barcode}/{name}")
                _require(_nonnegative_integer(qc["total_reads"]) and _nonnegative_integer(qc["technical_noise_reads"]), f"Invalid sequence call QC counts for {barcode}/{name}")

    asv_ids = set()
    asv_required = {
        "schema_version", "record_type", "run_id", "asv_id",
        "core_sequence", "final_surviving_cell_count",
    }
    for asv in asvs:
        _require_keys(asv, asv_required, set(), "ASV")
        _require(asv["schema_version"] == SCHEMA_VERSION and asv["record_type"] == "asv", "Invalid ASV record")
        _require(asv["run_id"] == run_id, "ASV run ID mismatch")
        asv_id = asv["asv_id"]
        _require(isinstance(asv_id, str) and asv_id.startswith("ASV_") and asv_id[4:].isdigit() and int(asv_id[4:]) > 0, f"Invalid ASV ID {asv_id}")
        _require(asv_id not in asv_ids, f"Duplicate ASV ID {asv_id}")
        asv_ids.add(asv_id)
        _require_keys(asv["core_sequence"], {"r1", "r2"}, set(), f"ASV sequence {asv_id}")
        _require(all(isinstance(value, str) for value in asv["core_sequence"].values()), f"Invalid ASV sequence {asv_id}")
        _require(asv["final_surviving_cell_count"] == asv_references[asv_id] > 0, f"ASV count mismatch for {asv_id}")
    _require(asv_ids == set(asv_references), "ASV foreign keys do not resolve exactly")

    sequence_ids = set()
    sequence_required = {
        "schema_version", "record_type", "run_id", "sequence_cluster_id",
        "parent_target", "core_sequence", "final_surviving_cell_count",
    }
    for sequence in target_sequences or []:
        _require_keys(sequence, sequence_required, {"phase_variation", "reference_matches"}, "target sequence")
        _require(sequence["schema_version"] == SCHEMA_VERSION and sequence["record_type"] == "target_sequence", "Invalid target-sequence record")
        _require(sequence["run_id"] == run_id, "Target-sequence run ID mismatch")
        sequence_id = sequence["sequence_cluster_id"]
        parent = sequence["parent_target"]
        prefix, separator, suffix = sequence_id.rpartition("_seq_")
        _require(separator and prefix == parent and suffix.isdigit() and int(suffix) > 0, f"Invalid target-sequence ID {sequence_id}")
        _require(sequence_id not in sequence_ids, f"Duplicate target-sequence ID {sequence_id}")
        sequence_ids.add(sequence_id)
        _require(parent in target_panel, f"Unknown parent target for {sequence_id}")
        _require_keys(sequence["core_sequence"], {"r1", "r2"}, set(), f"target sequence {sequence_id}")
        _require(
            all(isinstance(value, str) for value in sequence["core_sequence"].values()),
            f"Invalid target sequence {sequence_id}",
        )
        _require(sequence["final_surviving_cell_count"] == sequence_references[sequence_id] > 0, f"Target-sequence count mismatch for {sequence_id}")
        if "phase_variation" in sequence:
            phase = sequence["phase_variation"]
            _require_keys(
                phase,
                {"mode", "call", "motif", "repeat_count", "start", "end", "direct_similarity", "evidence"},
                set(), f"phase variation for {sequence_id}",
            )
            _require(all(isinstance(phase[key], str) for key in ("mode", "call", "evidence")), f"Invalid phase-variation text for {sequence_id}")
            _require(phase["motif"] is None or isinstance(phase["motif"], str), f"Invalid phase-variation motif for {sequence_id}")
            _require(all(phase[key] is None or isinstance(phase[key], int) for key in ("repeat_count", "start", "end")), f"Invalid phase-variation coordinates for {sequence_id}")
            _require(_number_or_null(phase["direct_similarity"]), f"Invalid phase-variation similarity for {sequence_id}")
        if "reference_matches" in sequence:
            _require(isinstance(sequence["reference_matches"], list), f"Invalid reference matches for {sequence_id}")
            reference_fields = {
                "subject", "r1_start", "r1_end", "r1_percent_identity",
                "r2_start", "r2_end", "r2_percent_identity",
            }
            for match in sequence["reference_matches"]:
                _require_keys(match, reference_fields, set(), f"reference match for {sequence_id}")
                _require(isinstance(match["subject"], str), f"Invalid reference subject for {sequence_id}")
                _require(all(isinstance(match[key], int) for key in ("r1_start", "r1_end", "r2_start", "r2_end")), f"Invalid reference coordinates for {sequence_id}")
                _require(all(_number_or_null(match[key]) for key in ("r1_percent_identity", "r2_percent_identity")), f"Invalid reference identity for {sequence_id}")
    _require(sequence_ids == set(sequence_references), "Target-sequence foreign keys do not resolve exactly")
    _require(run_summary["barcode_funnel"]["final_cells"] == len(cells), "Run final-cell count mismatch")
    _require(run_summary["asv_summary"]["final_asv_count"] == len(asvs), "Run ASV count mismatch")
    classification = run_summary["read_classification"]
    _require(
        run_summary["read_qc"]["passed_read_pairs"] == sum(classification.values()),
        "Passed-read and classification counts differ")
    _require(
        sum(run_summary["16s_r1_primer_starts"].values()) == classification["accepted_16s_reads"],
        "16S stagger and accepted-read counts differ")

def _write_jsonl(path: Path, records: list[dict]):
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")


def _write_json(path: Path, record: dict):
    with path.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, ensure_ascii=False, allow_nan=False, indent=2)
        handle.write("\n")


def write_canonical_results(args) -> dict:
    runtime_config = _read_json(args.runtime_config)
    run_id = runtime_config["run_id"]
    target_names = [item["target_name"] for item in runtime_config["target_panel"]]
    sequence_catalog = _load_sequence_catalog(args.sequence_list)
    cells = build_cells(
        run_id, args.raw_cell_table, args.filtered_counts, args.final_cell_table,
        target_names, sequence_catalog, args.sequence_qc)
    asvs = build_asvs(run_id, cells, args.global_asv)
    target_sequences = None
    if runtime_config["parameters"]["target_sequence_reconstruction"]["performed"]:
        target_sequences = build_target_sequences(
            run_id, cells, sequence_catalog, args.phase_variation,
            args.reference_matches)
    run_summary = build_run_summary(
        run_id, runtime_config, cells, asvs,
        args.read_qc_stats, args.barcode_cluster_stats,
        args.barcode_filter_stats, args.asv_stats, args.manifest,
        (args.packet_16s, args.packet_target, args.packet_unclassified),
        args.target_filter_stats, args.phase_variation)
    validate_canonical(run_id, cells, asvs, target_sequences, run_summary)

    staging = Path(args.staging_dir)
    staging.mkdir(parents=True, exist_ok=False)
    cell_tmp = staging / "cells.jsonl"
    asv_tmp = staging / "asvs.jsonl"
    run_tmp = staging / "run_summary.json"
    _write_jsonl(cell_tmp, cells)
    _write_jsonl(asv_tmp, asvs)
    if target_sequences is not None:
        _write_jsonl(staging / "target_sequences.jsonl", target_sequences)
    _write_json(run_tmp, run_summary)

    # Validate the serialized representation, not only the in-memory objects.
    reloaded_cells = [json.loads(line) for line in cell_tmp.read_text().splitlines()]
    reloaded_asvs = [json.loads(line) for line in asv_tmp.read_text().splitlines()]
    target_path = staging / "target_sequences.jsonl"
    reloaded_targets = (
        [json.loads(line) for line in target_path.read_text().splitlines()]
        if target_path.exists() else None)
    reloaded_summary = json.loads(run_tmp.read_text())
    validate_canonical(run_id, reloaded_cells, reloaded_asvs, reloaded_targets, reloaded_summary)

    for path in staging.iterdir():
        os.replace(path, Path(args.output_dir) / path.name)
    staging.rmdir()
    return {
        "cells": len(cells), "asvs": len(asvs),
        "target_sequences": len(target_sequences or []),
    }


def main():
    parser = argparse.ArgumentParser(description="Build v3.0.0 canonical results")
    for name in (
        "runtime_config", "raw_cell_table", "filtered_counts", "final_cell_table",
        "global_asv", "read_qc_stats", "barcode_cluster_stats",
        "barcode_filter_stats", "asv_stats", "manifest", "packet_16s",
        "packet_target", "packet_unclassified", "target_filter_stats",
        "staging_dir", "output_dir",
    ):
        parser.add_argument("--" + name.replace("_", "-"), dest=name, required=True)
    parser.add_argument("--sequence-list")
    parser.add_argument("--sequence-qc")
    parser.add_argument("--phase-variation")
    parser.add_argument("--reference-matches")
    args = parser.parse_args()
    counts = write_canonical_results(args)
    print(
        "Canonical v3.0.0 validated: "
        f"cells={counts['cells']} asvs={counts['asvs']} "
        f"target_sequences={counts['target_sequences']}", flush=True)


if __name__ == "__main__":
    main()
