"""Export the final DoTA-AMR cell-level result as JSON Lines."""

import argparse
import json
import os
from pathlib import Path

import pandas as pd

from helper_functions import get_arg_names


# 2026-08-10: Define stable public names for the taxonomy ranks in the TSV lineage.
# Reason: the final JSONL should expose structured taxa instead of requiring users to parse display text.
TAXONOMY_RANKS = {
    "R1": "root",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}


def _native_value(value):
    """Convert pandas/numpy values into JSON-compatible Python values."""
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


def _parse_taxonomy(lineage):
    """Return a rank-keyed taxonomy object from the displayed lineage."""
    if pd.isna(lineage):
        return {}

    ranks = {}
    for item in str(lineage).split(" | "):
        rank_code, separator, taxon = item.partition(" - ")
        if separator and rank_code in TAXONOMY_RANKS:
            ranks[TAXONOMY_RANKS[rank_code]] = None if taxon == "None" else taxon
    return ranks


def _is_detected(value):
    """Return whether an AMR assignment represents a detected target."""
    if pd.isna(value):
        return False
    return str(value).strip() not in {"", "0", "0.0", "None", "nan"}


def export_results(input_tsv, primers_file, output_jsonl):
    """Write one complete single-cell DoTA-AMR result per JSONL line."""
    dataframe = pd.read_csv(input_tsv, sep="\t", index_col="Barcode")
    gene_names = get_arg_names(primers_file)
    missing_genes = [gene for gene in gene_names if gene not in dataframe.columns]
    if missing_genes:
        raise ValueError(
            "Final result table is missing primer-defined AMR columns: "
            + ", ".join(missing_genes)
        )

    output_path = Path(output_jsonl)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 2026-08-10: Emit one nested record per cell and omit absent AMR targets.
    # Reason: JSONL is the primary machine-readable product and remains streamable for large experiments.
    with output_path.open("w", encoding="utf-8") as output_handle:
        for barcode, row in dataframe.iterrows():
            amr_assignments = [
                {"gene": gene, "assignment": _native_value(row[gene])}
                for gene in gene_names
                if _is_detected(row[gene])
            ]
            lineage = _native_value(row.get("Predicted taxonomy"))
            record = {
                "schema_version": "1.0",
                "cell_barcode": str(barcode),
                "taxonomy": {
                    "lineage": lineage,
                    "ranks": _parse_taxonomy(lineage),
                    "confidence": _native_value(row.get("Confidence")),
                    "contamination": _native_value(row.get("Contamination")),
                },
                "asv": {
                    "assignment": _native_value(row.get("Assigned_core_asv")),
                    "status": _native_value(row.get("Status")),
                    "reads": _native_value(row.get("Final_cell_asv_reads")),
                    "max_internal_distance": _native_value(row.get("Max_internal_distance")),
                },
                "read_qc": {
                    "total_16s_reads": _native_value(row.get("Total # of 16s reads")),
                    "technical_noise_reads": _native_value(row.get("Technical noise count")),
                    "reads_used_for_asv": _native_value(row.get("Reads_used_for_ASV")),
                },
                "amr": amr_assignments,
            }
            output_handle.write(json.dumps(record, ensure_ascii=False) + "\n")

    return len(dataframe)


def main():
    parser = argparse.ArgumentParser(
        description="Export the final DoTA-AMR single-cell result as JSONL."
    )
    parser.add_argument("--input_tsv", required=True)
    parser.add_argument("--primers_file", required=True)
    parser.add_argument("--output_jsonl", default="dota_amr_results.jsonl")
    args = parser.parse_args()

    for input_path in (args.input_tsv, args.primers_file):
        if not os.path.exists(input_path):
            parser.error(f"input file not found: {input_path}")

    cell_count = export_results(args.input_tsv, args.primers_file, args.output_jsonl)
    print(f"Wrote {cell_count} single-cell records to {args.output_jsonl}")


if __name__ == "__main__":
    main()
