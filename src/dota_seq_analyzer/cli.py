#!/usr/bin/env python3
"""Run the complete DoTA-Seq Analyzer workflow from one public command."""

# 2026-08-11: Rename the public entry point to DoTA-Seq Analyzer.
# Reason: users should run R1/R2 data with one command instead of invoking internal modules manually.

import argparse
import json
from importlib.metadata import PackageNotFoundError, version
import os
import shutil
import subprocess
import sys
import tomllib
import uuid
from pathlib import Path

from . import algorithm_config as config
from .helper_functions import get_target_modes
from .validate_inputs import check_reference_fasta


def _run_step(name: str, command: list[str], output_dir: Path) -> None:
    """Run one pipeline stage and stop immediately if it fails."""
    print(f"\n[DoTA-Seq Analyzer] {name}", flush=True)
    subprocess.run(command, cwd=output_dir, check=True)


def _find_project_database(relative_path: str) -> Path:
    """Find a database extracted beside either the current clone or source package."""
    source_root = Path(__file__).resolve().parents[2]
    for project_root in (Path.cwd(), source_root):
        candidate = project_root / relative_path
        if candidate.exists():
            return candidate.resolve()
    return (source_root / relative_path).resolve()


def _software_version(source_root: Path) -> str:
    """Return installed metadata, or the current source-tree project version."""
    try:
        return version("dota-seq-analyzer")
    except PackageNotFoundError:
        with (source_root / "pyproject.toml").open("rb") as handle:
            return tomllib.load(handle)["project"]["version"]


def _git_commit(source_root: Path):
    """Return the source checkout commit when Git metadata is available."""
    result = subprocess.run(
        ["git", "-C", str(source_root), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> None:
    source_dir = Path(__file__).resolve().parent
    source_root = source_dir.parents[1]
    parser = argparse.ArgumentParser(
        prog="dota-seq-analyzer",
        description="Profile targeted genes and phase variation in single-cell DoTA-Seq data.",
    )
    # 2026-09-10: Report the installed software version from the existing metadata source.
    # Reason: versioned releases should be identifiable without starting an analysis.
    parser.add_argument(
        "--version", action="version",
        version=f"%(prog)s {_software_version(source_root)}",
    )
    parser.add_argument("-1", "--r1", required=True, help="R1 FASTQ file")
    parser.add_argument("-2", "--r2", required=True, help="R2 FASTQ file")
    parser.add_argument("-p", "--primers", required=True, help="DoTA-Seq primer CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-@", "--threads", dest="analysis_workers", type=int, default=1, metavar="INT", help="Number of parallel workers/threads used by DoTA-seq analysis")
    parser.add_argument(
        "--keep-tmp", action="store_true",
        help="Preserve computational intermediate files after a successful run",
    )
    parser.add_argument("--taxonomy-db", help="Extracted Kraken2 taxonomy database directory")
    # 2026-08-28: Expose the Stage 2 threshold in the public CLI.
    # Reason: users can control low-count taxon filtering without a separate skip flag.
    parser.add_argument("--min-cells-per-taxon", type=int, default=10, help="Minimum cells required per taxon; 0 disables Stage 2 filtering (default: 10)")
    # 2026-08-10: Run reference annotation only when a FASTA is explicitly supplied.
    # Reason: DoTA-Seq Analyzer supports arbitrary targets that do not belong in the bundled AMR database.
    parser.add_argument("-r", "--reference", help="Optional reference FASTA for BLAST annotation")
    args = parser.parse_args()

    r1 = Path(args.r1).expanduser().resolve()
    r2 = Path(args.r2).expanduser().resolve()
    primers = Path(args.primers).expanduser().resolve()
    requested_output_dir = Path(args.output).expanduser().resolve()
    run_id = str(uuid.uuid4())
    if requested_output_dir.exists():
        parser.error(f"output directory already exists: {requested_output_dir}")
    output_dir = requested_output_dir.parent / f".{requested_output_dir.name}.incomplete.{run_id}"

    for label, path in (("R1", r1), ("R2", r2), ("primer CSV", primers)):
        if not path.is_file():
            parser.error(f"{label} file not found: {path}")
    if args.analysis_workers < 1:
        parser.error("--threads must be at least 1")
    if args.min_cells_per_taxon < 0:
        parser.error("--min-cells-per-taxon must be non-negative")

    taxonomy_db = (
        Path(args.taxonomy_db).expanduser().resolve()
        if args.taxonomy_db
        else _find_project_database("database/mnt/workspace2/jamie/ref/k2__gg2")
    )
    reference = Path(args.reference).expanduser().resolve() if args.reference else None
    if not taxonomy_db.is_dir():
        parser.error(
            "taxonomy database not found; extract database/dota-seq-analyzer-taxonomy-db.tar.gz "
            "or provide --taxonomy-db"
        )
    if reference is not None and not reference.is_file():
        parser.error(f"reference FASTA not found: {reference}")
    if reference is not None:
        reference_validation = check_reference_fasta(str(reference))
        if reference_validation != "Valid":
            parser.error(reference_validation)

    target_modes = get_target_modes(str(primers))
    # 2026-08-11: Trigger phase-variation analysis only for SSR primer targets.
    # Reason: SSR is the only supported phase-variation mode.
    pv_requested = any(mode == "ssr" for mode in target_modes.values())
    reconstruction_requested = pv_requested or reference is not None
    runtime_config = {
        "run_id": run_id,
        "software": {
            "name": "dota-seq-analyzer",
            "version": _software_version(source_root),
            "git_commit": _git_commit(source_root),
        },
        "inputs": {
            "r1": str(r1), "r2": str(r2), "primers": str(primers),
            "taxonomy_database": str(taxonomy_db),
            "reference_fasta": str(reference) if reference else None,
        },
        "target_panel": [
            {"target_name": target, "mode": mode or "detect"}
            for target, mode in target_modes.items()
        ],
        "parameters": {
            "analysis_workers": args.analysis_workers,
            "read_qc": {
                "minimum_read_length": config.MIN_READ_LENGTH,
                "minimum_mean_phred": config.MIN_MEAN_PHRED,
                "barcode_length": config.BARCODE_LENGTH,
                "minimum_barcode_q25_bases": config.MIN_BARCODE_Q25,
            },
            "primer_matching": {
                "maximum_shift": config.PRIMER_MAX_SHIFT,
                "maximum_mismatches": config.PRIMER_MAX_MISMATCHES,
                "r2_primer_start": config.R2_PRIMER_START,
                "valid_16s_r1_starts": list(config.VALID_16S_R1_STARTS),
            },
            "barcode_clustering": {
                "barcode_length": config.BARCODE_LENGTH,
                "maximum_shift": config.BARCODE_MAX_SHIFT,
            },
            "cell_taxonomy_filtering": {
                "minimum_16s_reads": config.MIN_16S_READS,
                "maximum_contamination": config.MAX_TAXONOMY_CONTAMINATION,
                "minimum_cells_per_taxon": args.min_cells_per_taxon,
            },
            "taxonomy_mle": {
                "p_match": config.MLE_P_MATCH, "p_none": config.MLE_P_NONE,
                "p_error": config.MLE_P_ERROR,
                "alpha_prior": config.MLE_ALPHA_PRIOR,
                "beta_prior": config.MLE_BETA_PRIOR,
                "minimum_confidence": config.MLE_MIN_CONFIDENCE,
                "minimum_noise_reads": config.MLE_MIN_NOISE_READS,
                "noise_cutoff_ratio": config.MLE_NOISE_CUTOFF_RATIO,
            },
            "asv": {
                "r1_start": config.ASV_R1_START, "r1_end": config.ASV_R1_END,
                "r2_start": config.ASV_R2_START, "r2_end": config.ASV_R2_END,
                "maximum_distance": config.ASV_MAX_DISTANCE,
                "maximum_shift": config.ASV_MAX_SHIFT,
                "minimum_reads": config.ASV_MIN_READS,
                "mixed_ratio_threshold": config.ASV_MIXED_RATIO_THRESHOLD,
                "filter_corrupted_single_asv": False,
                "taxonomy_conflict_minimum_cells": config.ASV_TAXONOMY_CONFLICT_MIN_CELLS,
                "taxonomy_conflict_dominant_phylum_fraction": config.ASV_TAXONOMY_CONFLICT_DOMINANCE,
            },
            "target_background_filtering": {"alpha": config.TARGET_BACKGROUND_ALPHA},
            "target_sequence_reconstruction": {
                "performed": reconstruction_requested,
                "maximum_shift": config.TARGET_SEQUENCE_MAX_SHIFT,
                "maximum_mismatches": config.TARGET_SEQUENCE_MAX_MISMATCHES,
                "alpha": config.TARGET_SEQUENCE_ALPHA,
                "r1_start": config.TARGET_R1_START, "r1_end": config.TARGET_R1_END,
                "r2_start": config.TARGET_R2_START, "r2_end": config.TARGET_R2_END,
                "include_all_targets": reference is not None,
                "mle": {
                    "p_match": config.MLE_P_MATCH, "p_none": config.MLE_P_NONE,
                    "p_error": config.MLE_P_ERROR,
                    "alpha_prior": config.MLE_ALPHA_PRIOR,
                    "beta_prior": config.MLE_BETA_PRIOR,
                    "minimum_confidence": config.MLE_MIN_CONFIDENCE,
                    "minimum_noise_reads": config.MLE_MIN_NOISE_READS,
                    "noise_cutoff_ratio": config.MLE_NOISE_CUTOFF_RATIO,
                },
            },
        },
    }

    for directory in (
        output_dir,
        output_dir / "tmp",
        output_dir / "figures",
    ):
        directory.mkdir(parents=True, exist_ok=True)

    with (output_dir / "tmp/runtime_config.json").open("w", encoding="utf-8") as handle:
        json.dump(runtime_config, handle, indent=2)
        handle.write("\n")

    python = sys.executable
    script = lambda name: str(source_dir / name)
    filtered_r1 = output_dir / "tmp/filtered_R1.fastq"
    filtered_r2 = output_dir / "tmp/filtered_R2.fastq"

    _run_step(
        "Filter reads",
        [
            python,
            script("filter_reads.py"),
            "--r1",
            str(r1),
            "--r2",
            str(r2),
            "--output-r1",
            str(filtered_r1),
            "--output-r2",
            str(filtered_r2),
            "--stats-json",
            "tmp/read_qc_stats.json",
        ],
        output_dir,
    )
    _run_step(
        "Validate inputs",
        [
            python,
            script("validate_inputs.py"),
            "--r1_fastq",
            str(filtered_r1),
            "--r2_fastq",
            str(filtered_r2),
            "--primers_file",
            str(primers),
        ],
        output_dir,
    )
    _run_step(
        "Extract targeted taxonomy reads",
        [
            python,
            script("extract_16s_reads.py"),
            "--r1_fastq",
            str(filtered_r1),
            "--r2_fastq",
            str(filtered_r2),
            "--primers_filename",
            str(primers),
            "--threads",
            str(args.analysis_workers),
            "--r1_16s_manifest",
            "tmp/16s_r1_primer_starts.tsv",
            "--max_shift_primer", str(config.PRIMER_MAX_SHIFT),
            "--max_mm_primer", str(config.PRIMER_MAX_MISMATCHES),
            "--primer_start_num", str(config.R2_PRIMER_START),
        ],
        output_dir,
    )
    _run_step(
        "Classify taxonomy reads",
        [
            "kraken2",
            "--db",
            str(taxonomy_db),
            "--threads",
            str(args.analysis_workers),
            "--paired",
            "tmp/kraken_R1.fastq",
            "tmp/kraken_R2.fastq",
            "--output",
            "tmp/kraken.output",
            "--report",
            "tmp/kraken.report",
        ],
        output_dir,
    )
    _run_step(
        "Assign read packets",
        [
            python,
            script("create_ID_packets.py"),
            "--r1_fastq",
            str(filtered_r1),
            "--r2_fastq",
            str(filtered_r2),
            "--primers_filename",
            str(primers),
            "--r1_16s_manifest",
            "tmp/16s_r1_primer_starts.tsv",
            "--kraken_output",
            "tmp/kraken.output",
            "--kraken_report",
            "tmp/kraken.report",
            "--max_shift_primer", str(config.PRIMER_MAX_SHIFT),
            "--max_mm_primer", str(config.PRIMER_MAX_MISMATCHES),
            "--primer_start_num", str(config.R2_PRIMER_START),
            "--barcode_len", str(config.BARCODE_LENGTH),
            "--threads",
            str(args.analysis_workers),
        ],
        output_dir,
    )
    _run_step(
        "Match reads to cell barcodes",
        [
            python,
            script("match_barcodes_to_IDs_revised.py"),
            "--r2_16s_fastq",
            "tmp/only_16s_R2.fastq",
            "--arg_r2_fastq",
            "tmp/arg_R2.fastq",
            "--unclassified_r2_fastq",
            "tmp/unclassified_R2.fastq",
            "--max_shift_barcode", str(config.BARCODE_MAX_SHIFT),
            "--barcode_len", str(config.BARCODE_LENGTH),
            "--stats_json", "tmp/barcode_cluster_stats.json",
        ],
        output_dir,
    )
    _run_step(
        "Build cell summaries",
        [
            python,
            script("barcode_summary.py"),
            "--b_with_ids_filename",
            "tmp/b_with_ids.txt",
            "--_16s_packet_filename",
            "tmp/packets_16s",
            "--arg_packet_filename",
            "tmp/packets_arg",
            "--primers_filename",
            str(primers),
            "--threads",
            str(args.analysis_workers),
            "--min_cells_per_taxon",
            str(args.min_cells_per_taxon),
            "--min_16s_reads", str(config.MIN_16S_READS),
            "--max_contam", str(config.MAX_TAXONOMY_CONTAMINATION),
            "--p_match", str(config.MLE_P_MATCH),
            "--p_none", str(config.MLE_P_NONE),
            "--p_error", str(config.MLE_P_ERROR),
            "--alpha_prior", str(config.MLE_ALPHA_PRIOR),
            "--beta_prior", str(config.MLE_BETA_PRIOR),
            "--min_confidence", str(config.MLE_MIN_CONFIDENCE),
            "--min_noise_reads", str(config.MLE_MIN_NOISE_READS),
            "--noise_cutoff_ratio", str(config.MLE_NOISE_CUTOFF_RATIO),
            "--stats_json", "tmp/barcode_filter_stats.json",
        ],
        output_dir,
    )
    _run_step(
        "Type and filter cell sequence variants",
        [
            python,
            script("asv_typing_revised.py"),
            "--barcode_summary_tsv",
            "tmp/barcode_summary.tsv",
            "--b_with_ids",
            "tmp/b_with_ids.txt",
            "--r1_16s_fastq",
            "tmp/only_16s_R1.fastq",
            "--r2_16s_fastq",
            "tmp/only_16s_R2.fastq",
            "--r1_16s_manifest",
            "tmp/16s_r1_primer_starts.tsv",
            "--primers_file",
            str(primers),
            "--filter_corrupted", "false",
            "--stats_json", "tmp/asv_stats.json",
        ],
        output_dir,
    )
    _run_step(
        "Filter target background",
        [
            python,
            script("filter_args.py"),
            "--input_arg_barcode_summary_tsv",
            "tmp/asv_barcode_summary.tsv",
            "--primers_file",
            str(primers),
            "--alpha", str(config.TARGET_BACKGROUND_ALPHA),
        ],
        output_dir,
    )
    # 2026-08-10: Reconstruct sequences only for requested PV or reference analyses.
    # Reason: blank Mode targets require only cell-level detection.
    if pv_requested or reference is not None:
        sequence_command = [
            python,
            script("sub_arg_database_revised.py"),
            "--filtered_counts_summary_arg_tsv",
            "tmp/filtered_counts_summary_arg.tsv",
            "--b_with_ids",
            "tmp/b_with_ids.txt",
            "--arg_packets",
            "tmp/packets_arg",
            "--r1_fastq",
            str(filtered_r1),
            "--r2_fastq",
            str(filtered_r2),
            "--primers_file",
            str(primers),
            "--filtered_sub_arg_barcode_summary_tsv",
            "tmp/cell_target_matrix.tsv",
            "--alpha", str(config.TARGET_SEQUENCE_ALPHA),
            "--max_shift_sub_arg", str(config.TARGET_SEQUENCE_MAX_SHIFT),
            "--max_mm_sub_arg", str(config.TARGET_SEQUENCE_MAX_MISMATCHES),
            "--p_match", str(config.MLE_P_MATCH),
            "--p_none", str(config.MLE_P_NONE),
            "--p_error", str(config.MLE_P_ERROR),
            "--alpha_prior", str(config.MLE_ALPHA_PRIOR),
            "--beta_prior", str(config.MLE_BETA_PRIOR),
            "--min_confidence", str(config.MLE_MIN_CONFIDENCE),
            "--min_noise_reads", str(config.MLE_MIN_NOISE_READS),
            "--noise_cutoff_ratio", str(config.MLE_NOISE_CUTOFF_RATIO),
        ]
        if reference is not None:
            sequence_command.append("--include_all_targets")
        _run_step("Reconstruct target sequences", sequence_command, output_dir)
    else:
        shutil.copyfile(
            output_dir / "tmp/filtered_counts_summary_arg.tsv",
            output_dir / "tmp/cell_target_matrix.tsv",
        )

    if pv_requested:
        pv_command = [
            python,
            script("phase_variation.py"),
            "--sequence_list",
            "tmp/sub_arg_seqs_list.txt",
            "--cell_matrix",
            "tmp/cell_target_matrix.tsv",
            "--primers_file",
            str(primers),
            "--output_tsv",
            "tmp/cell_phase_variation.tsv",
        ]
        _run_step(
            "Analyze phase variation",
            pv_command,
            output_dir,
        )

    if reference is not None:
        _run_step(
            "Match reconstructed sequences to reference",
            [
                python,
                script("blastn_sub_arg.py"),
                "--sub_arg_seqs_list",
                "tmp/sub_arg_seqs_list.txt",
                "--input_fasta",
                str(reference),
                "--blastn_sub_arg_tsv",
                "tmp/reference_matches.tsv",
                "--db",
                "tmp/blast_db/dota_seq_analyzer",
                "--final_barcode_summary_tsv",
                "tmp/cell_target_matrix.tsv",
                "--first_gene_column_num",
                "14",
            ],
            output_dir,
        )
    export_command = [
        python,
        script("export_results.py"),
        "--input_tsv",
        "tmp/cell_target_matrix.tsv",
        "--primers_file",
        str(primers),
        "--output_jsonl",
        "dota_seq_analyzer_results.jsonl",
    ]
    if pv_requested:
        export_command.extend([
            "--phase_variation_tsv", "tmp/cell_phase_variation.tsv"
        ])
    _run_step("Export JSONL", export_command, output_dir)

    canonical_command = [
        python, script("canonical_results.py"),
        "--runtime-config", "tmp/runtime_config.json",
        "--raw-cell-table", "tmp/asv_barcode_summary.tsv",
        "--filtered-counts", "tmp/filtered_counts_summary_arg.tsv",
        "--final-cell-table", "tmp/cell_target_matrix.tsv",
        "--global-asv", "tmp/global_asv.tsv",
        "--read-qc-stats", "tmp/read_qc_stats.json",
        "--barcode-cluster-stats", "tmp/barcode_cluster_stats.json",
        "--barcode-filter-stats", "tmp/barcode_filter_stats.json",
        "--asv-stats", "tmp/asv_stats.json",
        "--manifest", "tmp/16s_r1_primer_starts.tsv",
        "--packet-16s", "tmp/packets_16s",
        "--packet-target", "tmp/packets_arg",
        "--packet-unclassified", "tmp/packets_unclassified",
        "--target-filter-stats", "tmp/stats_filtering_summary_arg.tsv",
        "--staging-dir", f"tmp/canonical.{run_id}",
        "--output-dir", ".",
    ]
    if reconstruction_requested:
        canonical_command.extend([
            "--sequence-list", "tmp/sub_arg_seqs_list.txt",
            "--sequence-qc", "tmp/extra_mle_info_sub_arg.tsv",
        ])
    if pv_requested:
        canonical_command.extend([
            "--phase-variation", "tmp/cell_phase_variation.tsv",
        ])
    if reference is not None:
        canonical_command.extend([
            "--reference-matches", "tmp/reference_matches.tsv",
        ])
    _run_step("Build and validate canonical results", canonical_command, output_dir)

    # 2026-09-10: Promote canonical-derived reports after canonical validation.
    # Reason: reports/ contains user-facing views, while pre-canonical working TSVs stay in tmp/.
    _run_step(
        "Generate canonical-derived reports",
        [
            python,
            script("canonical_reports.py"),
            "--canonical-dir", ".",
            "--output-dir", "reports",
        ],
        output_dir,
    )

    # 2026-09-09: Generate final-result figures after canonical validation.
    # Reason: the ASV-target and primer-balance figures now read canonical v3 directly.
    _run_step(
        "Generate figures",
        [
            python,
            script("figures_program.py"),
            "--use_asvs_str",
            "yes",
            "--unfiltered_barcode_summary_tsv",
            "tmp/unfiltered_barcode_summary.tsv",
            "--primers_file",
            str(primers),
            "--b_with_ids",
            "tmp/b_with_ids.txt",
            "--canonical_dir",
            ".",
        ],
        output_dir,
    )

    # 2026-09-10: Remove intermediates only after every final artifact succeeds.
    # Reason: failed runs retain their incomplete staging directory for debugging.
    if not args.keep_tmp:
        shutil.rmtree(output_dir / "tmp")

    os.replace(output_dir, requested_output_dir)
    print(f"\nDoTA-Seq Analyzer complete: {requested_output_dir}", flush=True)


if __name__ == "__main__":
    main()
