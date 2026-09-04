#!/usr/bin/env python3
"""Run the complete DoTA-Seq Analyzer workflow from one public command."""

# 2026-08-11: Rename the public entry point to DoTA-Seq Analyzer.
# Reason: users should run R1/R2 data with one command instead of invoking internal modules manually.

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

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


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="dota-seq-analyzer",
        description="Profile targeted genes and phase variation in single-cell DoTA-Seq data.",
    )
    parser.add_argument("-1", "--r1", required=True, help="R1 FASTQ file")
    parser.add_argument("-2", "--r2", required=True, help="R2 FASTQ file")
    parser.add_argument("-p", "--primers", required=True, help="DoTA-Seq primer CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-@", "--threads", dest="analysis_workers", type=int, default=1, metavar="INT", help="Number of parallel workers/threads used by DoTA-seq analysis")
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
    output_dir = Path(args.output).expanduser().resolve()
    source_dir = Path(__file__).resolve().parent

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

    for directory in (
        output_dir,
        output_dir / "tmp",
        output_dir / "reports",
        output_dir / "figures",
    ):
        directory.mkdir(parents=True, exist_ok=True)

    python = sys.executable
    script = lambda name: str(source_dir / name)

    _run_step(
        "Validate inputs",
        [
            python,
            script("validate_inputs.py"),
            "--r1_fastq",
            str(r1),
            "--r2_fastq",
            str(r2),
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
            str(r1),
            "--r2_fastq",
            str(r2),
            "--primers_filename",
            str(primers),
            "--threads",
            str(args.analysis_workers),
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
            str(r1),
            "--r2_fastq",
            str(r2),
            "--primers_filename",
            str(primers),
            "--kraken_output",
            "tmp/kraken.output",
            "--kraken_report",
            "tmp/kraken.report",
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
            "--min_cells_per_taxon",
            str(args.min_cells_per_taxon),
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
            "--primers_file",
            str(primers),
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
            str(r1),
            "--r2_fastq",
            str(r2),
            "--primers_file",
            str(primers),
            "--filtered_sub_arg_barcode_summary_tsv",
            "reports/cell_target_matrix.tsv",
        ]
        if reference is not None:
            sequence_command.append("--include_all_targets")
        _run_step("Reconstruct target sequences", sequence_command, output_dir)
    else:
        shutil.copyfile(
            output_dir / "tmp/filtered_counts_summary_arg.tsv",
            output_dir / "reports/cell_target_matrix.tsv",
        )

    if pv_requested:
        pv_command = [
            python,
            script("phase_variation.py"),
            "--sequence_list",
            "tmp/sub_arg_seqs_list.txt",
            "--cell_matrix",
            "reports/cell_target_matrix.tsv",
            "--primers_file",
            str(primers),
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
                "reports/reference_matches.tsv",
                "--db",
                "tmp/blast_db/dota_seq_analyzer",
                "--final_barcode_summary_tsv",
                "reports/cell_target_matrix.tsv",
                "--first_gene_column_num",
                "14",
            ],
            output_dir,
        )
    _run_step(
        "Generate figures",
        [
            python,
            script("figures_program.py"),
            "--use_asvs_str",
            "yes",
            "--unfiltered_barcode_summary_tsv",
            "tmp/unfiltered_barcode_summary.tsv",
            "--final_asv_barcode_summary_tsv",
            "reports/cell_target_matrix.tsv",
            "--asv_barcode_summary_no_sub_args_tsv",
            "tmp/asv_barcode_summary.tsv",
            "--primers_file",
            str(primers),
            "--b_with_ids",
            "tmp/b_with_ids.txt",
            "--global_asv_tsv",
            "tmp/global_asv.tsv",
            "--first_gene_column_num",
            "14",
        ],
        output_dir,
    )
    export_command = [
        python,
        script("export_results.py"),
        "--input_tsv",
        "reports/cell_target_matrix.tsv",
        "--primers_file",
        str(primers),
        "--output_jsonl",
        "dota_seq_analyzer_results.jsonl",
    ]
    if pv_requested:
        export_command.extend([
            "--phase_variation_tsv", "reports/cell_phase_variation.tsv"
        ])
    _run_step("Export JSONL", export_command, output_dir)

    print(f"\nDoTA-Seq Analyzer complete: {output_dir}", flush=True)


if __name__ == "__main__":
    main()
