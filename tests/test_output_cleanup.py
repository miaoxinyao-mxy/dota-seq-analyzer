import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))

from dota_seq_analyzer import cli  # noqa: E402
import filter_args  # noqa: E402
import filter_sub_args  # noqa: E402
import sub_arg_database_revised  # noqa: E402
import phase_variation  # noqa: E402
import blastn_sub_arg  # noqa: E402


class OutputCleanupCLITests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.r1 = self.root / "reads_R1.fastq"
        self.r2 = self.root / "reads_R2.fastq"
        self.r1.write_text("", encoding="utf-8")
        self.r2.write_text("", encoding="utf-8")
        self.primers = self.root / "primers.csv"
        self.primers.write_text(
            "Primer,F,R,Mode\n"
            "16s,AAAA,TTTT,\n"
            "PV,CCCC,GGGG,ssr\n",
            encoding="utf-8",
        )
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">reference\nACGTACGT\n", encoding="utf-8")
        self.taxonomy = self.root / "taxonomy"
        self.taxonomy.mkdir()

    def tearDown(self):
        self.tempdir.cleanup()

    @staticmethod
    def _argument(command, option):
        return command[command.index(option) + 1]

    def _run_cli(self, output, keep_tmp=False, fail_at=None):
        calls = []

        def fake_step(name, command, staging):
            calls.append((name, list(command)))
            tmp = staging / "tmp"
            self.assertTrue(tmp.is_dir())
            (tmp / "lifecycle_sentinel.txt").write_text("present\n", encoding="utf-8")
            if name == fail_at:
                raise subprocess.CalledProcessError(1, command)
            if name == "Build and validate canonical results":
                (staging / "cells.jsonl").write_text("", encoding="utf-8")
                (staging / "asvs.jsonl").write_text("", encoding="utf-8")
                (staging / "target_sequences.jsonl").write_text("", encoding="utf-8")
                (staging / "run_summary.json").write_text(
                    json.dumps({"run_id": "fixture"}), encoding="utf-8")
            elif name == "Generate canonical-derived reports":
                reports = staging / "reports"
                reports.mkdir()
                for report in (
                    "summary.tsv", "cell_target_matrix.tsv", "asv_summary.tsv",
                    "taxonomy_summary.tsv", "target_summary.tsv",
                    "cell_phase_variation.tsv", "reference_matches.tsv",
                ):
                    (reports / report).write_text("fixture\n", encoding="utf-8")
            elif name == "Generate figures":
                self.assertTrue((tmp / "lifecycle_sentinel.txt").is_file())
                figures = staging / "figures"
                for figure in (
                    "taxa_target_table.png", "primer_balance_qc.png",
                    "barcode_group_size_qc.png",
                ):
                    (figures / figure).write_bytes(b"fixture")

        argv = [
            "dota-seq-analyzer",
            "-1", str(self.r1), "-2", str(self.r2),
            "-p", str(self.primers), "-o", str(output),
            "--taxonomy-db", str(self.taxonomy),
            "-r", str(self.reference),
        ]
        if keep_tmp:
            argv.append("--keep-tmp")
        with mock.patch.object(sys, "argv", argv), \
             mock.patch.object(cli, "_run_step", side_effect=fake_step), \
             mock.patch.object(cli, "check_reference_fasta", return_value="Valid"), \
             mock.patch.object(cli.uuid, "uuid4", return_value="phase2c-test-run"):
            if fail_at is None:
                cli.main()
            else:
                with self.assertRaises(subprocess.CalledProcessError):
                    cli.main()
        return calls

    def test_working_paths_reports_and_default_cleanup(self):
        output = self.root / "output"
        calls = self._run_cli(output)
        by_name = {name: command for name, command in calls}

        reconstruction = by_name["Reconstruct target sequences"]
        self.assertEqual(
            "tmp/cell_target_matrix.tsv",
            self._argument(reconstruction, "--filtered_sub_arg_barcode_summary_tsv"),
        )
        pv = by_name["Analyze phase variation"]
        self.assertEqual("tmp/cell_target_matrix.tsv",
                         self._argument(pv, "--cell_matrix"))
        self.assertEqual("tmp/cell_phase_variation.tsv",
                         self._argument(pv, "--output_tsv"))
        reference = by_name["Match reconstructed sequences to reference"]
        self.assertEqual("tmp/reference_matches.tsv",
                         self._argument(reference, "--blastn_sub_arg_tsv"))
        self.assertEqual("tmp/cell_target_matrix.tsv",
                         self._argument(reference, "--final_barcode_summary_tsv"))
        legacy = by_name["Export JSONL"]
        self.assertEqual("tmp/cell_target_matrix.tsv",
                         self._argument(legacy, "--input_tsv"))
        self.assertEqual("tmp/cell_phase_variation.tsv",
                         self._argument(legacy, "--phase_variation_tsv"))
        canonical = by_name["Build and validate canonical results"]
        self.assertEqual("tmp/cell_target_matrix.tsv",
                         self._argument(canonical, "--final-cell-table"))
        self.assertEqual("tmp/cell_phase_variation.tsv",
                         self._argument(canonical, "--phase-variation"))
        self.assertEqual("tmp/reference_matches.tsv",
                         self._argument(canonical, "--reference-matches"))

        report_command = by_name["Generate canonical-derived reports"]
        self.assertEqual(".", self._argument(report_command, "--canonical-dir"))
        self.assertEqual("reports", self._argument(report_command, "--output-dir"))
        self.assertNotIn(".phase2a_reports", " ".join(
            argument for _, command in calls for argument in command))
        self.assertLess(
            [name for name, _ in calls].index("Build and validate canonical results"),
            [name for name, _ in calls].index("Generate canonical-derived reports"),
        )
        self.assertLess(
            [name for name, _ in calls].index("Generate canonical-derived reports"),
            [name for name, _ in calls].index("Generate figures"),
        )

        figures = by_name["Generate figures"]
        self.assertEqual(".", self._argument(figures, "--canonical_dir"))
        self.assertEqual("tmp/unfiltered_barcode_summary.tsv",
                         self._argument(figures, "--unfiltered_barcode_summary_tsv"))
        self.assertEqual("tmp/b_with_ids.txt",
                         self._argument(figures, "--b_with_ids"))
        self.assertTrue(output.is_dir())
        self.assertFalse((output / "tmp").exists())
        self.assertTrue((output / "reports" / "summary.tsv").is_file())
        self.assertTrue((output / "figures" / "barcode_group_size_qc.png").is_file())

    def test_keep_tmp_preserves_complete_tmp(self):
        output = self.root / "kept"
        self._run_cli(output, keep_tmp=True)
        self.assertTrue((output / "tmp" / "lifecycle_sentinel.txt").is_file())

    def test_failure_preserves_incomplete_staging_and_tmp(self):
        output = self.root / "failed"
        self._run_cli(output, fail_at="Generate figures")
        self.assertFalse(output.exists())
        staging = self.root / ".failed.incomplete.phase2c-test-run"
        self.assertTrue((staging / "tmp" / "lifecycle_sentinel.txt").is_file())
        self.assertTrue((staging / "reports" / "summary.tsv").is_file())


class OptionalIntermediateTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_target_filter_binary_output_is_optional(self):
        primers = self.root / "primers.csv"
        primers.write_text(
            "Primer,F,R,Mode\n16s,AAAA,TTTT,\nTEM,CCCC,GGGG,\n",
            encoding="utf-8",
        )
        source = self.root / "asv.tsv"
        pd.DataFrame({
            "Barcode": ["A" * 20],
            "Total # of 16s reads": [10],
            "TEM": [1],
        }).set_index("Barcode").to_csv(source, sep="\t", index_label="Barcode")
        counts = self.root / "counts.tsv"
        stats = self.root / "stats.tsv"

        filter_args.filter_args(
            str(source), str(primers), str(counts), None, str(stats), 0.05)

        self.assertTrue(counts.is_file())
        self.assertTrue(stats.is_file())
        self.assertFalse((self.root / "filtered_binary_summary_arg.tsv").exists())

    def test_reconstruction_debug_outputs_are_optional(self):
        sequence_list = self.root / "sequences.tsv"
        sequence_list.write_text(
            "Sub-ARG_Arbitrary_Name\tCell_count\tCore_sequence\n"
            "TEM_seq_1\t5\tAAAA|TTTT\n",
            encoding="utf-8",
        )
        retained = filter_sub_args.run_sub_arg_denoising_pipeline(
            str(sequence_list), 0.05, None)
        self.assertEqual(["TEM_seq_1"], retained)

        df_arg = pd.DataFrame({"TEM": [1]}, index=["A" * 20])
        extra = self.root / "extra.tsv"
        generated = self.root / "generated.tsv"
        result = sub_arg_database_revised.write_sub_arg_barcode_summary(
            df_arg, [], {"A" * 20: {}}, ["TEM"],
            str(generated), None, str(extra),
            0.9, 0.05, 0.05, 1.0, 1.0, 0.95, 3, 0.1,
        )
        self.assertEqual("TEM", result.iloc[0]["TEM"])
        self.assertTrue(generated.is_file())
        self.assertTrue(extra.is_file())
        self.assertFalse((self.root / "sub_arg_barcode_summary.tsv").exists())


class OptionalWorkflowBehaviorTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def test_phase_variation_working_report(self):
        primers = self.root / "primers.csv"
        primers.write_text(
            "Primer,F,R,Mode\n16s,AAAA,TTTT,\nPV,CCCC,GGGG,ssr\n",
            encoding="utf-8",
        )
        sequences = self.root / "sequences.tsv"
        sequences.write_text(
            "Sub-ARG_Arbitrary_Name\tCell_count\tCore_sequence\n"
            "PV_seq_1\t1\tAAAA|TTT\n",
            encoding="utf-8",
        )
        matrix = self.root / "cell_matrix.tsv"
        pd.DataFrame(
            {"Barcode": ["A" * 20], "PV": ["PV_seq_1"]}
        ).set_index("Barcode").to_csv(matrix, sep="\t", index_label="Barcode")
        output = self.root / "cell_phase_variation.tsv"

        count = phase_variation.write_cell_calls(
            str(sequences), str(matrix), str(primers), str(output))

        report = pd.read_csv(output, sep="\t")
        self.assertEqual(1, count)
        self.assertEqual("PV_seq_1", report.iloc[0]["Sequence_assignment"])
        self.assertEqual("ssr_detected", report.iloc[0]["Call"])

    def test_reference_matching_working_report(self):
        sequences = self.root / "sequences.tsv"
        sequences.write_text(
            "Sub-ARG_Arbitrary_Name\tCell_count\tCore_sequence\n"
            "TEM_seq_1\t1\tAAAA|TTTT\n",
            encoding="utf-8",
        )
        matrix = self.root / "cell_matrix.tsv"
        pd.DataFrame(
            {"Barcode": ["A" * 20], "TEM": ["TEM_seq_1"]}
        ).set_index("Barcode").to_csv(matrix, sep="\t", index_label="Barcode")
        reference = self.root / "reference.fa"
        reference.write_text(">reference\nAAAATTTT\n", encoding="utf-8")
        output = self.root / "reference_matches.tsv"
        query = self.root / "query.fa"
        db = self.root / "blast_db" / "reference"
        blast_rows = (
            "read1\tref1\t100.000\t1\t4\t4\t10\t13\t0\tReference one\n"
            "read2\tref1\t100.000\t1\t4\t4\t20\t23\t0\tReference one\n"
        )

        def fake_blast(command, text):
            return blast_rows if command[0] == "blastn" else ""

        with mock.patch.object(
            blastn_sub_arg.subprocess, "check_output", side_effect=fake_blast
        ):
            blastn_sub_arg.make_summary_table(
                str(sequences), str(output), str(query), str(reference), str(db),
                str(matrix), 0,
            )

        report = pd.read_csv(output, sep="\t")
        self.assertEqual(["TEM_seq_1"], report["sub-ARG_name"].tolist())
        self.assertEqual(["Reference one"], report["subject"].tolist())
        self.assertEqual([10], report["read1_start"].tolist())
        self.assertEqual([23], report["read2_end"].tolist())


if __name__ == "__main__":
    unittest.main()
