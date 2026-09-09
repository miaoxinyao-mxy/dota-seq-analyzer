import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))
import canonical_reports  # noqa: E402

RUN_ID = "550e8400-e29b-41d4-a716-446655440000"


class CanonicalReportTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.canonical = self.root / "canonical"
        self.canonical.mkdir()
        self.output = self.root / "reports"
        self.run_summary = {
            "run_id": RUN_ID,
            "inputs": {"reference_fasta": "reference.fa"},
            "target_panel": [{"target_name": "TEM", "mode": "detect"},
                             {"target_name": "PV", "mode": "ssr"}],
            "read_qc": {"raw_read_pairs": 10, "passed_read_pairs": 8,
                        "filtered_read_pairs": 2},
            "read_classification": {"accepted_16s_reads": 4, "target_reads": 3,
                                    "unclassified_reads": 1},
            "barcode_funnel": {"raw_unique_barcodes": 2, "clustered_barcodes": 2,
                "after_stage1_taxonomy_filter": 1, "after_minimum_cells_per_taxon": 1,
                "before_asv_filter": 1, "after_asv_status_filter": 1,
                "asv_taxonomy_conflicts_removed": 0, "final_cells": 1},
            "asv_summary": {"final_asv_count": 1},
            "target_filtering_summary": [
                {"target_name": "TEM", "original_positive_cells": 1,
                 "filtered_out_cells": 1, "remaining_positive_cells": 0,
                 "retention_percent": 0.0},
                {"target_name": "PV", "original_positive_cells": 1,
                 "filtered_out_cells": 0, "remaining_positive_cells": 1,
                 "retention_percent": 100.0}],
            "phase_variation_summary": {"cell_calls": 1,
                                        "calls": {"ssr_detected": 1}},
        }
        self.cells = [{"run_id": RUN_ID, "cell_barcode": "A" * 20,
            "taxonomy": {"ranks": {"domain": "Bacteria", "phylum": "Bacillota",
                "class": "Bacilli", "order": None, "family": None,
                "genus": None, "species": None}, "confidence": 0.99,
                "contamination": 0.01},
            "taxonomy_evidence": {"total_16s_reads": 12,
                                  "technical_noise_reads": 1},
            "asv": {"asv_id": "ASV_1", "status": "single_ASV",
                "reads_used": 11, "raw_unique_core_sequences": 2,
                "dominant_raw_read_count": 10, "coexisting_2bp_reads": 1,
                "unauthorized_secondary_reads": 0, "final_cell_asv_reads": 11,
                "max_internal_distance": 1},
            "targets": [
                {"target_name": "TEM", "raw_read_count": 1,
                 "filtered_read_count": 0, "assignment_type": "filtered_out"},
                {"target_name": "PV", "raw_read_count": 3,
                 "filtered_read_count": 3, "assignment_type": "sequence_cluster",
                 "sequence_cluster_id": "PV_seq_1"}]}]
        self.asvs = [{"run_id": RUN_ID, "asv_id": "ASV_1",
            "core_sequence": {"r1": "AAAA", "r2": "TTTT"},
            "final_surviving_cell_count": 1}]
        self.sequences = [{"run_id": RUN_ID, "sequence_cluster_id": "PV_seq_1",
            "parent_target": "PV", "core_sequence": {"r1": "CCCC", "r2": "GGGG"},
            "final_surviving_cell_count": 1,
            "phase_variation": {"mode": "ssr", "call": "ssr_detected",
                "motif": "C", "repeat_count": 4, "start": 1, "end": 5,
                "direct_similarity": None, "evidence": "tandem_repeat_scan"},
            "reference_matches": [{"subject": "reference one", "r1_start": 1,
                "r1_end": 4, "r1_percent_identity": 100.0, "r2_start": 5,
                "r2_end": 8, "r2_percent_identity": 100.0}]}]
        (self.canonical / "run_summary.json").write_text(json.dumps(self.run_summary))
        for name, records in (("cells.jsonl", self.cells), ("asvs.jsonl", self.asvs),
                              ("target_sequences.jsonl", self.sequences)):
            (self.canonical / name).write_text(
                "".join(json.dumps(record) + "\n" for record in records))

    def tearDown(self):
        self.tempdir.cleanup()

    def test_generates_complete_reports_from_canonical_files(self):
        names = canonical_reports.generate_reports(self.canonical, self.output)
        self.assertEqual(names, ["summary.tsv", "cell_target_matrix.tsv",
            "asv_summary.tsv", "taxonomy_summary.tsv", "target_summary.tsv",
            "cell_phase_variation.tsv", "reference_matches.tsv"])
        matrix = pd.read_csv(self.output / "cell_target_matrix.tsv", sep="\t")
        self.assertEqual(matrix.loc[0, "TEM"], 0)
        self.assertEqual(matrix.loc[0, "PV"], "PV_seq_1")
        self.assertEqual(matrix.loc[0, "Assigned_core_asv"], "ASV_1")
        self.assertIn("P - Bacillota", matrix.loc[0, "Predicted taxonomy"])
        asvs = pd.read_csv(self.output / "asv_summary.tsv", sep="\t")
        self.assertEqual(asvs.loc[0, "core_sequence"], "AAAA|TTTT")
        pv = pd.read_csv(self.output / "cell_phase_variation.tsv", sep="\t")
        self.assertEqual(pv.loc[0, "Sequence_assignment"], "PV_seq_1")
        references = pd.read_csv(self.output / "reference_matches.tsv", sep="\t")
        self.assertEqual(references.loc[0, "read1_seq"], "CCCC")
        self.assertEqual(references.loc[0, "read1_percent_identity"], "100%")

    def test_optional_reports_are_absent_when_not_requested(self):
        self.run_summary["inputs"]["reference_fasta"] = None
        self.run_summary["target_panel"][1]["mode"] = "detect"
        (self.canonical / "run_summary.json").write_text(json.dumps(self.run_summary))
        (self.canonical / "target_sequences.jsonl").unlink()
        canonical_reports.generate_reports(self.canonical, self.output)
        self.assertFalse((self.output / "cell_phase_variation.tsv").exists())
        self.assertFalse((self.output / "reference_matches.tsv").exists())

    def test_empty_canonical_result_writes_header_only_reports(self):
        self.run_summary["inputs"]["reference_fasta"] = None
        self.run_summary["target_panel"][1]["mode"] = "detect"
        self.run_summary["barcode_funnel"]["final_cells"] = 0
        self.run_summary["asv_summary"]["final_asv_count"] = 0
        (self.canonical / "run_summary.json").write_text(json.dumps(self.run_summary))
        (self.canonical / "cells.jsonl").write_text("")
        (self.canonical / "asvs.jsonl").write_text("")
        (self.canonical / "target_sequences.jsonl").unlink()

        canonical_reports.generate_reports(self.canonical, self.output)

        matrix = pd.read_csv(self.output / "cell_target_matrix.tsv", sep="\t")
        asvs = pd.read_csv(self.output / "asv_summary.tsv", sep="\t")
        taxonomy = pd.read_csv(self.output / "taxonomy_summary.tsv", sep="\t")
        self.assertEqual(list(matrix.columns), [*canonical_reports.CELL_COLUMNS, "TEM", "PV"])
        self.assertTrue(matrix.empty)
        self.assertTrue(asvs.empty)
        self.assertTrue(taxonomy.empty)

    def test_phase_variation_report_requires_canonical_target_sequences(self):
        self.run_summary["inputs"]["reference_fasta"] = None
        (self.canonical / "run_summary.json").write_text(json.dumps(self.run_summary))
        (self.canonical / "target_sequences.jsonl").unlink()

        with self.assertRaisesRegex(ValueError, "phase-variation reports"):
            canonical_reports.generate_reports(self.canonical, self.output)

    def test_reference_report_requires_canonical_target_sequences(self):
        self.run_summary["target_panel"][1]["mode"] = "detect"
        (self.canonical / "run_summary.json").write_text(json.dumps(self.run_summary))
        (self.canonical / "target_sequences.jsonl").unlink()

        with self.assertRaisesRegex(ValueError, "reference-match reports"):
            canonical_reports.generate_reports(self.canonical, self.output)

    def test_mixed_run_ids_are_rejected(self):
        self.cells[0]["run_id"] = "different-run"
        (self.canonical / "cells.jsonl").write_text(json.dumps(self.cells[0]) + "\n")
        with self.assertRaisesRegex(ValueError, "cell run_id"):
            canonical_reports.generate_reports(self.canonical, self.output)


if __name__ == "__main__":
    unittest.main()
