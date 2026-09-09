import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))

import canonical_results as canonical  # noqa: E402


class CanonicalResultTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        base = {
            "Predicted taxonomy": [
                "R1 - Bacteria | P - Bacillota | C - Bacilli | O - None | F - None | G - None | S - None"
            ],
            "Confidence": [0.99], "Contamination": [0.01],
            "Total # of 16s reads": [12], "Technical noise count": [1],
            "Reads_used_for_ASV": [11], "Raw unique_core_sequences": [2],
            "Dominant_raw_read_count": [10], "Coexisting_2bp_reads": [1],
            "Unauthorized_secondary_reads": [0], "Final_cell_asv_reads": [11],
            "Max_internal_distance": [1], "Assigned_core_asv": ["ASV_1"],
            "Status": ["single_ASV"], "TEM": [4],
        }
        self.raw = pd.DataFrame(base, index=["A" * 20])
        self.filtered = self.raw.copy()
        self.filtered["TEM"] = 3
        self.final = self.filtered.copy()
        self.final["TEM"] = "TEM_seq_1"
        for name, frame in (("raw.tsv", self.raw), ("filtered.tsv", self.filtered), ("final.tsv", self.final)):
            frame.to_csv(self.root / name, sep="\t", index_label="Barcode")
        (self.root / "asvs.tsv").write_text(
            "Core_ASV_ID\tcell_count\tcore_sequence\nASV_1\t999\tAAAA|TTTT\n")
        (self.root / "sequences.tsv").write_text(
            "Sub-ARG_Arbitrary_Name\tCell_count\tCore_sequence\n"
            "TEM_seq_1\t999\tCCCC|GGGG\n")

    def tearDown(self):
        self.tempdir.cleanup()

    def valid_summary(self):
        mle = {
            "p_match": 0.9, "p_none": 0.09, "p_error": 0.01,
            "alpha_prior": 1.0, "beta_prior": 9.0,
            "minimum_confidence": 0.95, "minimum_noise_reads": 2,
            "noise_cutoff_ratio": 0.05,
        }
        return {
            "schema_version": "3.0.0", "record_type": "run_summary",
            "run_id": "550e8400-e29b-41d4-a716-446655440000",
            "status": "completed",
            "software": {"name": "dota-seq-analyzer", "version": "0.1.0", "git_commit": None},
            "inputs": {"r1": "r1.fastq", "r2": "r2.fastq", "primers": "primers.csv", "taxonomy_database": "db", "reference_fasta": None},
            "target_panel": [{"target_name": "TEM", "mode": "detect"}],
            "parameters": {
                "analysis_workers": 1,
                "read_qc": {"minimum_read_length": 130, "minimum_mean_phred": 25, "barcode_length": 20, "minimum_barcode_q25_bases": 15},
                "primer_matching": {"maximum_shift": 4, "maximum_mismatches": 4, "r2_primer_start": 42, "valid_16s_r1_starts": [0]},
                "barcode_clustering": {"barcode_length": 20, "maximum_shift": 1},
                "cell_taxonomy_filtering": {"minimum_16s_reads": 5, "maximum_contamination": 0.1, "minimum_cells_per_taxon": 10},
                "taxonomy_mle": mle,
                "asv": {"r1_start": 30, "r1_end": 120, "r2_start": 70, "r2_end": 120, "maximum_distance": 3, "maximum_shift": 3, "minimum_reads": 5, "mixed_ratio_threshold": 0.1, "filter_corrupted_single_asv": False, "taxonomy_conflict_minimum_cells": 10, "taxonomy_conflict_dominant_phylum_fraction": 0.99},
                "target_background_filtering": {"alpha": 0.05},
                "target_sequence_reconstruction": {"performed": True, "maximum_shift": 2, "maximum_mismatches": 0, "alpha": 0.05, "r1_start": 30, "r1_end": 120, "r2_start": 70, "r2_end": 120, "include_all_targets": False, "mle": mle.copy()},
            },
            "read_qc": {"raw_read_pairs": 0, "passed_read_pairs": 0, "filtered_read_pairs": 0, "rejected": {"read_length_or_quality_length": 0, "mean_phred": 0, "barcode_length": 0, "barcode_q25": 0}},
            "read_classification": {"accepted_16s_reads": 0, "target_reads": 0, "unclassified_reads": 0},
            "16s_r1_primer_starts": {"0": 0},
            "barcode_funnel": {"raw_unique_barcodes": 1, "clustered_barcodes": 1, "after_stage1_taxonomy_filter": 1, "after_minimum_cells_per_taxon": 1, "before_asv_filter": 1, "after_asv_status_filter": 1, "asv_taxonomy_conflicts_removed": 0, "final_cells": 1},
            "asv_summary": {"final_asv_count": 1},
            "target_filtering_summary": [{"target_name": "TEM", "original_positive_cells": 1, "filtered_out_cells": 0, "remaining_positive_cells": 1, "retention_percent": 100.0}],
        }

    def test_counts_are_recalculated_from_final_cell_references(self):
        catalog = canonical._load_sequence_catalog(self.root / "sequences.tsv")
        cells = canonical.build_cells(
            "550e8400-e29b-41d4-a716-446655440000", self.root / "raw.tsv", self.root / "filtered.tsv",
            self.root / "final.tsv", ["TEM"], catalog, None)
        asvs = canonical.build_asvs("550e8400-e29b-41d4-a716-446655440000", cells, self.root / "asvs.tsv")
        targets = canonical.build_target_sequences("550e8400-e29b-41d4-a716-446655440000", cells, catalog, None, None)
        self.assertEqual(asvs[0]["final_surviving_cell_count"], 1)
        self.assertEqual(targets[0]["final_surviving_cell_count"], 1)
        self.assertEqual(cells[0]["targets"][0]["assignment_type"], "sequence_cluster")

    def test_unresolved_and_filtered_assignment_names(self):
        catalog = canonical._load_sequence_catalog(self.root / "sequences.tsv")
        self.final["TEM"] = "TEM_parent"
        self.final.to_csv(self.root / "final.tsv", sep="\t", index_label="Barcode")
        cells = canonical.build_cells(
            "550e8400-e29b-41d4-a716-446655440000", self.root / "raw.tsv", self.root / "filtered.tsv",
            self.root / "final.tsv", ["TEM"], catalog, None)
        self.assertEqual(cells[0]["targets"][0]["assignment_type"], "sequence_unresolved")
        self.filtered["TEM"] = 0
        self.filtered.to_csv(self.root / "filtered.tsv", sep="\t", index_label="Barcode")
        cells = canonical.build_cells(
            "550e8400-e29b-41d4-a716-446655440000", self.root / "raw.tsv", self.root / "filtered.tsv",
            self.root / "final.tsv", ["TEM"], catalog, None)
        self.assertEqual(cells[0]["targets"][0]["assignment_type"], "filtered_out")

    def test_validator_rejects_foreign_key_and_count_mismatch(self):
        catalog = canonical._load_sequence_catalog(self.root / "sequences.tsv")
        cells = canonical.build_cells(
            "550e8400-e29b-41d4-a716-446655440000", self.root / "raw.tsv", self.root / "filtered.tsv",
            self.root / "final.tsv", ["TEM"], catalog, None)
        asvs = canonical.build_asvs("550e8400-e29b-41d4-a716-446655440000", cells, self.root / "asvs.tsv")
        targets = canonical.build_target_sequences("550e8400-e29b-41d4-a716-446655440000", cells, catalog, None, None)
        summary = self.valid_summary()
        canonical.validate_canonical("550e8400-e29b-41d4-a716-446655440000", cells, asvs, targets, summary)
        asvs[0]["final_surviving_cell_count"] = 2
        with self.assertRaisesRegex(ValueError, "ASV count mismatch"):
            canonical.validate_canonical("550e8400-e29b-41d4-a716-446655440000", cells, asvs, targets, summary)


    def test_validator_rejects_run_schema_violation(self):
        catalog = canonical._load_sequence_catalog(self.root / "sequences.tsv")
        cells = canonical.build_cells(
            "550e8400-e29b-41d4-a716-446655440000", self.root / "raw.tsv", self.root / "filtered.tsv",
            self.root / "final.tsv", ["TEM"], catalog, None)
        asvs = canonical.build_asvs("550e8400-e29b-41d4-a716-446655440000", cells, self.root / "asvs.tsv")
        targets = canonical.build_target_sequences("550e8400-e29b-41d4-a716-446655440000", cells, catalog, None, None)
        summary = self.valid_summary()
        del summary["parameters"]["read_qc"]["minimum_read_length"]
        with self.assertRaisesRegex(ValueError, "minimum_read_length"):
            canonical.validate_canonical("550e8400-e29b-41d4-a716-446655440000", cells, asvs, targets, summary)


if __name__ == "__main__":
    unittest.main()
