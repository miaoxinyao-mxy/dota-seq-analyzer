import csv
import json
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))

import figures_program  # noqa: E402


class PrimerBalanceFigureTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.primers = self.root / "primers.csv"
        self.primers.write_text(
            "Primer,F,R,Mode\n"
            "16s,AAAA,TTTT,\n"
            "target_1,CCCC,GGGG,\n"
            "target_2,ACAC,TGTG,\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.tempdir.cleanup()

    def test_empty_surviving_cell_table_writes_figure(self):
        summary = self.root / "empty.tsv"
        summary.write_text(
            "Barcode\tTotal # of 16s reads\ttarget_1\ttarget_2\n",
            encoding="utf-8",
        )
        output = self.root / "primer_balance.png"

        figures_program.make_primer_balance_figure(
            str(output), str(summary), str(self.primers), 72)

        self.assertTrue(output.is_file())
        self.assertGreater(output.stat().st_size, 0)

    def test_degenerate_barcode_group_sizes_do_not_emit_plot_warnings(self):
        import matplotlib.pyplot as plt

        for values in ([1, 1, 1], [0, 0], []):
            fig, ax = plt.subplots()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                figures_program.jackpottocurve(fig, ax, values, "All", 0)
            plt.close(fig)
            self.assertEqual([], [str(item.message) for item in caught])

    def test_all_zero_histogram_bins_do_not_emit_log_warning(self):
        summary = self.root / "below_plot_range.tsv"
        summary.write_text(
            "Barcode\tTotal # of 16s reads\ttarget_1\ttarget_2\n"
            f"{'A' * 20}\t19\t1\t1\n",
            encoding="utf-8",
        )
        output = self.root / "primer_balance.png"

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            figures_program.make_primer_balance_figure(
                str(output), str(summary), str(self.primers), 72)

        self.assertTrue(output.is_file())
        self.assertFalse(
            any("no positive values" in str(item.message).lower() for item in caught),
            [str(item.message) for item in caught],
        )


class CanonicalFigureDataTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.run_id = "figure-test-run"
        self.target_names = ["target_2", "target_1"]
        (self.root / "run_summary.json").write_text(json.dumps({
            "run_id": self.run_id,
            "target_panel": [
                {"target_name": name, "mode": "detect"}
                for name in self.target_names
            ],
        }), encoding="utf-8")
        # Physical JSONL order deliberately differs from numerical ASV order.
        asvs = [
            self._asv("ASV_3", 40),
            self._asv("ASV_1", 50),
            self._asv("ASV_2", 29),
        ]
        self._write_jsonl(self.root / "asvs.jsonl", asvs)

        cells = []
        for index in range(50):
            taxonomy = "Alpha" if index < 25 else "Beta"
            targets = []
            if index < 10:
                targets.append(self._target("target_2", 2, 1, "detected"))
            if index < 5:
                targets.append(self._target("target_1", 3, 0, "filtered_out"))
            cells.append(self._cell("ASV_1", index, taxonomy, targets))
        for index in range(40):
            targets = []
            if index < 20:
                targets.append(self._target("target_1", 1, 1, "detected"))
            cells.append(self._cell("ASV_3", index, "Gamma", targets))
        self.cells = cells
        self._write_jsonl(self.root / "cells.jsonl", cells)

        self.primers = self.root / "primers.csv"
        self.primers.write_text(
            "Primer,F,R,Mode\n"
            "16s,AAAA,TTTT,\n"
            "target_2,CCCC,GGGG,\n"
            "target_1,ACAC,TGTG,\n",
            encoding="utf-8",
        )
        self.legacy_matrix = self.root / "legacy_matrix.tsv"
        self.raw_counts = self.root / "raw_counts.tsv"
        self.global_asvs = self.root / "global_asv.tsv"
        self._write_legacy_inputs()

    def tearDown(self):
        self.tempdir.cleanup()

    def _asv(self, asv_id, count):
        return {
            "run_id": self.run_id,
            "asv_id": asv_id,
            "final_surviving_cell_count": count,
        }

    def _target(self, name, raw_count, filtered_count, assignment_type):
        return {
            "target_name": name,
            "raw_read_count": raw_count,
            "filtered_read_count": filtered_count,
            "assignment_type": assignment_type,
        }

    def _cell(self, asv_id, index, taxonomy, targets):
        return {
            "run_id": self.run_id,
            "cell_barcode": f"{asv_id}_{index}",
            "taxonomy": {"ranks": {
                "domain": "Bacteria", "phylum": taxonomy, "class": None,
                "order": None, "family": None, "genus": None,
                "species": None,
            }},
            "taxonomy_evidence": {"total_16s_reads": 7},
            "asv": {"asv_id": asv_id},
            "targets": targets,
        }

    def _write_jsonl(self, path, records):
        path.write_text(
            "".join(json.dumps(record) + "\n" for record in records),
            encoding="utf-8",
        )

    def _lineage(self, cell):
        return figures_program._canonical_taxonomy_lineage(
            cell["taxonomy"]["ranks"])

    def _write_legacy_inputs(self):
        with self.legacy_matrix.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([
                "Barcode", "Predicted taxonomy", "Assigned_core_asv",
                *self.target_names,
            ])
            for cell in self.cells:
                by_name = {target["target_name"]: target for target in cell["targets"]}
                values = []
                for name in self.target_names:
                    target = by_name.get(name)
                    values.append(
                        0 if target is None else target["filtered_read_count"])
                writer.writerow([
                    cell["cell_barcode"], self._lineage(cell),
                    cell["asv"]["asv_id"], *values,
                ])
        with self.raw_counts.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Barcode", "Total # of 16s reads", *self.target_names])
            for cell in self.cells:
                by_name = {target["target_name"]: target for target in cell["targets"]}
                writer.writerow([
                    cell["cell_barcode"], 7,
                    *(by_name.get(name, {}).get("raw_read_count", 0)
                      for name in self.target_names),
                ])
        self.global_asvs.write_text(
            "Core_ASV_ID\tcell_count\tcore_sequence\n"
            "ASV_1\t50\tA|T\n"
            "ASV_3\t40\tC|G\n"
            "ASV_2\t29\tG|C\n",
            encoding="utf-8",
        )

    def test_canonical_heatmap_data_matches_legacy_data(self):
        canonical, target_names = (
            figures_program.create_asv_target_matrix_from_canonical(
                str(self.root), min_cells_per_asv=30))
        legacy, legacy_target_names = figures_program.create_asv_arg_matrix(
            True, str(self.legacy_matrix), str(self.root / "legacy_plot.tsv"),
            30, 2, str(self.root / "unused_tax.tsv"), str(self.global_asvs))

        self.assertEqual(["ASV_1", "ASV_3"], canonical.index.to_list())
        self.assertEqual(self.target_names, target_names)
        self.assertEqual(legacy_target_names, target_names)
        self.assertEqual(legacy.index.to_list(), canonical.index.to_list())
        self.assertEqual(
            legacy["Predicted taxonomy"].to_list(),
            canonical["Predicted taxonomy"].to_list())
        self.assertTrue(
            np.array_equal(
                legacy[target_names].to_numpy(),
                canonical[target_names].to_numpy()))
        # Counter.most_common() retains the first-seen taxonomy in this tie.
        self.assertIn("P - Alpha", canonical.loc["ASV_1", "Predicted taxonomy"])

    def test_canonical_primer_balance_vectors_match_legacy_data(self):
        canonical = figures_program.get_target_to_16s_ratios_from_canonical(
            str(self.root))
        legacy = figures_program.get_ARG_to_16s_ratios(
            str(self.raw_counts), str(self.primers))

        self.assertEqual(self.target_names, list(canonical))
        self.assertEqual(legacy, canonical)
        self.assertEqual(10, len(canonical["target_2"]))
        # Includes five raw-positive/final-filtered cells plus twenty detected cells.
        self.assertEqual(25, len(canonical["target_1"]))

    def test_canonical_global_taxonomy_matches_legacy_bytes(self):
        legacy_output = self.root / "global_taxonomy_legacy.tsv"
        canonical_output = self.root / "global_taxonomy_canonical.tsv"
        figures_program.write_global_tax_classification_file(
            str(self.legacy_matrix), str(legacy_output))
        figures_program.write_global_tax_classification_from_canonical(
            str(self.root), str(canonical_output))

        self.assertEqual(legacy_output.read_bytes(), canonical_output.read_bytes())


if __name__ == "__main__":
    unittest.main()
