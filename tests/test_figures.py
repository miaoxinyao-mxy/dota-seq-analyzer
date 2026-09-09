import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import matplotlib

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


if __name__ == "__main__":
    unittest.main()
