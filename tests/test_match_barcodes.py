import random
import sys
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))

from match_barcodes_to_IDs_revised import (  # noqa: E402
    check_barcodes_match_revised,
    create_clustered_b_with_ids,
)


def old_create_clustered_b_with_ids(bcs_with_counts, all_b_with_ids, max_shift=1):
    sorted_bcs = [key for key, _ in bcs_with_counts.most_common()]
    clustered = {}
    for bc_s in sorted_bcs:
        for bc_d in clustered:
            if check_barcodes_match_revised(bc_s, bc_d, max_shift):
                for i in range(len(clustered[bc_d])):
                    clustered[bc_d][i].extend(all_b_with_ids[bc_s][i])
                break
        else:
            clustered[bc_s] = all_b_with_ids[bc_s]
    return clustered


def make_inputs(barcodes, counts=None):
    if counts is None:
        counts = [1] * len(barcodes)
    counter = Counter()
    all_ids = {}
    for number, (barcode, count) in enumerate(zip(barcodes, counts)):
        counter[barcode] = count
        all_ids[barcode] = [[f"{barcode}_{number}_{i}"] for i in range(3)]
    return counter, all_ids


class BarcodeClusteringRegressionTests(unittest.TestCase):
    def assert_same_result(self, barcodes, counts, max_shift):
        old_counts, old_ids = make_inputs(barcodes, counts)
        new_counts, new_ids = make_inputs(barcodes, counts)
        self.assertEqual(
            old_create_clustered_b_with_ids(old_counts, old_ids, max_shift),
            create_clustered_b_with_ids(new_counts, new_ids, max_shift),
        )

    def test_exact_and_shift_matches(self):
        self.assert_same_result(
            ["ACGTACGT", "ACGTACGT", "TACGTACG", "GTACGTAC", "CCCCCCCC"],
            [10, 5, 4, 3, 2],
            1,
        )

    def test_shift_zero_one_two(self):
        barcodes = ["ACGTACGTAC", "TACGTACGTA", "GTACGTACGT", "CCCCCCCCCC"]
        for max_shift in (0, 1, 2):
            self.assert_same_result(barcodes, [10, 4, 3, 2], max_shift)

    def test_multiple_candidates_and_ties(self):
        self.assert_same_result(
            ["ACGTACGT", "TACGTACG", "ACGTACGA", "GTACGTAC", "TTTTTTTT"],
            [5, 5, 5, 4, 1],
            1,
        )

    def test_unexpected_lengths_fall_back(self):
        self.assert_same_result(
            ["ACGTACGT", "TACGTACG", "ACGT", "CCCCCCCC"],
            [10, 4, 3, 2],
            1,
        )

    def test_random_synthetic_data(self):
        random.seed(20260904)
        alphabet = "ACGT"
        barcodes = ["".join(random.choice(alphabet) for _ in range(12)) for _ in range(80)]
        counts = [random.randint(1, 20) for _ in barcodes]
        for max_shift in (0, 1, 2):
            self.assert_same_result(barcodes, counts, max_shift)


if __name__ == "__main__":
    unittest.main()
