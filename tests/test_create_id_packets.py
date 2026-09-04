import copy
import multiprocessing
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))

import create_ID_packets as packets  # noqa: E402


def old_determine_gene(primers, f_seq, r_seq, max_shift, max_mm, primer_start_num, mapping):
    gene = [0] * max(0, len(primers) - 2)
    for fwd_primer, rev_primer in mapping:
        if rev_primer == r_seq[primer_start_num:primer_start_num + len(rev_primer)] and fwd_primer == f_seq[:len(fwd_primer)]:
            return mapping[(fwd_primer, rev_primer)]
    for i in range(len(primers[1:])):
        parsed = primers[i + 1].strip().split(",")
        if packets.check_primer_match_seq(r_seq, parsed[2], max_shift, max_mm, primer_start_num) and packets.check_primer_match_seq(f_seq, parsed[1], max_shift, max_mm):
            if i == 0:
                gene = "16s"
            else:
                gene[i - 1] = 1
            break
    return gene


class PrimerPrecomputationTests(unittest.TestCase):
    def setUp(self):
        self.primers = [
            "name,fwd,rev\n",
            "16s,ACGT,TTAA\n",
            "ARG1,GGGG,CCCC\n",
            "ARG2,TTTT,AAAA\n",
        ]
        self.mapping = packets.make_primers_to_genes_dict(self.primers)
        self.records = packets.make_primer_records(self.primers)

    def test_old_and_optimized_assignments_match(self):
        reads = [
            ("ACGT" + "A" * 20, "A" * 42 + "TTAA" + "G" * 10),
            ("GGGA" + "A" * 20, "A" * 42 + "CCCC" + "G" * 10),
            ("TTTC" + "A" * 20, "A" * 42 + "AAAA" + "G" * 10),
            ("NNNN" + "A" * 20, "A" * 42 + "NNNN" + "G" * 10),
            ("ACG" , "A" * 42 + "TTA"),
        ]
        for f_seq, r_seq in reads:
            expected = old_determine_gene(self.primers, f_seq, r_seq, 4, 4, 42, self.mapping)
            actual = packets.determine_gene_revised(self.primers, f_seq, r_seq, 4, 4, 42, self.mapping, self.records)
            self.assertEqual(expected, actual)

    def test_chunk_worker_matches_serial_classification(self):
        records = [
            (0, "@id0", "@id0", "id0", "ACGT" + "A" * 20, "A" * 42 + "TTAA" + "G" * 10, "I" * 56),
            (1, "@id1", "@id1", "id1", "GGGG" + "A" * 20, "A" * 42 + "CCCC" + "G" * 10, "I" * 56),
        ]
        expected = [
            packets.determine_gene_revised(self.primers, record[4], record[5], 4, 4, 42, self.mapping, self.records)
            for record in records
        ]
        with multiprocessing.Pool(
            processes=2,
            initializer=packets._initialize_primer_worker,
            initargs=(self.primers, self.mapping, self.records, 4, 4, 42),
        ) as pool:
            _, classified = pool.map(packets._classify_primer_chunk, [(0, records)])[0]
        self.assertEqual(expected, [gene for _, gene in classified])

    def test_first_match_order_and_record_reuse(self):
        f_seq = "ACGT" + "A" * 20
        r_seq = "A" * 42 + "TTAA" + "G" * 10
        self.assertEqual("16s", packets.determine_gene_revised(self.primers, f_seq, r_seq, 4, 4, 42, self.mapping, self.records))
        self.assertEqual(self.records, packets.make_primer_records(self.primers))


if __name__ == "__main__":
    unittest.main()
