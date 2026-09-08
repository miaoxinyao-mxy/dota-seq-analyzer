import io
import json
import multiprocessing
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "dota_seq_analyzer"))

import asv_typing_revised as asv  # noqa: E402
import create_ID_packets as packets  # noqa: E402
import extract_16s_reads as extract  # noqa: E402


FWD_16S = "CCTACGGGAGGCAGCAGT"
REV_16S = "GGACTACCAGGGTATCTAATCCTGTT"
FWD_ARG = "GCGATGTGCAGCACCAGTAA"
REV_ARG = "CACCGCTGCCGGTTTTATC"
INSERT = ("ACGTGCAATG" * 16)
R2_SUFFIX = "TGCATGCA" * 12
STAGGERS = {0: "", 4: "TGCA", 5: "AACGT", 8: "GATTACAG", 9: "ACGTACGTA"}


def write_fastq(path, records):
    with open(path, "w") as handle:
        for read_id, sequence in records:
            handle.write(f"@{read_id} 1:N:0:test\n{sequence}\n+\n{'I' * len(sequence)}\n")


def fastq_ids(path):
    with open(path) as handle:
        lines = handle.readlines()
    return [lines[i].split()[0].lstrip("@") for i in range(0, len(lines), 4)]


class StaggerAware16STests(unittest.TestCase):
    def test_valid_starts_trimming_and_asv_core(self):
        kraken_sequences = []
        asv_cores = []
        r2 = "A" * 150

        for expected_start, prefix in STAGGERS.items():
            r1 = prefix + FWD_16S + INSERT
            self.assertEqual(
                expected_start,
                extract.find_16s_r1_primer_start(r1, FWD_16S, max_mm=0),
            )

            outputs = [io.StringIO() for _ in range(5)]
            record = (
                expected_start, f"@read{expected_start}", f"@read{expected_start}",
                r1, r2, "I" * len(r1), "I" * len(r2),
            )
            extract._write_16s_record(
                record, expected_start, FWD_16S, REV_16S, 42,
                outputs[0], outputs[1], outputs[2], outputs[3], outputs[4],
            )
            kraken_sequences.append(outputs[2].getvalue().splitlines()[1])
            asv_cores.append(asv.extract_core(r1, r2, expected_start))

        self.assertEqual([INSERT] * len(STAGGERS), kraken_sequences)
        self.assertEqual(1, len(set(asv_cores)))

    def test_mismatch_tolerance_is_preserved_at_every_valid_start(self):
        for expected_start, prefix in STAGGERS.items():
            one_mismatch = list(FWD_16S)
            one_mismatch[6] = "A" if one_mismatch[6] != "A" else "C"
            one_mismatch = "".join(one_mismatch)
            self.assertEqual(
                expected_start,
                extract.find_16s_r1_primer_start(
                    prefix + one_mismatch + INSERT, FWD_16S, max_mm=1),
            )

            two_mismatches = list(FWD_16S)
            for position in (6, 13):
                two_mismatches[position] = (
                    "A" if two_mismatches[position] != "A" else "C")
            two_mismatches = "".join(two_mismatches)
            self.assertEqual(
                expected_start,
                extract.find_16s_r1_primer_start(
                    prefix + two_mismatches + INSERT, FWD_16S, max_mm=2),
            )

    def test_tied_start_uses_valid_start_order(self):
        primer = "AAAA"
        sequence = "A" * 20
        distances = [
            extract._primer_distance_at_start(sequence, primer, start, max_mm=0)
            for start in extract.VALID_16S_R1_STARTS
        ]
        self.assertEqual([0] * len(extract.VALID_16S_R1_STARTS), distances)
        self.assertEqual(
            0, extract.find_16s_r1_primer_start(sequence, primer, max_mm=0))

    def test_parallel_workers_preserve_stagger_and_target_assignments(self):
        r2_16s = "A" * 42 + REV_16S + R2_SUFFIX
        extraction_records = [
            (
                index, f"@read{start}", f"@read{start}",
                prefix + FWD_16S + INSERT, r2_16s,
                "I" * (len(prefix) + len(FWD_16S) + len(INSERT)),
                "I" * len(r2_16s),
            )
            for index, (start, prefix) in enumerate(STAGGERS.items())
        ]
        with multiprocessing.Pool(
            processes=2,
            initializer=extract._initialize_16s_worker,
            initargs=(FWD_16S, REV_16S, 4, 0, 42),
        ) as pool:
            _, starts = pool.map(extract._classify_16s_chunk, [extraction_records])[0]
        self.assertEqual(list(STAGGERS), starts)

        primers = [
            "Primer,F,R,Mode\n",
            f"16s,{FWD_16S},{REV_16S},single\n",
            f"ARG1,{FWD_ARG},{REV_ARG},single\n",
        ]
        mapping = packets.make_primers_to_genes_dict(primers)
        primer_records = packets.make_primer_records(primers)
        barcode = "ACGT" * 5
        arg_r2 = barcode + "A" * 22 + REV_ARG + R2_SUFFIX
        packet_records = [
            (0, "@read16s", "@read16s", "read16s", STAGGERS[8] + FWD_16S + INSERT, r2_16s, "I" * len(r2_16s), True),
            (1, "@read_arg", "@read_arg", "read_arg", FWD_ARG + INSERT, arg_r2, "I" * len(arg_r2), False),
        ]
        with multiprocessing.Pool(
            processes=2,
            initializer=packets._initialize_primer_worker,
            initargs=(primers, mapping, primer_records, 4, 0, 42),
        ) as pool:
            _, classified = pool.map(packets._classify_primer_chunk, [(0, packet_records)])[0]
        self.assertEqual(["16s", [1]], [gene for _, gene in classified])

    def test_manifest_is_packet_source_of_truth_and_arg_r2_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            primers = tmp / "primers.csv"
            primers.write_text(
                "Primer,F,R,Mode\n"
                f"16s,{FWD_16S},{REV_16S},single\n"
                f"ARG1,{FWD_ARG},{REV_ARG},single\n"
            )

            barcode = "ACGT" * 5
            r1_records = []
            r2_records = []
            for start, prefix in STAGGERS.items():
                read_id = f"read16s_{start}"
                r1_records.append((read_id, prefix + FWD_16S + INSERT))
                r2_records.append((read_id, barcode + "A" * 22 + REV_16S + R2_SUFFIX))
            arg_r1 = FWD_ARG + INSERT
            arg_r2 = barcode + "A" * 22 + REV_ARG + R2_SUFFIX
            r1_records.append(("read_arg", arg_r1))
            r2_records.append(("read_arg", arg_r2))
            r1_records.append(("read_none", "T" * 180))
            r2_records.append(("read_none", barcode + "T" * 160))

            r1 = tmp / "R1.fastq"
            r2 = tmp / "R2.fastq"
            write_fastq(r1, r1_records)
            write_fastq(r2, r2_records)

            only_r1 = tmp / "only_R1.fastq"
            only_r2 = tmp / "only_R2.fastq"
            kraken_r1 = tmp / "kraken_R1.fastq"
            kraken_r2 = tmp / "kraken_R2.fastq"
            manifest = tmp / "manifest.tsv"
            extract.create_16s_only_fastq(
                str(r1), str(r2), str(only_r1), str(only_r2),
                str(kraken_r1), str(kraken_r2), str(manifest), str(primers),
                max_shift=4, max_mm=0, primer_start_num=42, analysis_workers=1,
            )

            expected_16s_ids = [f"read16s_{start}" for start in STAGGERS]
            self.assertEqual(expected_16s_ids, fastq_ids(only_r1))
            manifest_lines = manifest.read_text().splitlines()[1:]
            self.assertEqual(
                expected_16s_ids,
                [line.split("\t")[1] for line in manifest_lines],
            )
            self.assertEqual(
                [str(start) for start in STAGGERS],
                [line.split("\t")[2] for line in manifest_lines],
            )
            with open(kraken_r2) as handle:
                kraken_r2_lines = handle.readlines()
            self.assertTrue(all(
                kraken_r2_lines[i + 1].strip() == R2_SUFFIX
                for i in range(0, len(kraken_r2_lines), 4)
            ))

            kraken_output = tmp / "kraken.output"
            kraken_output.write_text("".join(
                f"U\t{read_id}\t0\t0\tA:0\n" for read_id in expected_16s_ids
            ))
            kraken_report = tmp / "kraken.report"
            kraken_report.write_text("")
            packets_16s = tmp / "packets_16s"
            packets_arg = tmp / "packets_arg"
            packets_none = tmp / "packets_none"
            arg_r2_out = tmp / "arg_R2.fastq"
            none_r2_out = tmp / "none_R2.fastq"
            packets.generate_packets(
                str(r1), str(r2), str(primers), str(manifest),
                str(kraken_output), str(kraken_report), 20, 4, 0, 42,
                str(packets_16s), str(packets_arg), str(packets_none),
                str(arg_r2_out), str(none_r2_out), analysis_workers=1,
            )

            packet_ids = [json.loads(line)["ID"] for line in packets_16s.read_text().splitlines()]
            self.assertEqual(expected_16s_ids, packet_ids)
            self.assertEqual(["read_arg"], [json.loads(line)["ID"] for line in packets_arg.read_text().splitlines()])
            self.assertEqual(arg_r2, arg_r2_out.read_text().splitlines()[1])
            self.assertEqual(["read_none"], [json.loads(line)["ID"] for line in packets_none.read_text().splitlines()])


if __name__ == "__main__":
    unittest.main()
