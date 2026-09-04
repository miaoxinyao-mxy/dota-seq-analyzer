#!/usr/bin/env python3
"""Filter paired FASTQ reads before DoTA-Seq analysis."""

import argparse
import os
from typing import Dict

try:
    from .helper_functions import ensure_output_directories, open_maybe_gzip
except ImportError:
    from helper_functions import ensure_output_directories, open_maybe_gzip


MIN_READ_LENGTH = 130
MIN_MEAN_PHRED = 25
BARCODE_LENGTH = 20
MIN_BARCODE_Q25 = 15


def phred_scores(quality: str):
    """Convert Phred+33 FASTQ quality characters to integer scores."""
    return [ord(char) - 33 for char in quality]


def mean_phred(quality: str) -> float:
    """Return the mean Phred score for one FASTQ quality string."""
    if not quality:
        return 0.0
    return sum(phred_scores(quality)) / len(quality)


def passes_read_qc(
    r1_sequence: str,
    r1_quality: str,
    r2_sequence: str,
    r2_quality: str,
) -> tuple[bool, str]:
    """Apply the project read-length, read-quality, and barcode-quality rules."""
    if (
        len(r1_sequence) < MIN_READ_LENGTH
        or len(r2_sequence) < MIN_READ_LENGTH
        or len(r1_sequence) != len(r1_quality)
        or len(r2_sequence) != len(r2_quality)
    ):
        return False, "read_length_or_quality_length"

    if mean_phred(r1_quality) < MIN_MEAN_PHRED or mean_phred(r2_quality) < MIN_MEAN_PHRED:
        return False, "mean_phred"

    barcode_quality = phred_scores(r2_quality[:BARCODE_LENGTH])
    if len(r2_sequence) < BARCODE_LENGTH or len(barcode_quality) < BARCODE_LENGTH:
        return False, "barcode_length"
    if sum(score >= MIN_MEAN_PHRED for score in barcode_quality) < MIN_BARCODE_Q25:
        return False, "barcode_q25"

    return True, "passed"


def filter_paired_fastq(r1_input: str, r2_input: str, r1_output: str, r2_output: str) -> Dict[str, int]:
    """Filter paired FASTQ records while preserving order and original record text."""
    ensure_output_directories(r1_output, r2_output)
    counts = {
        "total_pairs": 0,
        "passed_pairs": 0,
        "read_length_or_quality_length": 0,
        "mean_phred": 0,
        "barcode_length": 0,
        "barcode_q25": 0,
    }
    with (
        open_maybe_gzip(r1_input, "r") as r1,
        open_maybe_gzip(r2_input, "r") as r2,
        open(r1_output, "w") as out1,
        open(r2_output, "w") as out2,
    ):
        while True:
            r1_header = r1.readline()
            r2_header = r2.readline()
            if r1_header == "" and r2_header == "":
                break
            if r1_header == "" or r2_header == "":
                raise ValueError("R1 and R2 FASTQ files contain different numbers of records")

            r1_sequence = r1.readline().rstrip("\n")
            r2_sequence = r2.readline().rstrip("\n")
            r1_plus = r1.readline().rstrip("\n")
            r2_plus = r2.readline().rstrip("\n")
            r1_quality = r1.readline().rstrip("\n")
            r2_quality = r2.readline().rstrip("\n")
            if not r1_plus.startswith("+") or not r2_plus.startswith("+"):
                raise ValueError("FASTQ record separator must begin with '+'")
            if r1_header.split()[0] != r2_header.split()[0]:
                raise ValueError(f"R1/R2 read ID mismatch: {r1_header.strip()} vs {r2_header.strip()}")

            counts["total_pairs"] += 1
            passed, reason = passes_read_qc(
                r1_sequence, r1_quality, r2_sequence, r2_quality
            )
            if passed:
                counts["passed_pairs"] += 1
                out1.write("".join((r1_header, r1_sequence + "\n", r1_plus + "\n", r1_quality + "\n")))
                out2.write("".join((r2_header, r2_sequence + "\n", r2_plus + "\n", r2_quality + "\n")))
            else:
                counts[reason] += 1

    counts["filtered_pairs"] = counts["total_pairs"] - counts["passed_pairs"]
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter paired reads before DoTA-Seq analysis.")
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--output-r1", required=True)
    parser.add_argument("--output-r2", required=True)
    args = parser.parse_args()
    counts = filter_paired_fastq(args.r1, args.r2, args.output_r1, args.output_r2)
    print("[filter-reads] " + " ".join(f"{key}={value}" for key, value in counts.items()), flush=True)


if __name__ == "__main__":
    main()
