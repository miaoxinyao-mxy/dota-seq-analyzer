#!/usr/bin/env python3
"""Generate the deterministic synthetic FASTQ files used by the public example.

The 16S templates are short segments of one sequence in the bundled
GreenGenes2-based database. Read identifiers, cell barcodes, target sequence,
and quality strings are synthetic. No internal test01/test02 reads are used.
"""

from pathlib import Path
import random


OUTPUT_DIR = Path(__file__).resolve().parent
NUM_CELLS = 12
SIXTEEN_S_READS_PER_CELL = 8
TARGET_POSITIVE_CELLS = 4
TARGET_READS_PER_POSITIVE_CELL = 8
READ_LENGTH = 150

SIXTEEN_S_R1_TEMPLATE = (
    "CCTACGGGAGGCAGCAGTAGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAAGGAA"
    "GAAGTATTTCGGTATGTAAACTTCTATCAGCAAGGAAGAAAATGACGGTACTTGACTAAGAAGCCCCGGCTAAA"
    "TACGTGCCAGCAGCCGCGGTAATACGTATGGGGC"
)
SIXTEEN_S_R2_TEMPLATE = (
    "GGACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAGCGTCAGTAATCGTCCAGTAAGC"
    "CGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCCCTCC"
)
TARGET_FWD = "ACGTTGCAACGTTGCAACGT"
TARGET_REV = "TGCACGATTGCACGATTGCA"
R2_PREFIX_AFTER_BARCODE = "ACGTACGTACGTACGTACGTAC"


def barcodes():
    """Return stable barcodes that do not match at the configured one-base shift."""
    rng = random.Random(101)
    result = []
    while len(result) < NUM_CELLS:
        candidate = "".join(rng.choice("ACGT") for _ in range(20))
        if all(
            candidate != existing
            and candidate[:-1] != existing[1:]
            and existing[:-1] != candidate[1:]
            for existing in result
        ):
            result.append(candidate)
    return result


def padded(prefix, motif):
    repeats = (READ_LENGTH - len(prefix) + len(motif) - 1) // len(motif)
    return (prefix + motif * repeats)[:READ_LENGTH]


def write_record(handle, read_id, mate, sequence):
    handle.write(
        f"@{read_id} {mate}:N:0:SYNTHETIC\n"
        f"{sequence}\n+\n{'I' * len(sequence)}\n"
    )


def main():
    r1_path = OUTPUT_DIR / "synthetic_R1.fastq"
    r2_path = OUTPUT_DIR / "synthetic_R2.fastq"
    with r1_path.open("w", encoding="ascii") as r1, r2_path.open(
        "w", encoding="ascii"
    ) as r2:
        for cell_number, barcode in enumerate(barcodes(), start=1):
            for read_number in range(1, SIXTEEN_S_READS_PER_CELL + 1):
                read_id = f"synthetic_cell{cell_number:02d}_16s_{read_number:02d}"
                write_record(r1, read_id, 1, SIXTEEN_S_R1_TEMPLATE[:READ_LENGTH])
                write_record(
                    r2,
                    read_id,
                    2,
                    (barcode + R2_PREFIX_AFTER_BARCODE + SIXTEEN_S_R2_TEMPLATE)[:READ_LENGTH],
                )
            target_reads = (
                TARGET_READS_PER_POSITIVE_CELL
                if cell_number <= TARGET_POSITIVE_CELLS
                else 0
            )
            for read_number in range(1, target_reads + 1):
                read_id = f"synthetic_cell{cell_number:02d}_target_{read_number:02d}"
                write_record(r1, read_id, 1, padded(TARGET_FWD, "GATTACA"))
                write_record(
                    r2,
                    read_id,
                    2,
                    padded(barcode + R2_PREFIX_AFTER_BARCODE + TARGET_REV, "TGCA"),
                )


if __name__ == "__main__":
    main()
