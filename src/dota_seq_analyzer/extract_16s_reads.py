#!/usr/bin/env python3
from typing import Optional, Tuple
import sys
# 2026-08-28: Exit non-zero when a required stage input is missing.
# Reason: the public CLI uses subprocess check=True to stop failed stages.
from helper_functions import open_maybe_gzip, ensure_output_directories
import argparse
import multiprocessing
import os

EXTRACTION_CHUNK_SIZE = 2048
_16S_WORKER_CONFIG = None

def determine_16s_primers(primers_filename: str) -> Tuple[str, str]:
    """Obtain the R1 & R2 primers for the 16s gene"""
    
    with open(primers_filename, 'r') as primers_file:
        for line in primers_file:
            line = line.strip("\n").split(",")
            if line[0] == "16s":
                return line[1], line[2]
    return None

def bounded_edit_distance(pattern: str, text: str, max_edits: int) -> int:
    """
    Calculate Levenshtein edit distance, stopping early once the distance
    must exceed `max_edits`.

    Optimized version:
    1. Keeps the same edit-distance logic.
    2. Avoids repeated min() calls.
    3. Avoids repeated list.append() calls by preallocating the current row.
    """
    plen = len(pattern)
    tlen = len(text)

    # If the length difference alone already exceeds the allowed edits,
    # the edit distance cannot pass.
    if abs(plen - tlen) > max_edits:
        return max_edits + 1

    previous = list(range(tlen + 1))

    for i, pattern_base in enumerate(pattern, start=1):
        current = [0] * (tlen + 1)
        current[0] = i
        row_min = i

        for j, text_base in enumerate(text, start=1):
            substitution_cost = 0 if pattern_base == text_base else 1

            deletion = previous[j] + 1
            insertion = current[j - 1] + 1
            substitution = previous[j - 1] + substitution_cost

            best = deletion
            if insertion < best:
                best = insertion
            if substitution < best:
                best = substitution

            current[j] = best

            if best < row_min:
                row_min = best

        if row_min > max_edits:
            return max_edits + 1

        previous = current

    return previous[tlen]

def check_primer_match_seq(
    seq: str, primer: str,
    max_shift: int, max_mm: int,
    centre: int = 0,
    max_indel: Optional[int] = None
) -> bool:
    """
    Check whether primer matches seq near the expected centre position.

    Find the best gapped match of `primer` in `seq` with:
      - primer aligned within `max_shift` positions of `centre`,
      - up to `max_mm` total edits (mismatches + insertions + deletions),
      - up to `max_indel` bases of length difference in the local window.

    Optimization:
    - if primer_fragment and seq_window have the same length,
      use mismatch counting instead of edit-distance DP.

    Returns True if there exists an alignment (where shift <= max shift), such that edit distance <= max_mm
    This indicates that the sequence does indeed contain the primer.
    If no alignment has edit distance <= max_mm, returns False.
    """
    seq = seq.strip().upper()
    primer = primer.strip().upper()

    slen = len(seq)
    plen = len(primer)

    if slen == 0 or plen == 0:
        return False

    if max_indel is None:
        max_indel = min(max_mm, 2)

    for shift in range(centre - max_shift, centre + max_shift + 1):

        if shift >= centre:
            seq_index_start = shift
            primer_index_start = 0
        else:
            seq_index_start = centre
            primer_index_start = abs(shift - centre)

        if seq_index_start < 0 or seq_index_start >= slen or primer_index_start >= plen:
            continue

        primer_fragment = primer[primer_index_start:]
        fragment_len = len(primer_fragment)

        min_window_len = max(1, fragment_len - max_indel)
        max_window_len = min(slen - seq_index_start, fragment_len + max_indel)

        if min_window_len > max_window_len:
            continue

        for window_len in range(min_window_len, max_window_len + 1):
            seq_window = seq[seq_index_start:seq_index_start + window_len]

            # Safety check: if Python slicing produced a shorter window,
            # skip it instead of comparing an unintended shorter sequence.
            if len(seq_window) != window_len:
                continue

            # Same length: Levenshtein distance equals mismatch count.
            # This avoids expensive DP but preserves the result.
            if fragment_len == window_len:
                mm = 0
                for a, b in zip(primer_fragment, seq_window):
                    if a != b:
                        mm += 1
                        if mm > max_mm:
                            break

                if mm <= max_mm:
                    return True

                continue

            # Different length: keep indel-aware edit distance.
            if bounded_edit_distance(primer_fragment, seq_window, max_mm) <= max_mm:
                return True

    return False


def _initialize_16s_worker(fwd_primer, rev_primer, max_shift, max_mm, primer_start_num):
    """Initialize invariant 16S matching data once per worker process."""
    global _16S_WORKER_CONFIG
    _16S_WORKER_CONFIG = (fwd_primer, rev_primer, max_shift, max_mm, primer_start_num)


def _classify_16s_chunk(records):
    """Classify one ordered chunk without writing output files."""
    fwd_primer, rev_primer, max_shift, max_mm, primer_start_num = _16S_WORKER_CONFIG
    decisions = [
        check_primer_match_seq(f_seq, fwd_primer, max_shift, max_mm)
        and check_primer_match_seq(r_seq, rev_primer, max_shift, max_mm, primer_start_num)
        for _, _, f_seq, r_seq, _, _ in records
    ]
    return records, decisions


def _read_16s_chunk(r1, r2, chunk_size):
    """Read contiguous paired FASTQ records while preserving input order."""
    records = []
    for _ in range(chunk_size):
        f_line = r1.readline()
        r_line = r2.readline()
        if f_line == "" or r_line == "":
            break

        id_f = f_line.strip()
        id_r = r_line.strip()
        assert id_f.split(" ")[0] == id_r.split(" ")[0]  # IDs should match
        f_seq = r1.readline().strip()
        r_seq = r2.readline().strip()
        for _ in range(2):
            f_quality = r1.readline().strip()
            r_quality = r2.readline().strip()
        records.append((id_f, id_r, f_seq, r_seq, f_quality, r_quality))
    return records


def _write_16s_record(record, is_16s, fwd_primer, rev_primer, primer_start_num,
                      only_16s_r1, only_16s_r2, kraken_only_16s_r1, kraken_only_16s_r2):
    """Write one already-classified record using the existing output format."""
    if not is_16s:
        return
    id_f, id_r, f_seq, r_seq, f_quality, r_quality = record
    only_16s_r1.write(f"{id_f}\n{f_seq}\n+\n{f_quality}\n")
    only_16s_r2.write(f"{id_r}\n{r_seq}\n+\n{r_quality}\n")
    trimmed_f_seq = f_seq[len(fwd_primer):]
    trimmed_r_seq = r_seq[(primer_start_num + len(rev_primer)):]
    trimmed_f_quality = f_quality[len(fwd_primer):]
    trimmed_r_quality = r_quality[(primer_start_num + len(rev_primer)):]
    kraken_only_16s_r1.write(f"{id_f}\n{trimmed_f_seq}\n+\n{trimmed_f_quality}\n")
    kraken_only_16s_r2.write(f"{id_r}\n{trimmed_r_seq}\n+\n{trimmed_r_quality}\n")


def _create_16s_only_fastq_parallel(
    r1_filename, r2_filename, only_16s_r1_filename, only_16s_r2_filename,
    kraken_only_16s_r1_filename, kraken_only_16s_r2_filename,
    fwd_primer, rev_primer, max_shift, max_mm, primer_start_num, analysis_workers
):
    """Classify read chunks in workers; keep all output writing in the parent."""
    with open_maybe_gzip(r1_filename, "rt") as r1, \
         open_maybe_gzip(r2_filename, "rt") as r2, \
         open(only_16s_r1_filename, "w") as only_16s_r1, \
         open(only_16s_r2_filename, "w") as only_16s_r2, \
         open(kraken_only_16s_r1_filename, "w") as kraken_only_16s_r1, \
         open(kraken_only_16s_r2_filename, "w") as kraken_only_16s_r2:
        def chunks():
            while True:
                records = _read_16s_chunk(r1, r2, EXTRACTION_CHUNK_SIZE)
                if not records:
                    return
                yield records

        with multiprocessing.Pool(
            processes=analysis_workers,
            initializer=_initialize_16s_worker,
            initargs=(fwd_primer, rev_primer, max_shift, max_mm, primer_start_num),
        ) as pool:
            processed = 0
            for records, decisions in pool.imap(_classify_16s_chunk, chunks()):
                for record, is_16s in zip(records, decisions):
                    _write_16s_record(
                        record, is_16s, fwd_primer, rev_primer, primer_start_num,
                        only_16s_r1, only_16s_r2, kraken_only_16s_r1, kraken_only_16s_r2,
                    )
                processed += len(records)
                if processed % 1000 == 0:
                    print(f"Processed {processed//1000},000 reads")


def create_16s_only_fastq(
    r1_filename: str, r2_filename: str,
    only_16s_r1_filename: str, only_16s_r2_filename: str,
    kraken_only_16s_r1_filename: str, kraken_only_16s_r2_filename: str,
    primers_filename: str,
    max_shift: int, max_mm: int,
    primer_start_num: int = 42,
    analysis_workers: int = 1,
):
    """
    Extract reads whose fwd & rev sequences both contain 16s primers.
    Write these 16s reads to fwd & rev fastq files, 
    with primers, barcodes, & overlap regions trimmed off - 
    this will be used for Kraken 16s taxonomic classification.
    Also write these 16s reads to other fwd & rev fastq files, without trimming - 
    this will be used for all other downstream processing, 
    where knowing the primers & barcodes is important.
    """

    # determine forward & reverse 16s primers, by parsing through primers file
    fwd_primer, rev_primer = determine_16s_primers(primers_filename)

    if analysis_workers < 1:
        raise ValueError("analysis_workers must be at least 1")
    # 2026-09-04: Parallelize only independent 16S primer decisions in chunks.
    # Reason: reuse the existing analysis worker setting without changing matching semantics.
    if analysis_workers > 1:
        return _create_16s_only_fastq_parallel(
            r1_filename, r2_filename, only_16s_r1_filename, only_16s_r2_filename,
            kraken_only_16s_r1_filename, kraken_only_16s_r2_filename,
            fwd_primer, rev_primer, max_shift, max_mm, primer_start_num, analysis_workers,
        )

    # open files, both original (to read from) and new (to write to)
    with open_maybe_gzip(r1_filename, "rt") as r1, \
    open_maybe_gzip(r2_filename, "rt") as r2, \
    open(only_16s_r1_filename, "w") as only_16s_r1, \
    open(only_16s_r2_filename, "w") as only_16s_r2, \
    open(kraken_only_16s_r1_filename, "w") as kraken_only_16s_r1, \
    open(kraken_only_16s_r2_filename, "w") as kraken_only_16s_r2:

        # read in the read_ID line
        i = 0
        f_line = r1.readline()
        r_line = r2.readline()

        while (f_line != "") and (r_line != ""):

            id_f = f_line.strip()
            id_r = r_line.strip()
            assert id_f.split(" ")[0] == id_r.split(" ")[0]  # IDs should match

            f_seq = r1.readline().strip()
            r_seq = r2.readline().strip()

            for _ in range(2):
                f_quality = r1.readline().strip()
                r_quality = r2.readline().strip()

            # note that 16s primers must be present in both fwd and rev reads
            if check_primer_match_seq(f_seq, fwd_primer, max_shift, max_mm) \
            and check_primer_match_seq(r_seq, rev_primer, max_shift, max_mm, primer_start_num):                    

                only_16s_r1.write(f"{id_f}\n{f_seq}\n+\n{f_quality}\n")
                only_16s_r2.write(f"{id_r}\n{r_seq}\n+\n{r_quality}\n")

                # trim sequences to prepare for Kraken 16s taxonomic classification
                # remove primers from both fwd & rev sequences
                # remove barcode & overlap from rev sequence
                # then add trimmed read to other 16s-only fastq files, to be used for Kraken input
                trimmed_f_seq = f_seq[len(fwd_primer):]
                trimmed_r_seq = r_seq[(primer_start_num + len(rev_primer)):]

            # 2026-08-10: Trim quality strings at the same coordinates as sequences.
                # Reason: FASTQ sequence and quality lines must have identical lengths.
                trimmed_f_quality = f_quality[len(fwd_primer):]
                trimmed_r_quality = r_quality[(primer_start_num + len(rev_primer)):]

                kraken_only_16s_r1.write(f"{id_f}\n{trimmed_f_seq}\n+\n{trimmed_f_quality}\n")
                kraken_only_16s_r2.write(f"{id_r}\n{trimmed_r_seq}\n+\n{trimmed_r_quality}\n")

            # prepare for next read
            f_line = r1.readline()
            r_line = r2.readline()

            i += 1
            if i % 1000 == 0:
                print(f"Processed {i//1000},000 reads")

            
def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    # 2026-08-10: Expose all paired-read CLI inputs as R1 and R2.
    # Reason: users provide sequencing reads as R1/R2 files.
    parser.add_argument("--r1_fastq", type=str, required=True)
    parser.add_argument("--r2_fastq", type=str, required=True)
    parser.add_argument("--primers_filename", type=str, required=True)
    # 2026-08-10: Route generated read files to the temporary-output directory by default.
    # Reason: intermediate FASTQ files should not clutter the result root.
    parser.add_argument("--r1_only_16s_fastq", type=str, default="tmp/only_16s_R1.fastq")
    parser.add_argument("--r2_only_16s_fastq", type=str, default="tmp/only_16s_R2.fastq")
    parser.add_argument("--kraken_r1_only_16s_fastq", type=str, default="tmp/kraken_R1.fastq")
    parser.add_argument("--kraken_r2_only_16s_fastq", type=str, default="tmp/kraken_R2.fastq")
    parser.add_argument("--max_shift_primer", type=int, default=4)
    parser.add_argument("--max_mm_primer", type=int, default=4)
    parser.add_argument("--primer_start_num", type=int, default=42)
    parser.add_argument("-@", "--threads", dest="analysis_workers", type=int, default=1, metavar="INT")

    args = parser.parse_args()

    # 2026-08-10: Materialize the organized output directories before writing files.
    # Reason: default tmp paths must work in a new result directory.
    ensure_output_directories(
        args.r1_only_16s_fastq, args.r2_only_16s_fastq,
        args.kraken_r1_only_16s_fastq, args.kraken_r2_only_16s_fastq)
    
    # make sure input file paths exist
    if not os.path.exists(args.r1_fastq):
        print(f"❌ Error: input file not found: {args.r1_fastq}")
        sys.exit(1)
    if not os.path.exists(args.r2_fastq):
        print(f"❌ Error: input file not found: {args.r2_fastq}")
        sys.exit(1)
    if not os.path.exists(args.primers_filename):
        print(f"❌ Error: input file not found: {args.primers_filename}")
        sys.exit(1)

    create_16s_only_fastq(
        args.r1_fastq, args.r2_fastq,
        args.r1_only_16s_fastq, args.r2_only_16s_fastq,
        args.kraken_r1_only_16s_fastq, args.kraken_r2_only_16s_fastq,
        args.primers_filename, args.max_shift_primer,
        args.max_mm_primer, args.primer_start_num, args.analysis_workers)

if __name__ == "__main__":
    main()
