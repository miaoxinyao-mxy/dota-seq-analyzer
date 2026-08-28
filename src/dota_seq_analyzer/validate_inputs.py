#!/usr/bin/env python3
import csv
import os
import sys
import argparse
try:
    from .helper_functions import get_arg_names, open_maybe_gzip
except ImportError:
    # 2026-08-28: Keep direct script execution compatible with package imports.
    # Reason: the public CLI imports this validator as a package module.
    from helper_functions import get_arg_names, open_maybe_gzip

def check_reference_fasta(reference_file: str):
    """Check that a reference file contains at least one structurally valid FASTA record."""
    saw_record = False
    saw_sequence = False
    with open_maybe_gzip(reference_file, "r") as reference:
        for line_number, raw_line in enumerate(reference, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if saw_record and not saw_sequence:
                    return f"<Reference FASTA> Record has no sequence before line {line_number}"
                if len(line) == 1:
                    return f"<Reference FASTA> Empty record ID on line {line_number}"
                saw_record = True
                saw_sequence = False
            else:
                if not saw_record:
                    return f"<Reference FASTA> Sequence appears before a header on line {line_number}"
                if any(char.isspace() for char in line):
                    return f"<Reference FASTA> Whitespace in sequence on line {line_number}"
                saw_sequence = True
        if not saw_record:
            return "<Reference FASTA> File contains no FASTA records"
        if not saw_sequence:
            return "<Reference FASTA> Last record has no sequence"
    return "Valid"

def check_input_files(fwd_fastq: str, rev_fastq: str, primers_file: str):
    """
    Check if input files are valid.
    Prints one error at a time, because some errors cannot be checked for unless a 
    pre-requisite component is functioning correctly (i.e. no error).
    Returns True if all input files pass all of their respective validation tests.
    Otherwise, returns False and prints the first observed error.
    """
    # For each input file, check if it's valid
    # If valid, continue to check next input file
    # If not valid, print out error message, and return from function

    # Primer File
    is_primer_file_valid = check_primers(primers_file)
    if is_primer_file_valid != "Valid":
        print(is_primer_file_valid)
        return False
    
    # FASTQ Files
    is_fastq_valid = check_paired_fastq(fwd_fastq, rev_fastq)
    if is_fastq_valid != "Valid":
        print(is_fastq_valid)
        return False

    # If all input files valid, return True
    return True

def check_paired_fastq(fwd_fastq: str, rev_fastq: str):
    """
    Preliminary check of validity of paired fastq input files.
    Specifically, check the following criteria:
    - For each pair of reads:
        - 1st line is ID: start with "@"; also, IDs for paired fwd & rev reads match
        - 2nd line is sequence: contains only A,T,C,G,N for this analysis
        - 3rd line is "+" sign
        - 4th line is quality: did not explicitly check this
    - Length of fwd & rev fastq files is the same (both have same # of reads)
    Returns either a string describing the first observed error (meaning it failed at least 1 test);
    or returns "Valid" if passed all tests.
    """
    with open_maybe_gzip(fwd_fastq, 'r') as fwd, \
    open_maybe_gzip(rev_fastq, 'r') as rev:

        id_f, id_r = fwd.readline().strip(), rev.readline().strip()

        i = 0
        while id_f != "" and id_r != "":
            # ID
            id_f = id_f.split(" ")[0]
            id_r = id_r.split(" ")[0] 
            if id_f != id_r: 
                return f"<FASTQ Files> Paired ID match error: IDs do not match on line {i*4 + 1}"
            if id_f[0] != "@": 
                return f"<FASTQ Files> Formatting error: ID does not start with '@' on line {i*4 + 1}"

            # Sequence
            f_seq = fwd.readline().strip()
            r_seq = rev.readline().strip()
            for char in f_seq:
                if char not in "ATCGN":
                    return f"<FASTQ Files> Formatting error: Forward sequence on line {i*4 + 2} contains letters other than A,T,C,G,N"
            for char in r_seq:
                if char not in "ATCGN":
                    return f"<FASTQ Files> Formatting error: Reverse sequence on line {i*4 + 2} contains letters other than A,T,C,G,N"
                
            # Line with "+" sign
            if not((fwd.readline().strip() == "+") and (rev.readline().strip() == "+")):
                return f"<FASTQ Files> Formatting error: Line separator should be '+' on line {i*4 + 3}"
            
            # Prepare for next read
            for _ in range(2):
                id_f, id_r = fwd.readline().strip(), rev.readline().strip()

            i += 1

        if id_f != id_r:
            return "<FASTQ Files> Paired file length error: R1 & R2 FASTQ files have different lengths (likely corresponds to different total # of reads)"
            # we expect f_line = r_line = "" at this point
            # if not true, then R1 & R2 files have different lengths - not good

        return "Valid" # input paired fastq files are valid, if all above conditions have been met
    
def check_primers(primers_file: str):
    """
    Check validity of input primer CSV file.
    Specifically, check the following criteria:
    - CSV format
    - 4 columns: Primer/Target, F, R, and Mode
    - 2nd & 3rd columns have headers "F" and "R", corresponding to forward & reverse primer sequences
    - No duplicate primer names
    - No duplicate paired primer sequences
    - Primer sequences only contain A,T,C,G
    - First primer is 16s, to optimize program run-time
    - Mode is blank or "ssr"; legacy single/family files remain accepted
    Returns either a string describing the first observed error (meaning it failed at least 1 test);
    or returns "Valid" if passed all tests. 
    """

    all_primer_names = []
    all_primer_seqs = []

    # 2026-08-11: Validate the new optional Mode column with CSV-aware parsing.
    # Reason: blank Mode values are valid and quoted target names must not break column parsing.
    with open(primers_file, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        rows = list(reader)
        if not rows or len(rows[0]) != 4:
            return "<Primer File> CSV error: expected 4 columns"

        header = rows[0]
        is_new_schema = header[0] in ["Primer", "Target"] and header[3] == "Mode"
        is_legacy_schema = header[0] == "Primer" and header[3] == "Sub-ARG"
        if not (is_new_schema or is_legacy_schema):
            return "<Primer File> Header error: expected Primer/Target,F,R,Mode"
        if not (header[1] == "F" and header[2] == "R"):
            return "<Primer File> Header name error: 2nd column should be labeled 'F', and 3rd column labeled 'R'"

        for i, line in enumerate(rows[1:], start=1):
            if len(line) != 4:
                return f"<Primer File> CSV error: row {i + 1} does not have 4 columns"

            primer = line[0].strip()
            seq = (line[1].strip(), line[2].strip())
            if primer in all_primer_names:
                return f"<Primer File> Duplicate primer name error: {primer} appears multiple times in primer CSV"
            elif seq in all_primer_seqs:
                return f"<Primer File> Duplicate primer sequence error: F & R sequence pair\n{seq}\nappears multiple times in primer CSV"
            else:
                all_primer_names.append(primer)
                all_primer_seqs.append(seq)

            for char in seq[0]:
                if char not in "ATCG":
                    return f"<Primer File> Sequence error: F primer {seq[0]} contains letters other than A,T,C,G"
            for char in seq[1]:
                if char not in "ATCG":
                    return f"<Primer File> Sequence error: R primer {seq[1]} contains letters other than A,T,C,G"

            mode = line[3].strip().casefold()
            if is_new_schema and mode not in ["", "ssr"]:
                return f"<Primer File> Mode error: {primer} should use blank or 'ssr'"
            if is_legacy_schema and mode not in ["single", "family"]:
                return f"<Primer File> Legacy Sub-ARG error: {primer} should use 'single' or 'family'"

            if i == 1 and primer.casefold() != "16s":
                return "<Primer File> Primer order error: First primer should be 16s, to optimize program run time"
        
    return "Valid" # input primer file is valid if it has met all the above conditions    


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    # 2026-08-10: Expose paired inputs as R1 and R2 to match sequencing conventions.
    # Reason: forward/reverse terminology is not the user-facing input model.
    parser.add_argument("--r1_fastq", type=str, required=True)
    parser.add_argument("--r2_fastq", type=str, required=True)
    parser.add_argument("--primers_file", type=str, required=True)

    args = parser.parse_args()
    
    # make sure input file paths exist
    if not os.path.exists(args.r1_fastq):
        print(f"❌ Error: input file not found: {args.r1_fastq}")
        sys.exit(1)
    if not os.path.exists(args.r2_fastq):
        print(f"❌ Error: input file not found: {args.r2_fastq}")
        sys.exit(1)
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        sys.exit(1)

    # exit from the program if not all input files are valid
    all_files_valid = check_input_files(args.r1_fastq, args.r2_fastq, args.primers_file)
    if not all_files_valid:
        sys.exit(1)

if __name__ == "__main__":
    main()
