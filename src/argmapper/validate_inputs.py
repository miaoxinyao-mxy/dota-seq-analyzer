import os
import sys
import argparse
from helper_functions import get_arg_names, open_maybe_gzip

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
            return "<FASTQ Files> Paired file length error: forward & reverse FASTQ files have different lengths (likely corresponds to different total # of reads)"
            # we expect f_line = r_line = "" at this point
            # if not true, then fwd & rev files have different lengths - not good

        return "Valid" # input paired fastq files are valid, if all above conditions have been met
    
def check_primers(primers_file: str):
    """
    Check validity of input primer CSV file.
    Specifically, check the following criteria:
    - CSV format
    - 4 columns
    - 2nd & 3rd columns have headers "F" and "R", corresponding to forward & reverse primer sequences
    - No duplicate primer names
    - No duplicate paired primer sequences
    - Primer sequences only contain A,T,C,G
    - First primer is 16s, to optimize program run-time
    - Sub-ARG classification is either "single" or "family"
    Returns either a string describing the first observed error (meaning it failed at least 1 test);
    or returns "Valid" if passed all tests. 
    """

    all_primer_names = []
    all_primer_seqs = []

    with open(primers_file, 'r') as f:
        i = 0
        line = f.readline()
        while line != "":
            line = line.split(",")
            if len(line) != 4:
                return "<Primer File> CSV error: either not in CSV format, or does not have 4 columns"
            
            if i == 0: # header line
                if not(line[1] == "F" and line[2] == "R"):
                    return "<Primer File> Header name error: 2nd column should be labeled 'F', and 3rd column labeled 'R'"
            
            elif i >= 1: # any primer line, including first primer

                primer = line[0]
                seq = (line[1], line[2])
                if primer in all_primer_names:
                    return f"<Primer File> Duplicate primer name error: {primer} appears multiple times in primer CSV"
                elif seq in all_primer_seqs:
                    return f"<Primer File> Duplicate primer sequence error: fwd & rev sequence pair\n{seq}\nappears multiple times in primer CSV"
                else:
                    all_primer_names.append(primer)
                    all_primer_seqs.append(seq)

                for char in line[1]:
                    if char not in "ATCG":
                        return f"<Primer File> Sequence error: Forward primer {line[1]} contains letters other than A,T,C,G"
                for char in line[2]:
                    if char not in "ATCG":
                        return f"<Primer File> Sequence error: Reverse primer {line[2]} contains letters other than A,T,C,G"
                    
                if line[3].strip().lower() not in ["single", "family"]:
                    return f"<Primer File> Sub-ARG classification error: sub-ARG status for {primer} should be either 'single' or 'family"

            if i == 1: # first primer line
                if not line[0] == "16s":
                    return "<Primer File> Primer order error: First primer should be 16s, to optimize program run time"
            
            line = f.readline()
            i += 1
        
    return "Valid" # input primer file is valid if it has met all the above conditions    


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    parser.add_argument("--fwd_fastq", type=str, required=True)
    parser.add_argument("--rev_fastq", type=str, required=True)
    parser.add_argument("--primers_file", type=str, required=True)

    args = parser.parse_args()
    
    # make sure input file paths exist
    if not os.path.exists(args.fwd_fastq):
        print(f"❌ Error: input file not found: {args.fwd_fastq}")
        return
    if not os.path.exists(args.rev_fastq):
        print(f"❌ Error: input file not found: {args.rev_fastq}")
        return
    if not os.path.exists(args.primers_file):
        print(f"❌ Error: input file not found: {args.primers_file}")
        return

    all_files_valid = check_input_files(args.fwd_fastq, args.rev_fastq, args.primers_file)
    if not all_files_valid:
        sys.exit(1)

if __name__ == "__main__":
    main()
