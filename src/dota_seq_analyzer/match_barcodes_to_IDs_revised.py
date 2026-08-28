#!/usr/bin/env python3
from collections import Counter
from typing import List, Dict, Tuple
import os
import argparse
import sys
from helper_functions import ensure_output_directories

def extract_all_b_with_ids(
    _16s_rev_fastq: str, arg_rev_fastq: str, unclassified_rev_fastq: str, 
    barcode_len: int = 20) -> Tuple[Counter, Dict]:
    """
    Match read IDs to their corresponding barcodes. Note that IDs are organized by their gene type (16s, ARG, or unclassified).
    Returns:
    - all_b_with_ids: a dictionary matching read IDs to barcodes
    - bcs_with_counts: a Counter tracking the size of each barcode group, 
      to be used to prioritize sequences in the upcoming clustering barcodes step
    """

    # preliminary steps
    print("Extracting barcodes from fastq...")
    gene_list_num_map = {"16s": 0, "arg": 1, "unclassified": 2}
    bcs_with_counts = Counter() # records the size of each barcode group (i.e. # of reads per barcode) 
                                # -> will be used later on to prioritize sequences, when clustering barcodes
    all_b_with_ids = {} # records which read IDs correspond to each barcode
                        # organizes these IDs by their gene type (16s, ARG, or unclassified)

    # match IDs to barcodes for reads of each type of gene (16s, ARG, unclassified)
    extract_b_with_ids_single(_16s_rev_fastq, "16s", gene_list_num_map, bcs_with_counts, all_b_with_ids, barcode_len)
    extract_b_with_ids_single(arg_rev_fastq, "arg", gene_list_num_map, bcs_with_counts, all_b_with_ids, barcode_len)
    extract_b_with_ids_single(unclassified_rev_fastq, "unclassified", gene_list_num_map, bcs_with_counts, all_b_with_ids, barcode_len)

    print("# of unique barcodes (before accounting for shift):", len(bcs_with_counts))
    return bcs_with_counts, all_b_with_ids

def extract_b_with_ids_single(
    fastq_filename: str, type_of_packet: str, gene_list_num_map: Dict[str, int], 
    bcs_with_counts: Counter, all_b_with_ids: Dict[str, List], barcode_len: int = 20):
    """Match IDs to their corresponding barcodes, for all reads of the given gene type (either 16s, ARG, or unclassified)"""
        
    # preliminary steps
    print(f"Extracting barcodes and matching {type_of_packet} IDs...")
    gene_list_num = gene_list_num_map[type_of_packet]

    # iterate through each read ID; for each ID:
    # - increment bcs_with_counts at the corresponding barcode
    # - add the read ID to its barcode-based position in the all_b_with_ids dictionary
    with open(fastq_filename, 'r') as f:
        line = f.readline()
        while line != "":
            # parse fastq text
            id = line.split(" ")[0].strip("@")
            bc = f.readline().strip("\n")[0:barcode_len]

            # increment bcs_with_counts
            bcs_with_counts[bc] += 1

            # append ID to corresponding position in all_b_with_IDs
            if bc not in all_b_with_ids:
                all_b_with_ids[bc] = [[], [], []]
            all_b_with_ids[bc][gene_list_num].append(id)
                
            for _ in range(3): # account for 4-line fastq format
                line = f.readline()

def create_clustered_b_with_ids(
    bcs_with_counts: Counter, all_b_with_ids: Dict, 
    max_shift: int = 1) -> Dict[str, List]:
    """
    Cluster similar barcodes, to account for sequencing errors.
    Return the modified dictionary matching now-clustered barcodes with all their corresponding read IDs.
    """

    print("Creating 'whitelist' of dominant barcodes, and hence forming barcode clusters.")
    print("Also adding IDs to their matching barcode clusters...")

    # preliminary steps
    sorted_bcs = [k for k, v in bcs_with_counts.most_common()] # sort barcodes by bargroup size
    clustered_b_with_ids = {} # same as all_b_with_ids, except that keys will only be dominant barcodes (not all extracted barcodes)

    # Iterate through all barcodes, where this barcode list has been sorted by bargroup size (i.e. from barcodes with the most # of reads, to least # of reads).
    # Reason: Dominant barcodes are identified as we go through this list - hence more dominant barcodes will be identified near the start of the list.
    #         Note that the dominant barcodes are the sequences representing their respective barcode clusters.
    #         We want the dominant barcodes to be the most common sequence in their respective clusters - combined with the fact that more dominant barcodes
    #         will be identified near the list's beginning, this is why we must sort the barcode list by size before iterating through it.
    for bc_s in sorted_bcs:
        in_dominant_bcs = False
        
        # iterate through all dominant barcodes identified so far
        for bc_d in clustered_b_with_ids: 
            
            # if current barcode can be clustered with one of the dominant barcodes, then merge their read IDs lists
            if check_barcodes_match_revised(bc_s, bc_d, max_shift):
                in_dominant_bcs = True
                # merge the 16s, arg, & unclassified ID lists for that barcode cluster
                for i in range(len(clustered_b_with_ids[bc_d])):
                    clustered_b_with_ids[bc_d][i].extend(all_b_with_ids[bc_s][i])
                bc_added = True
                break # similar barcode already in dominant barcodes list
                
        # if current barcode cannot be clustered with any dominant barcode, 
        # then add it to clustered_b_with_ids as a new dominant barcode
        if not in_dominant_bcs:
            clustered_b_with_ids[bc_s] = all_b_with_ids[bc_s]
            bc_added = True

        # barcode should have been added to clustered_b_with_ids already; if not, return error
        assert bc_added == True, f"Barcode {bc_s} not matched to any barcode in the dominant barcodes whitelist"

    return clustered_b_with_ids

def check_barcodes_match_revised(barcode1: str, barcode2: str, max_shift: int = 1) -> bool:
    """
    Checks if two barcodes match, given a max shift.
    Barcodes may be shifed by the max_shift (e.g. 1), but cannot have any mismatches, nor any internal shifts (gaps).
    Returns true if barcodes match; otherwise, returns false.
    """
    # iterate through each possible shift value
    for shift in range(-1*max_shift, 1):
        shift *= -1 # shift is in a negative range, and then converted back to positive
                    # reason: process non-zero shift condition (e.g. shift=1) before looking for exact match
                    # due to nature of this function's use
        len_to_check = len(barcode1) - shift  # assume barcodes 1 & 2 have same length
        if barcode1[0:len_to_check] == barcode2[(len(barcode1) - len_to_check):] \
        or barcode2[0:len_to_check] == barcode1[(len(barcode1) - len_to_check):]:
            return True
        
    return False

def format_list_to_str(list_: List) -> str:
    """Format a list as a string - namely, remove the square brackets [], and the single quotation marks"""
    return str(list_).replace("[", "").replace("]", "").replace("'", "")

def write_b_with_ids_to_file(clustered_b_with_ids: Dict, b_with_ids_filename: str):
    """Write each barcode with its corresponding IDs (organized by 16s vs ARG vs unclassified gene type), one per line, to a text file"""
    
    print("Writing barcodes with IDs to file...")

    with open(b_with_ids_filename, 'w') as f:
        for bc in clustered_b_with_ids:
            _16s_ids, arg_ids, unclassified_ids = tuple(clustered_b_with_ids[bc])
            line = f"{bc}: {format_list_to_str(_16s_ids)} | {format_list_to_str(arg_ids)} | {format_list_to_str(unclassified_ids)}\n"
            f.write(line)

def create_b_with_ids_file(
    _16s_rev_fastq: str, arg_rev_fastq: str, unclassified_rev_fastq: str,
    b_with_ids_filename: str, 
    barcode_shift: int = 1, barcode_len: int = 20):
    """
    Match all read IDs to their corresponding barcodes, and write to a file.
    This will be used in downstream per-cell analysis - specifically, in creating the barcode summary table.
    """
        
    bcs_with_counts, all_b_with_ids = extract_all_b_with_ids(_16s_rev_fastq, arg_rev_fastq, unclassified_rev_fastq, barcode_len) # match IDs to barcodes
    clustered_b_with_ids = create_clustered_b_with_ids(bcs_with_counts, all_b_with_ids, barcode_shift) # cluster similar barcodes, to account for sequencing errors
    write_b_with_ids_to_file(clustered_b_with_ids, b_with_ids_filename) # write data to a text file


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    # 2026-08-10: Expose barcode-matching FASTQ inputs with R2 terminology while retaining old aliases.
    # Reason: users provide sequencing files as R1/R2, and existing commands must remain compatible.
    parser.add_argument("--r2_16s_fastq", "--_16s_rev_fastq", dest="_16s_rev_fastq", type=str, required=True)
    parser.add_argument("--arg_r2_fastq", "--arg_rev_fastq", dest="arg_rev_fastq", type=str, required=True)
    parser.add_argument(
        "--unclassified_r2_fastq", "--unclassified_rev_fastq",
        dest="unclassified_rev_fastq", type=str, required=True,
    )
    parser.add_argument("--max_shift_barcode", type=int, default = 1)
    parser.add_argument("--barcode_len", type=int, default=20)
    # 2026-08-10: Route the barcode-to-read mapping to tmp by default.
    # Reason: this mapping is an internal input for downstream stages.
    parser.add_argument("--b_with_ids_filename", type=str, default = "tmp/b_with_ids.txt")

    args = parser.parse_args()

    # 2026-08-10: Materialize the temporary output directory before writing the mapping.
    # Reason: the default tmp path must work in a new result directory.
    ensure_output_directories(args.b_with_ids_filename)
    
    # make sure input file paths exist
    if not os.path.exists(args._16s_rev_fastq):
        print(f"❌ Error: input file not found: {args._16s_rev_fastq}")
        sys.exit(1)
    if not os.path.exists(args.arg_rev_fastq):
        print(f"❌ Error: input file not found: {args.arg_rev_fastq}")
        sys.exit(1)
    if not os.path.exists(args.unclassified_rev_fastq):
        print(f"❌ Error: input file not found: {args.unclassified_rev_fastq}")
        sys.exit(1)

    create_b_with_ids_file(
        args._16s_rev_fastq, args.arg_rev_fastq, args.unclassified_rev_fastq,
        args.b_with_ids_filename, args.max_shift_barcode, args.barcode_len)
    
if __name__ == "__main__":
    main()
