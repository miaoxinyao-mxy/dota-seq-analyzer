from collections import Counter
from typing import List, Dict, Tuple
import os
import argparse

def extract_all_b_with_ids(
    _16s_rev_fastq: str, arg_rev_fastq: str, unclassified_rev_fastq: str, 
    barcode_len: int = 20) -> Tuple[Counter, Dict]:

    print("Extracting barcodes from fastq...")
    gene_list_num_map = {"16s": 0, "arg": 1, "unclassified": 2}
    bcs_with_counts = Counter()
    all_b_with_ids = {}
    
    extract_b_with_ids_single(_16s_rev_fastq, "16s", gene_list_num_map, bcs_with_counts, all_b_with_ids, barcode_len)
    extract_b_with_ids_single(arg_rev_fastq, "arg", gene_list_num_map, bcs_with_counts, all_b_with_ids, barcode_len)
    extract_b_with_ids_single(unclassified_rev_fastq, "unclassified", gene_list_num_map, bcs_with_counts, all_b_with_ids, barcode_len)

    print("# of unique barcodes (before accounting for shift):", len(bcs_with_counts))
    return bcs_with_counts, all_b_with_ids

def extract_b_with_ids_single(
    fastq_filename: str, type_of_packet: str, gene_list_num_map: Dict[str, int], 
    bcs_with_counts: Counter, all_b_with_ids: Dict[str, List], barcode_len: int = 20):

    print(f"Extracting barcodes and matching {type_of_packet} IDs...")
    gene_list_num = gene_list_num_map[type_of_packet]

    with open(fastq_filename, 'r') as f:
        line = f.readline()
        while line != "":
            id = line.split(" ")[0].strip("@")
            bc = f.readline().strip("\n")[0:barcode_len]
            bcs_with_counts[bc] += 1

            if bc not in all_b_with_ids:
                all_b_with_ids[bc] = [[], [], []]
            all_b_with_ids[bc][gene_list_num].append(id)
                
            for _ in range(3): # account for 4-line fastq format
                line = f.readline()

def create_clustered_b_with_ids(
    bcs_with_counts: Counter, all_b_with_ids: Dict, 
    max_shift: int = 1) -> Dict[str, List]:

    print("Creating whitelist of dominant barcodes, and hence forming barcode clusters.")
    print("Also adding IDs to their matching barcode clusters...")

    sorted_bcs = [k for k, v in bcs_with_counts.most_common()]
    clustered_b_with_ids = {} # keys will only be dominant barcodes, not all extracted barcodes

    for bc_s in sorted_bcs:
        in_dominant_bcs = False
        for bc_d in clustered_b_with_ids: # iterate through all dominant barcodes identified so far
            if check_barcodes_match_revised(bc_s, bc_d, max_shift):
                in_dominant_bcs = True
                # merge the 16s, arg, & unclassified ID lists for that barcode cluster
                for i in range(len(clustered_b_with_ids[bc_d])):
                    clustered_b_with_ids[bc_d][i].extend(all_b_with_ids[bc_s][i])
                bc_added = True
                break # similar barcode already in dominant barcodes list
        if not in_dominant_bcs:
            clustered_b_with_ids[bc_s] = all_b_with_ids[bc_s]
            bc_added = True

        assert bc_added == True, f"Barcode {bc_s} not matched to any barcode in the dominant barcodes whitelist"

    return clustered_b_with_ids

def check_barcodes_match_revised(barcode1: str, barcode2: str, max_shift: int = 1) -> bool:
    """
    Checks if two barcodes match, given a max shift.
    Barcodes may be shifed by the max_shift (e.g. 1), but cannot have any mismatches, nor any internal shifts (gaps)
    Returns true if barcodes match
    """
    for shift in range(-1*max_shift, 1):
        shift *= -1 # shift is in a negative range, and then convered back to positive
                    # reason: process non-zero shift condition (e.g. shift=1) before looking for exact match
                    # due to nature of this function's use
        len_to_check = len(barcode1) - shift  # assume barcodes 1 & 2 have same length
        if barcode1[0:len_to_check] == barcode2[(len(barcode1) - len_to_check):] \
        or barcode2[0:len_to_check] == barcode1[(len(barcode1) - len_to_check):]:
            return True
        
    return False

def format_list_to_str(list_: List) -> str:
    return str(list_).replace("[", "").replace("]", "").replace("'", "")

def write_b_with_ids_to_file(clustered_b_with_ids: Dict, b_with_ids_filename: str):
    print("Writing barcodes with IDs to file...")

    with open(b_with_ids_filename, 'w') as f:
        for bc in clustered_b_with_ids:
            _16s_ids, arg_ids, unclassified_ids = tuple(clustered_b_with_ids[bc])
            line = f"{bc}: {format_list_to_str(_16s_ids)} | {format_list_to_str(arg_ids)} | {format_list_to_str(unclassified_ids)}\n"
            f.write(line)

def create_b_with_ids_file(
    _16s_rev_fastq: str, arg_rev_fastq: str, unclassified_rev_fastq: str,
    b_with_ids_filename: str, 
    _16s_packets: str, arg_packets: str, unclassified_packets: str, 
    barcode_shift: int = 1, barcode_len: int = 20):
    
    bcs_with_counts, all_b_with_ids = extract_all_b_with_ids(_16s_rev_fastq, arg_rev_fastq, unclassified_rev_fastq, barcode_len)
    clustered_b_with_ids = create_clustered_b_with_ids(bcs_with_counts, all_b_with_ids, barcode_shift)
    write_b_with_ids_to_file(clustered_b_with_ids, b_with_ids_filename)


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    parser.add_argument("--_16s_rev_fastq", type=str, required=True)
    parser.add_argument("--arg_rev_fastq", type=str, required=True)
    parser.add_argument("--unclassified_rev_fastq", type=str, required=True)
    parser.add_argument("--_16s_packet_filename", type=str, required=True)
    parser.add_argument("--arg_packet_filename", type=str, required=True)
    parser.add_argument("--unclassified_packet_filename", type=str, required=True)
    parser.add_argument("--max_shift_barcode", type=int, default = 1)
    parser.add_argument("--barcode_len", type=int, default=20)
    parser.add_argument("--b_with_ids_filename", type=str, default = "b_with_ids")

    args = parser.parse_args()
    
    # make sure input file paths exist
    if not os.path.exists(args._16s_rev_fastq):
        print(f"❌ Error: input file not found: {args._16s_rev_fastq}")
        return
    if not os.path.exists(args.arg_rev_fastq):
        print(f"❌ Error: input file not found: {args.arg_rev_fastq}")
        return
    if not os.path.exists(args.unclassified_rev_fastq):
        print(f"❌ Error: input file not found: {args.unclassified_rev_fastq}")
        return
    if not os.path.exists(args._16s_packet_filename):
        print(f"❌ Error: input file not found: {args._16s_packet_filename}")
        return
    if not os.path.exists(args.arg_packet_filename):
        print(f"❌ Error: input file not found: {args.arg_packet_filename}")
        return
    if not os.path.exists(args.unclassified_packet_filename):
        print(f"❌ Error: input file not found: {args.unclassified_packet_filename}")
        return

    create_b_with_ids_file(
        args._16s_rev_fastq, args.arg_rev_fastq, args.unclassified_rev_fastq,
        args.b_with_ids_filename,
        args._16s_packet_filename, args.arg_packet_filename, args.unclassified_packet_filename, 
        args.max_shift_barcode, args.barcode_len)
    
if __name__ == "__main__":
    main()