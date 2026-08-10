import json
import argparse
import os
from typing import List, Dict, Tuple
import pandas as pd
from mle_revised import parse_and_analyze_perfect_corrected_revised
from filter_barcodes import filter_barcodes_in_df
from helper_functions import get_arg_names, ensure_output_directories
    
def find_arg_data(packet_list: List[Dict], num_arg_genes: int):

    # 2026-08-10: Allocate counts from the primer panel, not a fixed 23-gene assumption.
    # Reason: the output vector must match the columns generated from the supplied primers file.
    gene_counts = [0] * num_arg_genes

    for packet in packet_list:
        gene = packet["gene"]
        # assume gene is a 1D matrix with 23 numbers
        # all numbers should be 0, except for one 1
        for i in range(len(gene)):
            if gene[i] == 1:
                gene_counts[i] += 1
                break # no other 1s to add, so break

    return gene_counts
    
# 2026-08-10: Added defaults to parameters that followed defaulted arguments. Reason: Python otherwise rejects this function before the pipeline can run.
def write_barcode_summary_to_tsv(b_with_ids_filename: str, 
    _16s_packet_filename: str, arg_packet_filename: str, 
    unfiltered_tsv_filename: str, tsv_filename: str, primers_filename: str,
    min_16s_reads: int = 5, max_contam: float = 0.1,
    p_match: float = 0.90, p_none: float = 0.09, p_error: float = 0.01,
    alpha_prior: float = 1.0, beta_prior: float = 9.0,
    min_confidence: float = 0.95, min_noise_reads: int = 2,
    noise_cutoff_ratio: float = 0.05):

    with open(b_with_ids_filename, 'r') as b_with_ids_file, \
    open(_16s_packet_filename, 'r') as _16s_packet_file, \
    open(arg_packet_filename, 'r') as arg_packet_file, \
    open(primers_filename, 'r') as primers_file: # csv file

        print("Writing to files") # status update to user
        df = pd.DataFrame(columns = find_column_names(primers_filename))
        _16s_packet_file_content = _16s_packet_file.readlines()
        arg_packet_file_content = arg_packet_file.readlines()

        i = 0
        for line in b_with_ids_file:
            barcode, _16s_packet_list, arg_packet_list = get_barcode_packet_lists(line, _16s_packet_file_content, arg_packet_file_content)
            total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination \
            = parse_and_analyze_perfect_corrected_revised(_16s_packet_list, 
                p_match, p_none, p_error, alpha_prior, beta_prior,
                min_confidence, min_noise_reads, noise_cutoff_ratio)
            # 2026-08-10: Pass the primer-derived ARG count into the summary writer.
            # Reason: empty ARG barcodes still need a row with the correct number of columns.
            data_arg = find_arg_data(arg_packet_list, len(get_arg_names(primers_filename)))
            write_content_tsv_row(barcode, total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination, data_arg, df, i) # content rows

            print(f"Processed {i} barcodes")
            i += 1

        df.to_csv(unfiltered_tsv_filename, sep = "\t", index_label = "Barcode")
        filter_barcodes_in_df(df, min_16s_reads, max_contam)
        df.to_csv(tsv_filename, sep = "\t", index_label = "Barcode")

def find_column_names(primers_filename):
    primer_names = get_arg_names(primers_filename)
    column_names = ["Predicted taxonomy", "Confidence", "Contamination", "Total # of 16s reads", "Technical noise count"] + primer_names
    return column_names

def write_content_tsv_row(barcode: str, total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination, data_arg, df, i):

    row_16s = [predicted_taxonomy, str(confidence), str(contamination), str(total_16s_reads), str(technical_noise_count)]
    row_arg = []
    for val in data_arg:
        row_arg.append(str(val))

    full_row = row_16s + row_arg

    df.loc[barcode] = full_row

def extract_id_packet(wanted_id: str, packet_file_content) -> Dict:

    for line in packet_file_content:
        if wanted_id in line:
            packet = json.loads(line.strip("\n"))
            assert packet["ID"] == wanted_id
            return packet
                
    # packet should have been found and returned; if not, give assertion error
    assert False, f"ID {wanted_id} could not be found in packet file"
        
def form_packet_list(id_list: List[str], packet_file_content) -> List[Dict]:
    packet_list = []
    for id in id_list:
        if id != "":
            packet_list.append(extract_id_packet(id, packet_file_content))
    return packet_list

def get_barcode_packet_lists(
    b_with_ids_line: str,
    _16s_packet_file_content,
    arg_packet_file_content) -> Tuple[str, List[Dict], List[Dict]]:

    # parse barcode line
    bc, ids = b_with_ids_line.strip("\n").split(": ")
    _16s_ids, arg_ids = ids.split(" | ")[0:2]
    _16s_ids = _16s_ids.split(", ")
    arg_ids = arg_ids.split(", ")

    _16s_packet_list = form_packet_list(_16s_ids, _16s_packet_file_content)
    arg_packet_list = form_packet_list(arg_ids, arg_packet_file_content)

    return bc, _16s_packet_list, arg_packet_list


def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    parser.add_argument("--b_with_ids_filename", type=str, required=True)
    parser.add_argument("--_16s_packet_filename", type=str, required=True)
    parser.add_argument("--arg_packet_filename", type=str, required=True)
    # 2026-08-10: Route barcode summary intermediates to tmp by default.
    # Reason: the final user-facing product is exported separately as JSONL.
    parser.add_argument("--unfiltered_tsv_filename", type=str, default = "tmp/unfiltered_barcode_summary.tsv")
    parser.add_argument("--tsv_filename", type=str, default = "tmp/barcode_summary.tsv")
    parser.add_argument("--primers_filename", type=str, required=True)
    parser.add_argument("--min_16s_reads", type=int, default = 5)
    parser.add_argument("--max_contam", type=float, default = 0.1)

    parser.add_argument("--p_match", type=float, default=0.90)
    parser.add_argument("--p_none", type=float, default=0.09)
    parser.add_argument("--p_error", type=float, default=0.01)
    parser.add_argument("--alpha_prior", type=float, default=1.0)
    parser.add_argument("--beta_prior", type=float, default=9.0)
    parser.add_argument("--min_confidence", type=float, default=0.95)
    parser.add_argument("--min_noise_reads", type=int, default=2)
    parser.add_argument("--noise_cutoff_ratio", type=float, default=0.05)

    args = parser.parse_args()

    # 2026-08-10: Materialize the temporary output directory before writing summaries.
    # Reason: the default tmp paths must work in a new result directory.
    ensure_output_directories(args.unfiltered_tsv_filename, args.tsv_filename)
    
    # make sure input file paths exist
    if not os.path.exists(args.b_with_ids_filename):
        print(f"❌ Error: input file not found: {args.b_with_ids_filename}")
        return
    if not os.path.exists(args._16s_packet_filename):
        print(f"❌ Error: input file not found: {args._16s_packet_filename}")
        return
    if not os.path.exists(args.arg_packet_filename):
        print(f"❌ Error: input file not found: {args.arg_packet_filename}")
        return
    if not os.path.exists(args.primers_filename):
        print(f"❌ Error: input file not found: {args.primers_filename}")
        return

    write_barcode_summary_to_tsv(
        args.b_with_ids_filename, 
        args._16s_packet_filename, args.arg_packet_filename, 
        args.unfiltered_tsv_filename, args.tsv_filename, args.primers_filename,
        args.min_16s_reads, args.max_contam, 
        args.p_match, args.p_none, args.p_error, args.alpha_prior, args.beta_prior,
        args.min_confidence, args.min_noise_reads, args.noise_cutoff_ratio)

if __name__ == "__main__":
    main()
