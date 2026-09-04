#!/usr/bin/env python3
import json
import argparse
import sys
import os
import multiprocessing
from typing import List, Dict, Tuple

BARCODE_CHUNK_SIZE = 2048
_BARCODE_WORKER_CONFIG = None
import pandas as pd
from mle_revised import parse_and_analyze_perfect_corrected_revised
from filter_barcodes import filter_barcodes_in_df
from helper_functions import get_arg_names, ensure_output_directories
    
def find_arg_data(packet_list: List[Dict], num_arg_genes: int):
    """Count the # of reads of each target gene, present in the given cell"""
    
    # 2026-08-10: Allocate counts from the primer panel, not a fixed 23-gene assumption.
    # Reason: the output vector must match the columns generated from the supplied primers file.
    gene_counts = [0] * num_arg_genes

    # Iterate through each packet (each representing one read) in the given cell, 
    # so as to count the # of reads of each target gene present
    for packet in packet_list:
        gene = packet["gene"]
        # assume gene is a 1D matrix, with a length corresponding to # of different possible ARGs
        # all numbers should be 0, except for one 1
        for i in range(len(gene)):
            if gene[i] == 1:
                gene_counts[i] += 1
                break # no other 1s to add, so break

    return gene_counts
    
# 2026-08-10: Added defaults to parameters that followed defaulted arguments. Reason: Python otherwise rejects this function before the pipeline can run.
# 2026-08-28: Added the configurable Stage 2 threshold at the end of the signature.
# Reason: preserve existing positional callers while exposing the new numeric control.

def _initialize_barcode_worker(packet_indexes, num_arg_genes, analysis_params):
    """Initialize read-only packet indexes and analysis settings once per worker."""
    global _BARCODE_WORKER_CONFIG
    _BARCODE_WORKER_CONFIG = (packet_indexes, num_arg_genes, analysis_params)


def _process_barcode_chunk(lines):
    """Analyze one ordered barcode chunk without writing output files."""
    (packet_indexes, num_arg_genes, analysis_params) = _BARCODE_WORKER_CONFIG
    _16s_packet_index, arg_packet_index = packet_indexes
    rows = {}
    for i, line in enumerate(lines):
        barcode, _16s_packet_list, arg_packet_list = get_barcode_packet_lists(
            line, _16s_packet_index, arg_packet_index)
        (p_match, p_none, p_error, alpha_prior, beta_prior,
         min_confidence, min_noise_reads, noise_cutoff_ratio) = analysis_params
        result = parse_and_analyze_perfect_corrected_revised(
            _16s_packet_list, p_match, p_none, p_error, alpha_prior, beta_prior,
            min_confidence, min_noise_reads, noise_cutoff_ratio)
        total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination = result
        data_arg = find_arg_data(arg_packet_list, num_arg_genes)
        write_content_tsv_row(
            barcode, total_16s_reads, technical_noise_count, predicted_taxonomy,
            confidence, contamination, data_arg, rows, i,
        )
    return list(rows.items())


def _iter_barcode_chunks(barcode_filename, chunk_size):
    """Yield contiguous barcode chunks in original file order."""
    with open(barcode_filename, 'r') as barcode_file:
        chunk = []
        for line in barcode_file:
            chunk.append(line)
            if len(chunk) == chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk


def _write_barcode_summary_parallel(
    b_with_ids_filename, _16s_packet_index, arg_packet_index,
    unfiltered_tsv_filename, tsv_filename, primers_filename,
    min_16s_reads, max_contam, p_match, p_none, p_error,
    alpha_prior, beta_prior, min_confidence, min_noise_reads,
    noise_cutoff_ratio, min_barcodes, num_arg_genes, analysis_workers,
):
    """Run barcode analysis in workers and finalize output in the parent."""
    analysis_params = (
        p_match, p_none, p_error, alpha_prior, beta_prior,
        min_confidence, min_noise_reads, noise_cutoff_ratio,
    )
    with multiprocessing.Pool(
        processes=analysis_workers,
        initializer=_initialize_barcode_worker,
        initargs=([_16s_packet_index, arg_packet_index], num_arg_genes, analysis_params),
    ) as pool:
        rows = {}
        processed = 0
        for chunk_rows in pool.imap(
            _process_barcode_chunk,
            _iter_barcode_chunks(b_with_ids_filename, BARCODE_CHUNK_SIZE),
        ):
            rows.update(chunk_rows)
            processed += len(chunk_rows)
            if processed % 10000 == 0:
                print(f"Processed {processed} barcodes")

    df = pd.DataFrame.from_dict(rows, orient="index", columns=find_column_names(primers_filename))
    df.to_csv(unfiltered_tsv_filename, sep="\t", index_label="Barcode")
    filter_barcodes_in_df(df, min_16s_reads, max_contam, min_barcodes=min_barcodes)
    df.to_csv(tsv_filename, sep="\t", index_label="Barcode")


def write_barcode_summary_to_tsv(b_with_ids_filename: str, 
    _16s_packet_filename: str, arg_packet_filename: str, 
    unfiltered_tsv_filename: str, tsv_filename: str, primers_filename: str,
    min_16s_reads: int = 5, max_contam: float = 0.1,
    p_match: float = 0.90, p_none: float = 0.09, p_error: float = 0.01,
    alpha_prior: float = 1.0, beta_prior: float = 9.0,
    min_confidence: float = 0.95, min_noise_reads: int = 2,
    noise_cutoff_ratio: float = 0.05, min_barcodes: int = 10, analysis_workers: int = 1):
    """Obtain and compile all per-cell taxonomic & target gene count data into a single 'barcode summary', and write to a TSV file"""
        
    with open(b_with_ids_filename, 'r') as b_with_ids_file, \
    open(_16s_packet_filename, 'r') as _16s_packet_file, \
    open(arg_packet_filename, 'r') as arg_packet_file, \
    open(primers_filename, 'r') as primers_file: # csv file

        print("Writing to files") # status update to user
        # 2026-08-27: collect rows in a dictionary and build the DataFrame once.
        # Reason: assigning one row at a time with df.loc caused quadratic runtime.
        rows = {}
        _16s_packet_file_content = _16s_packet_file.readlines()
        arg_packet_file_content = arg_packet_file.readlines()
        # 2026-08-27: index packets once so each barcode lookup is constant-time.
        # Reason: rescanning every packet file for every read caused quadratic runtime.
        _16s_packet_index = build_packet_index(_16s_packet_file_content)
        arg_packet_index = build_packet_index(arg_packet_file_content)
        num_arg_genes = len(get_arg_names(primers_filename))
        if analysis_workers < 1:
            raise ValueError("analysis_workers must be at least 1")
        # 2026-09-04: Analyze barcode chunks in workers while keeping output in the parent.
        # Reason: the per-barcode analysis is independent, but output order must remain stable.
        if analysis_workers > 1:
            _write_barcode_summary_parallel(
                b_with_ids_filename, _16s_packet_index, arg_packet_index,
                unfiltered_tsv_filename, tsv_filename, primers_filename,
                min_16s_reads, max_contam, p_match, p_none, p_error,
                alpha_prior, beta_prior, min_confidence, min_noise_reads,
                noise_cutoff_ratio, min_barcodes, num_arg_genes, analysis_workers,
            )
            return

        # iterate through barcodes, and obtain & compile information for each barcode
        i = 0
        for line in b_with_ids_file:
            barcode, _16s_packet_list, arg_packet_list = get_barcode_packet_lists(line, _16s_packet_index, arg_packet_index)
            total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination \
            = parse_and_analyze_perfect_corrected_revised(_16s_packet_list, 
                p_match, p_none, p_error, alpha_prior, beta_prior,
                min_confidence, min_noise_reads, noise_cutoff_ratio)
            # 2026-08-10: Pass the primer-derived ARG count into the summary writer.
            # Reason: empty ARG barcodes still need a row with the correct number of columns.
            data_arg = find_arg_data(arg_packet_list, num_arg_genes)
            write_content_tsv_row(barcode, total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination, data_arg, rows, i) # content rows

            i += 1
            if i % 10000 == 0:
                print(f"Processed {i} barcodes")

        df = pd.DataFrame.from_dict(rows, orient="index", columns=find_column_names(primers_filename))
        # filter barcodes, and write barcode summary to tsv
        df.to_csv(unfiltered_tsv_filename, sep = "\t", index_label = "Barcode")
        filter_barcodes_in_df(df, min_16s_reads, max_contam, min_barcodes=min_barcodes) # filter barcodes
        df.to_csv(tsv_filename, sep = "\t", index_label = "Barcode")

def find_column_names(primers_filename):
    """Obtain list of column names for the barcode summary"""
    primer_names = get_arg_names(primers_filename)
    column_names = ["Predicted taxonomy", "Confidence", "Contamination", "Total # of 16s reads", "Technical noise count"] + primer_names
    return column_names

def write_content_tsv_row(barcode: str, total_16s_reads, technical_noise_count, predicted_taxonomy, confidence, contamination, data_arg, rows, i):
    """Compile provided taxonomic & target gene count data, for the given barcode, and write it into a single row in the barcode summary dataframe"""
    row_16s = [predicted_taxonomy, str(confidence), str(contamination), str(total_16s_reads), str(technical_noise_count)]
    row_arg = []
    for val in data_arg:
        row_arg.append(str(val))

    full_row = row_16s + row_arg

    rows[barcode] = full_row

def build_packet_index(packet_file_content) -> Dict:
    """Build a lookup table from read ID to packet."""
    return {packet["ID"]: packet for packet in (json.loads(line) for line in packet_file_content)}

def extract_id_packet(wanted_id: str, packet_index: Dict) -> Dict:
    """Return the packet for a read ID from a pre-built index."""
    if wanted_id not in packet_index:
        raise KeyError(f"ID {wanted_id} could not be found in packet file")
    return packet_index[wanted_id]
        
def form_packet_list(id_list: List[str], packet_index) -> List[Dict]:
    """Compile a list of 'packets' corresponding to the IDs of each of the desired sequencing reads"""
    packet_list = []
    for id in id_list:
        if id != "": # skip over empty IDs
            packet_list.append(extract_id_packet(id, packet_index))
    return packet_list

def get_barcode_packet_lists(
    b_with_ids_line: str,
    _16s_packet_index,
    arg_packet_index) -> Tuple[str, List[Dict], List[Dict]]:
    """Obtain the lists of read 'packets' for the given cell's 16s reads, and target gene reads"""
        
    # parse barcode line
    bc, ids = b_with_ids_line.strip("\n").split(": ")
    _16s_ids, arg_ids = ids.split(" | ")[0:2]
    _16s_ids = _16s_ids.split(", ")
    arg_ids = arg_ids.split(", ")

    # obtain the packet lists
    _16s_packet_list = form_packet_list(_16s_ids, _16s_packet_index)
    arg_packet_list = form_packet_list(arg_ids, arg_packet_index)

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
    # 2026-08-28: Expose the Stage 2 threshold while retaining the current default.
    # Reason: users can control low-count taxon filtering without a separate skip flag.
    parser.add_argument("--min_cells_per_taxon", type=int, default=10)

    parser.add_argument("--p_match", type=float, default=0.90)
    parser.add_argument("--p_none", type=float, default=0.09)
    parser.add_argument("--p_error", type=float, default=0.01)
    parser.add_argument("--alpha_prior", type=float, default=1.0)
    parser.add_argument("--beta_prior", type=float, default=9.0)
    parser.add_argument("--min_confidence", type=float, default=0.95)
    parser.add_argument("--min_noise_reads", type=int, default=2)
    parser.add_argument("--noise_cutoff_ratio", type=float, default=0.05)
    parser.add_argument("-@", "--threads", dest="analysis_workers", type=int, default=1, metavar="INT")

    args = parser.parse_args()

    if args.min_cells_per_taxon < 0:
        parser.error("--min_cells_per_taxon must be non-negative")

    # 2026-08-10: Materialize the temporary output directory before writing summaries.
    # Reason: the default tmp paths must work in a new result directory.
    ensure_output_directories(args.unfiltered_tsv_filename, args.tsv_filename)
    
    # make sure input file paths exist
    if not os.path.exists(args.b_with_ids_filename):
        print(f"❌ Error: input file not found: {args.b_with_ids_filename}")
        sys.exit(1)
    if not os.path.exists(args._16s_packet_filename):
        print(f"❌ Error: input file not found: {args._16s_packet_filename}")
        sys.exit(1)
    if not os.path.exists(args.arg_packet_filename):
        print(f"❌ Error: input file not found: {args.arg_packet_filename}")
        sys.exit(1)
    if not os.path.exists(args.primers_filename):
        print(f"❌ Error: input file not found: {args.primers_filename}")
        sys.exit(1)

    write_barcode_summary_to_tsv(
        args.b_with_ids_filename, 
        args._16s_packet_filename, args.arg_packet_filename, 
        args.unfiltered_tsv_filename, args.tsv_filename, args.primers_filename,
        args.min_16s_reads, args.max_contam, 
        args.p_match, args.p_none, args.p_error, args.alpha_prior, args.beta_prior,
        args.min_confidence, args.min_noise_reads, args.noise_cutoff_ratio,
        args.min_cells_per_taxon, args.analysis_workers)

if __name__ == "__main__":
    main()
