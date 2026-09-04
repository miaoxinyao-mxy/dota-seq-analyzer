#!/usr/bin/env python3
from taxonomy_tree import create_taxonomy_tree
from extract_16s_reads import check_primer_match_seq
import json
import os
import argparse
import sys
import multiprocessing
from typing import List, Dict
from helper_functions import open_maybe_gzip, ensure_output_directories

PRIMER_CHUNK_SIZE = 2048
_PRIMER_WORKER_CONFIG = None

"""
This program converts raw data from fastq and kraken files into 
json-like objects containing the barcode, gene, and taxonomic classification for each read ID
"""

def determine_16s_taxonomy(t_line: str, tax_nodes_lists):
    """
    Determines taxonomic classification of a read ID at each applicable taxonomic rank 
    (e.g. classifiable, domain, ... genus, species)
    Iterates through corresponding list of nodes for lowest identified rank, then traverses up tree
    Returns a dictionary including either a taxonomic name (e.g. "E. Coli") or "None", for each taxonomic rank
    """

    taxonomy = {"classifiable": None, 'R1': None, 'P': None, 'C': None, 'O': None, 'F': None, 'G': None, 'S': None}

    # parsing
    t_line_parsed = t_line.strip().split("\t")        
    tax_lvl_and_name = t_line_parsed[2].split(" (taxid")[0]

    # 2026-08-10: Handle Kraken2 taxids for unclassified and root records before lookup.
    # Reason: standard Kraken2 writes taxid 0/1 instead of the custom names used downstream.
    if tax_lvl_and_name in {"0", "unclassified"}:
        taxonomy["classifiable"] = False
        return taxonomy
    if tax_lvl_and_name in {"1", "root"}:
        taxonomy["classifiable"] = True
        return taxonomy

    # 2026-08-10: Resolve standard Kraken2 taxid output through the report lookup.
    # Reason: Kraken2 writes a numeric taxid here, not a rank/name string.
    if "__" not in tax_lvl_and_name:
        taxid_node = tax_nodes_lists.get("_taxid", {}).get(t_line_parsed[2].strip())
        if taxid_node is not None:
            taxonomy["classifiable"] = True
            set_all_tax_vals(taxid_node, taxonomy)
            return taxonomy

    # deal with "unclassified" and "root" edge cases
    if tax_lvl_and_name == "unclassified":
        taxonomy["classifiable"] = False
        return taxonomy
    elif tax_lvl_and_name == "root":
        taxonomy["classifiable"] = True
        return taxonomy

    # must be classified since not assigned "unclassified"
    taxonomy["classifiable"] = True

    # determine taxonomy level and name
    tax_lvl_char, tax_name = tax_lvl_and_name.split("__")
    tax_lvl_char = tax_lvl_char.capitalize()
    # account for D vs R1 representation for domain
    if tax_lvl_char == "D":
        tax_lvl_char = "R1"

    # iterate through list of nodes at the given tax_lvl, to find matching node
    for node in tax_nodes_lists[tax_lvl_char]:
        if tax_name == node.tax_name:
            set_all_tax_vals(node, taxonomy)
            return taxonomy
        
def set_all_tax_vals(bottom_node, taxonomy):
    """Traverses up tree, assigning all taxonomy properties appropriately"""

    curr_node = bottom_node
    while curr_node.tax_name != "root":
        taxonomy[curr_node.tax_lvl_char] = curr_node.tax_name
        curr_node = curr_node.parent

def determine_gene_revised(
    primers: List[str], 
    f_seq: str, r_seq: str, 
    max_shift: int, max_mm: int, 
    primer_start_num: int, 
    primers_to_genes_dict: Dict, primer_records=None):

    """Determines the gene for a given read, using a primer matching algorithm; output is either "16s", or a 1D matrix"""

    # 2026-08-10: Size the gene vector from the primer panel instead of 23.
    # Reason: custom primer panels must not silently produce misaligned gene columns.
    gene = [0] * max(0, len(primers) - 2)

    # first check for exact primer match
    for fwd_primer, rev_primer in primers_to_genes_dict:
        # 2026-08-10: Offset the R2 slice end by the expected primer start.
        # Reason: using len(rev_primer) as an absolute endpoint produces an empty or truncated exact-match slice.
        if rev_primer == r_seq[primer_start_num : primer_start_num + len(rev_primer)] and fwd_primer == f_seq[0 : len(fwd_primer)]:
            return primers_to_genes_dict[(fwd_primer, rev_primer)]

    # if no exact match, then use primer match algorithm to determine gene

    if primer_records is None:
        primer_records = make_primer_records(primers)

    for i, (fwd_primer, rev_primer, _) in enumerate(primer_records):

        # if both R1 and R2 primers match, then assign this read's gene accordingly
        if check_primer_match_seq(r_seq, rev_primer, max_shift, max_mm, primer_start_num) \
        and check_primer_match_seq(f_seq, fwd_primer, max_shift, max_mm):
            
            if i == 0: # 16s gene
                gene = "16s"
            else:
                gene[i-1] = 1
            break # matching primer already found, so exit loop

    return gene

def make_primer_records(primers):
    """Parse primer rows once while preserving their original order."""
    primer_records = []
    for i, primer_line in enumerate(primers[1:]):
        parsed_primer_line = primer_line.strip().split(",")
        fwd_primer = parsed_primer_line[1]
        rev_primer = parsed_primer_line[2]
        gene = [0] * max(0, len(primers) - 2)
        if i == 0:
            gene = "16s"
        else:
            gene[i - 1] = 1
        primer_records.append((fwd_primer, rev_primer, gene))
    return primer_records

def make_primers_to_genes_dict(primers):
    """
    Create a dictionary matching primer sequences to their corresponding gene,
    where the gene is either "16s" or a 1D vector representing the given target gene.
    """
    primers_to_genes = {}

    # iterate through all primers
    for i in range(len(primers[1:])): # don't need to use header line
        # parse the primer line
        primer_line = primers[i+1]
        parsed_primer_line = primer_line.strip().split(",")
        fwd_primer = parsed_primer_line[1]
        rev_primer = parsed_primer_line[2]

        # determine the gene
        # 2026-08-10: Size the gene vector from the primer panel instead of 23.
        # Reason: custom primer panels must not silently produce misaligned gene columns.
        gene = [0] * max(0, len(primers) - 2)
        if i == 0: # 16s gene
            gene = "16s"
        else:
            gene[i-1] = 1

        # add to the primers_to_genes dictionary
        primers_to_genes[(fwd_primer, rev_primer)] = gene

    return primers_to_genes

def format_packet(id, barcode, gene, taxonomy):
    """Serialize one packet as valid JSON."""

    # 2026-09-02: Use the JSON encoder instead of manual string replacement.
    # Reason: IDs, barcodes, and taxonomy names may contain characters that need JSON escaping.
    return json.dumps({
        "ID": id,
        "barcode": barcode,
        "gene": gene,
        "taxonomy": taxonomy,
    })


def _initialize_primer_worker(primers, primers_to_genes_dict, primer_records, max_shift, max_mm, primer_start_num):
    global _PRIMER_WORKER_CONFIG
    _PRIMER_WORKER_CONFIG = (primers, primers_to_genes_dict, primer_records, max_shift, max_mm, primer_start_num)


def _classify_primer_chunk(index_and_records):
    chunk_index, records = index_and_records
    primers, primers_to_genes_dict, primer_records, max_shift, max_mm, primer_start_num = _PRIMER_WORKER_CONFIG
    classified = []
    for record in records:
        gene = determine_gene_revised(primers, record[4], record[5], max_shift, max_mm, primer_start_num, primers_to_genes_dict, primer_records)
        classified.append((record, gene))
    return chunk_index, classified


def _read_primer_chunks(fwd_file, rev_file, chunk_size=PRIMER_CHUNK_SIZE):
    i = 0
    f_line = fwd_file.readline()
    r_line = rev_file.readline()
    while f_line != "" and r_line != "":
        records = []
        while len(records) < chunk_size and f_line != "" and r_line != "":
            id_f_unparsed = f_line.strip()
            id_r_unparsed = r_line.strip()
            id_f = id_f_unparsed.split(" ")[0].strip("@")
            id_r = id_r_unparsed.split(" ")[0].strip("@")
            assert id_f == id_r, f"ID from forward and reverse fastq files do not match on line {i*4 + 1}"
            f_seq = fwd_file.readline().strip()
            r_seq = rev_file.readline().strip()
            fwd_file.readline()
            rev_file.readline()
            fwd_file.readline()
            r_quality = rev_file.readline().strip()
            records.append((i, id_f_unparsed, id_r_unparsed, id_f, f_seq, r_seq, r_quality))
            f_line = fwd_file.readline()
            r_line = rev_file.readline()
            i += 1
        yield records


def _generate_packets_parallel(
    fwd_fastq_filename, rev_fastq_filename, kraken_output_filename, tax_nodes_lists,
    primers, primers_to_genes_dict, primer_records, barcode_len, max_shift, max_mm,
    primer_start_num, out_16s_packet_filename, out_arg_packet_filename,
    out_unclassified_packet_filename, out_arg_rev_fastq, out_unclassified_rev_fastq,
    analysis_workers):
    unmatched_gene = [0] * max(0, len(primers) - 2)
    with open_maybe_gzip(fwd_fastq_filename, "r") as fwd_file, open_maybe_gzip(rev_fastq_filename, "r") as rev_file, open(kraken_output_filename, "r") as taxonomy_file, open(out_16s_packet_filename, "w") as out_16s_packet_file, open(out_unclassified_packet_filename, "w") as out_unclassified_packet_file, open(out_arg_packet_filename, "w") as out_arg_packet_file, open(out_arg_rev_fastq, "w") as out_arg_rev_fastq_file, open(out_unclassified_rev_fastq, "w") as out_unclassified_rev_fastq_file:
        t_line = taxonomy_file.readline()
        chunks = enumerate(_read_primer_chunks(fwd_file, rev_file))
        with multiprocessing.Pool(processes=analysis_workers, initializer=_initialize_primer_worker, initargs=(primers, primers_to_genes_dict, primer_records, max_shift, max_mm, primer_start_num)) as pool:
            for _, classified in pool.imap(_classify_primer_chunk, chunks):
                for record, gene in classified:
                    i, id_f_unparsed, id_r_unparsed, id_f, f_seq, r_seq, r_quality = record
                    barcode = r_seq[0:barcode_len]
                    taxonomy = None
                    if gene == "16s":
                        assert t_line != "", f"Taxonomy file does not have another line to match the 16s read on line {i*4 + 1} of fastq file"
                        id_t = t_line.split("\t")[1]
                        id_r = id_r_unparsed.split(" ")[0].strip("@")
                        assert id_r == id_t, f"ID from fastq file and kraken output file do not match on line {i*4 + 1} of kraken output file"
                        taxonomy = determine_16s_taxonomy(t_line, tax_nodes_lists)
                        t_line = taxonomy_file.readline()
                        out_16s_packet_file.write(format_packet(id_f, barcode, gene, taxonomy) + "\n")
                    elif gene == unmatched_gene:
                        out_unclassified_packet_file.write(format_packet(id_f, barcode, gene, taxonomy) + "\n")
                        out_unclassified_rev_fastq_file.write(f"{id_r_unparsed}\n{r_seq}\n+\n{r_quality}\n")
                    else:
                        out_arg_packet_file.write(format_packet(id_f, barcode, gene, taxonomy) + "\n")
                        out_arg_rev_fastq_file.write(f"{id_r_unparsed}\n{r_seq}\n+\n{r_quality}\n")


def generate_packets(
    fwd_fastq_filename: str, rev_fastq_filename: str, primers_filename: str, 
    kraken_output_filename: str, report_filename: str,
    barcode_len: int, max_shift: int, max_mm: int, primer_start_num: int,
    out_16s_packet_filename: str, out_arg_packet_filename: str, out_unclassified_packet_filename: str,
    out_arg_rev_fastq: str, out_unclassified_rev_fastq: str, analysis_workers: int = 1
    ):
    """
    Generates json-like "packets" of data (including barcode, gene, and taxonomy) for each read ID.
    Takes raw fastq, kraken, and primer files as input.
    Writes these packets to three files, based on type of gene: 16s, unclassified gene, and ARGs.
    Also write reverse fastq files, again based on type of gene; to be used in downstream processing.
    """

    # preliminary steps
    tax_nodes_lists = create_taxonomy_tree(report_filename)
    with open(primers_filename, 'r') as primers_file: # csv file
        primers = primers_file.readlines()
    # 2026-09-04: Parse primer rows once before processing reads.
    # Reason: primer strings and gene mappings are invariant across reads.
    primer_records = make_primer_records(primers)
    primers_to_genes_dict = make_primers_to_genes_dict(primers)
    unmatched_gene = [0] * max(0, len(primers) - 2)
    if analysis_workers < 1:
        raise ValueError("analysis_workers must be at least 1")
    if analysis_workers > 1:
        _generate_packets_parallel(
            fwd_fastq_filename, rev_fastq_filename, kraken_output_filename, tax_nodes_lists,
            primers, primers_to_genes_dict, primer_records, barcode_len, max_shift, max_mm,
            primer_start_num, out_16s_packet_filename, out_arg_packet_filename,
            out_unclassified_packet_filename, out_arg_rev_fastq, out_unclassified_rev_fastq,
            analysis_workers)
        return

        
    with open_maybe_gzip(fwd_fastq_filename, 'r') as fwd_file, \
    open_maybe_gzip(rev_fastq_filename, 'r') as rev_file, \
    open(kraken_output_filename, 'r') as taxonomy_file, \
    open(out_16s_packet_filename, 'w') as out_16s_packet_file, \
    open(out_unclassified_packet_filename, 'w') as out_unclassified_packet_file, \
    open(out_arg_packet_filename, 'w') as out_arg_packet_file, \
    open(out_arg_rev_fastq, 'w') as out_arg_rev_fastq_file, \
    open(out_unclassified_rev_fastq, 'w') as out_unclassified_rev_fastq_file:

        # read in the read_ID line
        i = 0
        f_line = fwd_file.readline()
        r_line = rev_file.readline()
        t_line = taxonomy_file.readline()

        # iterate through all lines in the fastq files
        while (f_line != "") and (r_line != ""):

            # parse lines to get ID from barcode file & taxonomy file
            id_f_unparsed = f_line.strip()
            id_r_unparsed = r_line.strip()
            id_f = id_f_unparsed.split(" ")[0].strip("@")
            id_r = id_r_unparsed.split(" ")[0].strip("@")

            # check that IDs for this line are equal - necessary for this analysis
            assert id_f == id_r, f"ID from forward and reverse fastq files do not match on line {i*4 + 1}"

            # parsing
            f_seq = fwd_file.readline().strip()
            r_seq = rev_file.readline().strip()
            for _ in range(2):
                fwd_file.readline()
                r_quality = rev_file.readline().strip()

            # assign values for read ID basic info
            barcode = r_seq[0:barcode_len]
            gene = determine_gene_revised(primers, f_seq, r_seq, max_shift, max_mm, primer_start_num, primers_to_genes_dict, primer_records)
            taxonomy = None
            
            # create packet and append to appropriate list, based on general group of gene
            if gene == "16s":
                # determine taxonomy for 16s reads
                # make sure there is another line of info in taxonomy file
                assert t_line != "", f"Taxonomy file does not have another line to match the 16s read on line {i*4 + 1} of fastq file"
                # check that IDs for this line are equal - necessary for this analysis
                id_t = t_line.split("\t")[1]
                assert id_r == id_t, f"ID from fastq file and kraken output file do not match on line {i*4 + 1} of kraken output file"
                # determine the taxonomy
                taxonomy = determine_16s_taxonomy(t_line, tax_nodes_lists)
                # prepare for next 16s read
                t_line = taxonomy_file.readline()
                
                # note that id_f, id_r, & id_t can be used interchangeably, because we already asserted that they're all equal
                out_16s_packet_file.write(format_packet(id_f, barcode, gene, taxonomy) + "\n")

            # 2026-08-10: Derive the unmatched vector from the input primer panel.
            # Reason: both vector content and length must match the supplied number of ARG targets.
            elif gene == unmatched_gene:
                out_unclassified_packet_file.write(format_packet(id_f, barcode, gene, taxonomy) + "\n")
                out_unclassified_rev_fastq_file.write(f"{id_r_unparsed}\n{r_seq}\n+\n{r_quality}\n")

            else: # arg genes
                out_arg_packet_file.write(format_packet(id_f, barcode, gene, taxonomy) + "\n")
                out_arg_rev_fastq_file.write(f"{id_r_unparsed}\n{r_seq}\n+\n{r_quality}\n")

            # prepare for reading in the next read_ID line
            f_line = fwd_file.readline()
            r_line = rev_file.readline()

            i += 1
            if i % 1000 == 0:
                print(f"Processed {i//1000},000 reads")
    
def main():
    parser = argparse.ArgumentParser()

    # take input parameters
    # 2026-08-10: Expose paired-read CLI inputs as R1 and R2.
    # Reason: users provide sequencing reads as R1/R2 files.
    parser.add_argument("--r1_fastq", type=str, required=True)
    parser.add_argument("--r2_fastq", type=str, required=True)
    parser.add_argument("--primers_filename", type=str, required=True)
    parser.add_argument("--kraken_output", type=str, required=True)
    parser.add_argument("--kraken_report", type=str, required=True)
    parser.add_argument("--max_shift_primer", type=int, default=4)
    parser.add_argument("--max_mm_primer", type=int, default=4)
    parser.add_argument("--primer_start_num", type=int, default=42)
    parser.add_argument("--barcode_len", type=int, default=20)
    parser.add_argument("-@", "--threads", dest="analysis_workers", type=int, default=1, metavar="INT", help="Number of parallel workers used by DoTA-seq analysis")
    # 2026-09-04: Reject invalid worker counts before starting the pipeline.
    # Reason: a non-positive process count cannot create a valid worker pool.
    # 2026-08-10: Route packet and derived R2 files to tmp by default.
    # Reason: these files are regenerable pipeline intermediates.
    parser.add_argument("--_16s_packet_filename", type=str, default="tmp/packets_16s")
    parser.add_argument("--arg_packet_filename", type=str, default="tmp/packets_arg")
    parser.add_argument("--unclassified_packet_filename", type=str, default="tmp/packets_unclassified")
    parser.add_argument("--arg_r2_fastq", type=str, default="tmp/arg_R2.fastq")
    parser.add_argument("--unclassified_r2_fastq", type=str, default="tmp/unclassified_R2.fastq")

    args = parser.parse_args()
    if args.analysis_workers < 1:
        parser.error("--threads must be at least 1")

    # 2026-08-10: Materialize the organized output directories before writing files.
    # Reason: default tmp paths must work in a new result directory.
    ensure_output_directories(
        args._16s_packet_filename, args.arg_packet_filename,
        args.unclassified_packet_filename, args.arg_r2_fastq,
        args.unclassified_r2_fastq)
    
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
    if not os.path.exists(args.kraken_output):
        print(f"❌ Error: input file not found: {args.kraken_output}")
        sys.exit(1)
    if not os.path.exists(args.kraken_report):
        print(f"❌ Error: input file not found: {args.kraken_report}")
        sys.exit(1)

    generate_packets(
        args.r1_fastq, args.r2_fastq, args.primers_filename,
        args.kraken_output, args.kraken_report,
        args.barcode_len, args.max_shift_primer, args.max_mm_primer, args.primer_start_num,
        args._16s_packet_filename, args.arg_packet_filename, args.unclassified_packet_filename,
        args.arg_r2_fastq, args.unclassified_r2_fastq, args.analysis_workers
    )

if __name__ == "__main__":
    main()
